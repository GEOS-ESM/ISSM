/*!\file:  esmfbinders.cpp
 * \brief: Binders developed for NASA's GEOS Earth System Model.
 *
 * Layout of every GEOS <-> ISSM exchange array (forcings, restarts, outputs):
 *
 *     array[ f*total_nodes + shift + lid ]
 *
 *   f           index into the field list (GCMForcingTerms / ISSMOutputTerms)
 *   total_nodes number of local vertices on this process summed over all models
 *   shift       number of local vertices of all preceding models (glaciers)
 *   lid         local id of the vertex within its own model on this process
 *
 * Every loop over models must therefore advance `shift` by the model's local
 * vertex count once it is done with that model.
 */

#include "./issm.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

/*
 * GEOS <-> ISSM FIELD CONTRACT  (keep in sync with the GEOS-side field definitions)
 *
 * The number AND ORDER of fields below must match the Fortran side exactly. Nothing
 * checks this at runtime: a count mismatch overruns the exchange arrays, and an order
 * mismatch silently hands GEOS the wrong field.
 *
 * Forcings, GEOS -> ISSM (gcm_forcings), field index f:
 *   0  SMBgcmEnum                surface mass balance
 *
 * Outputs / restarts, ISSM <-> GEOS (issm_outputs, gcm_restarts), field index f:
 *   0  SurfaceEnum               ice surface elevation
 *   1  ThicknessEnum             ice thickness
 *   2  VxEnum                    ice velocity, x
 *   3  VyEnum                    ice velocity, y
 *   4  MaskOceanLevelsetEnum     ocean mask
 *   5  MaskIceLevelsetEnum       ice mask
 *   6  FluxDivergenceEnum        div(H v)   (output only; ignored by InputFromRestarts)
 *
 * To add a field: append it at the END of the list, bump the count, and update the
 * Fortran side (GEOS_ISSMGridComp.F90). Never reorder or insert in the middle. =)
 */

/*GEOS 5 specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum};
const int ISSMOutputNumTerms = 7;
const int ISSMOutputTerms[ISSMOutputNumTerms]= {SurfaceEnum,ThicknessEnum,VxEnum,VyEnum,MaskOceanLevelsetEnum,MaskIceLevelsetEnum,FluxDivergenceEnum};

static int N = 0;                      // total number of input files (glaciers)
static int total_nodes = 0;            // total number of vertices across models on each process
static FemModel** femmodels = nullptr; // container for all FemModel objects

/* Copy one field of a GEOS array onto the vertices of an element. */
static void GetElementVertexValues(Element* element, const IssmDouble* field, int shift_local_vertices, std::vector<IssmDouble>& values){
    const int numvertices = element->GetNumberOfVertices();
    values.resize(numvertices);
    for (int iv = 0; iv < numvertices; ++iv){
        values[iv] = field[element->vertices[iv]->Lid() + shift_local_vertices];
    }
}

extern "C" {

    void FinalizeISSM(); // defined below; used by the re-initialization guard

    void InitializeISSM(const char* EXPDIR, int* ptotal_elements, int* ptotal_nodes, MPI_Fint* Fcomm){
        /*  Determine number of input (ISSM*.bin) files and create a FemModel */
        /*  for each one. Return total number of elements and nodes per       */
        /*  process across all models                                         */

        /* Guard against re-initialization without FinalizeISSM (would leak every FemModel) */
        if (femmodels) FinalizeISSM();

        int total_elements = 0; // #elements on process across all models
        total_nodes = 0;        // #vertices on process across all models (static)

        /* Convert Fortran MPI comm to C MPI comm */
        MPI_Comm Ccomm = MPI_Comm_f2c(*Fcomm);

        /* Scan directory for ISSM*.bin files (error_code overloads: report, don't terminate) */
        std::vector<std::string> binfiles;
        std::error_code ec;
        std::filesystem::directory_iterator it(EXPDIR, ec), end;
        if (ec) _error_("Cannot read experiment directory " << EXPDIR << ": " << ec.message() << "\n");

        for (; !ec && it != end; it.increment(ec)) {
            if (!it->is_regular_file(ec)) { ec.clear(); continue; }

            const std::filesystem::path path = it->path();
            const std::string stem = path.stem().string(); // file name without ".bin"

            // Match files like ISSM*.bin
            if (path.extension() == ".bin" && stem.rfind("ISSM", 0) == 0) {
                binfiles.push_back(stem);
            }
        }
        if (ec) _error_("Error while scanning experiment directory " << EXPDIR << ": " << ec.message() << "\n");

        std::sort(binfiles.begin(), binfiles.end()); // same model order on every process

        N = static_cast<int>(binfiles.size()); // total number of input files

        if (N == 0) {
            *ptotal_elements = 0;
            *ptotal_nodes    = 0;
            return;
        }

        /* Allocate FemModel pointer array (value-initialized to nullptr, so a   */
        /* failure part-way through the loop below still leaves a valid array)   */
        femmodels = new FemModel*[N]();

        /* Initialize each model and accumulate global sizes */
        for (int id = 0; id < N; ++id) {
            // set command line arguments:
            std::string solution = "TransientSolution";
            std::string expdir   = EXPDIR;
            std::string filename = binfiles[id];  // no .bin extension

            int argc = 4;
            char* argv[4];

            argv[0] = const_cast<char*>("issm");
            argv[1] = const_cast<char*>(solution.c_str());
            argv[2] = const_cast<char*>(expdir.c_str());
            argv[3] = const_cast<char*>(filename.c_str());

            // initialize model from command line arguments
            femmodels[id] = new FemModel(argc, argv, Ccomm);

            // add this model's local element/vertex counts to the totals
            total_elements += femmodels[id]->elements->Size();
            total_nodes    += femmodels[id]->vertices->Size();

            // specify that ISSM will be getting SMB from a GCM
            femmodels[id]->parameters->SetParam(SMBgcmEnum, SmbEnum);
            femmodels[id]->Restart();
        }

        /* Return (#elements, #nodes) across all models on each process */
        *ptotal_elements = total_elements;
        *ptotal_nodes    = total_nodes;
    }

    void InputFromRestarts(IssmDouble* gcm_restarts){
        /* Set inputs based on GEOS restarts. */

        // index of the surface elevation restart (needed to derive base elevation)
        int f0 = -1;
        for (int f = 0; f < ISSMOutputNumTerms; ++f) {
            if (ISSMOutputTerms[f] == SurfaceEnum) f0 = f;
        }
        _assert_(f0 >= 0);

        int shift_local_vertices = 0; // starting position of each model's vertices

        for (int id = 0; id < N; ++id) {
            /*Figure out number of elements and vertices local to process: */
            const int local_size_elements = femmodels[id]->elements->Size();
            const int local_size_vertices = femmodels[id]->vertices->Size();

            for (int f = 0; f < ISSMOutputNumTerms; ++f) {
                const int output_type = ISSMOutputTerms[f];
                const IssmDouble* field = gcm_restarts + f*total_nodes;

                for (int i = 0; i < local_size_elements; ++i) {
                    Element* element = xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));

                    std::vector<IssmDouble> values; // this field at the element's vertices
                    GetElementVertexValues(element, field, shift_local_vertices, values);

                    switch (output_type) {

                        case SurfaceEnum:
                            // ice surface elevation restart provided by GEOS
                            element->AddInput(SurfaceEnum,    values.data(), P1Enum);
                            element->AddInput(SurfaceOldEnum, values.data(), P1Enum);
                            break;

                        case VxEnum:
                            // ice velocity (x) restart provided by GEOS
                            element->AddInput(VxEnum,        values.data(), P1Enum);
                            element->AddInput(VxSurfaceEnum, values.data(), P1Enum);
                            element->AddInput(VxBaseEnum,    values.data(), P1Enum);
                            break;

                        case VyEnum:
                            // ice velocity (y) restart provided by GEOS
                            element->AddInput(VyEnum,        values.data(), P1Enum);
                            element->AddInput(VySurfaceEnum, values.data(), P1Enum);
                            element->AddInput(VyBaseEnum,    values.data(), P1Enum);
                            break;

                        case ThicknessEnum: {
                            // ice thickness restart provided by GEOS; also set base
                            // elevation and pressure restarts for consistency
                            const IssmDouble rho_ice = element->FindParam(MaterialsRhoIceEnum);
                            const IssmDouble g       = element->FindParam(ConstantsGEnum);

                            std::vector<IssmDouble> surface;
                            GetElementVertexValues(element, gcm_restarts + f0*total_nodes, shift_local_vertices, surface);

                            std::vector<IssmDouble> base(values.size()), pressure(values.size());
                            for (size_t iv = 0; iv < values.size(); ++iv) {
                                base[iv]     = surface[iv] - values[iv];
                                pressure[iv] = rho_ice*g*values[iv];
                            }

                            element->AddInput(ThicknessEnum,    values.data(),   P1Enum);
                            element->AddInput(ThicknessOldEnum, values.data(),   P1Enum);
                            element->AddInput(BaseEnum,         base.data(),     P1Enum);
                            element->AddInput(BaseOldEnum,      base.data(),     P1Enum);
                            element->AddInput(PressureEnum,     pressure.data(), P1Enum);
                            break;
                        }

                        case MaskOceanLevelsetEnum:
                            // ocean mask restart provided by GEOS
                            element->AddInput(MaskOceanLevelsetEnum, values.data(), P1Enum);
                            break;

                        case MaskIceLevelsetEnum:
                            // ice mask restart provided by GEOS
                            element->AddInput(MaskIceLevelsetEnum, values.data(), P1Enum);
                            break;

                        case FluxDivergenceEnum:
                            // output only: nothing to restore
                            break;

                        default:
                            _error_("Unknown restart type " << output_type << "\n");
                            break;
                    }
                }
            }

            // shift starting index for next model inputs
            shift_local_vertices += local_size_vertices;
        }
    }

    void RunISSM(IssmDouble dt, IssmDouble* gcm_forcings, IssmDouble* issm_outputs){
        /* Run ISSM for one time step given GEOS forcing. */
        /* Time step is set by GEOS.                      */
        /* Return a collection of ISSM outputs.           */

        int shift_local_vertices = 0; // starting position of each model's vertices

        for (int id = 0; id < N; ++id) {
            /*Figure out number of elements and vertices local to process: */
            const int local_size_vertices = femmodels[id]->vertices->Size();
            const int local_size_elements = femmodels[id]->elements->Size();

            /*Setup GCM forcings as element-wise input: */
            for (int f = 0; f < GCMForcingNumTerms; ++f) {
                const int forcing_type = GCMForcingTerms[f];
                const IssmDouble* field = gcm_forcings + f*total_nodes;

                for (int i = 0; i < local_size_elements; ++i) {
                    Element* element = xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));

                    switch (forcing_type) {
                        case SMBgcmEnum: {
                            std::vector<IssmDouble> smb_forcing;
                            GetElementVertexValues(element, field, shift_local_vertices, smb_forcing);
                            element->AddInput(SmbMassBalanceEnum, smb_forcing.data(), P1Enum);
                            break;
                        }
                        default:
                            _error_("Unknown forcing type " << forcing_type << "\n");
                            break;
                    }
                }
            }

            /*Before running, set the time interval to one time step: */
            IssmDouble start_time, final_time;
            femmodels[id]->parameters->FindParam(&start_time, TimeEnum);
            final_time = start_time + dt;
            femmodels[id]->parameters->SetParam(final_time, TimesteppingFinalTimeEnum);
            femmodels[id]->parameters->SetParam(dt, TimesteppingTimeStepEnum);

            /*Now, run: */
            femmodels[id]->Solve();

            /*Retrieve ISSM outputs and pass them back to the GCM: */
            for (int f = 0; f < ISSMOutputNumTerms; ++f) {
                const int output_type = ISSMOutputTerms[f];

                /*
                 * Fields that are continuous (P1) take the same value from every
                 * element sharing a vertex, so they are written directly. Flux
                 * divergence involves derivatives, which can differ between
                 * adjacent elements at a shared vertex, so it is accumulated per
                 * global vertex id inside the element loop and averaged afterwards.
                 */
                std::unique_ptr<Vector<IssmDouble>> fluxdiv_sum, fluxdiv_count;
                if (output_type == FluxDivergenceEnum) {
                    const int global_num_vertices = femmodels[id]->vertices->NumberOfVertices();
                    fluxdiv_sum.reset(new Vector<IssmDouble>(global_num_vertices));
                    fluxdiv_count.reset(new Vector<IssmDouble>(global_num_vertices));
                }

                for (int i = 0; i < local_size_elements; ++i) {
                    Element* element = xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));
                    const int numvertices = element->GetNumberOfVertices();

                    switch (output_type) {

                        case SurfaceEnum:
                        case ThicknessEnum:
                        case VxEnum:
                        case VyEnum:
                        case MaskOceanLevelsetEnum:
                        case MaskIceLevelsetEnum: {
                            Input* input = element->GetInput(output_type);
                            _assert_(input);

                            std::unique_ptr<Gauss> gauss(element->NewGauss());
                            for (int iv = 0; iv < numvertices; ++iv) {
                                IssmDouble value;
                                gauss->GaussVertex(iv);
                                input->GetInputValue(&value, gauss.get());

                                const int lid = element->vertices[iv]->Lid() + shift_local_vertices;
                                issm_outputs[lid + f*total_nodes] = value;
                            }
                            break;
                        }

                        case FluxDivergenceEnum: {
                            _assert_(numvertices == 3); // 2D (triangle) meshes only: dH, dvx, dvy hold 2 components

                            Input* H_input = element->GetInput(ThicknessEnum);
                            _assert_(H_input);
                            Input* vx_input = element->GetInput(VxEnum);
                            _assert_(vx_input);
                            Input* vy_input = element->GetInput(VyEnum);
                            _assert_(vy_input);

                            std::unique_ptr<Gauss> gauss(element->NewGauss());
                            IssmDouble* xyz_list = NULL;
                            element->GetVerticesCoordinates(&xyz_list);

                            for (int iv = 0; iv < numvertices; ++iv) {
                                IssmDouble H, vx, vy;
                                IssmDouble dH[2], dvx[2], dvy[2];
                                gauss->GaussVertex(iv);

                                H_input->GetInputValue(&H, gauss.get());
                                H_input->GetInputDerivativeValue(&dH[0], xyz_list, gauss.get());

                                vx_input->GetInputValue(&vx, gauss.get());
                                vy_input->GetInputValue(&vy, gauss.get());
                                vx_input->GetInputDerivativeValue(&dvx[0], xyz_list, gauss.get());
                                vy_input->GetInputDerivativeValue(&dvy[0], xyz_list, gauss.get());

                                // div(H v) = v . grad(H) + H div(v)
                                const IssmDouble fluxdiv = vx*dH[0] + vy*dH[1] + H*(dvx[0] + dvy[1]);

                                // Accumulate this element's contribution by global vertex id.
                                const int sid = element->vertices[iv]->Sid();
                                fluxdiv_sum->SetValue(sid, fluxdiv, ADD_VAL);
                                fluxdiv_count->SetValue(sid, 1., ADD_VAL);
                            }

                            xDelete<IssmDouble>(xyz_list);
                            break;
                        }

                        default:
                            _error_("Unknown output type " << output_type << "\n");
                            break;
                    }
                }

                /*
                 * Element loop done: assemble the per-vertex sums and counts across
                 * processes, average, and write once per local vertex.
                 */
                if (output_type == FluxDivergenceEnum) {
                    fluxdiv_sum->Assemble();
                    fluxdiv_count->Assemble();

                    IssmDouble* sums   = fluxdiv_sum->ToMPISerial();
                    IssmDouble* counts = fluxdiv_count->ToMPISerial();

                    for (int i = 0; i < local_size_vertices; ++i) {
                        Vertex* vertex = xDynamicCast<Vertex*>(femmodels[id]->vertices->GetObjectByOffset(i));

                        const int sid = vertex->Sid();
                        const int lid = vertex->Lid() + shift_local_vertices;

                        issm_outputs[lid + f*total_nodes] = (counts[sid] > 0.) ? sums[sid]/counts[sid] : 0.;
                    }

                    xDelete<IssmDouble>(sums);
                    xDelete<IssmDouble>(counts);
                }
            }

            /*For the next time around, save the final time as start time */
            femmodels[id]->parameters->SetParam(final_time, TimesteppingStartTimeEnum);

            // shift starting index for next model inputs/outputs
            shift_local_vertices += local_size_vertices;
        }
    }

    void FinalizeISSM(){
        // wrap up
        for (int i = 0; i < N; ++i) {
            delete femmodels[i];
            femmodels[i] = NULL;
        }
        delete[] femmodels;
        femmodels   = nullptr; // avoid dangling pointer / double free on repeated calls
        N           = 0;
        total_nodes = 0;
    }

    void GetNodesISSM(int* nodeIds, IssmDouble* nodeCoords){
        /*obtain nodes of mesh for creating ESMF version in Fortran interface */
        /*nodeIds are the global Id's of the nodes and nodeCoords are the     */
        /*(lon,lat) coordinates, as described in the ESMF reference document  */
        int shift_global_vertices = 0;
        int shift_local_vertices  = 0;

        for (int id = 0; id < N; ++id) {
            // #vertices local to process in this model
            const int local_size_vertices = femmodels[id]->vertices->Size();

            // #vertices across all processes in this model
            const int global_size_vertices = femmodels[id]->vertices->NumberOfVertices();

            for (int i = 0; i < local_size_vertices; ++i) {
                Vertex* vertex = xDynamicCast<Vertex*>(femmodels[id]->vertices->GetObjectByOffset(i));
                const int i0 = vertex->Lid() + shift_local_vertices;
                nodeIds[i0]          = vertex->Sid() + 1 + shift_global_vertices;
                nodeCoords[2*i0 + 0] = vertex->longitude;
                nodeCoords[2*i0 + 1] = vertex->latitude;
            }
            // shift for next model starting indices:
            shift_local_vertices  += local_size_vertices;
            shift_global_vertices += global_size_vertices;
        }
    }

    void GetElementsISSM(int* elementIds, int* elementConn, IssmDouble* elementCoords, int* glacIds){
        /*obtain elements of mesh for creating ESMF version in Fortran interface*/
        /*Element connectivity (elementConn) contains the indices of the nodes  */
        /*that form the element as described in the ESMF reference document     */
        int shift_local_elements  = 0;
        int shift_global_elements = 0;
        int shift_local_vertices  = 0;

        for (int id = 0; id < N; ++id) {
            // #elements local to process in this model
            const int local_size_elements = femmodels[id]->elements->Size();

            // #elements across all processes in this model
            const int global_size_elements = femmodels[id]->elements->NumberOfElements();

            // #vertices local to process in this model
            const int local_size_vertices = femmodels[id]->vertices->Size();

            for (int i = 0; i < local_size_elements; ++i) {
                Element* element = xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));
                _assert_(element->GetNumberOfVertices() == 3); // triangles only

                const int i0 = i + shift_local_elements;
                elementIds[i0] = element->Sid() + 1 + shift_global_elements;
                glacIds[i0]    = id;
                for (int j = 0; j < 3; ++j) {
                    elementConn[i0*3 + j] = element->vertices[j]->Lid() + 1 + shift_local_vertices;
                }

                // Triangle centroid in longitude/latitude: average on the unit sphere
                // (3D unit vectors), not lon/lat directly, so triangles that cross the
                // dateline / prime meridian or contain a pole are handled correctly.
                const IssmDouble d2r = PI/180.0;
                IssmDouble cx = 0, cy = 0, cz = 0;
                bool neg = false; // any negative longitude => mesh uses [-180,180]
                for (int j = 0; j < 3; ++j) {
                    const IssmDouble lo = element->vertices[j]->longitude*d2r;
                    const IssmDouble la = element->vertices[j]->latitude*d2r;
                    neg = neg || (lo < 0);
                    cx += cos(la)*cos(lo);
                    cy += cos(la)*sin(lo);
                    cz += sin(la);
                }
                IssmDouble centroid_lon = atan2(cy, cx)/d2r;
                IssmDouble centroid_lat = atan2(cz, sqrt(cx*cx + cy*cy))/d2r;
                if (!neg) centroid_lon = fmod(centroid_lon + 360.0, 360.0); // keep [0,360) convention

                elementCoords[2*i0 + 0] = centroid_lon;
                elementCoords[2*i0 + 1] = centroid_lat;
            }
            // shift for next model starting indices:
            shift_local_vertices  += local_size_vertices;
            shift_local_elements  += local_size_elements;
            shift_global_elements += global_size_elements;
        }
    }
}