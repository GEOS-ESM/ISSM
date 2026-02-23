/*!\file:  esmfbinder.cpp
 * \brief: Binders developed for NASA's GEOS Earth System Model.
 */ 

#include "./issm.h"
#include <filesystem>
#include <vector>       
#include <string>      
#include <algorithm>    

/*GEOS 5 specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum}; 
const int ISSMOutputNumTerms = 3;
const int ISSMOutputTerms[ISSMOutputNumTerms]= { SurfaceEnum, ThicknessEnum, VelEnum };

extern "C" {


	static int N = 0; // total number of input files (glaciers)
    static FemModel** femmodels = nullptr; //container for all FemModel objects

    void InitializeISSM(const char* EXPDIR,
                        int* ptotal_elements,
                        int* ptotal_nodes,
                        MPI_Fint* Fcomm)
    {
        /*  Determine number of input (.bin) files and create a FemModel */
        /*  for each one. Return total number of elements and nodes per  */
        /*  process across all models                                    */
        
        int total_elements = 0; // #elements on process across all models
        int total_nodes    = 0; // #nodes on process across all models
    
        /* Convert Fortran MPI comm to C MPI comm */
        MPI_Comm Ccomm = MPI_Comm_f2c(*Fcomm);
    
        /* 1) Scan directory for .bin files                   */
        std::vector<std::string> binfiles;
    
        for (auto& e : std::filesystem::directory_iterator(EXPDIR)) {
            if (e.is_regular_file() && e.path().extension() == ".bin") {
                binfiles.push_back(e.path().stem().string()); // remove ".bin"
            }
        }
    
        std::sort(binfiles.begin(), binfiles.end()); 
    
        N = static_cast<int>(binfiles.size()); // total number of input files
    
        if (N == 0) {
            *ptotal_elements = 0;
            *ptotal_nodes    = 0;
            return;
        }
    
        /* Allocate FemModel pointer array                 */
        femmodels = new FemModel*[N];
    
        /* Initialize each model and accumulate global sizes      */
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

            // get number of elements and nodes on this process for this model
            int local_elements = femmodels[id]->elements->Size();
            int local_nodes    = femmodels[id]->vertices->Size();

            // add to number of elements/nodes from previous models
            total_elements += local_elements;
            total_nodes    += local_nodes;

            // specify that ISSM will be getting SMB from a GCM
            femmodels[id]->parameters->SetParam(SMBgcmEnum, SmbEnum);
            femmodels[id]->Restart();
        }
    
        /* Return (#elements, #nodes) across all models on each process */
        *ptotal_elements = total_elements; // #elements on process across all models
        *ptotal_nodes    = total_nodes;    // #nodes on process across all models
    }

	void RunISSM(IssmDouble dt, IssmDouble* gcm_forcings, IssmDouble* issm_outputs){ /*{{{*/
		/* Run ISSM for one time step given GCM forcing. */
        /* Time step is set by GCM.                      */
        /* Return a collection of ISSM outputs.          */
        int local_size;
        int global_size;
		IssmDouble start_time,final_time;
        int shift;
        int i0;

        // get total number of elements across all models on this process
        global_size = 0;
        for (int id=0;id<N;id++){
            local_size=femmodels[id]->elements->Size();
            global_size += local_size;
        }

        shift = 0; // shift starting position of each model input
        for (int id=0;id<N;id++){
		
		/*Figure out number of elements local to process: */
		local_size=femmodels[id]->elements->Size();

		/*Setup GCM forcings as element-wise input: {{{ */
		for (int f=0;f<GCMForcingNumTerms;f++){

			int forcing_type=GCMForcingTerms[f];

			for (int i=0;i<local_size;i++){
				Element* element=dynamic_cast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));

                i0 = i + shift;
                
				switch(forcing_type){
					case SMBgcmEnum:
						/*{{{*/
						{

						/*Recover smb forcing from the gcm forcings*/
						IssmDouble smb_forcing=*(gcm_forcings+f*global_size+i0); 

						/*Add into the element as new forcing :*/
						element->AddInput(SmbMassBalanceEnum,&smb_forcing,P0Enum);
						}
						/*}}}*/
						break; 
					default: 
						{ _error_("Unknown forcing type " << forcing_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Before running, setup the time interval: */
		femmodels[id]->parameters->FindParam(&start_time,TimeEnum);
		final_time=start_time+dt;
		femmodels[id]->parameters->SetParam(final_time,TimesteppingFinalTimeEnum); //we are bypassing ISSM's initial final time!

		/*Now, run: */
		femmodels[id]->Solve();

		/*Retrieve ISSM outputs and pass them back to the GCM : {{{*/
		for (int f=0;f<ISSMOutputNumTerms;f++){

			int output_type=ISSMOutputTerms[f];

			for (int i=0;i<local_size;i++){
				Element* element=dynamic_cast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));
                i0 = i + shift;
				switch(output_type){
					case SurfaceEnum:
						/*{{{*/
						{

						IssmDouble surface;

						/*Recover surface from the ISSM element: */
						Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
						surface_input->GetInputAverage(&surface);

						*(issm_outputs+f*global_size+i0) = surface;

						}
						/*}}}*/
						break; 
				case ThicknessEnum:
						/*{{{*/
						{

						IssmDouble thickness;

						/*Recover thickness from the ISSM element: */
						Input* thickness_input = element->GetInput(ThicknessEnum); _assert_(thickness_input);
						thickness_input->GetInputAverage(&thickness);

						*(issm_outputs+f*global_size+i0) = thickness;

						}
						/*}}}*/
						break; 

				case VelEnum:
						/*{{{*/
						{

						IssmDouble vel;

						/*Recover flow speed from the ISSM element: */
						Input* vel_input = element->GetInput(VelEnum); _assert_(vel_input);
						vel_input->GetInputAverage(&vel);

						*(issm_outputs+f*global_size+i0) = vel;

						}
						/*}}}*/
						break; 

					default: 
						{ _error_("Unknown output type " << output_type << "\n"); }
						break;
				}
			}
		}

		/*For the next time around, save the final time as start time */
		femmodels[id]->parameters->SetParam(final_time,TimesteppingStartTimeEnum);


        shift += local_size;
        }
		
	} /*}}}*/

	void FinalizeISSM(){ /*{{{*/
        // wrap up
        for (int i = 0; i < N; ++i) {
		/*Output results: .outbin file for each model */
			OutputResultsx(femmodels[i]);	
			delete femmodels[i];
			femmodels[i]=NULL;
        }
        delete[] femmodels; 
	} /*}}}*/

    void GetNodesISSM(int* nodeIds, IssmDouble* nodeCoords){ 
        /*obtain nodes of mesh for creating ESMF version in Fortran interface */
        /*nodeIds are the global Id's of the nodes and nodeCoords are the     */
        /*(lon,lat) coordinates, as described in the ESMF reference document  */
        int shift_global=0;
        int shift_local=0;
        int i0;
        for (int id=0;id<N;id++){
            // #vertices local to process in this model
            int local_size = femmodels[id]->vertices->Size(); 

            // #vertices across all processes in this model
            int global_size = femmodels[id]->vertices->NumberOfVertices();
            
            for (int i=0;i<local_size;i++){
                Vertex* vertex = xDynamicCast<Vertex*>(femmodels[id]->vertices->GetObjectByOffset(i));
    			i0 = vertex->Lid() + shift_local;
                *(nodeIds+i0)     = vertex->Sid()+1+shift_global;
                *(nodeCoords+2*i0+0) = vertex->longitude;
                *(nodeCoords+2*i0+1) = vertex->latitude;
                
            }
            // shift for next model starting indices:
            shift_local += local_size;
            shift_global += global_size; 
        }    
    }

    void GetElementsISSM(int* elementIds,int* elementConn,IssmDouble* elementCoords,int* glacIds){
        /*obtain elements of mesh for creating ESMF version in Fortran interface*/
        /*Element connectivity (elementConn) contains the indices of the nodes  */
        /*that form the element as described in the ESMF reference document     */ 
        int shift_local_elements = 0;
        int shift_global_elements = 0;
        int shift_local_nodes = 0;
        int i0;
        for (int id=0;id<N;id++){
            // #elements local to process in this model
            int local_size_elements = femmodels[id]->elements->Size();

            // #elements across all processes in this model
            int global_size_elements = femmodels[id]->elements->NumberOfElements();

            // #nodes local to process in this model
            int local_size_nodes = femmodels[id]->vertices->Size();
            
            for(int i=0;i<local_size_elements;i++){
                Element* element=xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));
                i0 = i + shift_local_elements; 
                *(elementIds+i0)    = element->Sid()+1+shift_global_elements;
                *(glacIds+i0)    = id;
                *(elementConn + i0*3+0) = element->vertices[0]->Lid()+1+shift_local_nodes;
                *(elementConn + i0*3+1) = element->vertices[1]->Lid()+1+shift_local_nodes;
                *(elementConn + i0*3+2) = element->vertices[2]->Lid()+1+shift_local_nodes;
    
    			// Compute the triangle centroid in longitude/latitude
    			IssmDouble centroid_lon=0.0,centroid_lat=0.0;
    			for(int j=0;j<3;j++){
    			centroid_lon += element->vertices[j]->longitude / 3.0;
    			centroid_lat += element->vertices[j]->latitude / 3.0;
    			}
    		
    			*(elementCoords + 2*i0 + 0) = centroid_lon;
    			*(elementCoords + 2*i0 + 1) = centroid_lat;
    
            }
            // shift for next model starting indices:
            shift_local_nodes += local_size_nodes; 
            shift_local_elements += local_size_elements;
            shift_global_elements += global_size_elements;
            
        }
    }

}
