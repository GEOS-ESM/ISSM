/*!\file:  esmfbinder.cpp
 * \brief: ESMF binders for ISSM. Binders developed initially for the GEOS-5 framework.
 */ 

#include "./issm.h"
#include <filesystem>
#include <cstring>

/*GEOS 5 specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum}; 
const int ISSMOutputNumTerms = 3;
const int ISSMOutputTerms[ISSMOutputNumTerms]= { SurfaceEnum, ThicknessEnum, VelEnum };

extern "C" {


	static int N = 0;
	static FemModel** femmodels = nullptr;

	int init_models(const char* dir,
					char* files,
					int max_files,
					int len)
	{
		N = 0;

		for (auto& e : std::filesystem::directory_iterator(dir))
			if (e.is_regular_file() &&
				e.path().extension() == ".bin" &&
				N < max_files)
			{
				std::string name = e.path().filename().string();
				std::strncpy(files + N*len, name.c_str(), len-1);
				files[N*len + len-1] = '\0';
				++N;
			}

		// allocate the array of pointers for later use
		if (N > 0)
			femmodels = new FemModel*[N];

		return N;   
	}

	/*GEOS 5*/
      void InitializeISSM(int argc, char** argv, int* pnumberofelements, int* pnumberofnodes, MPI_Fint* Fcomm, int id){ /*{{{*/
		int numberofelements;
        int numberofnodes;
          
        /* convert Fortran MPI comm to C MPI comm */
        MPI_Comm Ccomm = MPI_Comm_f2c(*Fcomm);             
                
        /*Initialize femmodel from arguments provided command line: */
		femmodels[id] = new FemModel(argc,argv,Ccomm);

		/*Get number of nodes and elements local to each process: */
		numberofelements=femmodels[id]->elements->Size();
		numberofnodes=femmodels[id]->vertices->Size();

		/*Bypass SMB model, will be provided by GCM! */
		femmodels[id]->parameters->SetParam(SMBgcmEnum,SmbEnum); 

        /*Restart file: */
		femmodels[id]->Restart();

		/*Assign output pointers: */
		*pnumberofelements=numberofelements;
        *pnumberofnodes=numberofnodes;
	} /*}}}*/

	void RunISSM(IssmDouble dt, IssmDouble* gcm_forcings, IssmDouble* issm_outputs){ /*{{{*/
		int local_size;
        int global_size;
		IssmDouble start_time,final_time;
        int shift;
        int i0;

        // get total number of elements across all models
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

        for (int i = 0; i < N; ++i) {
		/*Output results: */
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
        int shift;
        shift = 0;
        int i0;
        for (int id=0;id<N;id++){
            int local_size = femmodels[id]->vertices->Size();
            for (int i=0;i<local_size;i++){
                Vertex* vertex = xDynamicCast<Vertex*>(femmodels[id]->vertices->GetObjectByOffset(i));
    			i0 = vertex->Lid() + shift;
                *(nodeIds+i0)     = vertex->Sid()+1+shift;
                *(nodeCoords+2*i0+0) = vertex->longitude;
                *(nodeCoords+2*i0+1) = vertex->latitude;
                
            }
            shift += local_size; 
        }    
    }

    void GetElementsISSM(int* elementIds,int* elementConn,IssmDouble* elementCoords,int* glacIds){
        /*obtain elements of mesh for creating ESMF version in Fortran interface*/
        /*Element connectivity (elementConn) contains the indices of the nodes  */
        /*that form the element as described in the ESMF reference document     */ 
        int shift_elements;
        int shift_nodes;
        shift_elements = 0;
        shift_nodes = 0;
        int i0;
        for (int id=0;id<N;id++){
            int local_size_elements = femmodels[id]->elements->Size();
            int local_size_nodes = femmodels[id]->vertices->Size();
            for(int i=0;i<local_size_elements;i++){
                Element* element=xDynamicCast<Element*>(femmodels[id]->elements->GetObjectByOffset(i));
                i0 = i + shift_elements; 
                *(elementIds+i0)    = element->Sid()+1+shift_elements;
                *(glacIds+i0)    = id;
                *(elementConn + i0*3+0) = element->vertices[0]->Lid()+1+shift_nodes;
                *(elementConn + i0*3+1) = element->vertices[1]->Lid()+1+shift_nodes;
                *(elementConn + i0*3+2) = element->vertices[2]->Lid()+1+shift_nodes;
    
    			// Compute the triangle centroid in longitude/latitude
    			IssmDouble centroid_lon=0.0,centroid_lat=0.0;
    			for(int j=0;j<3;j++){
    			centroid_lon += element->vertices[j]->longitude / 3.0;
    			centroid_lat += element->vertices[j]->latitude / 3.0;
    			}
    		
    			*(elementCoords + 2*i0 + 0) = centroid_lon;
    			*(elementCoords + 2*i0 + 1) = centroid_lat;
    
            }
            shift_elements += local_size_elements; 
            shift_nodes += local_size_nodes;
        }
    }

}
