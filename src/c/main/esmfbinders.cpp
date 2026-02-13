/*!\file:  esmfbinder.cpp
 * \brief: ESMF binders for ISSM. Binders developed initially for the GEOS-5 framework.
 */ 

#include "./issm.h"

/*GEOS 5 specific declarations:*/
const int GCMForcingNumTerms = 1;
const int GCMForcingTerms[GCMForcingNumTerms]= { SMBgcmEnum}; 
const int ISSMOutputNumTerms = 3;
const int ISSMOutputTerms[ISSMOutputNumTerms]= { SurfaceEnum, ThicknessEnum, VelEnum };

extern "C" {

	FemModel *femmodel;

	/*GEOS 5*/
      void InitializeISSM(int argc, char** argv, int* pnumberofelements, int* pnumberofnodes, MPI_Fint* Fcomm){ /*{{{*/
		int numberofelements;
        int numberofnodes;
          
        /* convert Fortran MPI comm to C MPI comm */
        MPI_Comm Ccomm = MPI_Comm_f2c(*Fcomm);             
                
        /*Initialize femmodel from arguments provided command line: */
		femmodel = new FemModel(argc,argv,Ccomm);

		/*Get number of nodes and elements local to each process: */
		numberofelements=femmodel->elements->Size();
		numberofnodes=femmodel->vertices->Size();

		/*Bypass SMB model, will be provided by GCM! */
		femmodel->parameters->SetParam(SMBgcmEnum,SmbEnum); 

        /*Restart file: */
		femmodel->Restart();

		/*Assign output pointers: */
		*pnumberofelements=numberofelements;
        *pnumberofnodes=numberofnodes;
	} /*}}}*/

	void RunISSM(IssmDouble dt, IssmDouble* gcmforcings, IssmDouble* issmoutputs){ /*{{{*/
		int numberofelements;
		IssmDouble start_time,final_time;
		
		/*Figure out number of elements local to process: */
		numberofelements=femmodel->elements->Size();

		/*Setup gcm forcings as element-wise input: {{{ */
		for (int f=0;f<GCMForcingNumTerms;f++){

			int forcing_type=GCMForcingTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(forcing_type){
					case SMBgcmEnum:
						/*{{{*/
						{

						/*Recover smb forcing from the gcm forcings*/
						IssmDouble smbforcing=*(gcmforcings+f*numberofelements+i); 

						/*Add into the element as new forcing :*/
						element->AddInput(SmbMassBalanceEnum,&smbforcing,P0Enum);
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

		/*Retrieve ISSM outputs and pass them back to the Gcm : {{{*/
		for (int f=0;f<ISSMOutputNumTerms;f++){

			int output_type=ISSMOutputTerms[f];

			for (int i=0;i<femmodel->elements->Size();i++){
				Element* element=dynamic_cast<Element*>(femmodel->elements->GetObjectByOffset(i));

				switch(output_type){
					case SurfaceEnum:
						/*{{{*/
						{

						IssmDouble surface;

						/*Recover surface from the ISSM element: */
						Input* surface_input = element->GetInput(SurfaceEnum); _assert_(surface_input);
						surface_input->GetInputAverage(&surface);

						*(issmoutputs+f*numberofelements+i) = surface;

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

						*(issmoutputs+f*numberofelements+i) = thickness;

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

						*(issmoutputs+f*numberofelements+i) = vel;

						}
						/*}}}*/
						break; 

					default: 
						{ _error_("Unknown output type " << output_type << "\n"); }
						break;
				}
			}
		}

		/*}}}*/

		/*Before running, setup the time interval: */
		femmodel->parameters->FindParam(&start_time,TimeEnum);
		final_time=start_time+dt;
		femmodel->parameters->SetParam(final_time,TimesteppingFinalTimeEnum); //we are bypassing ISSM's initial final time!

		/*Now, run: */
		femmodel->Solve();

		/*For the next time around, save the final time as start time */
		femmodel->parameters->SetParam(final_time,TimesteppingStartTimeEnum);
		
	} /*}}}*/

	void FinalizeISSM(){ /*{{{*/

		/*Output results: */
		OutputResultsx(femmodel);

		/*Wrap up: */
		delete femmodel; femmodel=NULL;
	} /*}}}*/

    void GetNodesISSM(int* nodeIds, IssmDouble* nodeCoords){ 
        /*obtain nodes of mesh for creating ESMF version in Fortran interface */
        /*nodeIds are the global Id's of the nodes and nodeCoords are the     */
        /*(lon,lat) coordinates, as described in the ESMF reference document  */
		int i0;
        for (int i=0;i<femmodel->vertices->Size();i++){
            Vertex* vertex = xDynamicCast<Vertex*>(femmodel->vertices->GetObjectByOffset(i));
			i0 = vertex->Lid();
            *(nodeIds+i0)     = vertex->Sid()+1;
            *(nodeCoords+2*i0+0) = vertex->longitude;
            *(nodeCoords+2*i0+1) = vertex->latitude;
        }
    }

    void GetElementsISSM(int* elementIds,int* elementConn,IssmDouble* elementCoords){
        /*obtain elements of mesh for creating ESMF version in Fortran interface*/
        /*Element connectivity (elementConn) contains the indices of the nodes  */
        /*that form the element as described in the ESMF reference document     */ 
        for(int i=0;i<femmodel->elements->Size();i++){
            Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
            *(elementIds + i)    = element->Sid()+1;
            *(elementConn + i*3+0) = element->vertices[0]->Lid()+1;
            *(elementConn + i*3+1) = element->vertices[1]->Lid()+1;
            *(elementConn + i*3+2) = element->vertices[2]->Lid()+1;

			// Compute the triangle centroid in longitude/latitude
			IssmDouble centroid_lon=0.0,centroid_lat=0.0;
			for(int j=0;j<3;j++){
			centroid_lon += element->vertices[j]->longitude / 3.0;
			centroid_lat += element->vertices[j]->latitude / 3.0;
			}
		
			*(elementCoords + 2*i + 0) = centroid_lon;
			*(elementCoords + 2*i + 1) = centroid_lat;

        }
    }

}
