! main_esmf.f90
! Minimal working example for calling ISSM Initialize-Run-Finalize methods,
! and creating an ESMF representation of the ISSM mesh
! (relevant ISSM functions are found in $ISSM_DIR/src/c/main/esmfbinders.cpp)

program main
    use iso_fortran_env, only: dp=>real64
    use iso_c_binding, only: c_ptr, c_double, c_f_pointer,c_null_char, c_loc, c_int, c_char
    use ESMF
    use netcdf
    implicit none

    ! Define the interface for the ISSM C++ functions
    interface
    subroutine InitializeISSM(expdir, num_elements, num_nodes, comm) bind(c, name="InitializeISSM")
        import :: c_char, c_int
        character(c_char), dimension(*) :: expdir
        integer(c_int)                  :: num_elements
        integer(c_int)                  :: num_nodes
        integer(c_int)                  :: comm
    end subroutine InitializeISSM
        
    subroutine RunISSM(dt, gcm_forcings, issm_outputs, elementConn) bind(C,NAME="RunISSM")
       import :: c_ptr, c_double, c_int
       real(c_double),   value   :: dt
       type(c_ptr),      value   :: gcm_forcings
       type(c_ptr),      value   :: issm_outputs
       type(c_ptr),      value   :: elementConn
    end subroutine RunISSM

    subroutine InputFromRestarts(gcm_restarts, elementConn) bind(C,NAME="InputFromRestarts")
       import :: c_ptr
       type(c_ptr),      value   :: gcm_restarts
       type(c_ptr),      value   :: elementConn
    end subroutine InputFromRestarts
   
    subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
       import :: c_ptr, c_int
       type(c_ptr),      value   :: nodeIds
       type(c_ptr),      value   :: nodeCoords
    end subroutine GetNodesISSM

    subroutine GetElementsISSM(elementIds, elementConn, elementCoords, glacIds) bind(C,NAME="GetElementsISSM")
       import :: c_ptr, c_int
       type(c_ptr),      value   :: elementIds
       type(c_ptr),      value   :: elementConn
       type(c_ptr),      value   :: elementCoords
       type(c_ptr),      value   :: glacIds
    end subroutine GetElementsISSM

    subroutine FinalizeISSM() bind(C,NAME="FinalizeISSM")
    end subroutine FinalizeISSM

    end interface
    
    ! declare ISSM-related variables
    integer(c_int)                 :: num_elements
    integer(c_int)                 :: num_nodes 
    integer(c_int)                 :: argc
    character(len=200), dimension(:), allocatable, target :: argv
    type(c_ptr), dimension(:), allocatable :: argv_ptr
    integer :: i,j
    real(dp) :: dt
    real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
    real(dp),    pointer, dimension(:)     :: SurfaceOnElements => null()
    real(dp),    pointer, dimension(:)     :: SurfaceOnNodes => null()
    real(dp),    pointer, dimension(:)     :: ICETHICK => null()
    real(dp),    pointer, dimension(:)     :: ICEVEL => null()
    real(dp),    pointer, dimension(:)     :: ICEVX => null()
    real(dp),    pointer, dimension(:)     :: ICEVY => null()
    real(dp),    pointer, dimension(:)     :: ICEEL => null()
    real(dp),    pointer, dimension(:)     :: issm_outputs => null()

    real(dp),    pointer, dimension(:)     :: SurfaceOwnNodes => null()
    

    ! filepath
    character(len=256) :: issm_path
    integer :: length, status

    ! for regridding test
    type(ESMF_RouteHandle) :: routehandle   ! routehandle for regridding nodes to elements
    type(ESMF_Field)       :: srcField      ! ice elevation on mesh
    type(ESMF_Field)       :: dstField      ! ice elevation on grid
    

    ! netcdf output
    integer :: ncid, dimid_nodes, dimid_elements, dimid_onodes
    integer :: varid_nlon,varid_nlat,varid_ecn1,varid_ecn2,varid_ecn3,varid_nid
    integer :: varid_eid,varid_elon,varid_elat,varid_esurf,varid_nsurf,varid_nowners
    integer :: varid_thick, varid_vel, varid_imask,varid_omask
    character(20) :: output_filename ! netcdf filename

    ! declare ESMF variables
    type(ESMF_VM)                  :: vm    
    integer(c_int)                 :: comm
    integer                        :: localPET
    integer                        :: petCount
    integer                        :: rc
    type(ESMF_Mesh)                :: mesh
    integer                        :: sdim
    integer, pointer, dimension(:) :: glacIds       => null()
    integer, pointer, dimension(:) :: elementIds    => null()
    integer, pointer, dimension(:) :: elementConn   => null()
    real(dp),pointer, dimension(:) :: elementCoords => null()
    real(dp),    pointer, dimension(:)     :: nodeCoords => null()
    integer,     pointer, dimension(:)     :: nodeIds => null()
    integer,     pointer, dimension(:)     :: nodeOwners => null()
    integer, allocatable  :: elementTypes(:)

    integer                        :: num_outputs = 6

    ! variables for masking the mesh seam (triangles that cross +/-180 longitude)
    real(dp)                               :: dlon,lon1,lon2,lon3
    integer, pointer, dimension(:)         :: elementMask => null()
    integer, pointer, dimension(:)         :: nodeMask    => null()
    integer                                :: n1,n2,n3

    ! test setting restarts from gcm
    real(dp),    pointer, dimension(:)     :: gcm_restarts       => null()
    real(dp),    pointer, dimension(:)     :: icemask_levelset   => null()
    real(dp),    pointer, dimension(:)     :: oceanmask_levelset => null()

    ! regridding onto nodes example
    real(dp),    pointer, dimension(:)     :: smb_elements => null()
    real(dp),    pointer, dimension(:)     :: smb_nodes => null()
    real(dp),    pointer, dimension(:)     :: smb_ownodes => null()
    type(ESMF_RouteHandle) :: routehandle_2 
    type(ESMF_Field)       :: srcField_2      
    type(ESMF_Field)       :: dstField_2 
    integer                        :: num_halo_nodes 
    integer                        :: num_owned_nodes 
    integer,     pointer, dimension(:)     :: haloSeqIndexList => null()
    type(ESMF_DistGrid)    :: nodalDistgrid
    type(ESMF_Array)       :: array0
    type(ESMF_Field)       :: field0
    type(ESMF_RouteHandle) :: haloHandle
    real(dp),    pointer, dimension(:)     :: smb_halos => null()
    
    ! stuff for file counting
    character(len=256) :: EXPDIR


    integer :: total_elements ! elements for all models (on each PET)
    integer :: total_nodes ! elements for all models (on each PET)
    

    ! initialize ESMF to and get vm info / comm for ISSM MPI
    call ESMF_Initialize(vm=vm, defaultlogfilename="VMDefaultBasicsEx.Log", logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

    call ESMF_VMGet(vm,mpiCommunicator=comm,localPET=localPET,petCount=petCount,rc=rc)    


    dt = 0.05   ! timestep in years

    ! Get the environment variable ISSM_DIR
    call get_environment_variable("ISSM_DIR", issm_path, length, status)

    EXPDIR = trim(issm_path)//"/examples/GEOSInput"//c_null_char

    call InitializeISSM(EXPDIR, num_elements, num_nodes, comm)

    print *, "rank: ", localPET, " num_elements: ", num_elements

    ! allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(nodeOwners(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(glacIds(num_elements))
    allocate(elementMask(num_elements))
    allocate(smb_elements(num_elements))
    allocate(smb_nodes(num_nodes))
   
    allocate(elementConn(3*num_elements))
    allocate(elementCoords(2*num_elements))
    
    ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
    allocate(SMBToISSM(num_nodes))
    allocate(SurfaceOnNodes(num_nodes))
    allocate(nodeMask(num_nodes))
    

    ! new for testing multiple outputs
    allocate(ICEEL(num_nodes))
    allocate(ICETHICK(num_nodes))
    allocate(ICEVEL(num_nodes))
    allocate(ICEVX(num_nodes))
    allocate(ICEVY(num_nodes))
    allocate(issm_outputs(num_outputs*num_nodes))
    allocate(gcm_restarts(num_outputs*num_nodes))
    allocate(icemask_levelset(num_nodes))
    allocate(oceanmask_levelset(num_nodes))
    
    
    

    ! create ESMF mesh corresponding to  ISSM mesh 
    ! get information about nodes and elements
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords))

    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords),c_loc(glacIds))

    print *, "max elementIds: ", maxval(elementIds)

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI
    call ESMF_VMBarrier(vm, rc=rc) 
    
    ! mask triangles that cross the seam (longitude +/- 180)
    elementMask(:) = 0
    nodeMask(:) = 0
    do j=1,num_elements
      n1 = elementConn(3*(j-1)+1)
      n2 = elementConn(3*(j-1)+2)
      n3 = elementConn(3*(j-1)+3)
      lon1 = nodeCoords(2*n1-1)
      lon2 = nodeCoords(2*n2-1)
      lon3 = nodeCoords(2*n3-1)
      dlon = maxval((/lon1,lon2,lon3/)) - minval((/lon1,lon2,lon3/))
      if ( dlon>180.0 ) then
        elementMask(j) = 1
        nodeMask(n1) = 1
        nodeMask(n2) = 1
        nodeMask(n3) = 1
      end if
    end do

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
            elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,elementMask=elementMask,&
            nodeMask=nodeMask,elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, rc=STATUS)
            
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error creating mesh'
        end if
    else if(rc == ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, "Succesfully created mesh"
        end if 
    end if

    !--------------------------------------------------------------
    ! need to set up mesh halo
    ! get nodeOwners, figure out how many owned nodes on this PET
    call ESMF_MeshGet(mesh,nodeOwners=nodeOwners)
    num_owned_nodes = count(nodeOwners == localPet)
    num_halo_nodes = num_nodes-num_owned_nodes

    allocate(haloSeqIndexList(num_halo_nodes))
    
    i=1
    do j=1,num_nodes
    if (nodeOwners(j)/= localPET) then
    haloSeqIndexList(i) = nodeIds(j)
    i = i+1
    end if 
    end do

    call ESMF_MeshGet(mesh, nodalDistgrid=nodalDistgrid, rc=rc)

    array0=ESMF_ArrayCreate(nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=haloSeqIndexList,rc=rc)

    field0=ESMF_FieldCreate(mesh, array=array0, meshLoc=ESMF_MESHLOC_NODE, rc=rc)

    call ESMF_FieldHaloStore(field0, routehandle=haloHandle, rc=rc)

    call ESMF_FieldGet(field0,farrayPtr=smb_halos,rc=rc)

    call ESMF_FieldHalo(field0, routehandle=haloHandle, rc=rc)


 !--------------------------------------------------------------

    allocate(SurfaceOwnNodes(num_owned_nodes))

    call ESMF_VMBarrier(vm, rc=rc)

    !-------------------------------------------------------------
    ! try regridding "smb" from elements to (owned) nodes
    smb_elements(:) = 1
    srcField_2 = ESMF_FieldCreate(mesh=mesh,farrayPtr=smb_elements,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)

    dstField_2 = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_NODE,rc=rc)


    call ESMF_FieldRegridStore(srcField=srcField_2, dstField=field0,routehandle=routehandle_2, &
    srcMaskValues=(/1/),unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_NEAREST_D,rc=rc)

    
    call ESMF_VMBarrier(vm, rc=rc)
    call ESMF_FieldRegrid(srcField_2, field0, routehandle_2, rc=rc)

    call ESMF_VMBarrier(vm, rc=rc)
    
    call ESMF_FieldHalo(field0, routehandle=haloHandle, rc=rc)
    call ESMF_FieldGet(field0,farrayPtr=smb_ownodes,rc=rc)


    if (localPET==1) then
        print *, "num nodes:", num_nodes
        print *, "size smb_ownodes:", smb_ownodes
    end if 

    !where (nodeOwners == localPET)
    !    smb_nodes = smb_ownodes
    !end where

    print *, "regrid with new field worked"

    STOP

    !print *, "smb regridded: ", smb_nodes

!-------------------------------------------------------------
    ! set smb and surface for test 
    SMBToISSM(:) = 0! smb_nodes(:)
    SurfaceOnElements(:) = 0 ! placeholder
    SurfaceOnNodes(:) = 0 ! placeholder    

    ICEEL(:) = 0 ! placeholder 
    ICETHICK(:) = 0 ! placeholder    
    ICEVEL(:) = 0 ! placeholder    
    issm_outputs(:) = 0 !placeholder

    call ESMF_VMBarrier(vm, rc=rc)
    !call the C++ routine for running a single time step
    call RunISSM(dt, c_loc(SMBToISSM), c_loc(issm_outputs), c_loc(elementConn))

    ! assign restarts to output for next time step, just as an example...
    gcm_restarts(:) = issm_outputs(:)

    ! try setting inputs...
    call InputFromRestarts(c_loc(gcm_restarts), c_loc(elementConn))
 
    ! try running again to see if we get correct result from restart inputs
    call ESMF_VMBarrier(vm, rc=rc)
    !call the C++ routine for running a single time step
    call RunISSM(dt, c_loc(SMBToISSM), c_loc(issm_outputs), c_loc(elementConn))

    SurfaceOnNodes(:) = issm_outputs(1:num_nodes)
    ICEEL(:) = issm_outputs(1:num_nodes)
    ICETHICK(:) = issm_outputs(num_nodes+1:2*num_nodes)
    ICEVX(:) = issm_outputs(2*num_nodes+1:3*num_nodes)
    ICEVY(:) = issm_outputs(3*num_nodes+1:4*num_nodes)
    oceanmask_levelset(:) = issm_outputs(4*num_nodes+1:5*num_nodes)
    icemask_levelset(:) = issm_outputs(5*num_nodes+1:6*num_nodes)
    ICEVEL = sqrt(ICEVX**2 + ICEVY**2)

    SurfaceOwnNodes(:) = pack(SurfaceOnNodes, nodeOwners == localPET)
 
    call ESMF_VMBarrier(vm, rc=rc)  
    srcField = ESMF_FieldCreate(mesh=mesh,farrayPtr=SurfaceOwnNodes,meshloc=ESMF_MESHLOC_NODE, & 
               datacopyflag=ESMF_DATACOPY_VALUE,rc=rc)
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error creating srcField'
        end if 
    else if (rc == ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Successfully created srcField'
        end if 
    end if

    ! create dstField (will be associated with SurfaceOnElements)
    dstField = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    
    call ESMF_FieldFill(dstField, dataFillScheme="const",const1=0.0_dp)
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error creating dstField'
        end if 
    else if (rc == ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Successfully created dstField'
        end if 
    end if

    !call ESMF_VMBarrier(vm, rc=rc)
    if (localPET == 0) then
        print *, "Creating routehandle...."
    end if 
    
    ! create routehandle for regridding fields from elements to vertices
    ! mask elements along mesh seam (longitudes crossing +/- 180)
    call ESMF_FieldRegridStore(srcField=srcField, dstField=dstField,routehandle=routehandle, &
    dstMaskValues=(/1/),unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_NONE,rc=rc)

    if (rc /= ESMF_SUCCESS) then
        print *, 'Error creating routehandles'
    else
        if (localPET == 0) then
        print *, "Created routehandles...."
        end if 
    end if
    
    if (localPET == 0) then
        print *, "trying to regrid...."
    end if 

    ! Perform regridding from elements to vertices
    call ESMF_FieldRegrid(srcField, dstField, routehandle, rc=rc)
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error regridding'
        end if 
    else
        if (localPET == 0) then
        print *, "Successfully regridded"
        end if 
    end if   

    ! associate regridded field with SurfaceOnElements pointer
    call ESMF_FieldGet(dstField,farrayPtr=SurfaceOnElements,rc=rc) 



    if (localPET == 0) then
    print *, "saving results to netcdf files..."
    end if 

    ! save a netcdf for each PET 
    do i = 0, petCount-1
        call ESMF_VMBarrier(vm, rc=rc) 
        if (i==localPET) then
        write(output_filename,'("data_rank",I0,".nc")') localPET
        ! definitions:
            rc = nf90_create(output_filename, NF90_CLOBBER, ncid)
            rc = nf90_def_dim(ncid,"num_nodes",num_nodes,dimid_nodes)
            rc = nf90_def_dim(ncid,"num_owned_nodes",num_owned_nodes,dimid_onodes)
            rc = nf90_def_dim(ncid,"num_elements",num_elements,dimid_elements)
            rc = nf90_def_var(ncid,"node_lon",NF90_REAL,(/dimid_nodes/),varid_nlon)
            rc = nf90_def_var(ncid,"node_lat",NF90_REAL,(/dimid_nodes/),varid_nlat)
            rc = nf90_def_var(ncid,"node_ids",NF90_REAL,(/dimid_nodes/),varid_nid)
            rc = nf90_def_var(ncid,"node_owners",NF90_REAL,(/dimid_nodes/),varid_nowners)
            rc = nf90_def_var(ncid,"element_lon",NF90_REAL,(/dimid_elements/),varid_elon)
            rc = nf90_def_var(ncid,"element_lat",NF90_REAL,(/dimid_elements/),varid_elat)
            rc = nf90_def_var(ncid,"element_ids",NF90_REAL,(/dimid_elements/),varid_eid)
            rc = nf90_def_var(ncid,"element_conn1",NF90_REAL,(/dimid_elements/),varid_ecn1)
            rc = nf90_def_var(ncid,"element_conn2",NF90_REAL,(/dimid_elements/),varid_ecn2)
            rc = nf90_def_var(ncid,"element_conn3",NF90_REAL,(/dimid_elements/),varid_ecn3)
            rc = nf90_def_var(ncid,"node_surf",NF90_REAL,(/dimid_nodes/),varid_nsurf)
            rc = nf90_def_var(ncid,"element_surf",NF90_REAL,(/dimid_elements/),varid_esurf)
            rc = nf90_def_var(ncid,"ice_thick",NF90_REAL,(/dimid_nodes/),varid_thick)
            rc = nf90_def_var(ncid,"ice_vel",NF90_REAL,(/dimid_nodes/),varid_vel)
            rc = nf90_def_var(ncid,"ice_mask",NF90_REAL,(/dimid_nodes/),varid_imask)
            rc = nf90_def_var(ncid,"ocean_mask",NF90_REAL,(/dimid_nodes/),varid_omask)
       
        ! set arrays:    
            rc = nf90_enddef(ncid)
            rc = nf90_put_var(ncid,varid_nlon,nodeCoords(1::2))
            rc = nf90_put_var(ncid,varid_nlat,nodeCoords(2::2))
            rc = nf90_put_var(ncid,varid_nid,nodeIds)
            rc = nf90_put_var(ncid,varid_nowners,nodeOwners)
            rc = nf90_put_var(ncid,varid_elon,elementCoords(1::2))
            rc = nf90_put_var(ncid,varid_elat,elementCoords(2::2))
            rc = nf90_put_var(ncid,varid_eid,elementIds(:))
            rc = nf90_put_var(ncid,varid_ecn1,elementConn(1::3))
            rc = nf90_put_var(ncid,varid_ecn2,elementConn(2::3))
            rc = nf90_put_var(ncid,varid_ecn3,elementConn(3::3))
            rc = nf90_put_var(ncid,varid_esurf,SurfaceOnElements(:))
            rc = nf90_put_var(ncid,varid_nsurf,SurfaceOnNodes(:))
            rc = nf90_put_var(ncid,varid_thick,ICETHICK(:))
            rc = nf90_put_var(ncid,varid_vel,ICEVEL(:))
            rc = nf90_put_var(ncid,varid_imask,icemask_levelset(:))
            rc = nf90_put_var(ncid,varid_omask,oceanmask_levelset(:))        
            rc = nf90_close(ncid)
        end if
        call ESMF_VMBarrier(vm, rc=rc) 
    end do
    
    call ESMF_VMBarrier(vm, rc=rc)    

    ! destroy fields
    call ESMF_FieldDestroy(srcField,rc=rc)
    call ESMF_FieldDestroy(dstField,rc=rc)
    call ESMF_FieldDestroy(srcField_2,rc=rc)
    call ESMF_FieldDestroy(dstField_2,rc=rc)

    ! deallocate pointers
    deallocate(nodeCoords)
    deallocate(nodeIds)
    deallocate(nodeOwners)
    deallocate(elementTypes)
    deallocate(elementIds)
    deallocate(elementConn)
    deallocate(elementCoords)
    deallocate(SMBToISSM)
    deallocate(ICEEL)
    deallocate(ICETHICK)
    deallocate(ICEVEL)
    deallocate(issm_outputs)
    deallocate(gcm_restarts)
    deallocate(icemask_levelset)
    deallocate(oceanmask_levelset)
    deallocate(ICEVX)
    deallocate(ICEVY)
    deallocate(smb_elements)
    deallocate(smb_nodes)
    
    
    ! call ISSM finalize (saves binary output .outbin file)
    call FinalizeISSM()
    
    ! finalize ESMF
    call ESMF_Finalize(rc=rc)
end program main
