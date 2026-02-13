! main_esmf.f90
! Minimal working example for calling ISSM Initialize-Run-Finalize methods,
! and creating an ESMF representation of the ISSM mesh
! (relevant ISSM functions are found in $ISSM_DIR/src/c/main/esmfbinders.cpp)

program main
    use iso_fortran_env, only: dp=>real64
    use iso_c_binding, only: c_ptr, c_double, c_f_pointer,c_null_char, c_loc, c_int
    use ESMF
    use netcdf
    implicit none

    ! Define the interface for the ISSM C++ functions
    interface
    subroutine InitializeISSM(argc, argv, num_elements, num_nodes, comm) bind(c, name="InitializeISSM")
        import :: c_ptr, c_int
        integer(c_int), value        :: argc
        type(c_ptr), dimension(argc) :: argv
        integer(c_int)               :: num_elements
        integer(c_int)               :: num_nodes
        integer(c_int)               :: comm
    end subroutine InitializeISSM
        
    subroutine RunISSM(dt, gcm_forcings, issm_outputs) bind(C,NAME="RunISSM")
       import :: c_ptr, c_double
       real(c_double),   value   :: dt
       type(c_ptr),      value   :: gcm_forcings
       type(c_ptr),      value   :: issm_outputs
    end subroutine RunISSM
   
    subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
       import :: c_ptr
       type(c_ptr),      value   :: nodeIds
       type(c_ptr),      value   :: nodeCoords
    end subroutine GetNodesISSM

    subroutine GetElementsISSM(elementIds, elementConn, elementCoords) bind(C,NAME="GetElementsISSM")
       import :: c_ptr
       type(c_ptr),      value   :: elementIds
       type(c_ptr),      value   :: elementConn
       type(c_ptr),      value   :: elementCoords
    end subroutine GetElementsISSM

    subroutine FinalizeISSM() bind(C,NAME="FinalizeISSM")
    end subroutine FinalizeISSM

    end interface
    
    ! declare ISSM-related variables
    integer(c_int)                 :: num_elements
    integer(c_int)                 :: num_nodes
    integer                        :: num_owned_nodes   
    integer(c_int)                 :: argc
    character(len=200), dimension(:), allocatable, target :: argv
    type(c_ptr), dimension(:), allocatable :: argv_ptr
    integer :: i
    real(dp) :: dt
    real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
    real(dp),    pointer, dimension(:)     :: SurfaceOnElements => null()
    real(dp),    pointer, dimension(:)     :: SurfaceOnNodes => null()
    real(dp),    pointer, dimension(:)     :: ICETHICK => null()
    real(dp),    pointer, dimension(:)     :: ICEVEL => null()
    real(dp),    pointer, dimension(:)     :: ICEEL => null()
    real(dp),    pointer, dimension(:)     :: issm_outputs => null()

    ! filepath
    character(len=256) :: issm_path
    integer :: length, status

    ! for regridding test
    type(ESMF_RouteHandle) :: routehandle ! routehandle for regridding
    type(ESMF_Field)       :: srcField    ! ice elevation on mesh
    type(ESMF_Field)       :: dstField    ! ice elevation on grid

    ! netcdf output
    integer :: ncid, dimid_nodes, dimid_elements, dimid_onodes
    integer :: varid_nlon,varid_nlat,varid_ecn1,varid_ecn2,varid_ecn3,varid_nid
    integer :: varid_eid,varid_elon,varid_elat,varid_esurf,varid_nsurf,varid_nowners
    character(20) :: output_filename ! netcdf filename

    ! declare ESMF variables
    type(ESMF_VM)                  :: vm    
    integer(c_int)                 :: comm
    integer                        :: localPET
    integer                        :: petCount
    integer                        :: rc
    type(ESMF_Mesh)                :: mesh
    integer                        :: sdim
    integer, pointer, dimension(:) :: elementIds    => null()
    integer, pointer, dimension(:) :: elementConn   => null()
    real(dp),pointer, dimension(:) :: elementCoords => null()
    real(dp),    pointer, dimension(:)     :: nodeCoords => null()
    integer,     pointer, dimension(:)     :: nodeIds => null()
    integer,     pointer, dimension(:)     :: nodeOwners => null()
    integer, allocatable  :: elementTypes(:)

    integer                        :: num_outputs = 3


    dt = 0.05   ! timestep in years

    ! Get the environment variable ISSM_DIR
    call get_environment_variable("ISSM_DIR", issm_path, length, status)

    ! Manually set argc and argv
    argc = 4  ! Example: 3 arguments
    allocate(argv(argc))
    argv(1) = "this arg does not matter"//c_null_char 
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = trim(issm_path)//"/examples/GEOSInput"//c_null_char
    argv(4) = "GreenlandGEOS"//c_null_char

    ! Convert Fortran strings to C pointers (argv)
    allocate(argv_ptr(argc))

    do i = 1, argc
        ! Ensure that we are only getting the memory address once per string
        argv_ptr(i) = c_loc(argv(i))
    end do

    ! initialize ESMF to and get vm info / comm for ISSM MPI
    call ESMF_Initialize(vm=vm, defaultlogfilename="VMDefaultBasicsEx.Log", logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

    call ESMF_VMGet(vm,mpiCommunicator=comm,localPET=localPET,petCount=petCount,rc=rc)    

    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call ESMF_VMBarrier(vm, rc=rc)
    call InitializeISSM(argc, argv_ptr,num_elements,num_nodes,comm)
    call ESMF_VMBarrier(vm, rc=rc)

    ! ! print out some info if desired:
    !print *, "number of elements on PET ", localPET, ": ", num_elements
    !print *, "number of nodes on PET ", localPET, ": ", num_nodes

    ! allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(nodeOwners(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(elementConn(3*num_elements))
    allocate(elementCoords(2*num_elements))
    
    ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
    allocate(SMBToISSM(num_elements))
    allocate(SurfaceOnElements(num_elements))

    ! new for testing multiple outputs
    allocate(ICEEL(num_elements))
    allocate(ICETHICK(num_elements))
    allocate(ICEVEL(num_elements))
    allocate(issm_outputs(num_outputs*num_elements))

    ! create ESMF mesh corresponding to  ISSM mesh 
    ! get information about nodes and elements
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords))
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI
    call ESMF_VMBarrier(vm, rc=rc) 
    
    ! create the ESMF mesh (later will be used for regridding)
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
           elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,& 
           elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error creating mesh'
        end if
    else if(rc == ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, "Succesfully created mesh"
        end if 
    end if

    call ESMF_VMBarrier(vm, rc=rc)

    call ESMF_MeshGet(mesh,nodeOwners=nodeOwners)

    num_owned_nodes = count(nodeOwners == localPet)

    ! set smb and surface for test 
    SMBToISSM(:) = 0
    SurfaceOnElements(:) = 0 ! placeholder
    SurfaceOnNodes(:) = 0 ! placeholder    

    ICEEL(:) = 0 ! placeholder 
    ICETHICK(:) = 0 ! placeholder    
    ICEVEL(:) = 0 ! placeholder    
    issm_outputs(:) = 0 !placeholder

    
    call ESMF_VMBarrier(vm, rc=rc)
    !call the C++ routine for running a single time step
    call RunISSM(dt, c_loc(SMBToISSM), c_loc(issm_outputs))


    SurfaceOnElements(:) = issm_outputs(1:num_elements)
    ICEEL(:) = big_array(1:num_elements)
    ICETHICK(:) = big_array(num_elements+1:2*num_elements)
    ICEVEL(:) = big_array(2*num_elements+1:3*num_elements)

 
    call ESMF_VMBarrier(vm, rc=rc)  
    srcField = ESMF_FieldCreate(mesh=mesh,farrayPtr=SurfaceOnElements,meshloc=ESMF_MESHLOC_ELEMENT, & 
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

    ! create dstField (will be associated with SurfaceOnNodes)
    dstField = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_NODE,rc=rc)
    
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
    call ESMF_FieldRegridStore(srcField=srcField, dstField=dstField,routehandle=routehandle, unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,rc=rc)
    
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error creating routehandle'
        end if 
    else
        if (localPET == 0) then
        print *, "Created routehandle...."
        end if 
    end if
    
    if (localPET == 0) then
        print *, "trying to regrid...."
    end if 

    ! Perform regridding from elements to vertices
    call ESMF_FieldRegrid(srcField, dstField, routeHandle, rc=rc)
    if (rc /= ESMF_SUCCESS) then
        if (localPET == 0) then
        print *, 'Error regridding'
        end if 
    else
        if (localPET == 0) then
        print *, "Successfully regridded"
        end if 
    end if   

    ! associate regridded field with SurfaceOnNodes pointer
    call ESMF_FieldGet(dstField,farrayPtr=SurfaceOnNodes,rc=rc) 

    if (localPET == 0) then
    print *, "saving results to netcdf files..."
    end if 

    ! save a netcdf for each PET 
    do i = 0, petCount-1
        call ESMF_VMBarrier(vm, rc=rc) 
        if (i==localPET) then
        write(output_filename,'("data_rank",I0,".nc")') localPET
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
            rc = nf90_def_var(ncid,"node_surf",NF90_REAL,(/dimid_onodes/),varid_nsurf)
            rc = nf90_def_var(ncid,"element_surf",NF90_REAL,(/dimid_elements/),varid_esurf)
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
            rc = nf90_close(ncid)
        end if
        call ESMF_VMBarrier(vm, rc=rc) 
    end do
    
    call ESMF_VMBarrier(vm, rc=rc)    

    ! destroy fields
    call ESMF_FieldDestroy(srcField,rc=rc)
    call ESMF_FieldDestroy(dstField,rc=rc)

    ! deallocate pointers
    deallocate(nodeCoords)
    deallocate(nodeIds)
    deallocate(nodeOwners)
    deallocate(elementTypes)
    deallocate(elementIds)
    deallocate(elementConn)
    deallocate(elementCoords)
    deallocate(SMBToISSM)
    deallocate(SurfaceOnElements)
    deallocate(ICEEL)
    deallocate(ICETHICK)
    deallocate(ICEVEL)
    deallocate(issm_outputs)
    
    ! call ISSM finalize (saves binary output .outbin file)
    call FinalizeISSM()
    
    ! finalize ESMF
    call ESMF_Finalize(rc=rc)
end program main
