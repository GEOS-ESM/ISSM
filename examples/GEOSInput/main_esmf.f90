! main_esmf.f90
! Minimal working example for calling ISSM Initialize-Run-Finalize methods,
! and creating an ESMF representation of the ISSM mesh
! (relevant ISSM functions are found in $ISSM_DIR/src/c/main/esmfbinders.cpp)

program main
    use iso_fortran_env, only: dp=>real64
    use iso_c_binding, only: c_ptr, c_double, c_f_pointer,c_null_char, c_loc, c_int
    use ESMF
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
        
    subroutine RunISSM(dt, gcmf, issmouts) bind(C,NAME="RunISSM")
       import :: c_ptr, c_double
       real(c_double),   value   :: dt
       type(c_ptr),      value   :: gcmf
       type(c_ptr),      value   :: issmouts
    end subroutine RunISSM
   
    subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
       import :: c_ptr
       type(c_ptr),      value   :: nodeIds
       type(c_ptr),      value   :: nodeCoords
    end subroutine GetNodesISSM

    subroutine GetElementsISSM(elementIds, elementConn) bind(C,NAME="GetElementsISSM")
       import :: c_ptr
       type(c_ptr),      value   :: elementIds
       type(c_ptr),      value   :: elementConn
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
    integer :: i
    real(dp) :: dt
    real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
    real(dp),    pointer, dimension(:)     :: SurfaceToGEOS => null()
    real(dp),    pointer, dimension(:)     :: nodeCoords => null()
    integer,     pointer, dimension(:)     :: nodeIds => null()

    ! declare ESMF variables
    type(ESMF_VM)                  :: vm    
    integer                        :: rc
    type(ESMF_Mesh)                :: mesh
    integer(c_int)                 :: comm
    integer                        :: sdim
    integer, pointer, dimension(:) :: elementIds
    integer, pointer, dimension(:) :: elementConn    
    integer, allocatable  :: elementTypes(:)
    type(ESMF_Field) :: field

    dt = 0.05   ! timestep in years

    ! Manually set argc and argv
    argc = 4  ! Example: 3 arguments
    allocate(argv(argc))
    argv(1) = "this arg does not matter" !"/discover/nobackup/projects/gmao/SIteam/ISSM/2025-09-02/ifort_2021.13.0-intelmpi_2021.13.0/ISSM/bin/issm.exe"//c_null_char
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = "/discover/nobackup/agstubbl/ISSM/projs/IRF-ISSM"//c_null_char
    argv(4) = "GreenlandGEOS"//c_null_char

    ! Convert Fortran strings to C pointers (argv)
    allocate(argv_ptr(argc))

    do i = 1, argc
        ! Ensure that we are only getting the memory address once per string
        argv_ptr(i) = c_loc(argv(i))
    end do

    ! initialize ESMF to and get vm info / comm for ISSM MPI
    call ESMF_Initialize(vm=vm, defaultlogfilename="VMDefaultBasicsEx.Log", logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

    call ESMF_VMGet(vm,mpiCommunicator=comm,rc=rc)    

    ! ! print the VM information if desired:
    !call ESMF_VMPrint(vm, rc=rc)
   
    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call ESMF_VMBarrier(vm, rc=rc)
    call InitializeISSM(argc, argv_ptr,num_elements,num_nodes,comm)
    call ESMF_VMBarrier(vm, rc=rc)

    print *, "number of elements: ", num_elements
    print *, "number of nodes: ", num_nodes

    ! allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(elementConn(3*num_elements))

    ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
    allocate(SMBToISSM(num_elements))
    allocate(SurfaceToGEOS(num_elements))

    ! create ESMF mesh corresponding to  ISSM mesh 
    ! get information about nodes and elements
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords))
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI

    !do i=2,size(nodeCoords)
    !if (modulo(i, 2) == 0) then
    !    print '(A, F6.2, A, F7.2)', 'lat: ', nodeCoords(i), ', lon: ', nodeCoords(i-1)
    !else
    !    continue
    !end if 
    !end do

    ! create the ESMF mesh (later will be used for regridding)
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
           elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)

    ! set smb and surface for test 
    SMBToISSM(:) = 0
    SurfaceToGEOS(:) = 0 ! placeholder    

    call ESMF_VMBarrier(vm, rc=rc)
    ! call the C++ routine for running a single time step
    call RunISSM(dt, c_loc(SMBToISSM), c_loc(SurfaceToGEOS))
    call ESMF_VMBarrier(vm, rc=rc)    

    ! call ISSM finalize (saves binary output .outbin file)
    call FinalizeISSM()
    call ESMF_VMBarrier(vm, rc=rc)     

    ! ! save mesh properties
    ! save mesh node ID's
    field = ESMF_FieldCreate(mesh=mesh,farrayPtr=nodeIds,meshloc=ESMF_MESHLOC_NODE)    
    call ESMF_FieldWrite(field,"mesh.nc",variableName="nodeIds")
    call ESMF_FieldDestroy(field,rc=rc)   
    
    ! save mesh node x-coord
    field = ESMF_FieldCreate(mesh=mesh,farray=nodeCoords(2*nodeIds-1),indexflag=ESMF_INDEX_DELOCAL,datacopyflag=ESMF_DATACOPY_VALUE,meshloc=ESMF_MESHLOC_NODE,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="nodeCoords_lon",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)


    ! save mesh node y-coord
    field = ESMF_FieldCreate(mesh=mesh,farray=nodeCoords(2*nodeIds),indexflag=ESMF_INDEX_DELOCAL,datacopyflag=ESMF_DATACOPY_VALUE,meshloc=ESMF_MESHLOC_NODE,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="nodeCoords_lat",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)   
    
    ! save element ID's
    field = ESMF_FieldCreate(mesh=mesh,farray=elementIds,indexflag=ESMF_INDEX_DELOCAL,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="elementIds",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)


    ! save element connectivity (node 1)
    field = ESMF_FieldCreate(mesh=mesh,farray=elementConn(3*elementIds-2),indexflag=ESMF_INDEX_DELOCAL,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="elementConn_n1",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)


    ! save element connectivity (node 2)
    field = ESMF_FieldCreate(mesh=mesh,farray=elementConn(3*elementIds-1),indexflag=ESMF_INDEX_DELOCAL,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="elementConn_n2",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)


    ! save element connectivity (node 3)
    field = ESMF_FieldCreate(mesh=mesh,farray=elementConn(3*elementIds),indexflag=ESMF_INDEX_DELOCAL,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="elementConn_n3",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)

    ! save ice-surface elevation
    field = ESMF_FieldCreate(mesh=mesh,farrayPtr=SurfaceToGEOS,meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
    call ESMF_FieldWrite(field,"mesh.nc",variableName="ice_elevation",rc=rc)
    call ESMF_FieldDestroy(field,rc=rc)


    ! finalize ESMF
    call ESMF_Finalize(rc=rc)
end program main
