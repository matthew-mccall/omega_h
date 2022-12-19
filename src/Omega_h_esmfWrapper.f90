subroutine esmfInit() bind(C, name='esmfInit')
  use iso_c_binding
  use ESMF
  implicit none
  type(ESMF_VM) :: vm
  integer :: finalrc, rc, result
  integer :: localPet, petCount

  write(*,*) 'Fortran fooF'
  finalrc = ESMF_SUCCESS
  call  ESMF_Initialize(vm=vm, defaultlogfilename="fooFortran.Log", &
                    logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  write(*,*) 'Fortran fooF done'
end

!  1.0   4 ------- 3
!        |  \   2  |
!        |    \    |
!        | 1    \  |
!  0.0   1 ------- 2
!
!        Node Id labels at corners
!       Element Id labels in centers
!       (Everything owned by PET 0)

subroutine esmfTestMeshGet(cNodeCoords) bind(C, name='esmfTestMeshGet')
  use iso_c_binding
  use ESMF
  implicit none
  type(c_ptr), value :: cNodeCoords
  real(c_double), pointer :: fNodeCoords(:)
  integer :: rc, localrc
  integer :: numNodes, numTriElems
  integer, allocatable :: nodeIds(:)
  real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
  integer, allocatable :: nodeOwners(:)
  integer, allocatable :: elemIds(:)
  integer, allocatable :: elemTypes(:)
  integer, allocatable :: elemConn(:)
  type(ESMF_Mesh) :: mesh
  write(*,*) 'Fortran esmfTestMeshGet'

  ! Set number of nodes
  numNodes=4

  ! Allocate and fill the node id array.
  allocate(nodeIds(numNodes))
  nodeIds=(/1,2,3,4/)

  ! Allocate and fill node coordinate array.
  allocate(nodeCoords(2*numNodes))
  nodeCoords=(/0.0,0.0, & ! node id 1
               1.0,0.0, & ! node id 2
               0.0,1.0, & ! node id 4
               1.0,1.0 /) ! node id 3

  ! Allocate and fill the node owner array.
  ! Since this Mesh is all on PET 0, it's just set ta all 0.
  allocate(nodeOwners(numNodes))
  nodeOwners=0 ! everything on PET 0

  ! Set the number of each type of element, plus the total number.
  numTriElems=2

  ! Allocate and fill the element id array.
  allocate(elemIds(numTriElems))
  elemIds=(/1,2/)

  ! Allocate and fill the element topology type array.
  allocate(elemTypes(numTriElems))
  elemTypes=(/ESMF_MESHELEMTYPE_TRI, & ! elem id 1
              ESMF_MESHELEMTYPE_TRI/)  ! elem id 2

  ! Allocate and fill the element connection type array.
  ! Note that entries in this array refer to the
  ! positions in the nodeIds, etc. arrays and that
  ! the order and number of entries for each element
  ! reflects that given in the Mesh options
  ! section for the corresponding entry
  ! in the elemTypes array. The number of
  ! entries in this elemConn array is the
  ! number of nodes in a quad. (4) times the
  ! number of quad. elements plus the number
  ! of nodes in a triangle (3) times the number
  ! of triangle elements.
  allocate(elemConn(3*numTriElems))
  elemConn=(/1,2,4, &  ! elem id 1
             2,3,4/)   ! elem id 2

  ! Create Mesh structure in 1 step
  mesh=ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         coordSys=ESMF_COORDSYS_CART, &
         nodeIds=nodeIds, nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, elementIds=elemIds,&
         elementTypes=elemTypes, elementConn=elemConn, &
         rc=localrc)
  if (localrc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call c_f_pointer(cNodeCoords,fNodeCoords,[8])

  call ESMF_MeshGet(mesh, &
    nodeIds=nodeIds, &
    nodeCoords=fNodeCoords, &
    elementConn=elemConn, &
    rc=localrc)
  if (localrc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! After the creation we are through with the arrays, so they may be
  ! deallocated.
  deallocate(nodeIds)
  deallocate(nodeCoords)
  deallocate(nodeOwners)
  deallocate(elemIds)
  deallocate(elemTypes)
  deallocate(elemConn)

  write(*,*) 'Fortran done esmfTestMeshGet'
end subroutine

subroutine esmfFinalize() bind(C, name='esmfFinalize')
  use iso_c_binding
  use ESMF
  implicit none
  integer :: rc
  write(*,*) 'Fortran esmfFinalize'
  call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc) !don't finalize mpi
  write(*,*) 'Fortran esmfFinalize done'
end
