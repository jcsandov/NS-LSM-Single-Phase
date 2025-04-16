module DataTypes

   use precision

   implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Derived types for better handling of tetrahedra data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   type pnode
   
      ! A pnode is a data structure which corresponds to a neighbour node
      ! to the interface. In addition, these nodes corresponds to the vertices
      ! of the tetrahedra obtained by the triangulation. These vertices are
      ! declared in the tetrahedron class as pointers to one of these pnodes
      ! for modifying SDistanceFreeSurface and phi_corrected
   
      integer          :: id = 0
      real(kind = rdf) :: x,y,z                  ! physical coordinates
      integer          :: i,j,k                  ! computational coordinates
      real(kind = rdf) :: phi_value              ! phi value        
      real(kind = rdf) :: SDistanceFreeSurface   ! signed distance to the fs 
      real(kind = rdf) :: WaterVolume       
      real(kind = rdf) :: phi_corrected          ! volume correction step
      
      integer          :: NumAsocTetrahedra = 0  ! number of tetrahedra with a 
                                                 ! vertex at this pnode
   
      ! we have an array with pointers to all the tetrahedra that have any 
      ! vertex at this pnode. The maximum value can be 24
   
      ! to avoid the problem of the circular dependency between this attribute
      ! becuase of its type (the typer tetrahedron is declared below), we take
      ! advantage of the Fortran 2003 feature which allows to declare an 
      ! allocatable variable of a type that hasn't been declared yet.

!      type(tetrahedron), dimension(:), allocatable :: AssociatedTetrahedra

      ! to simplify the references, I defined AssociatedTetrahedra just as
      ! an ID number

      integer, dimension(24) :: AssociatedTetrahedraID

      ! TO DO: I have to add, somehow, an atribute which is a list of
      ! pointer to tetrahedra, for the 1st and 2nd neighbours steps 
      
      real(kind = rdf) :: xiI = zero ! ξI - 3rd step reinitialisation algorithm

      !type(pnode), pointer :: NextPNode => null()

      !integer :: NodeType = 0 ! 1 = pnode
                               ! 2 = second-neighbour node
                               ! 3 = fully-submerged node


   end type pnode


   type triangle
      type(pnode) :: v1, v2, v3 ! vertices 
      integer :: TriangleType = 0 ! 1 : the triangle faces the free-surface
                                  ! 2 : the triangle is intersected by the free-surface
                                  ! 0 : otherwise
   end type triangle

   type tetrahedron
   
      integer                  :: id
      type(pnode),   pointer   :: v1,v2,v3,v4 ! every vertex is a 
                                              ! pointer pointing a pnode

      integer                  :: nisosurfaces = 0   ! # of reconstructed isosurfaces                                             
      type(triangle)           :: isosurface1
      type(triangle)           :: isosurface2

      real(kind = rdf) :: IsosurfaceArea = zero
      real(kind = rdf) :: WaterVolume
         
      real(kind = rdf) :: etaK ! piecewise constant function
   
      ! This is a logical variable that tells us if a tetrahedron
      ! belongs to the internal nodes region of a processor or not.
      ! It's gonna be important when we calculate the total volume
      ! enclosed by tetrahedrons in the global correction step of the
      ! geometric reinitialisation

      logical :: InternalTetrahedron

      !type(tetrahedron), pointer :: NextKTetrahedron => null()
      
      ! TO DO: check if this field is necessary
      !type(tetrahedron), pointer :: next => null()
      !type(tetrahedron), pointer :: NextAssociatedTetrahedron => null()

      !type(triangle), pointer :: Face1, Face2, Face3, Face4

      !integer :: TetrahedronType = 0 ! 1 = phase change tetrahedron
                                      ! 2 = fully-submerged tetrahedron
                                      ! 3 = air-phase tetrahedron

   end type tetrahedron

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! LINKED LISTS STRUCTURES
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

!   type nodes_llist
!      type(pnode), pointer :: head => null()      
!      type(pnode), pointer :: tail => null()      
!      type(pnode), pointer :: next => null()      
!
!      integer :: nnodes = 0
!   end type nodes_llist
!
!   type triangles_llist
!      type(triangle), pointer :: head => null()      
!      type(triangle), pointer :: tail => null()      
!      type(triangle), pointer :: next => null()      
!      integer :: ntriangles = 0
!   end type triangles_llist
!
!   type tetrahedra_llist
!      type(tetrahedron), pointer :: head => null()
!      type(tetrahedron), pointer :: tail => null()
!      type(tetrahedron), pointer :: next => null()
!      integer :: ntetrahedra = 0
!   end type tetrahedra_llist

!   contains

!   subroutine add_tetrahedron(list, new)
!   ! add a new item to the end of the list
!
!      type(tetrahedra_llist) :: list
!
!      ! dummy argument
!      type(tetrahedron), pointer :: new
!
!      list%ntetrahedra = list%ntetrahedra + 1
!
!      ! Check to see if the list is empty or not
!
!      if ( associated(list%head) ) then
!
!         ! the list IS NOT empty
!   
!         ! The last item of the list (tail) points its "next" attribute to 
!         ! the new element
!         list%tail%next => new
!   
!         ! The "next" attribute of the new element is nullified
!         nullify(new%next)
!
!         ! reset tail pointer to the new element
!
!         list%tail => new
!
!      else 
!      
!         ! the list IS empty
!
!         list%head => new
!         list%tail => new
!         nullify(list%tail%next)
!
!      end if 
!
!   end subroutine add_tetrahedron


end module DataTypes