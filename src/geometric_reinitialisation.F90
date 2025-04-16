subroutine geometric_reinitialisation( phizero , phicorrected, TotalWaterVolume, InterfaceNodesID , TriangulationCase )

use global_app
use global_mpi
use global_debug
!use global_lsm, only: phi_outputiter
use DataTypes
!use LinkedListModule

implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! input-output arguments
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real (kind = rdf), target, dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(inout) :: phicorrected
real (kind = rdf) , intent(in) :: TotalWaterVolume

! List with the ID of every pnode. The rest of the array is set to 0
integer, dimension(il:iu,jl:ju,kl:ku), intent(inout) :: InterfaceNodesID 

! debug iterators
integer :: i,j,k
integer, intent(in) , optional :: TriangulationCase

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! local  variables
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: nnodes, ntetrahedra, nvertices

! Array with all the nodes to be used for reinitialisation.It's got the 
! target atributte because elements from linked lists points to elements  
! of the Tetrahedra List
type(pnode), target , dimension(:), allocatable :: PNodesList

! Array with all the tetrahedra to be analised. It's got the target
! atributte because elements from linked lists points to elements of 
! the Tetrahedra List
type(tetrahedron), dimension(:), allocatable :: KTetrahedraList

nnodes      =  0
ntetrahedra =  0
nvertices   =  0

call GetTetrahedraList  ( phizero           ,                    &
                          nnodes            ,  ntetrahedra     , &
                          PNodesList        ,  KTetrahedraList , &
                          InterfaceNodesID                     , &
                          TriangulationCase                      &
                         )


!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
!
! I'll comment this writing while I'm implementing the parallelisation for the
! geometric redistancing
!
! if( mod( iteraciontiempo , phi_outputiter) == 0 ) then   
! 
!    call WriteTriangulation        ( nnodes     , ntetrahedra     ,  &
!                                     PNodesList , KTetrahedraList        )
! end if
!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


call GetIsosurfaceAndDistances ( nnodes            , ntetrahedra      , & 
                                 PNodesList        , KTetrahedraList  , &  
                                 InterfaceNodesID                         )

if ( GlobalMassCorrection ) then

    call PiecewiseConstantFunction ( ntetrahedra, KTetrahedraList )


    call OrthogonalProjection      (  nnodes      , ntetrahedra    , &
                                      PNodesList  , KTetrahedraList      )

end if

! If GlobalMassCorrection is false, then the Pnodes correction is just
! the signed distance
call NeighboursPhiCorrection  (  nnodes      , ntetrahedra     , & 
                                 PNodesList  , KTetrahedraList , &
                                 InterfaceNodesID              , &
                                 TotalWaterVolume                   )

! deallocate big arrays

deallocate ( PNodesList      )
deallocate ( KTetrahedraList )

end subroutine geometric_reinitialisation