subroutine OrthogonalProjection( nnodes          , &
                                 ntetrahedra     , &
                                 PNodesList      , &
                                 KTetrahedraList      )

   use DataTypes
   use precision

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! Array with all the nodes to be used for reinitialisation.It's got the 
   ! target atributte because elements from linked lists points to elements  
   ! of the Tetrahedra List

   integer, intent(in) :: nnodes, ntetrahedra
   type(pnode), dimension(:), allocatable, intent(inout) :: PNodesList
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
   
   ! This subroutine perform the 3rd step of the geometric reinitilisation
   ! Step

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer :: I, K, NI

   NI = 0

   do I = 1, nnodes

      if(PNodesList(I)%NumAsocTetrahedra > 0) then
         NI = PNodesList(I)%NumAsocTetrahedra
         do K = 1, NI ! loop over the simplices with a vertex at I

            PNodesList(I)%xiI =    PNodesList(I)%xiI & 
                                 + KTetrahedraList(PNodesList(I)%AssociatedTetrahedraID(K))%etaK/NI

         end do
      
      end if
   end do


end subroutine OrthogonalProjection