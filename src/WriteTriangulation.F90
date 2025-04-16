subroutine WriteTriangulation(nnodes, ntetrahedra, PNodesList, KTetrahedraList)
    
   use DataTypes

   implicit none
   integer, intent(in) :: nnodes, ntetrahedra

   type(pnode) , dimension(1:nnodes), intent(in) :: PNodesList
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList

   integer :: cont_nodes, cont_elemements
 
   if (myid == root) then
      open  (unit = 8, file = 'triangulation.dat', form = 'formatted')
 
      write (unit = 8, fmt = '(A)') 'TITLE = "TRIANGULATION"'
      write (unit = 8, fmt = '(A)') 'VARIABLES = "X", "Y", "Z"'
      write (unit = 8, fmt = '(A,I6,A,I6,A)') ' ZONE T="TRIANGULATION", DATAPACKING=POINT, NODES=', nnodes, &
                                              ', ELEMENTS=', ntetrahedra,', ZONETYPE=FETETRAHEDRON'
 
      ! Write the coordinates of every interface neighbour node
      do cont_nodes = 1, nnodes
         write (unit = 8, fmt = '(3(g13.6,1x))')  PNodesList(cont_nodes)%x, PNodesList(cont_nodes)%y, PNodesList(cont_nodes)%z 
      end do
 
      ! Write connectivity of the tetrahedra
      do cont_elemements = 1, ntetrahedra
         write (unit = 8, fmt = '(4(i6,1x))')  &
         KTetrahedraList(cont_elemements)%v1%id, KTetrahedraList(cont_elemements)%v2%id, &
         KTetrahedraList(cont_elemements)%v3%id, KTetrahedraList(cont_elemements)%v4%id 
      end do
 
      close (unit = 8)
   end if
 
 end subroutine WriteTriangulation