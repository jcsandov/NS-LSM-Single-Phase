subroutine write_triangulation(nnodes, ntetrahedra, NodesList, TetrahedraList)
    
   use DataTypes

   implicit none
   integer, intent(inout) :: nnodes, ntetrahedra

   type(pnode) , dimension(:), allocatable, intent(inout) :: NodesList
   type(tetrahedron), dimension(:), allocatable, intent(inout) :: TetrahedraList

   integer :: cont_nodes, cont_elemements
 
   if (myid == root) then
      open  (unit = 8, file = 'triangulation.dat', form = 'formatted')
 
      write (unit = 8, fmt = '(A)') 'TITLE = "TRIANGULATION"'
      write (unit = 8, fmt = '(A)') 'VARIABLES = "X", "Y", "Z"'
      write (unit = 8, fmt = '(A,I6,A,I6,A)') ' ZONE T="TRIANGULATION", DATAPACKING=POINT, NODES=',nnodes,', ELEMENTS=',ntetrahedra,', ZONETYPE=FETETRAHEDRON'
 
      ! Write the coordinates of every interface neighbour node
      do cont_nodes = 1, nnodes
         write (unit = 8, fmt = '(3(g13.6,1x))')  NodesList(cont_nodes)%x, NodesList(cont_nodes)%y, NodesList(cont_nodes)%z 
      end do
 
      ! Write connectivity of the tetrahedra
      do cont_elemements = 1, ntetrahedra
         write (unit = 8, fmt = '(4(i6,1x))')  TetrahedraList(cont_elemements)%v1%id, TetrahedraList(cont_elemements)%v2%id, &
                                                TetrahedraList(cont_elemements)%v3%id, TetrahedraList(cont_elemements)%v4%id 
      end do
 
      close (unit = 8)
   end if
 
 end subroutine write_triangulation