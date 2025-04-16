module LinkedListModule

   use DataTypes

   implicit none

   contains


   subroutine InitNodesList(list)
   ! initialise the empty list

      type(nodes_llist), intent(inout) :: list

      nullify(list%head)
      nullify(list%tail)
      list%nnodes = 0

   end subroutine InitNodesList


   subroutine InitTetrahedraList(list)
   ! initialise the empty list

      type(tetrahedra_llist), intent(inout) :: list

      nullify(list%head)
      nullify(list%tail)
      list%ntetrahedra = 0

   end subroutine InitTetrahedraList



   ! I need different subroutines to append linked lists, because
   ! the same entity (node, tetrahedron, etc) might point different
   ! entities depending on which list we are examining


   subroutine AddPNode(list, new)
   ! add a new item to the end of the list

      type(nodes_llist) :: list

      ! dummy argument
      type(pnode), pointer :: new

      list%nnodes = list%nnodes + 1

      ! Check to see if the list is empty or not

      if ( associated(list%head) ) then

         ! the list IS NOT empty
   
         ! The last item of the list (tail) points its "next" attribute to 
         ! the new element
         list%tail%NextPNode => new
   
         ! The "next" attribute of the new element is nullified
         nullify(new%NextPNode)

         ! reset tail pointer to the new element

         list%tail => new

      else 
      
         ! the list IS empty

         list%head => new
         list%tail => new
         nullify(list%tail%NextPNode)

      end if 

   end subroutine AddPNode


   
   subroutine AddKTetrahedron(list, new)
   ! add a new item to the end of the list

      type(tetrahedra_llist) :: list

      ! dummy argument
      type(tetrahedron), pointer :: new

      list%ntetrahedra = list%ntetrahedra + 1

      ! Check to see if the list is empty or not

      if ( associated(list%head) ) then

         ! the list IS NOT empty
   
         ! The last item of the K Tetrahedron list (tail) points its 
         ! "next" attribute to the new element

         list%tail%NextKTetrahedron => new
   
         ! The "next" attribute of the new element is nullified
         nullify(new%NextKTetrahedron)

         ! reset tail pointer to the new element

         list%tail => new

      else 
      
         ! the list IS empty

         list%head => new
         list%tail => new
         nullify(list%tail%NextKTetrahedron)

      end if 

   end subroutine AddKTetrahedron



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




   subroutine DeleteTetrahedron(list, first)
   
   ! Return a pointer to the first item in the linked list,
   ! and remove it from the list

      ! Dummy arguments
      type(tetrahedra_llist) :: list
      type(tetrahedron), pointer :: first

      ! Check if the list is empty
      if(associated(list%head)) then
      ! List is not empty. Now we check if there's more than
      ! one item in the list

         if ( associated(list%head%next) ) then
         ! More than one item in the list

            first     => list%head ! return pointer to the first item (dummy action)
            list%head => list%head%next ! Now the head is the next item
         else
         ! Only one item in the list

            first => list%head
            nullify(list%head, list%tail) ! List is now empty

         end if

      else ! list is empty

         nullify(first) ! return no element

      end if

   end subroutine DeleteTetrahedron


   subroutine ListTetrahedraID(list)

      ! dummy argument 
      type(tetrahedra_llist) :: list

      ! local variable
      type(tetrahedron), pointer :: ptr

      if (.not. associated(list%head) ) then
         ! list is empty

         !print *, 'List is empty ... wrrrrrr'

      else

         ptr => list%head

         ! loop to print tetrahedra id

         do
            
            !print *, 'Tetrahedron # ', ptr%id

            ! set pointer to the next item
            ptr => ptr%next

            ! Exit loop if there are no more items in the list
            if (.not. associated(ptr) ) exit

         end do
      end if

   end subroutine ListTetrahedraID


end module LinkedListModule