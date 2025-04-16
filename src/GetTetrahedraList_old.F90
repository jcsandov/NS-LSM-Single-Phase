subroutine GetTetrahedraList( phizero, nnodes, ntetrahedra, NodesList, TetrahedraList, InterfaceNodesID , &
                              TriangulationCase )

   use DataTypes
   use precision

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real (kind = rdf), target, dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero

   ! Array with all the nodes to be used for reinitialisation.It's got the 
   ! target atributte because elements from linked lists points to elements  
   ! of the Tetrahedra List
   type(pnode), target , dimension(:), allocatable, intent(inout) :: NodesList
   
   ! Array with all the tetrahedra to be analised. It's got the target
   ! atributte because elements from linked lists points to elements of 
   ! the Tetrahedra List
   type(tetrahedron), dimension(:), allocatable, intent(inout) :: TetrahedraList
   
   integer, intent(inout) :: nnodes, ntetrahedra

   ! This array contains the ID of every node that belongs to any phase-changing
   ! tetrahedron
   integer, dimension(il:iu,jl:ju,kl:ku), intent(inout):: InterfaceNodesID 

   integer, intent(in) , optional :: TriangulationCase

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer :: i,j,k ! main-loop indexes   
   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   integer :: isearch, jsearch, ksearch ! interface nodes searching range
   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   integer :: ContNodes, ContTetrahedra
        
   ! Phase change flags
   real(kind = rdf) :: SignPhi , phi_aux
   logical :: PhaseChange

   ! Tetrahedra setup variables

   integer, dimension(0:7,1:3) :: vertex_ijk_idx
   integer :: tet_cell_vertex
   integer :: iloc, jloc, kloc
   integer , dimension(4) :: ilocv, jlocv, klocv
   integer :: i1, i2, i3, i4, j1, j2, j3, j4, k1, k2, k3, k4

   integer, dimension(1:3,0:7) :: VertexOffset ! = 1 or 0 for the eight vertices
   integer, dimension(1:6,1:4) :: TetrahedronAux

   ! Tetrahedron loop
   integer :: CellVertexLoop, TetrahedronLoop, TetrahedronVertexLoop
   integer ::  iOffset, jOffset, kOffset ! local offset indexes

   ! Isosurface loops
   integer :: TriangleLoop, TriangleVertexLoop


   ! loop boundaries definition
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp
            
   ! processes on the domain boundaries
   
   if (myback  == mpi_proc_null)  i_mysta = il + igp + 1
   if (myleft  == mpi_proc_null)  j_mysta = jl + jgp + 1
   if (mydown  == mpi_proc_null)  k_mysta = kl + kgp + 1
   
   if (myfront == mpi_proc_null)  i_myend = iu - igp - 1
   if (myright == mpi_proc_null)  j_myend = ju - jgp - 1
   if (myup    == mpi_proc_null)  k_myend = ku - kgp - 1

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! CELL AND TETRAHEDRON NUMERATION
   !
   ! Cell vertices nomenclature for tetrahedron construction
   ! v0 = v000, v1 = v100, v2 = v110, v3 = v010
   ! v4 = v001, v5 = v101, v6 = v111, v7 = v011
   ! where v101 means v at (i+1,j,k+1)
   
   !                     vertex position
   !                     0 1 2 3 4 5 6 7
   VertexOffset(1,:) = (/0,1,1,0,0,1,1,0/) ! i - index offset
   VertexOffset(2,:) = (/0,0,1,1,0,0,1,1/) ! j - index offset
   VertexOffset(3,:) = (/0,0,0,0,1,1,1,1/) ! k - index offset
  
   ! Index notation for tetrahedra vertices (arbitrary convention)
   ! for example, the tetrahedron 1 is formed by the vertices 
   ! 0, 5, 1 and 6 of the cell:
                
   !    i,j+1,k+1 ----i+1,j+1,k+1    
   !     /|(7)           /|(6)                            
   !    / |             / |                                  
   ! i,j,k+1-------i+1,j,k+1                   
   !   |(4)           |(5)|                                   
   !   |  |           |   |                                
   !   |  |           |   |                                
   !   |  i,j+1,k-----|-i+1,j+1,k           
   !   | /(3)         |  /(2)                 
   !   |/             | /
   ! i,j,k-----------i+1,j,k
   !  (0)               (1)


   if ( present( TriangulationCase ) )  then

      print *, 'In TetrahedraList, CASE # ', TriangulationCase

      select case ( TriangulationCase )

         case(1) ! Main diagonal from v0 --> v6

            TetrahedronAux(1,:) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
            TetrahedronAux(2,:) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
            TetrahedronAux(3,:) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
            TetrahedronAux(4,:) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
            TetrahedronAux(5,:) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
            TetrahedronAux(6,:) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6

         case(2) ! Main diagonal from v1 --> v7

            TetrahedronAux(1,:) = (/1, 6, 2, 7/) ! cell vertices of the tetrahedron # 1
            TetrahedronAux(2,:) = (/1, 2, 3, 7/) ! cell vertices of the tetrahedron # 2
            TetrahedronAux(3,:) = (/1, 3, 0, 7/) ! cell vertices of the tetrahedron # 3
            TetrahedronAux(4,:) = (/1, 0, 4, 7/) ! cell vertices of the tetrahedron # 4
            TetrahedronAux(5,:) = (/1, 4, 5, 7/) ! cell vertices of the tetrahedron # 5
            TetrahedronAux(6,:) = (/1, 5, 6, 7/) ! cell vertices of the tetrahedron # 6

         case(3)! Main diagonal from v2 --> v4

            TetrahedronAux(1,:) = (/2, 7, 3, 4/) ! cell vertices of the tetrahedron # 1
            TetrahedronAux(2,:) = (/2, 3, 0, 4/) ! cell vertices of the tetrahedron # 2
            TetrahedronAux(3,:) = (/2, 0, 1, 4/) ! cell vertices of the tetrahedron # 3
            TetrahedronAux(4,:) = (/2, 1, 5, 4/) ! cell vertices of the tetrahedron # 4
            TetrahedronAux(5,:) = (/2, 5, 6, 4/) ! cell vertices of the tetrahedron # 5
            TetrahedronAux(6,:) = (/2, 6, 7, 4/) ! cell vertices of the tetrahedron # 6

         case(4) ! Main diagonal from v3 --> v5

            TetrahedronAux(1,:) = (/3, 4, 0, 5/) ! cell vertices of the tetrahedron # 1
            TetrahedronAux(2,:) = (/3, 0, 1, 5/) ! cell vertices of the tetrahedron # 2
            TetrahedronAux(3,:) = (/3, 1, 2, 5/) ! cell vertices of the tetrahedron # 3
            TetrahedronAux(4,:) = (/3, 2, 6, 5/) ! cell vertices of the tetrahedron # 4
            TetrahedronAux(5,:) = (/3, 6, 7, 5/) ! cell vertices of the tetrahedron # 5
            TetrahedronAux(6,:) = (/3, 7, 4, 5/) ! cell vertices of the tetrahedron # 6

      end select

   else

      TetrahedronAux(1,:) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
      TetrahedronAux(2,:) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
      TetrahedronAux(3,:) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
      TetrahedronAux(4,:) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
      TetrahedronAux(5,:) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
      TetrahedronAux(6,:) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6

   end if

   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! FIRST LOOP: Count total number of Pnodes and Tetrahedra and set Pnodes IDs
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   do k = kl, ku !k_mysta-1, k_myend+1
   do j = jl, ju !j_mysta-1, j_myend+1
   do i = il, iu !i_mysta-1, i_myend+1
               
      ! ------------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k is next to the phase change
      ! ------------------------------------------------------------------------------
      ! if there are changes in the signs of the phi functions in the neighbourhood
      ! of a node, then there is a phase change.
      ! ------------------------------------------------------------------------------
               
      kminus = max( k-1 , k_mysta-1 )
      kplus  = min( k+1 , k_myend+1 )
   
      jminus = max( j-1 , j_mysta-1 )
      jplus  = min( j+1 , j_myend+1 )
   
      iminus = max( i-1 , i_mysta-1 )
      iplus  = min( i+1 , i_myend+1 )
   
      SignPhi = sign( one , phizero(i,j,k) ) ! either +1.0 or -1.0
      
      PhaseChange = .false.
      
      if ( abs( phizero(i,j,k) ) < eps_sims ) PhaseChange = .true. 

      if (.not.(PhaseChange) ) then
         
         searchloop: do ksearch = kminus , kplus
                     do jsearch = jminus , jplus
                     do isearch = iminus , iplus
   
            if (  SignPhi * phizero (isearch,jsearch,ksearch) < zero       .or. &
                  abs( phizero (isearch,jsearch,ksearch) )    < eps_sims            ) then
   
               PhaseChange = .true.
               exit searchloop
   
            end if
      
         end do
         end do
         end do searchloop 

      end if

      ! Interface neighbourhood cells
      if( PhaseChange ) then
   
         ! Number of neighbour nodes to the free-surface
   
         ! Counting how many valid tetrahedron are inside every cell
         if ( i<=i_myend .and. j<=j_myend .and. k<=k_myend ) then

            do CellVertexLoop = 0,7
                     
               iOffset = VertexOffset(1,CellVertexLoop)
               jOffset = VertexOffset(2,CellVertexLoop)
               kOffset = VertexOffset(3,CellVertexLoop)
                     
               ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
               ! cell vertices (0 to 7) 
         
               vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
               vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
               vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1
         
            end do
         
            ! we check phase change over every tetrahedron contained in the cell
         
            do TetrahedronLoop = 1,6
         
               ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
               ! of the TetrahedronVertexLoop tetrahedron
         
               tet_cell_vertex = TetrahedronAux(TetrahedronLoop,1)
         
               ! get the i,j,k indexes of the vtx_tet vertex of the tet 
                           ! tetrahedron
         
               ilocv( 1 ) = vertex_ijk_idx(tet_cell_vertex,1)
               jlocv( 1 ) = vertex_ijk_idx(tet_cell_vertex,2)
               klocv( 1 ) = vertex_ijk_idx(tet_cell_vertex,3)
                           
               phi_aux = phizero( ilocv( 1 ) , jlocv( 1 ) , klocv( 1 ) )  

               SignPhi = sign( one , phi_aux )                        
   
               PhaseChange = .false.

               if ( abs( phi_aux ) < eps_sims ) PhaseChange = .true.
                        
               do TetrahedronVertexLoop = 2,4  ! four tetrahedron vertices
            
                  ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
                  ! of the TetrahedronVertexLoop tetrahedron
                     
                  tet_cell_vertex = TetrahedronAux(TetrahedronLoop,TetrahedronVertexLoop)
                     
                  ! get the i,j,k indexes of the vtx_tet vertex of the tet 
                  ! tetrahedron
                     
                  ilocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,1)
                  jlocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,2)
                  klocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,3)

                  phi_aux =  phizero ( ilocv( TetrahedronVertexLoop ) , & 
                                       jlocv( TetrahedronVertexLoop ) , &
                                       klocv( TetrahedronVertexLoop )      )
                                    
                  if ( SignPhi * phi_aux < zero .or. abs( phi_aux ) < eps_sims ) PhaseChange = .true.
               
               end do
         
               ! if there are changes in the signs of the phi functions in the neighbourhood
               ! of a node, then there is a phase change.
         
               if ( PhaseChange ) then
            
                  ntetrahedra = ntetrahedra + 1
                        
                  ! I tag the node i,j,k with a 1 value to sum all of them at the
                  ! end of the loop
                  
                  do TetrahedronVertexLoop = 1,4 ! four tetrahedron verticies
                              
                     ! I set again the tag of every changing phase node to 1
   
                     InterfaceNodesID( ilocv( TetrahedronVertexLoop ), &
                                       jlocv( TetrahedronVertexLoop ), &
                                       klocv( TetrahedronVertexLoop )    ) = 1
                           
                  end do
   
               end if ! Phase change within the tetrahedron
            end do ! TetrahedronLoop
         end if ! The cell construction is not outbounded
      end if !  Phase within the cell at i,j,k 
   end do
   end do
   end do

   ! the total amount of valid nodes is the sum of the 1 tags
   nnodes = sum( InterfaceNodesID )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! SECOND LOOP: Knowing the total number of Pnodes and Tetrahedra, I allocate
   ! memory for the Pnodes and Tetrahedra Lists and proceed to set all the Pnodes
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   allocate( NodesList(1:nnodes)             )
   allocate( TetrahedraList(1:ntetrahedra)   )

   ContNodes      = 1
   ContTetrahedra = 1

   do k = kl , ku ! k_mysta-1, k_myend+1
   do j = jl , ju ! j_mysta-1, j_myend+1
   do i = il , iu ! i_mysta-1, i_myend+1
                           
      if( InterfaceNodesID(i,j,k) > 0 ) then
   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
         ! Interface pnodes storing
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
         ! Set the attributes of the pnodes

         ! id
         InterfaceNodesID(i,j,k)  = ContNodes
         NodesList(ContNodes)%id  = ContNodes

         ! spatial position
         NodesList(ContNodes)%x  = x(i,j,k)
         NodesList(ContNodes)%y  = y(i,j,k)
         NodesList(ContNodes)%z  = z(i,j,k)

         ! computational indexes
         NodesList(ContNodes)%i  = i
         NodesList(ContNodes)%j  = j
         NodesList(ContNodes)%k  = k

         ! Phi value. Its a pointer to modify those specific locations
         ! only
         NodesList(ContNodes)%phi_value  => phizero(i,j,k)

         ! Distance to the free surface for first step of reinitialistion
         ! method. It's initialised as a high number to be update a few
         ! steps ahead

         NodesList(ContNodes)%SDistanceFreeSurface  = 99.9_rdf

         ! phih* value initialisation
         NodesList(ContNodes)%phi_corrected  = zero

         ! # of associated tetrahedrons
         NodesList(ContNodes)%NumAsocTetrahedra = 0

         ! update the nodes counter
         ContNodes = ContNodes + 1

      end if
   end do
   end do
   end do

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! THIRD LOOP: Having set all the Pnodes, now I set all the tetrahedra  
   ! and set the pointers of their vertices to the corresponding Pnodes
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   do k = k_mysta-1, k_myend
   do j = j_mysta-1, j_myend
   do i = i_mysta-1, i_myend
                            
      !if( InterfaceNodesID(i,j,k) > 0 ) then
               
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Tetrahedra identification and values setup
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      do CellVertexLoop = 0,7
               
         iOffset = VertexOffset(1,CellVertexLoop)
         jOffset = VertexOffset(2,CellVertexLoop)
         kOffset = VertexOffset(3,CellVertexLoop)
               
         ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
         ! cell vertices (0 to 7) 
   
         vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
         vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
         vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1
   
      end do
   
      ! we check phase change over every tetrahedron contained in the cell
   
      do TetrahedronLoop = 1,6
   
         ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
         ! of the TetrahedronVertexLoop tetrahedron
         
         tet_cell_vertex = TetrahedronAux(TetrahedronLoop,1)
         
         ! get the i,j,k indexes of the vtx_tet vertex of the tet tetrahedron
         
         iloc = vertex_ijk_idx(tet_cell_vertex,1)
         jloc = vertex_ijk_idx(tet_cell_vertex,2)
         kloc = vertex_ijk_idx(tet_cell_vertex,3)
                           
         SignPhi = sign( one , phizero( iloc , jloc , kloc ) )                        
   
         PhaseChange = .false.

         if ( abs( phizero( iloc , jloc , kloc ) ) < eps_sims ) PhaseChange = .true.

         do TetrahedronVertexLoop = 2,4  ! four tetrahedron verticies
   
            ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
            ! of the TetrahedronVertexLoop tetrahedron
   
            tet_cell_vertex = TetrahedronAux(TetrahedronLoop,TetrahedronVertexLoop)
   
            ! get the i,j,k indexes of the vtx_tet vertex of the tet 
            ! tetrahedron
   
            iloc = vertex_ijk_idx(tet_cell_vertex,1)
            jloc = vertex_ijk_idx(tet_cell_vertex,2)
            kloc = vertex_ijk_idx(tet_cell_vertex,3)
                     
            if (  SignPhi * phizero (iloc , jloc , kloc) < zero .or. &
                  abs( phizero (iloc , jloc , kloc) )    < eps_sims     ) then 
               
               PhaseChange = .true.
            
            end if
            
         end do 
   
         ! if there are changes in the signs of the phi functions in the neighbourhood
         ! of a node, then there is a phase change.
   
         if ( PhaseChange ) then

            ! Four tetrahedron vertices indexes

            i1 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,1),1)
            j1 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,1),2)
            k1 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,1),3)
   
            i2 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,2),1)
            j2 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,2),2)
            k2 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,2),3)
   
            i3 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,3),1)
            j3 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,3),2)
            k3 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,3),3)
   
            i4 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,4),1)
            j4 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,4),2)
            k4 = vertex_ijk_idx(TetrahedronAux(TetrahedronLoop,4),3)

            ! Setting tetrahedron ID

            TetrahedraList(ContTetrahedra)%id = ContTetrahedra

            ! Setting tetrahedron vertices pointing to pnodes
            TetrahedraList( ContTetrahedra )%v1 => NodesList( InterfaceNodesID( i1,j1,k1 ) )
            TetrahedraList( ContTetrahedra )%v2 => NodesList( InterfaceNodesID( i2,j2,k2 ) )
            TetrahedraList( ContTetrahedra )%v3 => NodesList( InterfaceNodesID( i3,j3,k3 ) )
            TetrahedraList( ContTetrahedra )%v4 => NodesList( InterfaceNodesID( i4,j4,k4 ) )

            ! Update the number of associated tetrahedra for all the pnodes
            ! that make up the vertices and set the pointers to the 
            ! corresponding tetrahedron

            ! vertex 1
            NodesList(InterfaceNodesID(i1,j1,k1))%NumAsocTetrahedra = &
            NodesList(InterfaceNodesID(i1,j1,k1))%NumAsocTetrahedra + 1

            NodesList(InterfaceNodesID(i1,j1,k1))%AssociatedTetrahedraID(  &
            NodesList(InterfaceNodesID(i1,j1,k1))%NumAsocTetrahedra ) = TetrahedraList(ContTetrahedra)%id

            ! vertex 2

            NodesList(InterfaceNodesID(i2,j2,k2))%NumAsocTetrahedra = &
            NodesList(InterfaceNodesID(i2,j2,k2))%NumAsocTetrahedra + 1

            NodesList(InterfaceNodesID(i2,j2,k2))%AssociatedTetrahedraID(  &
            NodesList(InterfaceNodesID(i2,j2,k2))%NumAsocTetrahedra ) = TetrahedraList(ContTetrahedra)%id

            ! vertex 3
            NodesList(InterfaceNodesID(i3,j3,k3))%NumAsocTetrahedra = &
            NodesList(InterfaceNodesID(i3,j3,k3))%NumAsocTetrahedra + 1

            NodesList(InterfaceNodesID(i3,j3,k3))%AssociatedTetrahedraID(  &
            NodesList(InterfaceNodesID(i3,j3,k3))%NumAsocTetrahedra ) = TetrahedraList(ContTetrahedra)%id

            ! vertex 4
            NodesList(InterfaceNodesID(i4,j4,k4))%NumAsocTetrahedra = &
            NodesList(InterfaceNodesID(i4,j4,k4))%NumAsocTetrahedra + 1

            NodesList(InterfaceNodesID(i4,j4,k4))%AssociatedTetrahedraID(  &
            NodesList(InterfaceNodesID(i4,j4,k4))%NumAsocTetrahedra ) = TetrahedraList(ContTetrahedra)%id

            ! Update the tetrahedra counter
            ContTetrahedra  = ContTetrahedra + 1

            end if ! Phase change tetrahedron

         end do ! TetrahedronLoop
      !end if ! InterfaceNodesID(i,j,k) > 0
   end do ! i
   end do ! j
   end do ! k

end subroutine GetTetrahedraList