subroutine GetTetrahedraList( phizero           , &
                              nnodes            , &
                              ntetrahedra       , &
                              PNodesList        , &
                              KTetrahedraList   , &
                              InterfaceNodesID  , &
                              TriangulationCase )

   use DataTypes
   use precision
   use TetrahedronMethods

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero

   ! Array with all the nodes to be used for reinitialisation.It's got the 
   ! target atributte because elements from linked lists points to elements  
   ! of the Tetrahedra List
   type(pnode), target , dimension(:), allocatable, intent(inout) :: PNodesList
   
   ! Array with all the tetrahedra to be analised. It's got the target
   ! atributte because elements from linked lists points to elements of 
   ! the Tetrahedra List
   type(tetrahedron), dimension(:), allocatable, intent(inout) :: KTetrahedraList
   
   integer, intent(inout) :: nnodes, ntetrahedra

   ! This array contains the ID of every node that belongs to any phase-changing
   ! tetrahedron
   integer, dimension(il:iu,jl:ju,kl:ku), intent(inout):: InterfaceNodesID 

   integer, intent(in) , optional :: TriangulationCase

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   integer, dimension(il:iu,jl:ju,kl:ku) :: KTetrahedronsWithinCellList 
   real(kind = rdf), dimension(0:7) :: PhiVertices , NodesListID
   integer :: NodesIDAux1, NodesIDAux2, NodesIDAux3, NodesIDAux4
   integer :: nPhaseChangeTetrahedrons , nCellVerticesPnodes



   integer :: i,j,k ! main-loop indexes   
   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   integer :: ista, iend, jsta, jend, ksta, kend

   integer :: isearch, jsearch, ksearch ! interface nodes searching range
   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   integer :: ContNodes, ContTetrahedra
        
   ! Blanking Flag
   logical :: BlankingFlag

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
   integer, dimension(1:6,1:4) :: TetrahedronCellVertices

   ! Tetrahedron loop
   integer :: CellVertexLoop, TetrahedronLoop, TetrahedronVertexLoop


   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 
   
   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 
   
   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp


   ! loop boundaries definition
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp
            
   ! processes on the domain boundaries
   ! I only keep these correction, because of the way I build the cells. From
   ! i,j,k to i+1,j+1,k+1   

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
   !  
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

      !print *, 'In GetTetrahedraList, CASE # ', TriangulationCase

      select case ( TriangulationCase )
      
         case(1) ! Main diagonal from v0 --> v6
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6
      
         case(2) ! Main diagonal from v1 --> v7
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/1, 6, 2, 7/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/1, 2, 3, 7/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/1, 3, 0, 7/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/1, 0, 4, 7/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/1, 4, 5, 7/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/1, 5, 6, 7/) ! cell vertices of the tetrahedron # 6
      
         case(3)! Main diagonal from v2 --> v4
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/2, 7, 3, 4/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/2, 3, 0, 4/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/2, 0, 1, 4/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/2, 1, 5, 4/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/2, 5, 6, 4/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/2, 6, 7, 4/) ! cell vertices of the tetrahedron # 6
      
         case(4) ! Main diagonal from v3 --> v5
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/3, 4, 0, 5/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/3, 0, 1, 5/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/3, 1, 2, 5/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/3, 2, 6, 5/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/3, 6, 7, 5/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/3, 7, 4, 5/) ! cell vertices of the tetrahedron # 6
      
      end select
      
   else
      
      TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
      TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
      TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
      TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
      TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
      TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6

   end if

   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! FIRST LOOP: Count total number of Pnodes and Tetrahedra and set Pnodes IDs
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   KTetrahedronsWithinCellList = 0

   do k = ksta, kend-1 !kl, ku-1 !k_mysta-1, k_myend+1
   do j = jsta, jend-1 !jl, ju-1 !j_mysta-1, j_myend+1
   do i = ista, iend-1 !il, iu-1 !i_mysta-1, i_myend+1

      BlankingFlag = .false.

      if ( nblke/=0 ) then
      
         do nb = 1,nblke
            
            if ( i >= le_blk_ia(1,nb) .and. i < le_blk_ib(1,nb) .and. &
                 j >= le_blk_ja(1,nb) .and. j < le_blk_jb(1,nb) .and. &
                 k >= le_blk_ka(1,nb) .and. k < le_blk_kb(1,nb) ) then

               BlankingFlag = .true.
         
            end if
         
         end do
      
      end if

      if ( BlankingFlag ) cycle

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Phase change checking
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Phi values at those vertices
      PhiVertices(0) = phizero( i   , j   , k   )
      PhiVertices(1) = phizero( i+1 , j   , k   )
      PhiVertices(2) = phizero( i+1 , j+1 , k   )
      PhiVertices(3) = phizero( i   , j+1 , k   )

      PhiVertices(4) = phizero( i   , j   , k+1 )
      PhiVertices(5) = phizero( i+1 , j   , k+1 )
      PhiVertices(6) = phizero( i+1 , j+1 , k+1 )
      PhiVertices(7) = phizero( i   , j+1 , k+1 )

      ! FULL WATER CELL OR FULL AIR CELL
      if ( all ( PhiVertices  > eps_sims ) .or. all ( PhiVertices  < -eps_sims ) ) then
      
         cycle

      ! PHASE CHANGE CELL CANDIDATE
      else

         nPhaseChangeTetrahedrons = 0 ! all the bits are set to 0
         nCellVerticesPnodes      = 0 ! all the bits are set to 0
         
         call CellPnodeTetrahedronIdentification (  PhiVertices                    , &
                                                    TriangulationCase              , &
                                                    nCellVerticesPnodes            , &
                                                    nPhaseChangeTetrahedrons         &
                                                 )

         ! Update the Pnodes distribution over the domain 
         if( btest ( nCellVerticesPnodes , 0 ) )  InterfaceNodesID( i   , j   , k   ) = 1
         if( btest ( nCellVerticesPnodes , 1 ) )  InterfaceNodesID( i+1 , j   , k   ) = 1
         if( btest ( nCellVerticesPnodes , 2 ) )  InterfaceNodesID( i+1 , j+1 , k   ) = 1
         if( btest ( nCellVerticesPnodes , 3 ) )  InterfaceNodesID( i   , j+1 , k   ) = 1
         if( btest ( nCellVerticesPnodes , 4 ) )  InterfaceNodesID( i   , j   , k+1 ) = 1
         if( btest ( nCellVerticesPnodes , 5 ) )  InterfaceNodesID( i+1 , j   , k+1 ) = 1
         if( btest ( nCellVerticesPnodes , 6 ) )  InterfaceNodesID( i+1 , j+1 , k+1 ) = 1
         if( btest ( nCellVerticesPnodes , 7 ) )  InterfaceNodesID( i   , j+1 , k+1 ) = 1

         ! nPhaseChangeTetrahedrons is an integer that has its bits set in the postion
         ! correspoing to the number of a changing-phase tetrahedron. If it's 0, the cell
         ! doesn't have changing-phase tetrahedrons inside
         KTetrahedronsWithinCellList(i,j,k) = nPhaseChangeTetrahedrons

         ! Update the number of KTetrahedrons ( popcnt returns the number of bits 
         ! set (’1’ bits) in the binary representation of nPhaseChangeTetrahedrons )

         ntetrahedra = ntetrahedra + popcnt ( nPhaseChangeTetrahedrons )
      
      end if                     

   end do
   end do
   end do

   ! the total amount of valid nodes is the sum of the 1 tags
   nnodes = sum( InterfaceNodesID )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! SECOND LOOP: Knowing the total number of Pnodes and Tetrahedra, I allocate
   ! memory for the Pnodes and Tetrahedra Lists and proceed to set all the Pnodes
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   allocate( PNodesList      ( 1:nnodes      ) )
   allocate( KTetrahedraList ( 1:ntetrahedra ) )

   ContNodes      = 1
   ContTetrahedra = 1

   do k = ksta, kend !kl , ku ! k_mysta-1, k_myend+1
   do j = jsta, jend !jl , ju ! j_mysta-1, j_myend+1
   do i = ista, iend !il , iu ! i_mysta-1, i_myend+1
                           
      if( InterfaceNodesID(i,j,k) > 0 ) then
   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
         ! Interface pnodes storing
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
         ! Set the attributes of the pnodes

         ! id
         InterfaceNodesID(i,j,k)   = ContNodes
         PNodesList(ContNodes)%id  = ContNodes

         ! spatial position
         PNodesList(ContNodes)%x  = x(i,j,k)
         PNodesList(ContNodes)%y  = y(i,j,k)
         PNodesList(ContNodes)%z  = z(i,j,k)

         ! computational indexes within the processor 
         PNodesList(ContNodes)%i  = i
         PNodesList(ContNodes)%j  = j
         PNodesList(ContNodes)%k  = k

         ! Initial phi value
         PNodesList(ContNodes)%phi_value = phizero(i,j,k)

         ! Distance to the free surface for first step of reinitialistion
         ! method. It's initialised as a high number to be update a few
         ! steps ahead
         PNodesList(ContNodes)%SDistanceFreeSurface  = 99.9_rdf

         ! phih* value initialisation
         PNodesList(ContNodes)%phi_corrected  = zero

         ! # of associated tetrahedrons
         PNodesList(ContNodes)%NumAsocTetrahedra = 0

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

   do k = ksta, kend-1
   do j = jsta, jend-1
   do i = ista, iend-1
                            
      if ( KTetrahedronsWithinCellList(i,j,k) > 0 ) then
      
         ! Pnodes IDs. They might be 0 or something > 0 if (i,j,k) node is a Pnode
         NodesListID(0) = InterfaceNodesID ( i   , j   , k   )
         NodesListID(1) = InterfaceNodesID ( i+1 , j   , k   )
         NodesListID(2) = InterfaceNodesID ( i+1 , j+1 , k   )
         NodesListID(3) = InterfaceNodesID ( i   , j+1 , k   )
   
         NodesListID(4) = InterfaceNodesID ( i   , j   , k+1 )
         NodesListID(5) = InterfaceNodesID ( i+1 , j   , k+1 )
         NodesListID(6) = InterfaceNodesID ( i+1 , j+1 , k+1 )
         NodesListID(7) = InterfaceNodesID ( i   , j+1 , k+1 )

         ! I check the bits of KTetrahedronsWithinCellList(i,j,k) to see which tetrahedrons
         ! within the cell are effectively KTetrahedrons

         do TetrahedronLoop = 1,6

            if( btest ( KTetrahedronsWithinCellList(i,j,k) , TetrahedronLoop ) ) then

               ! Setting tetrahedron ID

               KTetrahedraList(ContTetrahedra)%id = ContTetrahedra

               ! Set the flag variable to know if the tetrahedron is within the internal 
               ! nodes of the procs or it's in the ghost nodes region

               KTetrahedraList( ContTetrahedra )%InternalTetrahedron = .true.

               if ( i < i_mysta .or. i > i_myend .or. &
                    j < j_mysta .or. j > j_myend .or. &
                    k < k_mysta .or. k > k_myend ) then

                  KTetrahedraList( ContTetrahedra )%InternalTetrahedron = .false.

               end if

               ! ID of the four nodes of the next tetrahedron to be indexed
               NodesIDAux1 = NodesListID( TetrahedronCellVertices( TetrahedronLoop , 1) )
               NodesIDAux2 = NodesListID( TetrahedronCellVertices( TetrahedronLoop , 2) )
               NodesIDAux3 = NodesListID( TetrahedronCellVertices( TetrahedronLoop , 3) )
               NodesIDAux4 = NodesListID( TetrahedronCellVertices( TetrahedronLoop , 4) )

               ! Point the Tetrahedron vertices to the respective Pnodes in PNodesList
               KTetrahedraList( ContTetrahedra )%v1 => PNodesList( NodesIDAux1 )
               KTetrahedraList( ContTetrahedra )%v2 => PNodesList( NodesIDAux2 )
               KTetrahedraList( ContTetrahedra )%v3 => PNodesList( NodesIDAux3 )
               KTetrahedraList( ContTetrahedra )%v4 => PNodesList( NodesIDAux4 )

               ! Update the number of associated tetrahedra for all the pnodes
               ! that make up the vertices and set the pointers to the 
               ! corresponding tetrahedron
   
               ! VERTEX 1 # of associated tetrahedrons update
               PNodesList( NodesIDAux1 )%NumAsocTetrahedra = PNodesList( NodesIDAux1 )%NumAsocTetrahedra + 1

               ! VERTEX 1 : we append the ID of the new indexed tetrahedron to the list of tetrahedros that
               !            have the node NodesIDAux1 as a vertex
               PNodesList( NodesIDAux1 )%AssociatedTetrahedraID( PNodesList( NodesIDAux1 )%NumAsocTetrahedra ) = &
               KTetrahedraList(ContTetrahedra)%id
   
               ! VERTEX 2
               PNodesList( NodesIDAux2 )%NumAsocTetrahedra = PNodesList( NodesIDAux2 )%NumAsocTetrahedra + 1

               PNodesList( NodesIDAux2 )%AssociatedTetrahedraID( PNodesList( NodesIDAux2 )%NumAsocTetrahedra ) = &
               KTetrahedraList(ContTetrahedra)%id

               ! VERTEX 3
               PNodesList( NodesIDAux3 )%NumAsocTetrahedra = PNodesList( NodesIDAux3 )%NumAsocTetrahedra + 1

               PNodesList( NodesIDAux3 )%AssociatedTetrahedraID( PNodesList( NodesIDAux3 )%NumAsocTetrahedra ) = &
               KTetrahedraList(ContTetrahedra)%id

               ! VERTEX 4
               PNodesList( NodesIDAux4 )%NumAsocTetrahedra = PNodesList( NodesIDAux4 )%NumAsocTetrahedra + 1

               PNodesList( NodesIDAux4 )%AssociatedTetrahedraID( PNodesList( NodesIDAux4 )%NumAsocTetrahedra ) = &
               KTetrahedraList(ContTetrahedra)%id

               ! Update the tetrahedra counter
               ContTetrahedra  = ContTetrahedra + 1
            
            end if

         end do

      end if                     

   end do ! i
   end do ! j
   end do ! k

end subroutine GetTetrahedraList