subroutine TotalWaterVolume( phizero           , &
                             VolumeKTetrahedra , VolumeSecondTetrahedra , VolumeBulkCells , &
                             AdvectionNodes    , &
                             WriteGlobalVolume , &
                             TriangulationCase                                                 )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This subroutine computes the total volume of fluid given a spatial phi distribution
   ! To do that, we compute the volume of every water cell. For the interface cells, the
   ! volume is obtained by the application of Marching Tetrahedron Method.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   use precision
   use DataTypes
   use TetrahedronMethods

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! INPUT arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real(kind = rdf) , target, dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero
   logical, optional :: WriteGlobalVolume 
   integer, intent(in) , optional :: TriangulationCase
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! OUTPUT arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   ! TotalWaterVolume = VolumeKTetrahedra + VolumeSecondTetrahedra + VolumeBulkCells
   real(kind = rdf) , intent(out) :: VolumeKTetrahedra, VolumeSecondTetrahedra, VolumeBulkCells
   integer, dimension(il:iu,jl:ju,kl:ku), intent(out), optional :: AdvectionNodes 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   logical :: FileExist
   ! Phase change flags
   real(kind = rdf) :: SignPhi
   logical :: PhaseChangeNode
   real(kind = rdf), dimension(0:7,1:3) :: VerticesCoordinates
   real(kind = rdf), dimension(0:7)     :: PhiVertices
   real(kind = rdf) :: CellWaterVolume , SecondTetrahedronsVolume , PhaseChangeTetrahedronsVolume
   integer :: TriangulationCaseAux

   logical :: BlankingFlag

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Auxiliary counters 
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer :: i, j, k
   integer :: isearch, jsearch, ksearch ! interface nodes searching range
   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   
   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   integer :: ista, iend, jsta, jend, ksta, kend

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

   if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1
   if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1
   if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1

   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp

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
          
   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- 

   ! Initialisation of the id of the pnodes
   if( present( AdvectionNodes ) ) AdvectionNodes = 0

   ! Check if there is an impossed triangulation, or we just use the first one   
   TriangulationCaseAux = 1
   if ( present( TriangulationCase ) ) TriangulationCaseAux = TriangulationCase

   VolumeKTetrahedra       = zero
   VolumeSecondTetrahedra  = zero
   VolumeBulkCells         = zero

   ! Loop over the whole domain
   do k =  k_mysta , k_myend ! kl, ku !k_mysta-1, k_myend+1
   do j =  j_mysta , j_myend ! jl, ju !j_mysta-1, j_myend+1
   do i =  i_mysta , i_myend ! il, iu !i_mysta-1, i_myend+1

      ! Flag to identify if I'm within the blanking region
      BlankingFlag = .false.

      if (nblk /= 0) then
         
         do nb = 1, nblk
            
            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. &     
                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then

               BlankingFlag = .true.
         
            end if

         end do

      end if

      ! This keep AdvectionNodes = 0 within (boundaries included) the blank region
      if ( BlankingFlag ) cycle

      ! ------------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k is next to the phase change for advection
      ! ------------------------------------------------------------------------------
      ! if there are changes in the signs of the phi functions in the neighbourhood
      ! of a node, then there is a phase change.
      ! ------------------------------------------------------------------------------
               
      kminus = max( k-1 , ksta ) ! kl )
      kplus  = min( k+1 , kend ) ! ku )
   
      jminus = max( j-1 , jsta ) ! jl )
      jplus  = min( j+1 , jend ) ! ju )
   
      iminus = max( i-1 , ista ) ! il )
      iplus  = min( i+1 , iend ) ! iu )
      
      ! Sign of the base node (i,j,k)
      SignPhi = sign( one , phizero(i,j,k) ) ! either +1.0 or -1.0
      
      PhaseChangeNode = .false.
      
      if ( abs( phizero(i,j,k) ) < eps_sims ) PhaseChangeNode = .true. 

      if ( .not.( PhaseChangeNode ) ) then
         
         searchloop1: &
         do ksearch = kminus , kplus
         do jsearch = jminus , jplus
         do isearch = iminus , iplus
   
            BlankingFlag = .false.

            if (nblke /= 0) then
               do nb = 1, nblke
                  
                  if ( isearch > le_blk_ia(1,nb) .and. isearch < le_blk_ib(1,nb) .and. &     
                       jsearch > le_blk_ja(1,nb) .and. jsearch < le_blk_jb(1,nb) .and. &
                       ksearch > le_blk_ka(1,nb) .and. ksearch < le_blk_kb(1,nb) ) then
      
                     BlankingFlag = .true.
               
                  end if
      
               end do
      
            end if
      
            if ( BlankingFlag ) cycle 

            ! Check if there are sign changes of ϕ within any cell around i,j,k
            if (  SignPhi * phizero (isearch,jsearch,ksearch)   < zero      .or. &
                       abs( phizero (isearch,jsearch,ksearch) ) < eps_sims ) then
   
               PhaseChangeNode = .true.
               exit searchloop1
   
            end if
      
         end do
         end do
         end do searchloop1 

      end if

      ! Interface neighbourhood cells
      if( PhaseChangeNode ) then

         ! Any node somehow surrounded by the free-surface is gonna be advected         
         !if( present( AdvectionNodes ) ) AdvectionNodes(i,j,k) = 1
         if( present( AdvectionNodes ) ) then
            AdvectionNodes( iminus : iplus , jminus : jplus , kminus : kplus ) = 1
         end if

      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Water volume components computation
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
      ! Coordinates of the eight vertices of the cell
      VerticesCoordinates(0,:) = (/ x( i   , j   , k   ) , y( i   , j   , k   ) , z( i   , j   , k   ) /)
      VerticesCoordinates(1,:) = (/ x( i+1 , j   , k   ) , y( i+1 , j   , k   ) , z( i+1 , j   , k   ) /)
      VerticesCoordinates(2,:) = (/ x( i+1 , j+1 , k   ) , y( i+1 , j+1 , k   ) , z( i+1 , j+1 , k   ) /)
      VerticesCoordinates(3,:) = (/ x( i   , j+1 , k   ) , y( i   , j+1 , k   ) , z( i   , j+1 , k   ) /)
   
      VerticesCoordinates(4,:) = (/ x( i   , j   , k+1 ) , y( i   , j   , k+1 ) , z( i   , j   , k+1 ) /)
      VerticesCoordinates(5,:) = (/ x( i+1 , j   , k+1 ) , y( i+1 , j   , k+1 ) , z( i+1 , j   , k+1 ) /)
      VerticesCoordinates(6,:) = (/ x( i+1 , j+1 , k+1 ) , y( i+1 , j+1 , k+1 ) , z( i+1 , j+1 , k+1 ) /)
      VerticesCoordinates(7,:) = (/ x( i   , j+1 , k+1 ) , y( i   , j+1 , k+1 ) , z( i   , j+1 , k+1 ) /)
   
      ! Phi values at those vertices
      PhiVertices(0) = phizero( i   , j   , k   )
      PhiVertices(1) = phizero( i+1 , j   , k   )
      PhiVertices(2) = phizero( i+1 , j+1 , k   )
      PhiVertices(3) = phizero( i   , j+1 , k   )
   
      PhiVertices(4) = phizero( i   , j   , k+1 )
      PhiVertices(5) = phizero( i+1 , j   , k+1 )
      PhiVertices(6) = phizero( i+1 , j+1 , k+1 )
      PhiVertices(7) = phizero( i   , j+1 , k+1 )
   
      CellWaterVolume               = zero
      SecondTetrahedronsVolume      = zero
      PhaseChangeTetrahedronsVolume = zero
   
      ! FULL WATER CELL
      if ( all ( PhiVertices  > eps_sims ) ) then
   
         CellWaterVolume = HexahedronVolume( VerticesCoordinates )
         
      ! FULL AIR CELL
      else if ( all ( PhiVertices  < -eps_sims ) ) then
   
         cycle
   
      ! PHASE CHANGE CELL CANDIDATE
      else
   
         call CellAnalysis( PhiVertices                     , &
                            VerticesCoordinates             , &
                            TriangulationCaseAux            , &
                            PhaseChangeTetrahedronsVolume   , &
                            SecondTetrahedronsVolume          &
                           )
      end if
   
      VolumeKTetrahedra       = VolumeKTetrahedra      + PhaseChangeTetrahedronsVolume
      VolumeSecondTetrahedra  = VolumeSecondTetrahedra + SecondTetrahedronsVolume
      VolumeBulkCells         = VolumeBulkCells        + CellWaterVolume

   end do
   end do
   end do


end subroutine TotalWaterVolume
