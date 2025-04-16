subroutine TotalWaterVolumeSussman( phizero , VolumeKTetrahedra, VolumeSecondTetrahedra, VolumeBulkCells )

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
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real(kind = rdf) , target, dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero
   
   ! TotalWaterVolume = VolumeKTetrahedra + VolumeSecondTetrahedra + VolumeBulkCells
   real(kind = rdf) , intent(out) :: VolumeKTetrahedra, VolumeSecondTetrahedra, VolumeBulkCells


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   logical :: FileExist
   real(kind = rdf) :: VolumeAux   
   real(kind = rdf), dimension(3) :: vertA, vertB, vertC, vertD, vertE, vertF, vertG, vertH
   real(kind = rdf), dimension(3) :: AC, AF, AH, BE, CH, DB, ED, FC, GA, GB, GC, GD, GE, GF, GH, HF
   real(kind = rdf), dimension(8) :: xvertices, yvertices, zvertices

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Marching Tetrahedron Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   type(pnode), target, dimension(:), allocatable :: NodesListAux
   
   ! Array with all the tetrahedra to be analised. It's got the target
   ! atributte because elements from linked lists points to elements of 
   ! the Tetrahedra List
   type(tetrahedron), dimension(:), allocatable :: TetrahedraListAux
   
   integer :: nnodes, ntetrahedra

   ! This array contains the ID of every node that belongs to any phase-changing
   ! tetrahedron
   integer, dimension(il:iu,jl:ju,kl:ku) :: InterfaceNodesIDAux 

   integer :: isearch, jsearch, ksearch ! interface nodes searching range
   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   integer :: ContNodes, ContTetrahedra

   ! Phase change flags
   real(kind = rdf) :: SignPhi, phi_aux, sum_phi_aux, VolumeSecondTetrahedraAux
   logical :: PhaseChangeNode, PhaseChangeCell , PhaseChangeTetrahedron
                                                           
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

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   real(kind = rdf), dimension(4,3)    :: VerticesCoordinates
   real(kind = rdf), dimension (4)     :: phi_tetrahedron
   real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   real(kind = rdf)                    :: IsosurfaceArea
   real(kind = rdf)                    :: WaterVolumeAux
   real(kind = rdf), dimension(4)      :: VerticesDistances
   real(kind = rdf), dimension(4,3)    :: vertices_coordinates_aux

   logical :: SecondTetrahedron ! Single phase tetrahedron within a phase-change cell  
   logical :: WaterVA, WaterVB, WaterVC, WaterVD, WaterVE, WaterVF, WaterVG, WaterVH 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Auxiliary counters for distance computation
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer :: i, j, k, ii, jj, kk
   integer :: Ktet

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
          
   TetrahedronAux(1,:) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
   TetrahedronAux(2,:) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
   TetrahedronAux(3,:) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
   TetrahedronAux(4,:) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
   TetrahedronAux(5,:) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
   TetrahedronAux(6,:) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6

   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- 

   ! Initialisation of the id of the pnodes
   nnodes = 0
   ntetrahedra = 0
   InterfaceNodesIDAux  = 0
   VolumeSecondTetrahedra = zero

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
      
      PhaseChangeNode = .false.
      
      if ( abs( phizero(i,j,k) ) < eps_sims ) PhaseChangeNode = .true. 

      if (.not.(PhaseChangeNode) ) then
         
         searchloop: do ksearch = kminus , kplus
                     do jsearch = jminus , jplus
                     do isearch = iminus , iplus
   
            if (  SignPhi * phizero (isearch,jsearch,ksearch)   < zero      .or. &
                       abs( phizero (isearch,jsearch,ksearch) ) < eps_sims         ) then
   
               PhaseChangeNode = .true.
               exit searchloop
   
            end if
      
         end do
         end do
         end do searchloop 

      end if

      ! Interface neighbourhood cells
      if( PhaseChangeNode ) then

         if ( i<=i_myend .and. j<=j_myend .and. k<=k_myend ) then

            CellVertexLoop = 0

            iOffset = VertexOffset(1,CellVertexLoop)
            jOffset = VertexOffset(2,CellVertexLoop)
            kOffset = VertexOffset(3,CellVertexLoop)
                  
            ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
            ! cell vertices (0 to 7) 
         
            vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
            vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
            vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1

            SignPhi = sign( one , phizero( i + iOffset, j + jOffset , k + kOffset ) )

            PhaseChangeCell = .false.

            if ( abs( phizero( i+iOffset, j+jOffset, k+kOffset ) ) < eps_sims ) PhaseChangeCell = .true.

            do CellVertexLoop = 1,7
                     
               iOffset = VertexOffset(1,CellVertexLoop)
               jOffset = VertexOffset(2,CellVertexLoop)
               kOffset = VertexOffset(3,CellVertexLoop)
                     
               ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
               ! cell vertices (0 to 7) 
         
               vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
               vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
               vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1

               if ( SignPhi * phizero( i+iOffset, j+jOffset, k+kOffset )   < zero  .or. & 
                         abs( phizero( i+iOffset, j+jOffset, k+kOffset ) ) < eps_sims     ) then
                   
                  PhaseChangeCell = .true.
               
               end if

            end do

            if ( PhaseChangeCell ) then
         
            ! we check phase change over every tetrahedron contained in the cell

               do TetrahedronLoop = 1,6

                  PhaseChangeTetrahedron = .false. 
                  SecondTetrahedron = .true.     
                  TetrahedronVertexLoop = 1

            
                  ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
                  ! of the TetrahedronVertexLoop tetrahedron
            
                  tet_cell_vertex = TetrahedronAux( TetrahedronLoop , TetrahedronVertexLoop )
            
                  ! get the i,j,k indexes of the vtx_tet vertex of the tet 
                              ! tetrahedron
            
                  ilocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,1)
                  jlocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,2)
                  klocv( TetrahedronVertexLoop ) = vertex_ijk_idx(tet_cell_vertex,3)
                              
                  phi_aux =  phizero ( ilocv( TetrahedronVertexLoop ) , & 
                                       jlocv( TetrahedronVertexLoop ) , &
                                       klocv( TetrahedronVertexLoop )      )   
   
                  SignPhi = sign( one , phi_aux )                        

                  ! This condition is to ensure the Second Tetrahedron candidate
                  !  holds POSITIVE values of ϕ only
                  if ( phi_aux <= eps_sims ) SecondTetrahedron = .false.

                  ! If one of the vertices lies on the free surface, it's a phase
                  ! changing tetrahedron      
                  if ( abs( phi_aux ) < eps_sims ) PhaseChangeTetrahedron = .true.
                                 
                  do TetrahedronVertexLoop = 2,4 ! four tetrahedron vertices
               
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
                                       
                     if ( phi_aux <= eps_sims ) SecondTetrahedron = .false.

                     ! If there's sign change, it's a phase changing tetrahedron
                     if ( SignPhi * phi_aux < zero .or. abs( phi_aux ) < eps_sims ) PhaseChangeTetrahedron = .true.
                  
                  end do
            
                  ! if there are changes in the signs of the phi functions in the neighbourhood
                  ! of a node, then there is a phase change.            

                  if ( SecondTetrahedron .and. .not.(PhaseChangeTetrahedron) ) then

                     do TetrahedronVertexLoop = 1,4

                        tet_cell_vertex = TetrahedronAux(TetrahedronLoop,TetrahedronVertexLoop)
                           
                        ! get the i,j,k indexes of the vtx_tet vertex of the tet 
                        ! tetrahedron
                           
                        ii = vertex_ijk_idx(tet_cell_vertex,1)
                        jj = vertex_ijk_idx(tet_cell_vertex,2)
                        kk = vertex_ijk_idx(tet_cell_vertex,3)

                        vertices_coordinates_aux(TetrahedronVertexLoop,:) = (/x(ii,jj,kk), y(ii,jj,kk), z(ii,jj,kk)/)

                     end do 

                     VolumeSecondTetrahedraAux = GetTetrahedronVolume( vertices_coordinates_aux )                     

                     VolumeSecondTetrahedra = VolumeSecondTetrahedra + VolumeSecondTetrahedraAux

                  end if


                  if ( PhaseChangeTetrahedron ) then
               
                     ntetrahedra = ntetrahedra + 1
                           
                     ! I tag the node i,j,k with a 1 value to sum all of them at the
                     ! end of the loop
                     
                     do TetrahedronVertexLoop = 1,4 ! four tetrahedron verticies
                                 
                        ! I set again the tag of every non-changing phase node to 1
      
                        InterfaceNodesIDAux( ilocv( TetrahedronVertexLoop ), &
                                             jlocv( TetrahedronVertexLoop ), &
                                             klocv( TetrahedronVertexLoop )    ) = 1
                     end do
      
                  end if ! Phase change within the tetrahedron
               end do ! TetrahedronLoop
            end if ! Phase change at the cell with its first corner at i,j,k
         end if ! The cell construction is not outbounded
      end if !  Phase next to the node i,j,k 
   end do
   end do
   end do

   ! the total amount of valid nodes is the sum of the 1 tags
   nnodes = sum(InterfaceNodesIDAux)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! SECOND LOOP: Knowing the total number of Pnodes and Tetrahedra, I allocate
   ! memory for the Pnodes and Tetrahedra Lists and proceed to set all the Pnodes
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   allocate(NodesListAux(1:nnodes))
   allocate(TetrahedraListAux(1:ntetrahedra))

   ContNodes      = 1
   ContTetrahedra = 1

   do k = kl , ku ! k_mysta-1, k_myend+1
   do j = jl , ju ! j_mysta-1, j_myend+1
   do i = il , iu ! i_mysta-1, i_myend+1
                           
      if( InterfaceNodesIDAux(i,j,k) > 0 ) then
   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
         ! Interface pnodes storing
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
         ! Set the attributes of the pnodes

         ! id
         InterfaceNodesIDAux(i,j,k)  = ContNodes
         NodesListAux(ContNodes)%id  = ContNodes

         ! spatial position
         NodesListAux(ContNodes)%x  = x(i,j,k)
         NodesListAux(ContNodes)%y  = y(i,j,k)
         NodesListAux(ContNodes)%z  = z(i,j,k)

         ! computational indexes
         NodesListAux(ContNodes)%i  = i
         NodesListAux(ContNodes)%j  = j
         NodesListAux(ContNodes)%k  = k

         ! Phi value. Its a pointer to modify those specific locations
         ! only
         NodesListAux(ContNodes)%phi_value  => phizero(i,j,k)

         ! # of associated tetrahedrons
         NodesListAux(ContNodes)%NumAsocTetrahedra = 0

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
                            
      if( InterfaceNodesIDAux(i,j,k) > 0 ) then
               
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         ! Tetrahedra identification and values setup
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

         CellVertexLoop = 0

         iOffset = VertexOffset(1,CellVertexLoop)
         jOffset = VertexOffset(2,CellVertexLoop)
         kOffset = VertexOffset(3,CellVertexLoop)
               
         ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
         ! cell vertices (0 to 7) 
         
         vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
         vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
         vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1

         SignPhi = sign( one , phizero( i + iOffset, j + jOffset , k + kOffset ) )
         PhaseChangeCell = .false.

         do CellVertexLoop = 1,7
                  
            iOffset = VertexOffset(1,CellVertexLoop)
            jOffset = VertexOffset(2,CellVertexLoop)
            kOffset = VertexOffset(3,CellVertexLoop)
                  
            ! vertex_ijk_idx contians the i,j,k indexes for each of the 8
            ! cell vertices (0 to 7) 
         
            vertex_ijk_idx(CellVertexLoop,1) = i + iOffset ! i or i+1
            vertex_ijk_idx(CellVertexLoop,2) = j + jOffset ! j or j+1
            vertex_ijk_idx(CellVertexLoop,3) = k + kOffset ! k or k+1

            if ( SignPhi * phizero( i+iOffset, j+jOffset, k+kOffset )   < zero  .or. & 
                      abs( phizero( i+iOffset, j+jOffset, k+kOffset ) ) < eps_sims     ) then
                
               PhaseChangeCell = .true.
            
            end if

         end do

      
         if ( PhaseChangeCell ) then

         ! we check phase change over every tetrahedron contained in the cell
      
            do TetrahedronLoop = 1,6
      
               sum_phi_aux = zero
         
               ! tet_cell_vertex is the local index (0 to 7) of the vtx_tet vertex  
               ! of the TetrahedronVertexLoop tetrahedron
               
               tet_cell_vertex = TetrahedronAux(TetrahedronLoop,1)
               
               ! get the i,j,k indexes of the vtx_tet vertex of the tet tetrahedron
               
               iloc = vertex_ijk_idx(tet_cell_vertex,1)
               jloc = vertex_ijk_idx(tet_cell_vertex,2)
               kloc = vertex_ijk_idx(tet_cell_vertex,3)
                                 
               SignPhi = sign( one , phizero( iloc , jloc , kloc ) )                        
                        
               PhaseChangeTetrahedron = .false.
      
               if ( abs( phizero( iloc , jloc , kloc ) ) < eps_sims ) PhaseChangeTetrahedron = .true.
      
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
                     
                     PhaseChangeTetrahedron = .true.
                  
                  end if
                        
               end do 
         
               ! if there are changes in the signs of the phi functions in the neighbourhood
               ! of a node, then there is a phase change.
         
               !if ( sum_phi_aux > eps_sims ) PhaseChangeTetrahedron = .true.
      
               if ( PhaseChangeTetrahedron ) then
      
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
      
                  TetrahedraListAux(ContTetrahedra)%id = ContTetrahedra
      
                  ! Setting tetrahedron vertices pointing to pnodes
                  TetrahedraListAux( ContTetrahedra )%v1 => NodesListAux( InterfaceNodesIDAux( i1,j1,k1 ) )
                  TetrahedraListAux( ContTetrahedra )%v2 => NodesListAux( InterfaceNodesIDAux( i2,j2,k2 ) )
                  TetrahedraListAux( ContTetrahedra )%v3 => NodesListAux( InterfaceNodesIDAux( i3,j3,k3 ) )
                  TetrahedraListAux( ContTetrahedra )%v4 => NodesListAux( InterfaceNodesIDAux( i4,j4,k4 ) )
      
                  ! Update the number of associated tetrahedra for all the pnodes
                  ! that make up the vertices and set the pointers to the 
                  ! corresponding tetrahedron
      
                  ! vertex 1
                  NodesListAux(InterfaceNodesIDAux(i1,j1,k1))%NumAsocTetrahedra = &
                  NodesListAux(InterfaceNodesIDAux(i1,j1,k1))%NumAsocTetrahedra + 1
      
                  NodesListAux(InterfaceNodesIDAux(i1,j1,k1))%AssociatedTetrahedraID(  &
                  NodesListAux(InterfaceNodesIDAux(i1,j1,k1))%NumAsocTetrahedra ) = TetrahedraListAux(ContTetrahedra)%id
      
                  ! vertex 2
      
                  NodesListAux(InterfaceNodesIDAux(i2,j2,k2))%NumAsocTetrahedra = &
                  NodesListAux(InterfaceNodesIDAux(i2,j2,k2))%NumAsocTetrahedra + 1
      
                  NodesListAux(InterfaceNodesIDAux(i2,j2,k2))%AssociatedTetrahedraID(  &
                  NodesListAux(InterfaceNodesIDAux(i2,j2,k2))%NumAsocTetrahedra ) = TetrahedraListAux(ContTetrahedra)%id
      
                  ! vertex 3
                  NodesListAux(InterfaceNodesIDAux(i3,j3,k3))%NumAsocTetrahedra = &
                  NodesListAux(InterfaceNodesIDAux(i3,j3,k3))%NumAsocTetrahedra + 1
      
                  NodesListAux(InterfaceNodesIDAux(i3,j3,k3))%AssociatedTetrahedraID(  &
                  NodesListAux(InterfaceNodesIDAux(i3,j3,k3))%NumAsocTetrahedra ) = TetrahedraListAux(ContTetrahedra)%id
      
                  ! vertex 4
                  NodesListAux(InterfaceNodesIDAux(i4,j4,k4))%NumAsocTetrahedra = &
                  NodesListAux(InterfaceNodesIDAux(i4,j4,k4))%NumAsocTetrahedra + 1
      
                  NodesListAux(InterfaceNodesIDAux(i4,j4,k4))%AssociatedTetrahedraID(  &
                  NodesListAux(InterfaceNodesIDAux(i4,j4,k4))%NumAsocTetrahedra ) = TetrahedraListAux(ContTetrahedra)%id
      
                  ! Update the tetrahedra counter
                  ContTetrahedra  = ContTetrahedra + 1
      
               end if ! Phase change tetrahedron
            end do ! TetrahedronLoop
         end if
      end if ! InterfaceNodesIDAux(i,j,k) > 0
   end do ! i
   end do ! j
   end do ! k

   ! Marching Tetrahedron Step
   do Ktet = 1, ntetrahedra

      VerticesCoordinates(1,:) = (/ TetrahedraListAux(Ktet)%v1%x, &
                                    TetrahedraListAux(Ktet)%v1%y, &
                                    TetrahedraListAux(Ktet)%v1%z /)  

      VerticesCoordinates(2,:) = (/ TetrahedraListAux(Ktet)%v2%x, &
                                    TetrahedraListAux(Ktet)%v2%y, &
                                    TetrahedraListAux(Ktet)%v2%z /)  
      
      VerticesCoordinates(3,:) = (/ TetrahedraListAux(Ktet)%v3%x, &
                                    TetrahedraListAux(Ktet)%v3%y, &
                                    TetrahedraListAux(Ktet)%v3%z /)  

      VerticesCoordinates(4,:) = (/ TetrahedraListAux(Ktet)%v4%x, &
                                    TetrahedraListAux(Ktet)%v4%y, &
                                    TetrahedraListAux(Ktet)%v4%z /)  


      phi_tetrahedron(1) = TetrahedraListAux(Ktet)%v1%phi_value
      phi_tetrahedron(2) = TetrahedraListAux(Ktet)%v2%phi_value
      phi_tetrahedron(3) = TetrahedraListAux(Ktet)%v3%phi_value
      phi_tetrahedron(4) = TetrahedraListAux(Ktet)%v4%phi_value

      ! TO DO: MarchingTetrahedron could receive pointers to nodes
      ! instead of vertices coordinates and phi_tetrahedron values
      ! (just doing the previous steps internally)

      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 phi_tetrahedron        , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 IsosurfaceArea         , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )

      TetrahedraListAux(Ktet)%WaterVolume = WaterVolumeAux

   end do   

   ! I look over the whole domain. If ϕ>0 and I'm on a pnode, then the volume is the one obtained by TetrahedronMethods
   ! if ϕ(i,j,k)>0 and all my neighbours ϕ(i±1, j±1, k±1)>0, then the volume is the cell volume
   
!   TotalWaterVolumeVar = zero
 
   VolumeKTetrahedra      = zero
   VolumeBulkCells        = zero   

   do i = i_mysta-1, i_myend
   do j = j_mysta-1, j_myend
   do k = k_mysta-1, k_myend
   
      ! we check if I'm in a full-water cell
   
      WaterVA = phizero( i   , j   , k   ) > eps_sims 
      WaterVB = phizero( i+1 , j   , k   ) > eps_sims 
      WaterVC = phizero( i+1 , j+1 , k   ) > eps_sims 
      WaterVD = phizero( i   , j+1 , k   ) > eps_sims 
      WaterVE = phizero( i   , j   , k+1 ) > eps_sims 
      WaterVF = phizero( i+1 , j   , k+1 ) > eps_sims 
      WaterVG = phizero( i+1 , j+1 , k+1 ) > eps_sims 
      WaterVH = phizero( i   , j+1 , k+1 ) > eps_sims 

      if ( WaterVA .and. WaterVB .and. WaterVC .and. WaterVD .and. WaterVE .and. WaterVF .and. WaterVG .and. WaterVH  ) then

         !   NODE DISTRIBUTION OF THE COMPUTATIONAL CELL
         !
         !    i,j+1,k+1 ----i+1,j+1,k+1    
         !     /|(H)           /|(G)                            
         !    / |             / |                                  
         ! i,j,k+1-------i+1,j,k+1                   
         !   |(E)           |(F)|                                   
         !   |  |           |   |                                
         !   |  |           |   |                                
         !   |  i,j+1,k-----|-i+1,j+1,k           
         !   | /(D)         |  /(C)                 
         !   |/             | /
         ! i,j,k-----------i+1,j,k
         !  (A)               (B)

         ! 8 vertices coordinates
         vertA(:) = (/ x( i   , j   , k   ), y( i   , j   , k   ), z( i   , j   , k   ) /)
         vertB(:) = (/ x( i+1 , j   , k   ), y( i+1 , j   , k   ), z( i+1 , j   , k   ) /)
         vertC(:) = (/ x( i+1 , j+1 , k   ), y( i+1 , j+1 , k   ), z( i+1 , j+1 , k   ) /)
         vertD(:) = (/ x( i   , j+1 , k   ), y( i   , j+1 , k   ), z( i   , j+1 , k   ) /)
         vertE(:) = (/ x( i   , j   , k+1 ), y( i   , j   , k+1 ), z( i   , j   , k+1 ) /)
         vertF(:) = (/ x( i+1 , j   , k+1 ), y( i+1 , j   , k+1 ), z( i+1 , j   , k+1 ) /)
         vertG(:) = (/ x( i+1 , j+1 , k+1 ), y( i+1 , j+1 , k+1 ), z( i+1 , j+1 , k+1 ) /)
         vertH(:) = (/ x( i   , j+1 , k+1 ), y( i   , j+1 , k+1 ), z( i   , j+1 , k+1 ) /)
   
         ! edges vectors
         AC = vertC - vertA
         AF = vertF - vertA
         AH = vertH - vertA
   
         BE = vertE - vertB
   
         CH = vertH - vertC
   
         DB = vertB - vertD
   
         ED = vertD - vertE
   
         FC = vertC - vertF
   
         GA = vertA - vertG
         GB = vertB - vertG
         GC = vertC - vertG
         GD = vertD - vertG
         GE = vertE - vertG
         GF = vertF - vertG
         GH = vertH - vertG
   
         HF = vertF - vertH
   
         ! volume formula for a general hexahedronn obtained from
         ! Davies, D. E., & Salmond, D. J. (1985). Calculation of the volume of a general  
         ! hexahedron for flow predictions. AIAA journal, 23(6), 954-956.
   
         VolumeAux =   -one / twelve * &
                        (    Dot_Product( (GA+GB) , ( CrossProduct(DB,AC) ) ) &
                          +  Dot_Product( (GA+GE) , ( CrossProduct(BE,AF) ) ) &
                          +  Dot_Product( (GA+GD) , ( CrossProduct(ED,AH) ) ) &
                          +  Dot_Product( GF      , ( CrossProduct(HF,GE) ) ) &
                          +  Dot_Product( GH      , ( CrossProduct(CH,GD) ) ) &
                          +  Dot_Product( GC      , ( CrossProduct(FC,GB) ) ) &
                        )
   
         VolumeBulkCells = VolumeBulkCells + VolumeAux 
   
      end if
   
   end do
   end do
   end do
   
   ! All the rest of the fluid volume is contributed by the one contained within tetrahedra
   ! at the interface
   
   do Ktet = 1, ntetrahedra
      VolumeKTetrahedra = VolumeKTetrahedra + TetrahedraListAux(Ktet)%WaterVolume
   end do
   
   ! Deallocating this local version of these big arrays
   deallocate(NodesListAux, TetrahedraListAux)
   
   ! Global mass file writing
   
   inquire(file = "GlobalMassSussman.txt", exist = FileExist)
   
   if (FileExist) then
     open( 12, file = "GlobalMassSussman.txt", status = "old", &
               position = "append", action = "write"       )
   else
     open( 12, file = "GlobalMassSussman.txt", status = "new", & 
               action = "write")
   end if
   
   write(12, *) VolumeKTetrahedra + VolumeSecondTetrahedra + VolumeBulkCells
   close(12)

end subroutine TotalWaterVolumeSussman
