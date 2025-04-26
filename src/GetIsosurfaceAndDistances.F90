subroutine GetIsosurfaceAndDistances( nnodes          , &
                                      ntetrahedra     , &
                                      PNodesList      , &
                                      KTetrahedraList , &
                                      InterfaceNodesID     )

   ! This subroutines get the isosurface reconstruction on every tetrahedron

   use DataTypes
   use precision
   use TetrahedronMethods
   use AdvectionMethods

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer, intent(in) :: nnodes, ntetrahedra
   type(pnode), dimension(1:nnodes), intent(inout) :: PNodesList
   type(tetrahedron), dimension(1:ntetrahedra), intent(inout) :: KTetrahedraList
   integer, dimension(il:iu,jl:ju,kl:ku), intent(inout) :: InterfaceNodesID 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   real(kind = rdf), dimension (4,3)   :: VerticesCoordinates
   real(kind = rdf), dimension (4)     :: phi_tetrahedron
   real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   real(kind = rdf)                    :: IsosurfaceAreaAux
   real(kind = rdf)                    :: WaterVolumeAux
   real(kind = rdf), dimension(4)      :: VerticesDistances

   ! Auxiliary variables for distance computation
   real(kind = rdf) :: exsign, exsign2, DistanceAux, DistancePoint
   real(kind = rdf), dimension(3) :: point
   real(kind = rdf), dimension(3,3) :: VerticesTriangle

   ! Auxiliary counters for distance computation
   integer :: iidx, jidx, kidx 
   integer :: isearch, jsearch, ksearch 
   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   integer :: K, I, Ilocal, NKI

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   integer :: ista, iend, jsta, jend, ksta, kend

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
            
   ! Marching Tetrahedron Step
   do K = 1, ntetrahedra

      VerticesCoordinates(1,:) = (/ KTetrahedraList(K)%v1%x, &
                                    KTetrahedraList(K)%v1%y, &
                                    KTetrahedraList(K)%v1%z /)  

      VerticesCoordinates(2,:) = (/ KTetrahedraList(K)%v2%x, &
                                    KTetrahedraList(K)%v2%y, &
                                    KTetrahedraList(K)%v2%z /)  
      
      VerticesCoordinates(3,:) = (/ KTetrahedraList(K)%v3%x, &
                                    KTetrahedraList(K)%v3%y, &
                                    KTetrahedraList(K)%v3%z /)  

      VerticesCoordinates(4,:) = (/ KTetrahedraList(K)%v4%x, &
                                    KTetrahedraList(K)%v4%y, &
                                    KTetrahedraList(K)%v4%z /)  

      phi_tetrahedron(1) = KTetrahedraList(K)%v1%phi_value
      phi_tetrahedron(2) = KTetrahedraList(K)%v2%phi_value
      phi_tetrahedron(3) = KTetrahedraList(K)%v3%phi_value
      phi_tetrahedron(4) = KTetrahedraList(K)%v4%phi_value

      ! TO DO: MarchingTetrahedron could receive pointers to nodes
      ! instead of vertices coordinates and phi_tetrahedron values
      ! (just doing the previous steps internally)

      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 phi_tetrahedron        , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 IsosurfaceAreaAux      , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )

      KTetrahedraList(K)%IsosurfaceArea = IsosurfaceAreaAux
      KTetrahedraList(K)%WaterVolume    = WaterVolumeAux

      ! if there's an isosurface inside the tetrahedron
      if (ntriangles > 0) then

         ! vertex 1 of the first isosurface triangle
         KTetrahedraList(K)%isosurface1%v1%x = VerticesIsosurfaces(1,1,1)
         KTetrahedraList(K)%isosurface1%v1%y = VerticesIsosurfaces(1,1,2)
         KTetrahedraList(K)%isosurface1%v1%z = VerticesIsosurfaces(1,1,3)
   
         ! vertex 2 of the first isosurface triangle
         KTetrahedraList(K)%isosurface1%v2%x = VerticesIsosurfaces(1,2,1)
         KTetrahedraList(K)%isosurface1%v2%y = VerticesIsosurfaces(1,2,2)
         KTetrahedraList(K)%isosurface1%v2%z = VerticesIsosurfaces(1,2,3)
   
         ! vertex 3 of the first isosurface triangle
         KTetrahedraList(K)%isosurface1%v3%x = VerticesIsosurfaces(1,3,1)
         KTetrahedraList(K)%isosurface1%v3%y = VerticesIsosurfaces(1,3,2)
         KTetrahedraList(K)%isosurface1%v3%z = VerticesIsosurfaces(1,3,3)
   
         KTetrahedraList(K)%nisosurfaces = ntriangles

      end if

      ! if there's a second isosurface inside the tetrahedron
      if (ntriangles > 1) then

         ! vertex 1 of the second isosurface triangle
         KTetrahedraList(K)%isosurface2%v1%x = VerticesIsosurfaces(2,1,1)
         KTetrahedraList(K)%isosurface2%v1%y = VerticesIsosurfaces(2,1,2)
         KTetrahedraList(K)%isosurface2%v1%z = VerticesIsosurfaces(2,1,3)
   
         ! vertex 2 of the second isosurface triangle
         KTetrahedraList(K)%isosurface2%v2%x = VerticesIsosurfaces(2,2,1)
         KTetrahedraList(K)%isosurface2%v2%y = VerticesIsosurfaces(2,2,2)
         KTetrahedraList(K)%isosurface2%v2%z = VerticesIsosurfaces(2,2,3)
   
         ! vertex 3 of the second isosurface triangle
         KTetrahedraList(K)%isosurface2%v3%x = VerticesIsosurfaces(2,3,1)
         KTetrahedraList(K)%isosurface2%v3%y = VerticesIsosurfaces(2,3,2)
         KTetrahedraList(K)%isosurface2%v3%z = VerticesIsosurfaces(2,3,3)

      end if

   end do   

   ! Computing the distance to the vertices of the tetrahedra based on
   ! an extended search radius (i,j,k ± 3 search radius)

   do I = 1, nnodes

      ! local node location
      iidx = PNodesList(I)%i
      jidx = PNodesList(I)%j
      kidx = PNodesList(I)%k

      ! geometrical location of the point of interest
      point(1) = PNodesList(I)%x
      point(2) = PNodesList(I)%y
      point(3) = PNodesList(I)%z

      ! initialisation of the distance to the point as a huge number 
      DistancePoint = ten

      ! search radius to compute the distance from the point
      ! to the free-surface. This loop evaluate the distance to the
      ! adjacent reconstructed isosurfaces represented by triangles inside
      ! tetrahedrons

      !kminus = maxval( (/ kidx-3 , ksta /) )
      !kplus  = minval( (/ kidx+3 , kend /) )

      !jminus = maxval( (/ jidx-3 , jsta /) )
      !jplus  = minval( (/ jidx+3 , jend /) )

      !iminus = maxval( (/ iidx-3 , ista /) )
      !iplus  = minval( (/ iidx+3 , iend /) )

      kminus = max( kidx-4 , ksta )
      kplus  = min( kidx+4 , kend )

      jminus = max( jidx-4 , jsta )
      jplus  = min( jidx+4 , jend )

      iminus = max( iidx-4 , ista )
      iplus  = min( iidx+4 , iend )

      do ksearch = kminus , kplus
      do jsearch = jminus , jplus
      do isearch = iminus , iplus

         exsign  = zero 
         exsign2 = zero
         
         ! check if the neighbour is a pnode (including the node itself) 
         if ( InterfaceNodesID( isearch, jsearch, ksearch ) > 0 ) then

            ! ID of the node I'm asking for its tetrahedrons
            Ilocal = InterfaceNodesID( isearch, jsearch, ksearch )

            ! check how many phase-changing tetrahedra are associated 
            ! to the node
            if ( PNodesList(Ilocal)%NumAsocTetrahedra > 0 ) then
                     
               NKI = PNodesList(Ilocal)%NumAsocTetrahedra
                        
               do K = 1, NKI ! loop over the simplices with a vertex at Ilocal
                     
                  ! check how many isosurfaces are associated to every associated tetrahedron
                  if (       KTetrahedraList( PNodesList(Ilocal)%AssociatedTetrahedraID(K) )%nisosurfaces > 0 &
                       .and. KTetrahedraList( PNodesList(Ilocal)%AssociatedTetrahedraID(K) )%IsosurfaceArea > zero ) then
                           
                     ! vertex 1 of the isosurface 1                                 
                     VerticesTriangle(1,1) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v1%x
                     VerticesTriangle(1,2) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v1%y
                     VerticesTriangle(1,3) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v1%z

                     ! vertex 2 of the isosurface 1                                 
                     VerticesTriangle(2,1) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v2%x
                     VerticesTriangle(2,2) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v2%y
                     VerticesTriangle(2,3) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v2%z

                     ! vertex 3 of the isosurface 1                                 
                     VerticesTriangle(3,1) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v3%x
                     VerticesTriangle(3,2) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v3%y
                     VerticesTriangle(3,3) = KTetrahedraList( &
                                             PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface1%v3%z

                     DistanceAux = DistancePointTriangle( VerticesTriangle, point )

                     ! if DistanceAux < DistancePoint --> exsign = 1. Otherwise, exsign = 0
                     exsign = ( sign( one , DistancePoint - DistanceAux ) + one ) / two

                     ! if DistanceAux < DistancePoint, DistancePoint is replaced by DistanceAux
                     DistancePoint =  DistancePoint + exsign * ( DistanceAux - DistancePoint )

                     if ( KTetrahedraList( PNodesList(Ilocal)%AssociatedTetrahedraID(K) )%nisosurfaces &
                        > 1) then

                        ! vertex 1 of the isosurface 2                                 
                        VerticesTriangle(1,1) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v1%x
                        VerticesTriangle(1,2) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v1%y
                        VerticesTriangle(1,3) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v1%z
   
                        ! vertex 2 of the isosurface 2                                 
                        VerticesTriangle(2,1) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v2%x
                        VerticesTriangle(2,2) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v2%y
                        VerticesTriangle(2,3) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v2%z
   
                        ! vertex 3 of the isosurface 2                                 
                        VerticesTriangle(3,1) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v3%x
                        VerticesTriangle(3,2) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v3%y
                        VerticesTriangle(3,3) = KTetrahedraList( &
                                                PNodesList(Ilocal)%AssociatedTetrahedraID(K))%isosurface2%v3%z

                        DistanceAux = DistancePointTriangle( VerticesTriangle, point )
   
                        ! if DistanceAux < DistancePoint --> exsign = 1. Otherwise, exsign = 0
                        exsign = (sign(one , DistancePoint - DistanceAux) + one)/two
                        
                        ! if DistanceAux < DistancePoint, DistancePoint is replaced by DistanceAux
                        DistancePoint =  DistancePoint + exsign * ( DistanceAux - DistancePoint )

                     end if

                  end if

                  ! if abs(PNodesList(I)%SDistanceFreeSurface) < DistancePoint --> exsign2 = 1. 
                  ! Otherwise, exsign2 = 0
                  exsign2 = (sign(one , abs(PNodesList(I)%SDistanceFreeSurface) - DistancePoint) + one)/two

                  ! We update the geometrical distance to the free-surface
                  PNodesList(I)%SDistanceFreeSurface =    abs( PNodesList(I)%SDistanceFreeSurface ) &
                                                      +  exsign2 * (   DistancePoint &
                                                                     - abs(PNodesList(I)%SDistanceFreeSurface) ) 

                  ! we add the sign depending on the phase
                  PNodesList(I)%SDistanceFreeSurface =   sign(one, PNodesList(I)%phi_value) &
                                                      * PNodesList(I)%SDistanceFreeSurface

               end do ! do K = 1, NKI

            else
               print *, 'Pnode with no associated Tetrahedra'                  
            end if ! if(PNodesList(Ilocal)%NumAsocTetrahedra > 0)

         else

            ! I set these flags to compute the geometric distance to the nodes surrounding the pnodes
            ! all the nodes arround the pnodes to -1, except those which are pnodes too
         
            InterfaceNodesID( isearch , jsearch , ksearch ) = -1 
         
         end if ! if ( InterfaceNodesID( isearch, jsearch, ksearch ) > 0 ) then

      end do ! do isearch = iminus , iplus
      end do ! do jsearch = jminus , jplus
      end do ! do ksearch = kminus , kplus

   end do ! do I = 1, nnodes



end subroutine GetIsosurfaceAndDistances
