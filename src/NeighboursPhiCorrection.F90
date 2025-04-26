subroutine NeighboursPhiCorrection ( nnodes            , & 
                                     ntetrahedra       , &
                                     PNodesList        , & 
                                     KTetrahedraList   , &
                                     InterfaceNodesID  , &
                                     TotalWaterVolume       )

   use DataTypes
   use precision
   use TetrahedronMethods
   use NewtonRaphsonSolver

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! Array with all the nodes to be used for reinitialisation.It's got the 
   ! target atributte because elements from linked lists points to elements  
   ! of the Tetrahedra List
   type(pnode), target , dimension(:), allocatable, intent(inout) :: PNodesList
   
   ! Array with all the tetrahedra to be analised. It's got the target
   ! atributte because elements from linked lists points to elements of 
   ! the Tetrahedra List
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
   
   integer, intent(in) :: nnodes, ntetrahedra
   integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: InterfaceNodesID 
   real (kind = rdf) , intent(in) :: TotalWaterVolume

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! This routine corresponds to the step 4 of the geometric
   ! redistancing algorithm

   real(kind = rdf) :: CCorrected, C_InitialGuest, C_i, PhiCorrectionAux
   integer :: counter
   integer :: Inode
   integer :: i,j,k, isearch, jsearch, ksearch


   ! Surrounding nodes reinitialisation

   integer :: iminus, iplus, jminus, jplus, kminus, kplus
   integer :: iminusDist, iplusDist, jminusDist, jplusDist, kminusDist, kplusDist
   integer :: isearchDist, jsearchDist, ksearchDist
   integer :: igc, jgc, kgc
   integer :: Ilocal , NKI
   real (kind = rdf), dimension(3) :: point
   real (kind = rdf):: DistancePoint , DistanceAux, exsign, exsign2, PhiCorrectedAux
   real (kind = rdf), dimension(3,3) :: VerticesTriangle
   logical :: PnodeAround

   ! dummy variables for Marching Tetrahedra

   REAL(kind = rdf), dimension(4,3)    :: VerticesCoordinates
   REAL(kind = rdf), dimension (4)     :: PhiAux
   REAL(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   REAL(kind = rdf)                    :: IsosurfaceArea
   REAL(kind = rdf)                    :: WaterVolumeAux
   REAL(kind = rdf), dimension(4)      :: VerticesDistances
   real(kind = rdf) :: TotalVolumeAux

   real(kind = rdf) :: PhaseChangeFactor

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   integer :: ista , iend
   integer :: jsta , jend
   integer :: ksta , kend

   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp


   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 
   
   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 
   
   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! TotalTetrahedraVolume is the summation of all the water volume contained in the 
   ! triangulation USING THE ORIGINAL PHI DISTRIBUTION
   ! TotalTetrahedraVolume  = zero 
   
   ! TotalTetrahedraVolume is the summation of all the water volume contained in the 
   ! triangulation using the MODIFIED PHI DISTRIBUTION
   ! TotalVolumeAux = zero

   ! Initial guest for C_i before the iterative solver
   C_i = zero !0.1_rdf !0.01_rdf ! zero !0.1_rdf
   C_InitialGuest = C_i
   CCorrected = zero

   if ( GlobalMassCorrection ) then

      call NewtonRaphson ( ntetrahedra         , &  
                           KTetrahedraList     , & 
                           C_InitialGuest      , & 
                           C_i                 , & 
                           TotalWaterVolume      &
                          )
      
      print *, 'C_i AFTER Newton Raphson = ', C_i

      CCorrected = C_i

   end if

   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
   ! ϕ CORRECTION OVER PNODES
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          

   ! ϕ  = ϕh* + φh
   ! φh = C*ξh

   ! ϕh* : Signed Distance to the free surface (SDistanceFreeSurface)
   ! ξh  : orthogonal projection of etah (xiI)
   ! C   : constant that globally preserves volume (CCorrected)

   ! I change the value of the ϕ array
   do Inode = 1 , nnodes


      ! There is no phase change --> PhaseChangeFactor = 1
      ! There is    phase change --> PhaseChangeFactor = 0
      
      PhaseChangeFactor = abs( sign( one , PNodesList(Inode)%SDistanceFreeSurface ) + & 
                               sign( one , PNodesList(Inode)%SDistanceFreeSurface   + &
                                           CCorrected * PNodesList(Inode)%xiI     )      ) / two


      PhiCorrectionAux = PhaseChangeFactor * CCorrected * PNodesList(Inode)%xiI


      PNodesList(Inode)%phi_value = PNodesList(Inode)%SDistanceFreeSurface + PhiCorrectionAux

      ! Now I change the ϕ values in phicorrected

      i = PNodesList(Inode)%i
      j = PNodesList(Inode)%j
      k = PNodesList(Inode)%k

      ! Update the the ϕ value of interior nodes only
      if ( i >= i_mysta .and. i <= i_myend  .and. &
           j >= j_mysta .and. j <= j_myend  .and. & 
           k >= k_mysta .and. k <= k_myend          ) then

            phicorrected(i,j,k) = PNodesList(Inode)%phi_value
      
      end if

   end do

   ! igc = i - geometric correction

   do kgc = k_mysta , k_myend ! ksta , kend
   do jgc = j_mysta , j_myend ! jsta , jend
   do igc = i_mysta , i_myend ! ista , iend

      ! InterfaceNodesID(igc,jgc,kgc)  = -1 --> <= 3 nodes away from the pnodes
      ! InterfaceNodesID(igc,jgc,kgc)  =  0 --> >  3 nodes away from the pnodes
      ! InterfaceNodesID(igc,jgc,kgc) >=  1 -->    I'm at a pnode (I don't want to change these ones)

      if ( InterfaceNodesID(igc,jgc,kgc) >= 0 ) cycle

      !kminus = max( kgc - 3 , ksta )
      !kplus  = min( kgc + 3 , kend )

      !jminus = max( jgc - 3 , jsta )
      !jplus  = min( jgc + 3 , jend )

      !iminus = max( igc - 3 , ista )
      !iplus  = min( igc + 3 , iend )

      kminus = max( kgc - 4 , ksta )
      kplus  = min( kgc + 4 , kend )

      jminus = max( jgc - 4 , jsta )
      jplus  = min( jgc + 4 , jend )

      iminus = max( igc - 4 , ista )
      iplus  = min( igc + 4 , iend )

      PnodeAround = any( InterfaceNodesID( iminus:iplus , jminus:jplus , kminus:kplus ) > 0 )

      if (PnodeAround) then

         point(1) = x( igc , jgc , kgc )
         point(2) = y( igc , jgc , kgc )
         point(3) = z( igc , jgc , kgc )

         do ksearchDist = kminus , kplus
         do jsearchDist = jminus , jplus
         do isearchDist = iminus , iplus

            ! check if the neighbour is a pnode
            if ( InterfaceNodesID( isearchDist, jsearchDist, ksearchDist ) > 0 ) then

               exsign  = zero 
               exsign2 = zero
               DistancePoint = ten

               ! ID of the node I'm asking for its tetrahedrons
               Ilocal = InterfaceNodesID( isearchDist, jsearchDist, ksearchDist )

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

                        ! if DistanceAux <= DistancePoint --> exsign = 1 --> Replace DistancePoint by DistanceAux
                        ! if DistanceAux  > DistancePoint --> exsign = 0 --> Keep current DistancePoint value
                        exsign = ( sign( one , DistancePoint - DistanceAux ) + one )/two

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
   
                           ! if DistanceAux <= DistancePoint --> exsign = 1 --> Replace DistancePoint by DistanceAux
                           ! if DistanceAux  > DistancePoint --> exsign = 0 --> Keep current DistancePoint value
                           exsign = (sign(one , DistancePoint - DistanceAux) + one)/two
                        
                           ! if DistanceAux < DistancePoint, DistancePoint is replaced by DistanceAux
                           DistancePoint =  DistancePoint + exsign * ( DistanceAux - DistancePoint )

                        end if

                     end if

                     ! if phi( igc , jgc , kgc ) < DistancePoint --> exsign2 = 1. 
                     ! Otherwise, exsign2 = 0
                     exsign2 = ( sign( one , ( abs( phicorrected(igc, jgc, kgc) )  - DistancePoint ) ) + one )/two

                     ! We update the geometrical distance to the free-surface
                     PhiCorrectedAux =  abs( phicorrected(igc, jgc, kgc) ) +  &
                                             exsign2 * ( DistancePoint - abs( phicorrected(igc, jgc, kgc) ) ) 

                     ! we add the sign depending on the phase
                     phicorrected( igc , jgc , kgc ) =   sign(one, phicorrected(igc, jgc, kgc)) * &
                                                         PhiCorrectedAux

                  end do ! do K = 1, NKI

               else
                  print *, 'Pnode without associated Tetrahedra'                  
               end if ! if(PNodesList(Ilocal)%NumAsocTetrahedra > 0)

            end if ! if ( InterfaceNodesID( isearch, jsearch, ksearch ) > 0 ) then

         end do
         end do
         end do

      end if

   end do
   end do
   end do


end subroutine NeighboursPhiCorrection