subroutine PiecewiseConstantFunction(ntetrahedra, KTetrahedraList)

   use DataTypes
   use precision
   use TetrahedronMethods
   use global_lsm, only: ConvergenceToleranceGeomReini

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer, intent(in) :: ntetrahedra
   type(tetrahedron), dimension(1:ntetrahedra), intent(inout) :: KTetrahedraList
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real(kind = rdf), dimension(4,3)    :: VerticesCoordinates
   real(kind = rdf), dimension(4)      :: phihast
   real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   real(kind = rdf), dimension(4)      :: VerticesDistances

   real(kind = rdf) :: DeltaVolume , DeltaVolumeIni
   real(kind = rdf) :: IsosurfaceArea, WaterVolumePhiH, WaterVolumeAux
   real(kind = rdf) :: VolumeError
   real(kind = rdf) :: etaKlocal, SDAux , etaKAux

   integer :: K
   integer :: ConvergenceCounter, MaxNumberOfIterations
   logical :: IsConverged

   DeltaVolume = ten
   VolumeError = ConvergenceToleranceGeomReini!/100.0_rdf
   MaxNumberOfIterations = 40

   do K = 1, ntetrahedra

      ! WaterVolumePhiH = V(ϕh) : Water Volume obtained with the original 
      ! linear reconstruction of ϕ, ϕh

      WaterVolumePhiH = KTetrahedraList(K)%WaterVolume

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

      ! phihast = ϕh* (see the paper, table I), which is the 
      ! reconstruction of ϕ using the signed distance function
      ! computed from nodal points  

      phihast(1) = KTetrahedraList(K)%v1%SDistanceFreeSurface
      phihast(2) = KTetrahedraList(K)%v2%SDistanceFreeSurface
      phihast(3) = KTetrahedraList(K)%v3%SDistanceFreeSurface
      phihast(4) = KTetrahedraList(K)%v4%SDistanceFreeSurface

      ! I perform MarchingTetrahedron to compute the water volume
      ! enclosed by ϕh* to compare the volume with the one enclosed
      ! by ϕh

      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 phihast                , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 IsosurfaceArea         , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )

      ! ΔV(ϕh, ϕh*) = VK(ϕh) - VK(ϕh*)
      DeltaVolumeIni = WaterVolumePhiH - WaterVolumeAux     
      DeltaVolume    = DeltaVolumeIni


      ! If V(ϕh) or V(ϕh*) < VolumeError, I skip the convergence loop
      ! by setting DeltaVolume to 0

      if (       WaterVolumePhiH < eps_sims  &
            .or. WaterVolumeAux  < eps_sims    ) then

         DeltaVolume = zero

      end if

      ! counter and correction term initialisation
      ConvergenceCounter = 0
      etaKlocal = zero
      IsConverged = .false.

      ! I look for a local correction for ϕh*, ηK, to minimise ΔV(ϕh, ϕh* + ηK) inside
      ! each tetrahedron
   
!      ConvergenceLoop: do while( abs(DeltaVolume)   >  VolumeError .and.      &
!                                 ConvergenceCounter <= MaxNumberOfIterations  )

      ConvergenceLoop: do while( (.not. IsConverged ) .and. ConvergenceCounter <= MaxNumberOfIterations )

         ! if Isosurface ≈ 0, ηK is undefined. If so, we drop the local mass 
         ! correction (setting ηK = 0) and we keep the geometric distance for 
         ! the reinitialisation

         if ( IsosurfaceArea  < eps_sims ) then
            etaKlocal = zero
            exit ConvergenceLoop
         end if

        ! Originally in the paper ηK = - DeltaVolume / IsosurfaceArea, but we checked 
        ! that it diverges with the - sign so we tried without it for the Zalesak's disk 
        ! and it converged in a few iterations

         etaKAux = DeltaVolume / IsosurfaceArea

         ! if etaKlocal produces phase changes, I exit the convergence loop

         !if( sign( one , KTetrahedraList(K)%v1%SDistanceFreeSurface ) * sign( one , phihast(1) + etaKAux ) < zero .or. & 
         !    sign( one , KTetrahedraList(K)%v2%SDistanceFreeSurface ) * sign( one , phihast(2) + etaKAux ) < zero .or. &
         !    sign( one , KTetrahedraList(K)%v3%SDistanceFreeSurface ) * sign( one , phihast(3) + etaKAux ) < zero .or. & 
         !    sign( one , KTetrahedraList(K)%v4%SDistanceFreeSurface ) * sign( one , phihast(4) + etaKAux ) < zero          ) then
         !
         !   exit ConvergenceLoop
         !   
         !end if

         ! correction of ϕh* --> ϕh* + ηK
         phihast(1) = phihast(1) + PhaseChangeLimiterCorrection ( phihast(1), etaKAux, PCR_PiecewiseFunction )
         phihast(2) = phihast(2) + PhaseChangeLimiterCorrection ( phihast(2), etaKAux, PCR_PiecewiseFunction )
         phihast(3) = phihast(3) + PhaseChangeLimiterCorrection ( phihast(3), etaKAux, PCR_PiecewiseFunction )
         phihast(4) = phihast(4) + PhaseChangeLimiterCorrection ( phihast(4), etaKAux, PCR_PiecewiseFunction )
   
         ! reconstruction of the free surface using (ϕh* + ηK) to compute the
         ! volume enclosed by the corrected ϕh* distribution over the 
         ! tetrahedron
         
         call MarchingTetrahedron(  VerticesCoordinates    , &
                                    phihast                , &
                                    VerticesIsosurfaces    , &
                                    ntriangles             , &
                                    IsosurfaceArea         , &
                                    WaterVolumeAux         , &
                                    VerticesDistances           )
   
         ! I check ΔV(ϕh, ϕh* + ηK) for the next while iteration
         DeltaVolume = WaterVolumePhiH - WaterVolumeAux

         if ( WaterVolumeAux     < VolumeError           ) DeltaVolume = zero

         if ( abs( DeltaVolume ) < abs( DeltaVolumeIni ) ) etaKlocal   = etaKAux
         if ( abs( DeltaVolume ) <      VolumeError      ) IsConverged = .true.

         ! I update the counter of the convergence loop
         ConvergenceCounter = ConvergenceCounter + 1

      end do ConvergenceLoop

      !print *, ' '

      ! if after all the iterations, ηK didn't help to reduce ΔV(ϕh, ϕh* + ηK), 
      ! I just set it to 0, which means I keep the signed distance estimation as 
      ! the correction for ϕ 
      
      if ( abs( DeltaVolume ) > abs( DeltaVolumeIni ) ) etaKlocal = zero

      ! Store the converged value of etaK for every tetrahedron
      KTetrahedraList(K)%etaK = etaKlocal

   end do ! do K = 1, ntetrahedra

end subroutine PiecewiseConstantFunction