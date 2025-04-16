subroutine PiecewiseConstantFunction ( ntetrahedra , KTetrahedraList )

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
   real(kind = rdf), dimension(4)      :: PhiHast
   real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   real(kind = rdf), dimension(4)      :: VerticesDistances

   real(kind = rdf) :: DeltaVolume , DeltaVolumeIni
   real(kind = rdf) :: AreaIsosurfacePhiH, WaterVolumePhiH
   real(kind = rdf) :: VolumeErrorTolerance
   real(kind = rdf) :: etaKlocal, SDAux , etaKAux

   integer :: K
   integer :: ConvergenceCounter, MaxNumberOfIterations, ConvergenceRate, ConvergenceAverage
   logical :: IsConverged

   DeltaVolume = ten
   VolumeErrorTolerance = ConvergenceToleranceGeomReini/1000.0_rdf
   MaxNumberOfIterations = 40

   ConvergenceRate    = 0
   ConvergenceAverage = 0

   do K = 1, ntetrahedra

      ! WaterVolumePhiH = V(ϕh) : Water Volume obtained with the original 
      ! linear reconstruction of ϕ, ϕh just after the advection step

      WaterVolumePhiH    = KTetrahedraList(K)%WaterVolume
      AreaIsosurfacePhiH = KTetrahedraList(K)%IsosurfaceArea


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

      ! PhiHast = ϕh* (see the paper, table I), which is the 
      ! reconstruction of the zero level-set of ϕ using the 
      ! signed distance function computed from nodal points  

      PhiHast(1) = KTetrahedraList(K)%v1%SDistanceFreeSurface
      PhiHast(2) = KTetrahedraList(K)%v2%SDistanceFreeSurface
      PhiHast(3) = KTetrahedraList(K)%v3%SDistanceFreeSurface
      PhiHast(4) = KTetrahedraList(K)%v4%SDistanceFreeSurface

      etaKlocal = zero

      call SecantMethodLocalCorrection ( VerticesCoordinates      , &
                                         WaterVolumePhiH          , &
                                         AreaIsosurfacePhiH       , &
                                         PhiHast                  , &
                                         VolumeErrorTolerance     , &
                                         MaxNumberOfIterations    , &
                                         etaKAux                  , &
                                         IsConverged              , &
                                         ConvergenceAverage         &    
                                       )

      if ( IsConverged ) then
         
         etaKlocal          = etaKAux

         ConvergenceRate    = ConvergenceRate    + 1
         ConvergenceAverage = ConvergenceAverage + 1
      
      end if

      ! Store the converged value of etaK for every tetrahedron
      KTetrahedraList(K)%etaK = etaKlocal

   end do ! do K = 1, ntetrahedra


end subroutine PiecewiseConstantFunction