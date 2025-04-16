subroutine NeighboursPhiCorrection( nnodes, ntetrahedra, NodesList, TetrahedraList )

   use DataTypes
   use precision
   use TetrahedronMethods
   use global_lsm, only: ConvergenceToleranceGeomReini

   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! Array with all the nodes to be used for reinitialisation.It's got the 
   ! target atributte because elements from linked lists points to elements  
   ! of the Tetrahedra List
   type(pnode), target , dimension(:), allocatable, intent(inout) :: NodesList
   
   ! Array with all the tetrahedra to be analised. It's got the target
   ! atributte because elements from linked lists points to elements of 
   ! the Tetrahedra List
!   type(tetrahedron), target, dimension(:), allocatable, intent(in) :: TetrahedraList
    type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: TetrahedraList
  
   integer, intent(inout) :: nnodes, ntetrahedra

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! This routine corresponds to the step 4 of the geometric
   ! redistancing algorithm

   real(kind = rdf) :: TotalTetrahedraVolume, TotalVolumeAux
   real(kind = rdf) :: DeltaVolume_i, DeltaVolume_i1, DeltaVolume_i2
   real(kind = rdf) :: CCorrected, C_i, C_i1, C_i2
   real(kind = rdf) :: m_i
   integer :: counter

   ! ***** NEWTON - RHAPSON *****

   real(kind = rdf) :: Der, epsNR

   ! dummy variables for Marching Tetrahedra

   real(kind = rdf), dimension(4,3)    :: VerticesCoordinates
   real(kind = rdf), dimension (4)     :: PhiAux
   real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   real(kind = rdf)                    :: IsosurfaceArea
   real(kind = rdf)                    :: WaterVolumeAux
   real(kind = rdf), dimension(4)      :: VerticesDistances
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! TotalTetrahedraVolume is the summation of all the water volume contained in the 
   ! triangulation USING THE ORIGINA PHI DISTRIBUTION
   TotalTetrahedraVolume  = zero 
   
   ! TotalTetrahedraVolume is the summation of all the water volume contained in the 
   ! triangulation using the MODIFIED PHI DISTRIBUTION
   TotalVolumeAux = zero

   ! Delta Volume (function to be solved: ΔV( ϕh, ϕh* + ηK) = 0 )
   ! δV^i
   DeltaVolume_i = one
   ! δV^{i-1}
   DeltaVolume_i1 = one/two ! zero !one  !one/three ! 1/3
   ! δV^{i-2}
   DeltaVolume_i2 = 0.0000001_rdf !one/two  ! 1/2

   ! Global mass correction factor
   C_i  = one
   C_i1 = one/two !zero !one/five !zero!one/two  !one/five ! 1/5
   C_i2 = 0.000001_rdf !zero !one/four !zero !one/four ! 1/4


   counter = 1

   ! ***** NEWTON - RHAPSON *****

   Der = one
   epsNR = 0.000001_rdf


   ! This is the total water volume for the original phi distribution
   ! over phase-change tetrahedra set

   do K = 1, ntetrahedra

      TotalTetrahedraVolume =    TotalTetrahedraVolume  &
                              +  TetrahedraList(K)%WaterVolume
   end do

   !print *, ' '
   !print *, ' ConvergenceToleranceGeomReini = ', ConvergenceToleranceGeomReini
   !print *, ' '

   ! Convergence iteration for determining C (CCorrected)
   do while ( abs(DeltaVolume_i) > 0.000001_rdf)! ConvergenceToleranceGeomReini )

      m_i = (C_i1 - C_i2) / (   DeltaVolume_i1 - DeltaVolume_i2  )

      C_i = C_i2 - m_i * DeltaVolume_i2

      !print *, 'counter = ', counter
      counter = counter + 1
      !print *, 'C_i before do Tetrahedra = ', C_i

!      C_i = C_i1 - m_i * DeltaVolume_i1

!      ! LET'S TRY NEWTON-RHAPSON

!      C_i = C_i1 - DeltaVolume_i1/Der

      ! Compute the whole volume with the iteratively corrected
      ! C value

      TotalVolumeAux = zero

      ! loop over all the tetrahedra set
      do K = 1, ntetrahedra

         VerticesCoordinates(1,:) = (/ TetrahedraList(K)%v1%x, &
                                       TetrahedraList(K)%v1%y, &
                                       TetrahedraList(K)%v1%z /)  
   
         VerticesCoordinates(2,:) = (/ TetrahedraList(K)%v2%x, &
                                       TetrahedraList(K)%v2%y, &
                                       TetrahedraList(K)%v2%z /)  
         
         VerticesCoordinates(3,:) = (/ TetrahedraList(K)%v3%x, &
                                       TetrahedraList(K)%v3%y, &
                                       TetrahedraList(K)%v3%z /)  
   
         VerticesCoordinates(4,:) = (/ TetrahedraList(K)%v4%x, &
                                       TetrahedraList(K)%v4%y, &
                                       TetrahedraList(K)%v4%z /)  
   
         ! PhiAux = phih* + C * ξh 
   
         PhiAux(1) =   TetrahedraList(K)%v1%SDistanceFreeSurface &
                     + C_i * TetrahedraList(K)%v1%xiI

!         ! Newton-Rhapson attempt
!         PhiAux(1) =   TetrahedraList(K)%v1%SDistanceFreeSurface &
!                     + (C_i + epsNR) * TetrahedraList(K)%v1%xiI

         PhiAux(2) =   TetrahedraList(K)%v2%SDistanceFreeSurface &
                     + C_i * TetrahedraList(K)%v2%xiI

!         PhiAux(2) =   TetrahedraList(K)%v2%SDistanceFreeSurface &
!                     + (C_i + epsNR) * TetrahedraList(K)%v2%xiI

         PhiAux(3) =   TetrahedraList(K)%v3%SDistanceFreeSurface &
                     + C_i * TetrahedraList(K)%v3%xiI

!         PhiAux(3) =   TetrahedraList(K)%v3%SDistanceFreeSurface &
!                     + (C_i + epsNR) * TetrahedraList(K)%v3%xiI

         PhiAux(4) =   TetrahedraList(K)%v4%SDistanceFreeSurface &
                     + C_i * TetrahedraList(K)%v4%xiI
   
!         PhiAux(4) =   TetrahedraList(K)%v4%SDistanceFreeSurface &
!                     + (C_i + epsNR) * TetrahedraList(K)%v4%xiI

         call MarchingTetrahedron(  VerticesCoordinates    , &
                                    PhiAux                 , &
                                    VerticesIsosurfaces    , &
                                    ntriangles             , &
                                    IsosurfaceArea         , &
                                    WaterVolumeAux         , &
                                    VerticesDistances           )
   
         ! We update the auxiliary total volume (the one computed
         ! using phiAux)
         TotalVolumeAux = TotalVolumeAux + WaterVolumeAux

         if (abs(PhiAux(1))> one) then
            !print *, 'i,j,k = ', TetrahedraList(K)%v1%i, &
            !                     TetrahedraList(K)%v1%j, &
            !                     TetrahedraList(K)%v1%k, &
            !                    ', PhiAux(1) = ', PhiAux(1)

         end if

         if (abs(PhiAux(2))> one) then
            !print *, 'i,j,k = ', TetrahedraList(K)%v2%i, &
            !                     TetrahedraList(K)%v2%j, &
            !                     TetrahedraList(K)%v2%k, &
            !                    ', PhiAux(2) = ', PhiAux(2)
         end if

         if (abs(PhiAux(3))> one) then
            !print *, 'i,j,k = ', TetrahedraList(K)%v3%i, &
            !                     TetrahedraList(K)%v3%j, &
            !                     TetrahedraList(K)%v3%k, &
            !                    ', PhiAux(3) = ', PhiAux(3)
         end if


         if (abs(PhiAux(4))> one) then
            !print *, 'i,j,k = ', TetrahedraList(K)%v4%i, &
            !                     TetrahedraList(K)%v4%j, &
            !                     TetrahedraList(K)%v4%k, &
            !                    ', PhiAux(4) = ', PhiAux(4)
         end if

      end do ! loop over all the tetrahedra set

      !Der = (DeltaVolume_i - DeltaVolume_i1)/epsNR

      DeltaVolume_i2 = DeltaVolume_i1
      DeltaVolume_i1 = DeltaVolume_i

      DeltaVolume_i = TotalTetrahedraVolume - TotalVolumeAux

      ! Iteration values update

      C_i2  = C_i1 ! C^{i-2} = C^{i-1}
      C_i1  = C_i  ! C^{i-1} = C^{i}

      !print *, 'TotalTetrahedraVolume = ', TotalTetrahedraVolume
      !print *, 'TotalVolumeAux = ', TotalVolumeAux
      !print *, 'DeltaVolume_i = ', DeltaVolume_i
      !print *, 'C_i = ', C_i 

   end do ! do while 

   CCorrected = C_i

   ! phi correction over pnodes

   ! Φ  = Φh* + φh
   ! φh = C*ξh

   ! Φh* : Signed Distance to the free surface (SDistanceFreeSurface)
   ! ξh  : orthogonal projection of etah (xiI)
   ! C   : constant that globally preserves volume (CCorrected)

   do I = 1 , nnodes

      ! I change the value of the phi array

      NodesList(I)%phi_value =    NodesList(I)%SDistanceFreeSurface &
                                 + CCorrected * NodesList(I)%xiI

      ! Now I change the phi values in phicorrected

      phicorrected( NodesList(I)%i , &
                    NodesList(I)%j , &
                    NodesList(I)%k     ) = NodesList(I)%phi_value

   end do

end subroutine NeighboursPhiCorrection