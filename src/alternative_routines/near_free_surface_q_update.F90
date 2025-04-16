
subroutine near_free_surface_q_update( PressureInterpolation, VelocityInterpolation )

   use InterpolationMethods

   implicit none

   logical, intent(in) :: PressureInterpolation, VelocityInterpolation

   integer, parameter :: nvars = 15

   real ( kind = rdf ), dimension(3) :: nvec
   real ( kind = rdf ), dimension(3) :: rijk , FreeSurfacePosition

   real ( kind = rdf ), dimension(8,3) :: CellVerticesCoordinates 
   real ( kind = rdf ), dimension(6)   :: DistanceToFaces


   real ( kind = rdf ), dimension(8)   :: rsignVertices    , &
                                          u_CellArray      , &
                                          v_CellArray      , &
                                          w_CellArray      , &
                                          dudx_CellArray   , &
                                          dudy_CellArray   , &
                                          dudz_CellArray   , &
                                          dvdx_CellArray   , &
                                          dvdy_CellArray   , &
                                          dvdz_CellArray   , &
                                          dwdx_CellArray   , &
                                          dwdy_CellArray   , &
                                          dwdz_CellArray   , &
                                          dphidx_CellArray , &
                                          dphidy_CellArray , &
                                          dphidz_CellArray 

   real ( kind = rdf ), dimension(3,3) :: VelocityGradientV1 , &
                                          VelocityGradientV2 , &
                                          VelocityGradientV3 , &
                                          VelocityGradientV4 , &
                                          VelocityGradientV5 , &
                                          VelocityGradientV6 , &
                                          VelocityGradientV7 , &
                                          VelocityGradientV8 

   real ( kind = rdf ), dimension(3) ::   PhiGradientV1 , &
                                          PhiGradientV2 , &
                                          PhiGradientV3 , &
                                          PhiGradientV4 , &
                                          PhiGradientV5 , &
                                          PhiGradientV6 , &
                                          PhiGradientV7 , &
                                          PhiGradientV8 


   real ( kind = rdf ), dimension(nvars,8) :: VarsToInterpolate_CellArray 
   real ( kind = rdf ), dimension(nvars)   :: VarsInterpolatedfs 

   real ( kind = rdf ) :: pfs , ufs , vfs , wfs , dui_dxj_ni_nj_fs 
   real ( kind = rdf ), dimension(3,3) :: VelocityGradient_fs 
   real ( kind = rdf ), dimension(3)   :: nvec_fs 
 
   real ( kind = rdf ), dimension(3) :: r_ip1, r_jp1, r_kp1 

   integer :: OctantID 
   integer, dimension (3,4,3) :: FacesOffset

   real ( kind = rdf ), dimension(4,3) :: VerticesCoordinates
   real ( kind = rdf ), dimension(4,4) :: VarsToInterpolate_CuadrilateralArray
   real ( kind = rdf ), dimension(4)   :: VarsInterpolatedProjectedFace
   real ( kind = rdf ), dimension(3)   :: FaceLineIntersectionPoint
   real ( kind = rdf ) :: ppf , upf , vpf , wpf  

   real ( kind = rdf ) :: TotalDistance , IntpCoeff1 , IntpCoeff2

   logical :: AirNodeAround , PointWithinCell, FaceIntersectionFlag
   integer :: i,j,k
   integer :: iv1,jv1,kv1
   integer :: ivertex , jvertex , kvertex
   integer :: isum , jsum

   integer :: VertexLoop

   ! loop over the whole local domain looking for nodes to be interpolated

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , i_myend
            
      ! ----------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
      ! ----------------------------------------------------------------------------
      ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
      ! phase (rsign = 1), then this if is true. This condition identifies air nodes
      ! next to the free-surface
      ! ----------------------------------------------------------------------------
   

      ! The first if is to reduce as much as I can the area of the domain where I'm gonna
      ! the adjacency to the free-surface. This is because the search over the 
      ! rsign(i-1:i+1,j-1:j+1,k-1:k+1) neighbourhood is computationally costly.

      if ( rsign(i,j,k) > one_half .and. phi(i,j,k) < ten  ) then

         ! I check the i±1,j±1,k±1 directions. If any of them lie within the air phase,
         ! then AirNodeAround < 1/2 == .true. and the interpolation is performed

         AirNodeAround =   rsign( i-1 , j   , k   ) * rsign( i+1 , j   , k   ) &
                         * rsign( i   , j-1 , k   ) * rsign( i   , j+1 , k   ) &
                         * rsign( i   , j   , k-1 ) * rsign( i   , j   , k+1 )  < one_half

         if ( AirNodeAround ) then

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 1. Normal vector computation at i,j,k: nvec = -∇ϕ/|∇ϕ|
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! phi gradient
            !
            !  ∂ϕ     ∂ϕ     ∂ξ^k    
            ! ---- = ---- * ------ 
            ! ∂x_j   ∂ξ^k    ∂x_j  
            !

            ! With this call, I´m assigning nvec = ∇ϕ
            call PhiGradientVector(i,j,k,nvec)

            ! Reorient and normalise: nvec = -∇ϕ / |∇ϕ| (nvec points towards the air-phase
            ! now)
            nvec = -nvec / norm2(nvec)

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 2. Location computation and cell identification of the free-surface 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! position of the target node
            rijk = (/x(i,j,k) , y(i,j,k) , z(i,j,k)/)

            ! free surface position: rnode + ϕ * n = rfs
            FreeSurfacePosition = rijk + phi(i,j,k) * nvec

            ! Identification of the cell where the free surface lies on

            PointWithinCell = .false.

            ! I look for within which octant, the free surface is to perform the 
            ! three-linear interpolation of the velocity and pressure at the fs.

            searchloop : do iv1 = i-1 , i 
                         do jv1 = j-1 , j 
                         do kv1 = k-1 , k 

               ! 
               !     i,j+1,k+1 ----i+1,j+1,k+1    
               !      /|(8)           /|(7)                            
               !     / |             / |                                  
               !  i,j,k+1-------i+1,j,k+1                   
               !    |(5)           |(6)|                                   
               !    |  |           |   |                                
               !    |  |           |   |                                
               !    |  i,j+1,k-----|-i+1,j+1,k           
               !    | /(4)         |  /(3)                 
               !    |/             | /
               !  i,j,k-----------i+1,j,k
               !   (1)               (2)
               ! 

               ! Phi values at those vertices
               rsignVertices(1) = rsign( iv1   , jv1   , kv1   )
               rsignVertices(2) = rsign( iv1+1 , jv1   , kv1   )
               rsignVertices(3) = rsign( iv1+1 , jv1+1 , kv1   )
               rsignVertices(4) = rsign( iv1   , jv1+1 , kv1   )
               rsignVertices(5) = rsign( iv1   , jv1   , kv1+1 )
               rsignVertices(6) = rsign( iv1+1 , jv1   , kv1+1 )
               rsignVertices(7) = rsign( iv1+1 , jv1+1 , kv1+1 )
               rsignVertices(8) = rsign( iv1   , jv1+1 , kv1+1 )

               ! single phase cell
               if ( all( rsignVertices > one_half ) .or. &
                    all( rsignVertices < one_half ) ) then

                  cycle

               ! changing phase cell
               else

                  ! Coordinates of the eight vertices of the cell
                  CellVerticesCoordinates(1,:) = (/ x( iv1   , jv1   , kv1   ) , &
                                                    y( iv1   , jv1   , kv1   ) , &
                                                    z( iv1   , jv1   , kv1   ) /)

                  CellVerticesCoordinates(2,:) = (/ x( iv1+1 , jv1   , kv1   ) , &
                                                    y( iv1+1 , jv1   , kv1   ) , &
                                                    z( iv1+1 , jv1   , kv1   ) /)

                  CellVerticesCoordinates(3,:) = (/ x( iv1+1 , jv1+1 , kv1   ) , &
                                                    y( iv1+1 , jv1+1 , kv1   ) , &
                                                    z( iv1+1 , jv1+1 , kv1   ) /)

                  CellVerticesCoordinates(4,:) = (/ x( iv1   , jv1+1 , kv1   ) , &
                                                    y( iv1   , jv1+1 , kv1   ) , &
                                                    z( iv1   , jv1+1 , kv1   ) /)
   
                  CellVerticesCoordinates(5,:) = (/ x( iv1   , jv1   , kv1+1 ) , &
                                                    y( iv1   , jv1   , kv1+1 ) , &
                                                    z( iv1   , jv1   , kv1+1 ) /)

                  CellVerticesCoordinates(6,:) = (/ x( iv1+1 , jv1   , kv1+1 ) , &
                                                    y( iv1+1 , jv1   , kv1+1 ) , &
                                                    z( iv1+1 , jv1   , kv1+1 ) /)

                  CellVerticesCoordinates(7,:) = (/ x( iv1+1 , jv1+1 , kv1+1 ) , &
                                                    y( iv1+1 , jv1+1 , kv1+1 ) , &
                                                    z( iv1+1 , jv1+1 , kv1+1 ) /)

                  CellVerticesCoordinates(8,:) = (/ x( iv1   , jv1+1 , kv1+1 ) , &
                                                    y( iv1   , jv1+1 , kv1+1 ) , &
                                                    z( iv1   , jv1+1 , kv1+1 ) /)

                  ! verification if the fs lies on the current cell

                  call PointWithinCellCheck ( CellVerticesCoordinates , &
                                              FreeSurfacePosition     , &
                                              PointWithinCell         , &
                                              DistanceToFaces           &
                                            )

                  if ( PointWithinCell ) exit searchloop

               end if

            end do
            end do
            end do searchloop

            if ( .not. PointWithinCell ) then

               !print *, ' '
               !write(*, '(A, A, I0, A, I0, A, I0, A)') " Error at " , "i = ", i, " , j = ", j, ", k = ", k
               !write(*, '(A)') "Free-surface position not found"
               !print *, ' '

               stop

            end if

            ! I came out this searchloop with iv1, jv1 and kv1 which tells me the v1 node
            ! of the cell where the free surface is located. 

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 3. Velocities, velocity gradient and phi gradient at the free surface by 
            !    trilinear interpolation
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! velocities at the eight vertices constructed from iv1, jv1, kv1

            u_CellArray(1) = q( 2 , iv1   , jv1   , kv1   )
            u_CellArray(2) = q( 2 , iv1+1 , jv1   , kv1   )
            u_CellArray(3) = q( 2 , iv1+1 , jv1+1 , kv1   )
            u_CellArray(4) = q( 2 , iv1   , jv1+1 , kv1   )
            u_CellArray(5) = q( 2 , iv1   , jv1   , kv1+1 )
            u_CellArray(6) = q( 2 , iv1+1 , jv1   , kv1+1 )
            u_CellArray(7) = q( 2 , iv1+1 , jv1+1 , kv1+1 )
            u_CellArray(8) = q( 2 , iv1   , jv1+1 , kv1+1 )

            v_CellArray(1) = q( 3 , iv1   , jv1   , kv1   )
            v_CellArray(2) = q( 3 , iv1+1 , jv1   , kv1   )
            v_CellArray(3) = q( 3 , iv1+1 , jv1+1 , kv1   )
            v_CellArray(4) = q( 3 , iv1   , jv1+1 , kv1   )
            v_CellArray(5) = q( 3 , iv1   , jv1   , kv1+1 )
            v_CellArray(6) = q( 3 , iv1+1 , jv1   , kv1+1 )
            v_CellArray(7) = q( 3 , iv1+1 , jv1+1 , kv1+1 )
            v_CellArray(8) = q( 3 , iv1   , jv1+1 , kv1+1 )

            w_CellArray(1) = q( 4 , iv1   , jv1   , kv1   )
            w_CellArray(2) = q( 4 , iv1+1 , jv1   , kv1   )
            w_CellArray(3) = q( 4 , iv1+1 , jv1+1 , kv1   )
            w_CellArray(4) = q( 4 , iv1   , jv1+1 , kv1   )
            w_CellArray(5) = q( 4 , iv1   , jv1   , kv1+1 )
            w_CellArray(6) = q( 4 , iv1+1 , jv1   , kv1+1 )
            w_CellArray(7) = q( 4 , iv1+1 , jv1+1 , kv1+1 )
            w_CellArray(8) = q( 4 , iv1   , jv1+1 , kv1+1 )

            ! velocity gradient the eight vertices constructed from iv0, jv0, kv0

            call VelocityGradientTensor( iv1   , jv1   , kv1   , VelocityGradientV1 )
            call VelocityGradientTensor( iv1+1 , jv1   , kv1   , VelocityGradientV2 )
            call VelocityGradientTensor( iv1+1 , jv1+1 , kv1   , VelocityGradientV3 )
            call VelocityGradientTensor( iv1   , jv1+1 , kv1   , VelocityGradientV4 )
            call VelocityGradientTensor( iv1   , jv1   , kv1+1 , VelocityGradientV5 )
            call VelocityGradientTensor( iv1+1 , jv1   , kv1+1 , VelocityGradientV6 )
            call VelocityGradientTensor( iv1+1 , jv1+1 , kv1+1 , VelocityGradientV7 )
            call VelocityGradientTensor( iv1   , jv1+1 , kv1+1 , VelocityGradientV8 )

            !------------------------------------------
            ! ∂u/∂xj
            !------------------------------------------

            dudx_CellArray(1) = VelocityGradientV1(1,1)
            dudx_CellArray(2) = VelocityGradientV2(1,1)
            dudx_CellArray(3) = VelocityGradientV3(1,1)
            dudx_CellArray(4) = VelocityGradientV4(1,1)
            dudx_CellArray(5) = VelocityGradientV5(1,1)
            dudx_CellArray(6) = VelocityGradientV6(1,1)
            dudx_CellArray(7) = VelocityGradientV7(1,1)
            dudx_CellArray(8) = VelocityGradientV8(1,1)

            dudy_CellArray(1) = VelocityGradientV1(1,2)
            dudy_CellArray(2) = VelocityGradientV2(1,2)
            dudy_CellArray(3) = VelocityGradientV3(1,2)
            dudy_CellArray(4) = VelocityGradientV4(1,2)
            dudy_CellArray(5) = VelocityGradientV5(1,2)
            dudy_CellArray(6) = VelocityGradientV6(1,2)
            dudy_CellArray(7) = VelocityGradientV7(1,2)
            dudy_CellArray(8) = VelocityGradientV8(1,2)

            dudz_CellArray(1) = VelocityGradientV1(1,3)
            dudz_CellArray(2) = VelocityGradientV2(1,3)
            dudz_CellArray(3) = VelocityGradientV3(1,3)
            dudz_CellArray(4) = VelocityGradientV4(1,3)
            dudz_CellArray(5) = VelocityGradientV5(1,3)
            dudz_CellArray(6) = VelocityGradientV6(1,3)
            dudz_CellArray(7) = VelocityGradientV7(1,3)
            dudz_CellArray(8) = VelocityGradientV8(1,3)

            !------------------------------------------
            ! ∂v/∂xj
            !------------------------------------------

            dvdx_CellArray(1) = VelocityGradientV1(2,1)
            dvdx_CellArray(2) = VelocityGradientV2(2,1)
            dvdx_CellArray(3) = VelocityGradientV3(2,1)
            dvdx_CellArray(4) = VelocityGradientV4(2,1)
            dvdx_CellArray(5) = VelocityGradientV5(2,1)
            dvdx_CellArray(6) = VelocityGradientV6(2,1)
            dvdx_CellArray(7) = VelocityGradientV7(2,1)
            dvdx_CellArray(8) = VelocityGradientV8(2,1)
            
            dvdy_CellArray(1) = VelocityGradientV1(2,2)
            dvdy_CellArray(2) = VelocityGradientV2(2,2)
            dvdy_CellArray(3) = VelocityGradientV3(2,2)
            dvdy_CellArray(4) = VelocityGradientV4(2,2)
            dvdy_CellArray(5) = VelocityGradientV5(2,2)
            dvdy_CellArray(6) = VelocityGradientV6(2,2)
            dvdy_CellArray(7) = VelocityGradientV7(2,2)
            dvdy_CellArray(8) = VelocityGradientV8(2,2)
            
            dvdz_CellArray(1) = VelocityGradientV1(2,3)
            dvdz_CellArray(2) = VelocityGradientV2(2,3)
            dvdz_CellArray(3) = VelocityGradientV3(2,3)
            dvdz_CellArray(4) = VelocityGradientV4(2,3)
            dvdz_CellArray(5) = VelocityGradientV5(2,3)
            dvdz_CellArray(6) = VelocityGradientV6(2,3)
            dvdz_CellArray(7) = VelocityGradientV7(2,3)
            dvdz_CellArray(8) = VelocityGradientV8(2,3)

            !------------------------------------------
            ! ∂w/∂xj
            !------------------------------------------

            dwdx_CellArray(1) = VelocityGradientV1(3,1)
            dwdx_CellArray(2) = VelocityGradientV2(3,1)
            dwdx_CellArray(3) = VelocityGradientV3(3,1)
            dwdx_CellArray(4) = VelocityGradientV4(3,1)
            dwdx_CellArray(5) = VelocityGradientV5(3,1)
            dwdx_CellArray(6) = VelocityGradientV6(3,1)
            dwdx_CellArray(7) = VelocityGradientV7(3,1)
            dwdx_CellArray(8) = VelocityGradientV8(3,1)
            
            dwdy_CellArray(1) = VelocityGradientV1(3,2)
            dwdy_CellArray(2) = VelocityGradientV2(3,2)
            dwdy_CellArray(3) = VelocityGradientV3(3,2)
            dwdy_CellArray(4) = VelocityGradientV4(3,2)
            dwdy_CellArray(5) = VelocityGradientV5(3,2)
            dwdy_CellArray(6) = VelocityGradientV6(3,2)
            dwdy_CellArray(7) = VelocityGradientV7(3,2)
            dwdy_CellArray(8) = VelocityGradientV8(3,2)
            
            dwdz_CellArray(1) = VelocityGradientV1(3,3)
            dwdz_CellArray(2) = VelocityGradientV2(3,3)
            dwdz_CellArray(3) = VelocityGradientV3(3,3)
            dwdz_CellArray(4) = VelocityGradientV4(3,3)
            dwdz_CellArray(5) = VelocityGradientV5(3,3)
            dwdz_CellArray(6) = VelocityGradientV6(3,3)
            dwdz_CellArray(7) = VelocityGradientV7(3,3)
            dwdz_CellArray(8) = VelocityGradientV8(3,3)


            ! normal vector at the eight vertices constructed from iv1, jv1, kv1

            call PhiGradientVector( iv1   , jv1   , kv1   , PhiGradientV1 )
            call PhiGradientVector( iv1+1 , jv1   , kv1   , PhiGradientV2 )
            call PhiGradientVector( iv1+1 , jv1+1 , kv1   , PhiGradientV3 )
            call PhiGradientVector( iv1   , jv1+1 , kv1   , PhiGradientV4 )
            call PhiGradientVector( iv1   , jv1   , kv1+1 , PhiGradientV5 )
            call PhiGradientVector( iv1+1 , jv1   , kv1+1 , PhiGradientV6 )
            call PhiGradientVector( iv1+1 , jv1+1 , kv1+1 , PhiGradientV7 )
            call PhiGradientVector( iv1   , jv1+1 , kv1+1 , PhiGradientV8 )

            !------------------------------------------
            ! ∂ϕ/∂x
            !------------------------------------------

            dphidx_CellArray(1) = PhiGradientV1(1)
            dphidx_CellArray(2) = PhiGradientV2(1)
            dphidx_CellArray(3) = PhiGradientV3(1)
            dphidx_CellArray(4) = PhiGradientV4(1)
            dphidx_CellArray(5) = PhiGradientV5(1)
            dphidx_CellArray(6) = PhiGradientV6(1)
            dphidx_CellArray(7) = PhiGradientV7(1)
            dphidx_CellArray(8) = PhiGradientV8(1)

            !------------------------------------------
            ! ∂ϕ/∂y
            !------------------------------------------

            dphidy_CellArray(1) = PhiGradientV1(2)
            dphidy_CellArray(2) = PhiGradientV2(2)
            dphidy_CellArray(3) = PhiGradientV3(2)
            dphidy_CellArray(4) = PhiGradientV4(2)
            dphidy_CellArray(5) = PhiGradientV5(2)
            dphidy_CellArray(6) = PhiGradientV6(2)
            dphidy_CellArray(7) = PhiGradientV7(2)
            dphidy_CellArray(8) = PhiGradientV8(2)

            !------------------------------------------
            ! ∂ϕ/∂z
            !------------------------------------------

            dphidz_CellArray(1) = PhiGradientV1(3)
            dphidz_CellArray(2) = PhiGradientV2(3)
            dphidz_CellArray(3) = PhiGradientV3(3)
            dphidz_CellArray(4) = PhiGradientV4(3)
            dphidz_CellArray(5) = PhiGradientV5(3)
            dphidz_CellArray(6) = PhiGradientV6(3)
            dphidz_CellArray(7) = PhiGradientV7(3)
            dphidz_CellArray(8) = PhiGradientV8(3)

            ! Trilinear interpolation

            ! nvars = 3 velocities + 9 velocity gradients + 3 ϕ gradients = 15 variables to
            ! interpolate at the free surface

            VarsToInterpolate_CellArray( 1  , 1:8 ) =      u_CellArray 
            VarsToInterpolate_CellArray( 2  , 1:8 ) =      v_CellArray 
            VarsToInterpolate_CellArray( 3  , 1:8 ) =      w_CellArray 
            VarsToInterpolate_CellArray( 4  , 1:8 ) =   dudx_CellArray 
            VarsToInterpolate_CellArray( 5  , 1:8 ) =   dudy_CellArray 
            VarsToInterpolate_CellArray( 6  , 1:8 ) =   dudz_CellArray 
            VarsToInterpolate_CellArray( 7  , 1:8 ) =   dvdx_CellArray 
            VarsToInterpolate_CellArray( 8  , 1:8 ) =   dvdy_CellArray 
            VarsToInterpolate_CellArray( 9  , 1:8 ) =   dvdz_CellArray 
            VarsToInterpolate_CellArray( 10 , 1:8 ) =   dwdx_CellArray 
            VarsToInterpolate_CellArray( 11 , 1:8 ) =   dwdy_CellArray 
            VarsToInterpolate_CellArray( 12 , 1:8 ) =   dwdz_CellArray 
            VarsToInterpolate_CellArray( 13 , 1:8 ) = dphidx_CellArray 
            VarsToInterpolate_CellArray( 14 , 1:8 ) = dphidy_CellArray 
            VarsToInterpolate_CellArray( 15 , 1:8 ) = dphidz_CellArray 
            
            VarsInterpolatedfs = zero

            call TrilinearInterpolation( CellVerticesCoordinates     , &
                                         FreeSurfacePosition         , &
                                         DistanceToFaces             , &
                                         nvars                       , &
                                         VarsToInterpolate_CellArray , &
                                         VarsInterpolatedfs            &
                                       )

            ! Velocities at the free surface
            ufs = VarsInterpolatedfs(1)
            vfs = VarsInterpolatedfs(2)
            wfs = VarsInterpolatedfs(3)

            !print *, ' '
            !print *, ' ufs, vfs, wfs = ', ufs, vfs, wfs

            ! Velocity gradient at the free surface for the NDBC

            ! ∂u/∂xj
            VelocityGradient_fs(1,1) = VarsInterpolatedfs(4)
            VelocityGradient_fs(1,2) = VarsInterpolatedfs(5)
            VelocityGradient_fs(1,3) = VarsInterpolatedfs(6)            
            
            ! ∂v/∂xj
            VelocityGradient_fs(2,1) = VarsInterpolatedfs(7)
            VelocityGradient_fs(2,2) = VarsInterpolatedfs(8)
            VelocityGradient_fs(2,3) = VarsInterpolatedfs(9)            

            ! ∂w/∂xj
            VelocityGradient_fs(3,1) = VarsInterpolatedfs(10)
            VelocityGradient_fs(3,2) = VarsInterpolatedfs(11)
            VelocityGradient_fs(3,3) = VarsInterpolatedfs(12)            

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 4.  Normal Dynamic Boundary Condition at the free surface
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !
            !
            !         /  2      ∂ui            \    
            ! pfs  + (  ---- * ----- * ni * nj  )    = 0 
            !         \  Re     ∂xj            / fs
            !
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            nvec_fs =  VarsInterpolatedfs(13:15) ! ∂ϕ/∂xj
            nvec_fs = -nvec_fs / norm2(nvec_fs) ! unitary normal vector

            dui_dxj_ni_nj_fs = zero

            do jsum = 1,3
            do isum = 1,3
                  
               dui_dxj_ni_nj_fs =   dui_dxj_ni_nj_fs &
                                  + VelocityGradient_fs(isum,jsum) * nvec_fs(isum) * nvec_fs(jsum) 
      
            end do
            end do

            !pfs = - ( two / ren ) * dui_dxj_ni_nj_fs

            pfs = zero

            !print *, 'pfs = ', pfs

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 5. Projection to the adjacent face
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! I build the local vector base as ei = ri+1 - rijk (base vector i direction)
            !                                  ej = rj+1 - rijk (base vector j direction)
            !                                  ek = rk+1 - rijk (base vector k direction)

            ! Local vector base
            r_ip1 = (/ x(i+1,j,k) , y(i+1,j,k) , z(i+1,j,k) /)
            r_jp1 = (/ x(i,j+1,k) , y(i,j+1,k) , z(i,j+1,k) /)
            r_kp1 = (/ x(i,j,k+1) , y(i,j,k+1) , z(i,j,k+1) /)

            OctantID    = OctantIdentification( rijk , r_ip1 , r_jp1 , r_kp1 , -nvec )
            FacesOffset = GetOctantFacesIndexes (OctantID)

            FaceIntersectionFlag      = .false.
            FaceLineIntersectionPoint = zero
            
            ! Face 1

            ! I build up the array with the position of the vertices of the quadrilateral
            do VertexLoop = 1,4

               ivertex = i + FacesOffset(1,VertexLoop,1)
               jvertex = j + FacesOffset(1,VertexLoop,2)
               kvertex = k + FacesOffset(1,VertexLoop,3) 

               VerticesCoordinates(VertexLoop,1) = x( ivertex , jvertex , kvertex )  
               VerticesCoordinates(VertexLoop,2) = y( ivertex , jvertex , kvertex )
               VerticesCoordinates(VertexLoop,3) = z( ivertex , jvertex , kvertex )

            end do 

            ! Check if the ray from rijk intersects the plane defined by
            ! VerticesCoordinates

            call LineQuadrilateralIntersectionTest (  VerticesCoordinates        , &
                                                      rijk                       , &
                                                     -nvec                       , &
                                                      FaceIntersectionFlag       , &
                                                      FaceLineIntersectionPoint    &  
                                                   ) 

            if ( FaceIntersectionFlag ) then

               !print *, ' '
               !print *, 'FaceIntersectionFlag 1st face '
               !write(*, '(A, 3(F9.5, A))') "FaceLineIntersectionPoint = ", FaceLineIntersectionPoint(1), " ", FaceLineIntersectionPoint(2), " ", FaceLineIntersectionPoint(3)


               do VertexLoop = 1,4
   
                  ivertex = i + FacesOffset(1,VertexLoop,1)
                  jvertex = j + FacesOffset(1,VertexLoop,2)
                  kvertex = k + FacesOffset(1,VertexLoop,3) 
   
                  VarsToInterpolate_CuadrilateralArray( 1 , VertexLoop ) = q( 1 , ivertex , jvertex , kvertex )
                  VarsToInterpolate_CuadrilateralArray( 2 , VertexLoop ) = q( 2 , ivertex , jvertex , kvertex )
                  VarsToInterpolate_CuadrilateralArray( 3 , VertexLoop ) = q( 3 , ivertex , jvertex , kvertex ) 
                  VarsToInterpolate_CuadrilateralArray( 4 , VertexLoop ) = q( 4 , ivertex , jvertex , kvertex )
   
               end do 

            end if


            ! Face 2
            if ( .not. FaceIntersectionFlag ) then

               do VertexLoop = 1,4
   
                  ivertex = i + FacesOffset(2,VertexLoop,1)
                  jvertex = j + FacesOffset(2,VertexLoop,2)
                  kvertex = k + FacesOffset(2,VertexLoop,3) 
   
                  VerticesCoordinates(VertexLoop,1) = x( ivertex , jvertex , kvertex )  
                  VerticesCoordinates(VertexLoop,2) = y( ivertex , jvertex , kvertex )
                  VerticesCoordinates(VertexLoop,3) = z( ivertex , jvertex , kvertex )
   
               end do 

               call LineQuadrilateralIntersectionTest (  VerticesCoordinates        , &
                                                         rijk                       , &
                                                        -nvec                       , &
                                                         FaceIntersectionFlag       , &
                                                         FaceLineIntersectionPoint    &  
                                                      ) 

               if ( FaceIntersectionFlag ) then
   
                  do VertexLoop = 1,4
      
                     !print *, ' '
                     !print *, 'FaceIntersectionFlag 2nd face '
                     !write(*, '(A, 3(F9.5, A))') "FaceLineIntersectionPoint = ", FaceLineIntersectionPoint(1), " ", FaceLineIntersectionPoint(2), " ", FaceLineIntersectionPoint(3)

                     ivertex = i + FacesOffset(2,VertexLoop,1)
                     jvertex = j + FacesOffset(2,VertexLoop,2)
                     kvertex = k + FacesOffset(2,VertexLoop,3) 
      
                     VarsToInterpolate_CuadrilateralArray( 1 , VertexLoop ) = q( 1 , ivertex , jvertex , kvertex )
                     VarsToInterpolate_CuadrilateralArray( 2 , VertexLoop ) = q( 2 , ivertex , jvertex , kvertex )
                     VarsToInterpolate_CuadrilateralArray( 3 , VertexLoop ) = q( 3 , ivertex , jvertex , kvertex ) 
                     VarsToInterpolate_CuadrilateralArray( 4 , VertexLoop ) = q( 4 , ivertex , jvertex , kvertex )
      
                  end do 
   
               end if

            end if

            ! Face 3

            if ( .not. FaceIntersectionFlag ) then

               do VertexLoop = 1,4
   
                  ivertex = i + FacesOffset(3,VertexLoop,1)
                  jvertex = j + FacesOffset(3,VertexLoop,2)
                  kvertex = k + FacesOffset(3,VertexLoop,3) 
   
                  VerticesCoordinates(VertexLoop,1) = x( ivertex , jvertex , kvertex )  
                  VerticesCoordinates(VertexLoop,2) = y( ivertex , jvertex , kvertex )
                  VerticesCoordinates(VertexLoop,3) = z( ivertex , jvertex , kvertex )
   
               end do 

               call LineQuadrilateralIntersectionTest (  VerticesCoordinates        , &
                                                         rijk                       , &
                                                        -nvec                       , &
                                                         FaceIntersectionFlag       , &
                                                         FaceLineIntersectionPoint    &  
                                                      ) 

               if ( FaceIntersectionFlag ) then
   
                     !print *, ' '
                     !print *, 'FaceIntersectionFlag 3rd face '
                     !write(*, '(A, 3(F9.5, A))') "FaceLineIntersectionPoint = ", FaceLineIntersectionPoint(1), " ", FaceLineIntersectionPoint(2), " ", FaceLineIntersectionPoint(3)

                  do VertexLoop = 1,4
      
                     ivertex = i + FacesOffset(3,VertexLoop,1)
                     jvertex = j + FacesOffset(3,VertexLoop,2)
                     kvertex = k + FacesOffset(3,VertexLoop,3) 
      
                     VarsToInterpolate_CuadrilateralArray( 1 , VertexLoop ) = q( 1 , ivertex , jvertex , kvertex )
                     VarsToInterpolate_CuadrilateralArray( 2 , VertexLoop ) = q( 2 , ivertex , jvertex , kvertex )
                     VarsToInterpolate_CuadrilateralArray( 3 , VertexLoop ) = q( 3 , ivertex , jvertex , kvertex ) 
                     VarsToInterpolate_CuadrilateralArray( 4 , VertexLoop ) = q( 4 , ivertex , jvertex , kvertex )
      
                  end do 
   

               else
   
                  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
                  !print *, ' '
                  !print *, 'Error in the projection to the adjacent face at node'
                  !write(*, '(A, I0, A, I0, A, I0, A)') "i = ", i, " , j = ", j, ", k = ", k
                  !print *, ' '
                  !write(*, '(A, 3(F9.5, A))') "rijk: ", rijk(1), " ", rijk(2), " ", rijk(3)
                  !write(*, '(A, 3(F9.5, A))') "-nvec: ", -nvec(1), " ", -nvec(2), " ", -nvec(3)
                  !write(*, '(A, I0)') "OctantID = ", OctantID
                  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


                  stop
   
               end if

            end if

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 6. Bi-cuadratic interpolation
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! NOTE: after VarsInterpolatedProjectedFace, we can add an optional boolean 
            ! input variable to print the results of the Bilinear Interpolation

            call BilinearInterpolation( VerticesCoordinates                  , & 
                                        FaceLineIntersectionPoint            , &
                                        4                                    , &
                                        VarsToInterpolate_CuadrilateralArray , &
                                        VarsInterpolatedProjectedFace          &
                                       )

            ! pf : projection face

            ppf = VarsInterpolatedProjectedFace(1)
            upf = VarsInterpolatedProjectedFace(2)
            vpf = VarsInterpolatedProjectedFace(3)
            wpf = VarsInterpolatedProjectedFace(4)

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! 7. Flow field update next to the free surface
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ! Total Distance from the free surface to the adjacent face along the normal
            ! vector line

            TotalDistance =   norm2( rijk - FreeSurfacePosition       )  &
                            + norm2( rijk - FaceLineIntersectionPoint ) 

            IntpCoeff1 = one - norm2( rijk - FreeSurfacePosition       ) / TotalDistance
            IntpCoeff2 = one - norm2( rijk - FaceLineIntersectionPoint ) / TotalDistance

            ! Final update

            if ( PressureInterpolation ) then
            
               q(1,i,j,k) = IntpCoeff1 * pfs + IntpCoeff2 * ppf

            end if

            if ( VelocityInterpolation ) then
            
               q(2,i,j,k) = IntpCoeff1 * ufs + IntpCoeff2 * upf
               q(3,i,j,k) = IntpCoeff1 * vfs + IntpCoeff2 * vpf
               q(4,i,j,k) = IntpCoeff1 * wfs + IntpCoeff2 * wpf

            end if

            !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
            
            !print *, ' '
            !write(*, '(A, I0, A, I0, A, I0, A)') "i = ", i, " , j = ", j, ", k = ", k
            !write(*, '(A, 3(F15.12, A))') "rijk position: ", rijk(1), " ", rijk(2), " ", rijk(3)
            !write(*, '(A, 3(F15.12, A))') "Free Surface Position: ", FreeSurfacePosition(1), " ", FreeSurfacePosition(2), " ", FreeSurfacePosition(3)
            !write(*, '(A, 3(F15.12, A))') "Face Line Intersection Point: ", FaceLineIntersectionPoint(1), " ", FaceLineIntersectionPoint(2), " ", FaceLineIntersectionPoint(3)
            !write(*, '(A, (F15.12, A))') "InterpCoeff1: ", IntpCoeff1
            !write(*, '(A, (F15.12, A))') "InterpCoeff2: ", IntpCoeff2
            !print *, ' '
            !print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '

            !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

         end if

      end if

   end do
   end do
   end do

 
end subroutine near_free_surface_q_update


