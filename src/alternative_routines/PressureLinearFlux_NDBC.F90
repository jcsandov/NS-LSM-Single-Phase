function PressureLinearFlux_NDBC(  pL       , pC       , pR        , &
                                   JL       , JC       , JR        , &
                                   metricsL , metricsC , metricsR  , &
                                   dcsi                            , & 
                                   phiL     , phiC     , phiR      , & 
                                   GradPhiL , GradPhiC , GradPhiR  , &  
                                   GradVelL , GradVelC , GradVelR  , &
                                   PrintResults                      &
                                 )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! This function computes the flux related to the pressure in the linear flux
   ! term of Navier-Stokes equations for nodes next to free surface. It uses the 
   ! normal dynamic boundary condition to compute the pressure at the free surface.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   !                                                          these derivatives, may need to be adjusted  
   !                 _                   _     _        _     depending on the phase of the L,R sides.                                             
   !                | ∂/∂ξ (p/J * ∂ξ/∂x)  |   | ∂f(1)/∂ξ |    To compute them we use the "irregular stars"
   ! PressureFlux = | ∂/∂ξ (p/J * ∂ξ/∂y)  | = | ∂f(2)/∂ξ | -> algorithm described in "A Computer Study of                                
   !                | ∂/∂ξ (p/J * ∂ξ/∂z)  |   | ∂f(3)/∂ξ |    Finite-Amplitude Water Waves" by Chan &  
   !                •-                   -•   •-        -•    Street (JCP, 1970).                    
   !
   !
   !                    ϕ = 0 
   !                     /  
   !        air         /      water
   !       (ϕ<0)       /       (ϕ>0)
   !                  /
   !     ----o-------x--------o------------------o---  --> ξ
   !        i-1     /fs       i                 i+1
   !        (L)    /         (C)                (R)
   !         |    /           |                  |
   !         |   /            |                  | 
   !         |  /             |                  |
   !         | |              |                  |
   !         | |<--- ΔξL  --->|<-----  ΔξR ----->|
   !         
   !         
   !  If it is necessary to modify the stencil, the general expression for the derivative becomes
   !  ( Ferziger, Computational Methods for Fluid Dynamics, 4th edition, section 3.4 )  
   !
   !  ∂f/∂ξ = a * fL + b * fC + c * fR 
   !  
   !  Where 
   !  
   !  a = (       - ΔξR^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !  b = ( ΔξR^2 - ΔξL^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !  c = (         ΔξL^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !
   ! As fL or fR may lie on the free surface, the pressure needs to be computed using the normal 
   ! dynamic boundary condition (NDBC) in non-dimensional form
   ! 
   !
   !
   !         2      ∂ui            |    
   ! pfs  + ---- * ----- * ni * nj |    = 0 
   !         Re     ∂xj            |fs
   !
   !
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Input variables
   real ( kind = rdf ), intent(in) :: pL, pC, pR ! pressure
   real ( kind = rdf ), intent(in) :: JL, JC, JR ! Jacobian
   real ( kind = rdf ), dimension(3)   , intent(in) :: metricsL, metricsC, metricsR ! ∂ξ/∂xj
   real ( kind = rdf ), intent(in) :: dcsi ! Δξ
   real ( kind = rdf ), intent(in) :: phiL , phiC, phiR ! ϕ
   real ( kind = rdf ), dimension(3)   , intent(in) :: GradPhiL , GradPhiC , GradPhiR ! ∇ϕ
   real ( kind = rdf ), dimension(3,3) , intent(in) :: GradVelL , GradVelC , GradVelR ! ∇u
   logical, optional :: PrintResults

   ! Local variables

   integer :: i,j
   real ( kind = rdf ) :: dcsi_dx_L ,  dcsi_dy_L , dcsi_dz_L ! ∂ξ/∂xj_L 
   real ( kind = rdf ) :: dcsi_dx_C ,  dcsi_dy_C , dcsi_dz_C ! ∂ξ/∂xj_C
   real ( kind = rdf ) :: dcsi_dx_R ,  dcsi_dy_R , dcsi_dz_R ! ∂ξ/∂xj_R
   real ( kind = rdf ), dimension(3) :: fL, fC, fR ! p/J * ∂ξ/∂xj
   real ( kind = rdf ) :: InterpCoeff ! Interpolation coefficient (0 < InterpCoeff < 1)
   real ( kind = rdf ) :: dcsiL , dcsiR
   real ( kind = rdf ) :: tol

   ! Variables at the free surface

   real ( kind = rdf ) :: J_fs , dcsi_dx_fs , dcsi_dy_fs , dcsi_dz_fs ! Jacobian and metrics at the fs
   real ( kind = rdf ), dimension(3)   :: GradPhi_fs , nvec_fs ! ∇ϕ_fs & normal vector at the fs
   real ( kind = rdf ), dimension(3,3) :: GradVel_fs ! ∇u_fs
   real ( kind = rdf ) :: dui_dxj_ni_nj ! ∂ui/∂xj * ni * nj
   real ( kind = rdf ) :: pfs ! Pressure at the free surface obtained by NDBC
   real ( kind = rdf ) :: a , b, c ! Weights for the stencil of general derivative

   ! Output variable
   real ( kind = rdf ), dimension(3) :: PressureLinearFlux_NDBC ! ∂/∂ξ (p/J * ∂ξ/∂xj)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Set the tolerance

   tol = ten**(-eight)

   ! Metrics ∂ξ/∂x, ∂ξ/∂y & ∂ξ/∂z at L,C & R nodes

   dcsi_dx_L = metricsL(1)
   dcsi_dy_L = metricsL(2)
   dcsi_dz_L = metricsL(3)

   dcsi_dx_C = metricsC(1)
   dcsi_dy_C = metricsC(2)
   dcsi_dz_C = metricsC(3)

   dcsi_dx_R = metricsR(1)
   dcsi_dy_R = metricsR(2)
   dcsi_dz_R = metricsR(3)

   ! Initialise ΔξL and ΔξR as the mesh grid size
   dcsiL = dcsi
   dcsiR = dcsi

   ! Initialise fL, fc and fR without considering phase change
   
   ! p/J * ∂ξ/∂x
   fL(1) = pL / JL * dcsi_dx_L 
   fC(1) = pC / JC * dcsi_dx_C 
   fR(1) = pR / JR * dcsi_dx_R 

   ! p/J * ∂ξ/∂y
   fL(2) = pL / JL * dcsi_dy_L 
   fC(2) = pC / JC * dcsi_dy_C 
   fR(2) = pR / JR * dcsi_dy_R 

   ! p/J * ∂ξ/∂z
   fL(3) = pL / JL * dcsi_dz_L 
   fC(3) = pC / JC * dcsi_dz_C 
   fR(3) = pR / JR * dcsi_dz_R 

!   if( phiL < -eps_sims ) then
   if( phiL < zero ) then
   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  LEFT NODE WITHIN THE AIR PHASE                     
      !                       
      !        air       ϕ = 0    water
      !       (ϕ<0)       /       (ϕ>0)
      !                  /
      !     ----o-------x--------o------------------o---  --> ξ
      !        i-1     /fs       i                 i+1
      !        (L)    /         (C)                (R)
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
      ! Interpolation coefficient:
      ! Free surface approaches node C ==>  InterpCoeff --> 0
      ! Free surface approaches node L ==>  InterpCoeff --> 1

      InterpCoeff = abs( phiC / ( phiC - phiL + eps_sims ) )
      !InterpCoeff = abs( phiC / abs( phiC - phiL ) )

      ! I modify dcsiL interpolating Δξ to the free surface
      dcsiL = InterpCoeff * dcsi 
   
      ! I compute J, ∂ξ/∂xj, ∇ϕ (normal vector), and ∇u at the free surface as an weighted combination
      ! of the values at nodes L and C
   
      J_fs        = InterpCoeff * JL        + ( one - InterpCoeff ) * JC       
      
      dcsi_dx_fs  = InterpCoeff * dcsi_dx_L + ( one - InterpCoeff ) * dcsi_dx_C  
      dcsi_dy_fs  = InterpCoeff * dcsi_dy_L + ( one - InterpCoeff ) * dcsi_dy_C  
      dcsi_dz_fs  = InterpCoeff * dcsi_dz_L + ( one - InterpCoeff ) * dcsi_dz_C  
      
      GradPhi_fs  = InterpCoeff * GradPhiL  + ( one - InterpCoeff ) * GradPhiC 
      GradVel_fs  = InterpCoeff * GradVelL  + ( one - InterpCoeff ) * GradVelC 

      ! unitary normal vector at the free surface n = -∇ϕ / |∇ϕ|
      !nvec_fs = -GradPhi_fs / ( norm2( GradPhi_fs ) + eps_sims )
      nvec_fs = -GradPhi_fs / ( norm2( GradPhi_fs ) )

      ! Computation of ∂ui/∂xj * ni * nj

      dui_dxj_ni_nj = zero

      do j = 1,3
         do i = 1,3
            
            dui_dxj_ni_nj = dui_dxj_ni_nj + GradVel_fs(i,j) * nvec_fs(i) * nvec_fs(j) 

         end do
      end do

      ! Normal Dynamic Boundary Condition p_fs = -2/Re * ( ∂ui / ∂xj * ni * nj )_fs
      !pfs = -two / ren * dui_dxj_ni_nj
      pfs = zero

      ! Update of fL using the values at the free surface
      fL(1) = pfs / J_fs * dcsi_dx_fs 
      fL(2) = pfs / J_fs * dcsi_dy_fs 
      fL(3) = pfs / J_fs * dcsi_dz_fs 

   end if


!   if ( phiR < -eps_sims ) then
   if ( phiR < zero ) then

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  RIGHT NODE WITHIN THE AIR PHASE                     
      !                       
      !                water            ϕ = 0      air
      !                (ϕ>0)               \      (ϕ<0)
      !                                     \   
      !     ----o----------------o-----------x-------o---  --> ξ
      !        i-1               i         fs \     i+1
      !        (L)              (C)            \    (R)
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Interpolation coefficient:
      ! Free surface approaches node C ==>  InterpCoeff --> 0
      ! Free surface approaches node R ==>  InterpCoeff --> 1

      InterpCoeff = abs( phiC / ( phiC - phiR + eps_sims ) )
      !InterpCoeff = abs( phiC / ( phiC - phiR ) )

      ! I modify dcsiR interpolating Δξ to the free surface
      dcsiR = InterpCoeff * dcsi 
   
      ! I compute J, ∂ξ/∂xj, ∇ϕ (normal vector), and ∇u at the free surface as an weighted combination
      ! of the values at nodes L and C
   
      J_fs        = InterpCoeff * JR        + ( one - InterpCoeff ) * JC

      dcsi_dx_fs  = InterpCoeff * dcsi_dx_R + ( one - InterpCoeff ) * dcsi_dx_C  
      dcsi_dy_fs  = InterpCoeff * dcsi_dy_R + ( one - InterpCoeff ) * dcsi_dy_C  
      dcsi_dz_fs  = InterpCoeff * dcsi_dz_R + ( one - InterpCoeff ) * dcsi_dz_C  
      
      GradPhi_fs  = InterpCoeff * GradPhiR  + ( one - InterpCoeff ) * GradPhiC 
      GradVel_fs  = InterpCoeff * GradVelR  + ( one - InterpCoeff ) * GradVelC 

      ! unitary normal vector at the free surface n = -∇ϕ / |∇ϕ|
      nvec_fs = -GradPhi_fs / norm2( GradPhi_fs )

      ! Computation of ∂ui/∂xj * ni * nj

      dui_dxj_ni_nj = zero

      do j = 1,3
         do i = 1,3
            
            dui_dxj_ni_nj = dui_dxj_ni_nj + GradVel_fs(i,j) * nvec_fs(i) * nvec_fs(j) 

         end do
      end do

      ! Normal Dynamic Boundary Condition p_fs = -2/Re * ( ∂ui / ∂xj * ni * nj )_fs
      !pfs = -two / ren * dui_dxj_ni_nj
      pfs = zero

      ! Update of fR using the values at the free surface
      fR(1) = pfs / J_fs * dcsi_dx_fs 
      fR(2) = pfs / J_fs * dcsi_dy_fs 
      fR(3) = pfs / J_fs * dcsi_dz_fs 
   
   end if

   ! Stencil coefficients for ∂f/∂ξ = a * fL + b * fC + c * fR (Ferziger)
   
   ! I think there's a typo error with this coefficients, as a should go with
   ! dcsiL instead of dcsiR
   !a =   -dcsiR**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
   !b = (  dcsiR**2 - dcsiL**2 ) / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
   !c =    dcsiL**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  


   if ( tol * dcsiR > dcsiL ) then

      PressureLinearFlux_NDBC = (fR - fC) / dcsi

   else if ( tol * dcsiL > dcsiR ) then

      PressureLinearFlux_NDBC = (fC - fL) / dcsi

   else

      a =   -dcsiR**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
      b = (  dcsiR**2 - dcsiL**2 ) / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
      c =    dcsiL**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  

      ! Finally, I compute PressureLinearFlux_NDBC = ∂f/∂ξ 
      PressureLinearFlux_NDBC = (a * fL) + (b * fC) + (c * fR)

   end if


   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   if ( present ( PrintResults ) ) then

      if ( PrintResults ) then

         print *, ' '
         
         write(*,'(A,F15.12)') 'dcsiL = ', dcsiL
         write(*,'(A,F15.12)') 'dcsiR = ', dcsiR

         if ( phiR < zero ) write(*,'(A,F15.12)') 'R InterpCoeff = ', InterpCoeff
         if ( phiL < zero ) write(*,'(A,F15.12)') 'L InterpCoeff = ', InterpCoeff

         print *, ' '

         write(*,'(A,F15.12)') 'a = ', a
         write(*,'(A,F15.12)') 'b = ', b
         write(*,'(A,F15.12)') 'c = ', c

         print *, ' '

         write(*,'(A,F15.12)') 'fL(3) = ', fL(3)
         write(*,'(A,F15.12)') 'fC(3) = ', fC(3)
         write(*,'(A,F15.12)') 'fR(3) = ', fR(3)

         print *, ' '

         write(*,'(A,F15.12)') 'PressureLinearFlux_NDBC(3) = ', PressureLinearFlux_NDBC(3)

      end if

   end if

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


end function PressureLinearFlux_NDBC