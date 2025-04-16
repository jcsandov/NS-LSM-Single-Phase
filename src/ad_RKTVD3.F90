subroutine ad_RKTVD3( phiAD_n1, phiAD_n , AdvectionNodes )
!---------------------------------------------------------
!calcula el siguiente valor de de la superficie libre       
!phi en n+1 (tiempo siguiente)                              
!---------------------------------------------------------

   implicit none

   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: rh_AD
   real (kind = rdf),dimension(3,3) :: coef_RKTVD_AD
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(out) :: phiAD_n1
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: phiAD_n
   integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: AdvectionNodes 
   character(len = 256) :: debugname

   !-------------------------------------------------------------------------------------

   !Nodos incluyendo el borde
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! RK coefficients (obtained from: Cox, C., Liang, C., & Plesniak, M. W. (2016). 
   ! A high-order solver for unsteady incompressible Navier–Stokes equations using 
   ! the flux reconstruction method on unstructured grids with implicit dual time 
   ! stepping. Journal of Computational Physics, 314, 414-435.
   ! https://doi.org/10.1016/j.jcp.2016.03.016 )
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! RK stage 1
   coef_RKTVD_AD(1,1) = one  ! 1.0_rdf ! one 
   coef_RKTVD_AD(2,1) = zero ! 0.0_rdf ! zero
   coef_RKTVD_AD(3,1) = one  ! 1.0_rdf ! one 
  
   ! RK stage 2
   coef_RKTVD_AD(1,2) =  three / four ! 3.0_rdf/4.0_rdf ! three / four
   coef_RKTVD_AD(2,2) =  one   / four ! 1.0_rdf/4.0_rdf ! one   / four
   coef_RKTVD_AD(3,2) =  one   / four ! 1.0_rdf/4.0_rdf ! one   / four

   ! RK stage 3
   coef_RKTVD_AD(1,3) = one    / three ! 1.0_rdf/3.0_rdf ! one / three
   coef_RKTVD_AD(2,3) = two    / three ! 2.0_rdf/3.0_rdf ! two / three
   coef_RKTVD_AD(3,3) = two    / three ! 2.0_rdf/3.0_rdf ! two / three

   ! Initialising ϕ^n+1 = ϕ^n
   phiAD_n1 = phiAD_n
   
   do RK = 1,3
    
      rh_AD = zero

      if (orderAD == 0) call calc_RH_AD      ( phiAD_n1 , rh_AD , AdvectionNodes )    ! Calculo el RHS de la eqn de LS usando WENO 3
      if (orderAD == 1) call calc_RH_AD_WENO3( phiAD_n1 , rh_AD , AdvectionNodes )    ! Calculo el RHS de la eqn de LS usando WENO 3

      !if(is_obstacle == 1) call bcond_lsm_obstacle(phiAD_n1) ! extrapolacion upwind de phi al obstaculo
    
      do k = k_mysta , k_myend
      do j = j_mysta , j_myend
      do i = i_mysta , i_myend
          
         if ( AdvectionNodes(i,j,k) > 0 ) then
            
         !if(act_obstacle_ad(i,j,k) == 1) then

            phiAD_n1(i,j,k) = coef_RKTVD_AD(1,RK) * phiAD_n(i,j,k)                        + &
                              coef_RKTVD_AD(2,RK) * phiAD_n1(i,j,k)                       - &
                              coef_RKTVD_AD(3,RK) * aj(i,j,k) * delti_lsm * rh_AD(i,j,k)


         !end if
            
         end if

      end do
      end do
      end do
  
      call bcond_lsm( phiAD_n1 )    
      call rhs_exchng3_3d (phiAD_n1)

   end do ! RK

end subroutine ad_RKTVD3
