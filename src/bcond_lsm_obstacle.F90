
subroutine bcond_lsm_obstacle (phiobstacle)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! Actualizacion de la presion y velocidades 
  ! (non-slip condition) en el borde de obstaculo
  ! 
  ! 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use global_app
  use global_mpi
  use global_obstacle

  implicit none
  
  ! solution vector
  ! 
  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: phiobstacle
  real (kind = rdf) :: coef_1, coef_2

  integer :: i, j, k, l_o
  integer :: i1,i2,i3,i4, j1,j2,j3,j4
  coef_1 = four_third
  coef_2 = -one_third	
	
  if(is_obstacle == 1) then
    do l_o = 1,8 
      if(obstacle_zone(l_o) == 1) then
        do k = k_la_o(l_o), k_lb_o(l_o)
          do j = j_la_o(l_o), j_lb_o(l_o)
            do i = i_la_o(l_o), i_lb_o(l_o)
              i1 = i + i_idx_o(1,l_o)
              i2 = i + i_idx_o(2,l_o)
              i3 = i + i_idx_o(3,l_o)
              i4 = i + i_idx_o(4,l_o)
              j1 = j + j_idx_o(1,l_o)
              j2 = j + j_idx_o(2,l_o)
              j3 = j + j_idx_o(3,l_o)
              j4 = j + j_idx_o(4,l_o)
              
              phiobstacle(i,j,k) = one/two * ( coef_1 * phiobstacle(i1,j1,k) + coef_2 * phiobstacle(i2,j2,k) + &
                                       coef_1 * phiobstacle(i3,j3,k) + coef_2 * phiobstacle(i4,j4,k)    )
                                      
            end do
          end do
        end do		
      end if
    end do
  end if
  
  

end subroutine bcond_lsm_obstacle

