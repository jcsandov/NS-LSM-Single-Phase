
subroutine bcond_obstacle (il, iu, jl, ju, kl, ku, q)

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
  
  ! extents
  ! 
  integer, intent(in) :: il
  integer, intent(in) :: jl
  integer, intent(in) :: kl
  integer, intent(in) :: iu
  integer, intent(in) :: ju
  integer, intent(in) :: ku


  ! solution vector
  ! 
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: q
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
              
              q(1,i,j,k) = one/two * ( coef_1 * q(1,i1,j1,k) + coef_2 * q(1,i2,j2,k) + &
                                       coef_1 * q(1,i3,j3,k) + coef_2 * q(1,i4,j4,k)    )
                                      
              if(l_o == 8 .OR. l_o == 4) then !Paredes laterales direccion x
                q(2,i,j,k) = 0.0
                q(3,i,j,k) = 0.0!*one/two * ( coef_1 * q(3,i1,j1,k) + coef_2 * q(3,i2,j2,k) + &
                                         !coef_1 * q(3,i3,j3,k) + coef_2 * q(3,i4,j4,k)    )
                q(4,i,j,k) = 0.0!*one/two * ( coef_1 * q(4,i1,j1,k) + coef_2 * q(4,i2,j2,k) + &
                                        ! coef_1 * q(4,i3,j3,k) + coef_2 * q(4,i4,j4,k)    )
              end if
              
              if(l_o == 2 .OR. l_o == 6) then !Paredes laterales direccion y
                q(2,i,j,k) = 0.0!*one/two * ( coef_1 * q(2,i1,j1,k) + coef_2 * q(2,i2,j2,k) + &
                                         !coef_1 * q(2,i3,j3,k) + coef_2 * q(2,i4,j4,k)    )
                q(4,i,j,k) = 0.0!*one/two * ( coef_1 * q(4,i1,j1,k) + coef_2 * q(4,i2,j2,k) + &
                                         !coef_1 * q(4,i3,j3,k) + coef_2 * q(4,i4,j4,k)    )
                q(3,i,j,k) = 0.0
              end if
              

              if(l_o == 1 .OR. l_o == 3 .OR. l_o == 5 .OR. l_o == 7) then  !esquinas
                q(2:4,i,j,k) = 0.0!*one/two * ( coef_1 * q(2:4,i1,j1,k) + coef_2 * q(2:4,i2,j2,k) + &
                                           !coef_1 * q(2:4,i3,j3,k) + coef_2 * q(2:4,i4,j4,k)    )
              end if
              
            end do
          end do
        end do		
      end if
    end do
  end if
  
  

end subroutine bcond_obstacle

