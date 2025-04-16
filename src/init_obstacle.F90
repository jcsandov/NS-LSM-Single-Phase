subroutine init_obstacle ()

  use global
  use global_param
  use global_app
  use global_mpi
  use global_osg
  use checksum

  use global_obstacle
  use global_debug
  use global_lsm, only: nobstacles

implicit none

integer :: i,j,k,l,nz
integer :: i_ini_o, i_fin_o, j_ini_o, j_fin_o

if(nobstacles.ge.1) then

	call read_nodos() 
	call encontrar_nodos()
	call convec_quick_nodos()
	call obstacle_lsm()

else

  	allocate(act_obstacle_ad(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)))
	allocate(act_obstacle_rn(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)))
	allocate(csi_obs(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)),   &
			 eta_obs(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1))     )

	is_obstacle		= 0
	act_obstacle_ad = 1
	act_obstacle_rn = 1
	csi_obs 		= 0
	eta_obs			= 0

end if

!call outputD3_real(real(act_obstacle_ad(:,:,:)),'act')
!call outputD3_real(real(act_obstacle_rn(:,:,:)),'act_rn')
!call outputD3_real(real(boundary_obstacle(:,:,:,1)),'ob1')
!call outputD3_real(real(boundary_obstacle(:,:,:,2)),'ob2')
!call outputD3_real(real(boundary_obstacle(:,:,:,3)),'ob3')
!call outputD3_real(real(boundary_obstacle(:,:,:,4)),'ob4')

!debug prints

!print *, i_ini_o, i_fin_o
!print *, j_ini_o, j_fin_o

!l=8
!call mpi_barrier(mpi_comm_world, ierr)
!do nz = 0,nproc
!call mpi_barrier(mpi_comm_world, ierr)
!	if(myid == nz) then
!		print *, k_la_o(l) , k_lb_o(l) ,myid
!		print *, j_la_o(l) , j_lb_o(l) ,myid
!		print *, i_la_o(l) , i_lb_o(l), myid
!		print *, is_obstacle, obstacle_zone(l),myid
!		print *, " "
!	end if
!end do
!------------------------------------------------------------------------------------------

contains


subroutine read_nodos()

!leo las posiciones i,j de los nodos donde empieza el obstaculo
	implicit none

	open(60, file = "nodos_obstaculo", form = "formatted")
		read(60,*) i_ini_o, i_fin_o
		read(60,*) j_ini_o, j_fin_o
	close(60)
                
! TO DO: si se quisieran tener varios obstaculos, lo que habria que hacer es 
! transformar estas variables (i_ini_o...) en arreglos y hacer loops para irlos
! leyendo

end subroutine read_nodos


subroutine encontrar_nodos()
	!cada procesador encuentra sus nodos correspondientes de obstaculo(inicio fin,i,j,k) y si tiene
	!siempre pensado que en el plano k no se dividen los procesadores
	!8 regiones de posibles tipos de nodos obstaculo(archivo excel)
	
	implicit none
	
	integer :: I_l, I_u
	
	is_obstacle = 0
	obstacle_zone = 0
	k_la_o = 0 ; k_lb_o = 0;
	j_la_o = 0 ; j_lb_o = 0;
	i_la_o = 0 ; i_lb_o = 0;
	
	!Zona 1
	
	if( gi_ia(1) .LE. i_ini_o .AND. i_ini_o .LE. gi_ib(1) .AND. &
		gi_ja(1) .LE. j_ini_o .AND. j_ini_o .LE. gi_jb(1)        ) then
			obstacle_zone(1) = 1
			k_la_o(1) = li_ka(1)			   ; k_lb_o(1) = li_kb(1)
			j_la_o(1) = j_ini_o - gi_ja(1) + 1 ; j_lb_o(1) = j_ini_o - gi_ja(1) + 1 
			i_la_o(1) = i_ini_o - gi_ia(1) + 1 ; i_lb_o(1) = i_ini_o - gi_ia(1) + 1
			
			i_idx_o(1,1) = -1 ; i_idx_o(2,1) = -2 ; i_idx_o(3,1) =  0 ; i_idx_o(4,1) =  0
			j_idx_o(1,1) =  0 ; j_idx_o(2,1) =  0 ; j_idx_o(3,1) = -1 ; j_idx_o(4,1) = -2
	end if
	
	!Zona 2
	
	I_l = max(gi_ia(1), i_ini_o + 1)  !interseccion de intervalos (limite inferior)
	I_u = min(gi_ib(1), i_fin_o - 1)  !limite superior
	if( I_l .LE. I_u  .AND. &
		gi_ja(1) .LE. j_ini_o .AND. j_ini_o .LE. gi_jb(1)        ) then
			obstacle_zone(2) = 1
			k_la_o(2) = li_ka(1)			   		 ; k_lb_o(2) = li_kb(1)
			j_la_o(2) = j_ini_o - gi_ja(1) + 1 		 ; j_lb_o(2) = j_ini_o - gi_ja(1) + 1 
			i_la_o(2) = I_l - gi_ia(1) + 1			 ; i_lb_o(2) = I_u - gi_ia(1) + 1
			
			i_idx_o(1,2) =  0 ; i_idx_o(2,2) =  0 ; i_idx_o(3,2) =  0 ; i_idx_o(4,2) =  0
			j_idx_o(1,2) = -1 ; j_idx_o(2,2) = -2 ; j_idx_o(3,2) = -1 ; j_idx_o(4,2) = -2			
	end if
	
	!Zona 3
	
	if( gi_ia(1) .LE. i_fin_o .AND. i_fin_o .LE. gi_ib(1) .AND. &
		gi_ja(1) .LE. j_ini_o .AND. j_ini_o .LE. gi_jb(1)        ) then
			obstacle_zone(3) = 1
			k_la_o(3) = li_ka(1)			   ; k_lb_o(3) = li_kb(1)
			j_la_o(3) = j_ini_o - gi_ja(1) + 1 ; j_lb_o(3) = j_ini_o - gi_ja(1) + 1 
			i_la_o(3) = i_fin_o - gi_ia(1) + 1 ; i_lb_o(3) = i_fin_o - gi_ia(1) + 1
			
			i_idx_o(1,3) =  1 ; i_idx_o(2,3) =  2 ; i_idx_o(3,3) =  0 ; i_idx_o(4,3) =  0
			j_idx_o(1,3) =  0 ; j_idx_o(2,3) =  0 ; j_idx_o(3,3) = -1 ; j_idx_o(4,3) = -2			
	end if
	
	!Zona 4
	
	I_l = max(gi_ja(1), j_ini_o + 1)  !interseccion de intervalos (limite inferior)
	I_u = min(gi_jb(1), j_fin_o - 1)  !limite superior
	if( I_l .LE. I_u  .AND. &
		gi_ia(1) .LE. i_fin_o .AND. i_fin_o .LE. gi_ib(1)        ) then
			obstacle_zone(4) = 1
			k_la_o(4) = li_ka(1)			   		 ; k_lb_o(4) = li_kb(1)
			j_la_o(4) = I_l - gi_ja(1) + 1 		 	 ; j_lb_o(4) = I_u - gi_ja(1) + 1 
			i_la_o(4) = i_fin_o - gi_ia(1) + 1		 ; i_lb_o(4) = i_fin_o - gi_ia(1) + 1
			
			i_idx_o(1,4) =  1 ; i_idx_o(2,4) =  2 ; i_idx_o(3,4) =  1 ; i_idx_o(4,4) =  2
			j_idx_o(1,4) =  0 ; j_idx_o(2,4) =  0 ; j_idx_o(3,4) =  0 ; j_idx_o(4,4) =  0				
	end if
	
	!Zona 5
	
	if( gi_ia(1) .LE. i_fin_o .AND. i_fin_o .LE. gi_ib(1) .AND. &
		gi_ja(1) .LE. j_fin_o .AND. j_fin_o .LE. gi_jb(1)        ) then
			obstacle_zone(5) = 1
			k_la_o(5) = li_ka(1)			   ; k_lb_o(5) = li_kb(1)
			j_la_o(5) = j_fin_o - gi_ja(1) + 1 ; j_lb_o(5) = j_fin_o - gi_ja(1) + 1 
			i_la_o(5) = i_fin_o - gi_ia(1) + 1 ; i_lb_o(5) = i_fin_o - gi_ia(1) + 1
			
			i_idx_o(1,5) =  1 ; i_idx_o(2,5) =  2 ; i_idx_o(3,5) =  0 ; i_idx_o(4,5) =  0
			j_idx_o(1,5) =  0 ; j_idx_o(2,5) =  0 ; j_idx_o(3,5) =  1 ; j_idx_o(4,5) =  2			
	end if

	!Zona 6
	
	I_l = max(gi_ia(1), i_ini_o + 1)  !interseccion de intervalos (limite inferior)
	I_u = min(gi_ib(1), i_fin_o - 1)  !limite superior
	if( I_l .LE. I_u  .AND. &
		gi_ja(1) .LE. j_fin_o .AND. j_fin_o .LE. gi_jb(1)        ) then
			obstacle_zone(6) = 1
			k_la_o(6) = li_ka(1)			   		 ; k_lb_o(6) = li_kb(1)
			j_la_o(6) = j_fin_o - gi_ja(1) + 1 		 ; j_lb_o(6) = j_fin_o - gi_ja(1) + 1 
			i_la_o(6) = I_l - gi_ia(1) + 1			 ; i_lb_o(6) = I_u - gi_ia(1) + 1
			
			i_idx_o(1,6) =  0 ; i_idx_o(2,6) =  0 ; i_idx_o(3,6) =  0 ; i_idx_o(4,6) =  0
			j_idx_o(1,6) =  1 ; j_idx_o(2,6) =  2 ; j_idx_o(3,6) =  1 ; j_idx_o(4,6) =  2					
	end if			

	!Zona 7
	
	if( gi_ia(1) .LE. i_ini_o .AND. i_ini_o .LE. gi_ib(1) .AND. &
		gi_ja(1) .LE. j_fin_o .AND. j_fin_o .LE. gi_jb(1)        ) then
			obstacle_zone(7) = 1
			k_la_o(7) = li_ka(1)			   ; k_lb_o(7) = li_kb(1)
			j_la_o(7) = j_fin_o - gi_ja(1) + 1 ; j_lb_o(7) = j_fin_o - gi_ja(1) + 1 
			i_la_o(7) = i_ini_o - gi_ia(1) + 1 ; i_lb_o(7) = i_ini_o - gi_ia(1) + 1
			
			i_idx_o(1,7) = -1 ; i_idx_o(2,7) = -2 ; i_idx_o(3,7) =  0 ; i_idx_o(4,7) =  0
			j_idx_o(1,7) =  0 ; j_idx_o(2,7) =  0 ; j_idx_o(3,7) =  1 ; j_idx_o(4,7) =  2			
	end if
	
	!Zona 8
	
	I_l = max(gi_ja(1), j_ini_o + 1)  !interseccion de intervalos (limite inferior)
	I_u = min(gi_jb(1), j_fin_o - 1)  !limite superior
	if( I_l .LE. I_u  .AND. &
		gi_ia(1) .LE. i_ini_o .AND. i_ini_o .LE. gi_ib(1)        ) then
			obstacle_zone(8) = 1
			k_la_o(8) = li_ka(1)			   		 ; k_lb_o(8) = li_kb(1)
			j_la_o(8) = I_l - gi_ja(1) + 1 		 	 ; j_lb_o(8) = I_u - gi_ja(1) + 1 
			i_la_o(8) = i_ini_o - gi_ia(1) + 1		 ; i_lb_o(8) = i_ini_o - gi_ia(1) + 1
			
			i_idx_o(1,8) = -1 ; i_idx_o(2,8) = -2 ; i_idx_o(3,8) = -1 ; i_idx_o(4,8) = -2
			j_idx_o(1,8) =  0 ; j_idx_o(2,8) =  0 ; j_idx_o(3,8) =  0 ; j_idx_o(4,8) =  0	
	end if

	! Que el procesador este completamente contenido dentro del obstaculo

	if( i_ini_o .LE. gi_ia(1) .AND. gi_ib(1) .LE. i_fin_o .AND. &
	j_ini_o .LE. gi_ja(1) .AND. gi_jb(1) .LE. j_fin_o        ) then	
	
		is_obstacle = 1	

	end if

	if(sum(obstacle_zone) .GT. 0.OR.is_obstacle==1) then

	! If is_obstacle==1 antes de entrar en este if es porque el procesador esta completamente
	! dentro del obstaculo. De todas maneras se fija is_obstacle=1 por si se entra al if por
	! la primera condicion

		is_obstacle = 1	
	
		! Save the local indexes of each processor that is inside (including the borders) 
		! the obstacle to set velocities to zero in solver_daf (J. Sandoval, 2020)

		li_obs_ia=max(1,i_ini_o-gi_ia(1)+1)
		li_obs_ib=min(i_fin_o-gi_ia(1)+1,gi_ib(1)-gi_ia(1)+1)
	
		li_obs_ja=max(1,j_ini_o-gi_ja(1)+1)
		li_obs_jb=min(j_fin_o-gi_ja(1)+1,gi_jb(1)-gi_ja(1)+1)

		li_obs_ka=li_ka(1)
		li_obs_kb=li_kb(1)

    end if


end subroutine encontrar_nodos

subroutine convec_quick_nodos()

!identifica nodos para aplicar solo upwind en rhs_convec_quick
	Implicit none
	integer, dimension(6) :: l_o_csi, l_o_eta
	integer :: l_o, l_aux, i1,j1,i3,j3
	
	allocate(csi_obs(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)),   &
			 eta_obs(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1))     )

	l_o_csi = (/1, 8, 7, 3, 4, 5/)
	l_o_eta = (/1, 2, 3, 7, 6, 5/)	
	
	csi_obs = 0
	eta_obs = 0
	
	if(is_obstacle == 1) then
	  do l_aux = 1,6
      l_o = l_o_csi(l_aux) 
        if(obstacle_zone(l_o) == 1) then
          do k = k_la_o(l_o), k_lb_o(l_o)
            do j = j_la_o(l_o), j_lb_o(l_o)
              do i = i_la_o(l_o), i_lb_o(l_o)

                 i1 = i + i_idx_o(1,l_o)
                 j1 = j + j_idx_o(1,l_o)
                 csi_obs(i1, j1, k) = 1

              end do
            end do
          end do		
        end if
      end do
	end if
	
	if(is_obstacle == 1) then
	  do l_aux = 1,6
      l_o = l_o_eta(l_aux) 
        if(obstacle_zone(l_o) == 1) then
          do k = k_la_o(l_o), k_lb_o(l_o)
            do j = j_la_o(l_o), j_lb_o(l_o)
              do i = i_la_o(l_o), i_lb_o(l_o)

                 i3 = i + i_idx_o(3,l_o)
                 j3 = j + j_idx_o(3,l_o)
                 eta_obs(i3, j3, k) = 1

              end do
            end do
          end do		
        end if
      end do
	end if		 		 
			 
	
end subroutine

subroutine obstacle_lsm()

  implicit none
  integer :: i_loc, j_loc, k_loc
	allocate(obstacle_lsm_ad(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1),4))
	allocate(boundary_obstacle(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1),4))
  allocate(act_obstacle_ad(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)))
  allocate(act_obstacle_rn(le_ia(1):le_ib(1), le_ja(1):le_jb(1), le_ka(1):le_kb(1)))
  
  obstacle_lsm_ad = 0
  boundary_obstacle = 0
  act_obstacle_ad = 1
  
  do k = gi_ka(1), gi_kb(1)
    do j = gi_ja(1), gi_jb(1)
      do i = gi_ia(1), gi_ib(1)
      
        i_loc = i - gi_ia(1) + 1
        j_loc = j - gi_ja(1) + 1
        k_loc = k - gi_ka(1) + 1
        

        ! Aqui identifico en que sector del anillo de nodos
        ! estoy

        if(i_ini_o -1 .le. i .AND. i .le. i_fin_o +1 .AND. &
           j_ini_o -1 .le. j .AND. j .le. j_fin_o +1) then


        	if(i == i_ini_o-1) then
	          obstacle_lsm_ad(i_loc,j_loc,k_loc,1) = 1
        	end if
        	
        	if(i == i_fin_o+1) then
	          obstacle_lsm_ad(i_loc,j_loc,k_loc,3) = 1
        	end if

        	if(j == j_ini_o-1) then
	          obstacle_lsm_ad(i_loc,j_loc,k_loc,2) = 1
        	end if
        	
        	if(j == j_fin_o+1) then
	          obstacle_lsm_ad(i_loc,j_loc,k_loc,4) = 1
        	end if

        end if
        
        ! Veo los nodos locales que corresponden al borde

        if(i == i_ini_o .AND. &
           j_ini_o .le. j .AND. j .le. j_fin_o ) then
           
          boundary_obstacle(i_loc,j_loc,k_loc,1) = 1
        end if
        
        if(j == j_ini_o .AND. &
           i_ini_o .le. i .AND. i .le. i_fin_o ) then
           
          boundary_obstacle(i_loc,j_loc,k_loc,2) = 1
        end if
        
        if(i == i_fin_o .AND. &
           j_ini_o .le. j .AND. j .le. j_fin_o ) then
           
          boundary_obstacle(i_loc,j_loc,k_loc,3) = 1
        end if

        if(j == j_fin_o .AND. &
           i_ini_o .le. i .AND. i .le. i_fin_o ) then
           
          boundary_obstacle(i_loc,j_loc,k_loc,4) = 1
        end if

        ! Le indico que no aplique adveccion en los nodos que comprenden el 
        ! obstaculo (incluido el borde)
        
        if(i_ini_o .le. i .AND. i .le. i_fin_o .AND. &
           j_ini_o.le. j .AND. j .le. j_fin_o ) then
          act_obstacle_ad(i_loc,j_loc,k_loc) = 0
        end if   
        
                             
      end do
    end do
  end do
  
  act_obstacle_rn = act_obstacle_ad  !en este caso corresponde al mismo valor   
  
end subroutine obstacle_lsm



end subroutine init_obstacle
 