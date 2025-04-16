subroutine calc_RH_AD_obstacle(phi_AD,rightH_AD)

!calcula lado derecho usando esquema upwind
!en la zona alrededor de donde se encuentra el obstaculo

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD


!local

real (kind = rdf) :: ucon
real (kind = rdf), dimension(:,:,:,:), allocatable :: up
real (kind = rdf), dimension(:,:,:,:), allocatable :: um

!index
integer :: iRHAD_sta,iRHAD_end
integer :: jRHAD_sta,jRHAD_end
integer :: kRHAD_sta,kRHAD_end


!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


allocate (up(1:3,il:iu,jl:ju,kl:ku), &
            um(1:3,il:iu,jl:ju,kl:ku))


do k=kl,ku
do j=jl,ju
do i=il,iu

  !csi direction
  ucon=one_half*ucn_j(1,i,j,k)
  up(1,i,j,k)=ucon+abs(ucon)
  um(1,i,j,k)=ucon-abs(ucon)
  if(boundary_obstacle(i,j,k,1) == 1 .AND. ucn_j(1,i,j,k) .le. zero) then  !fuerzo a que calcule hacia adentro del dominio en el borde
    up(1,i,j,k) = um(1,i,j,k)
    um(1,i,j,k) = zero
  end if
  
  if(boundary_obstacle(i,j,k,3) == 1 .AND. ucn_j(1,i,j,k) .ge. zero) then  !idem
    um(1,i,j,k) = up(1,i,j,k)
    up(1,i,j,k) = zero
  end if
  
  !eta direction
  ucon=one_half*ucn_j(2,i,j,k)
  up(2,i,j,k)=ucon+abs(ucon)
  um(2,i,j,k)=ucon-abs(ucon)
  
  if(boundary_obstacle(i,j,k,2) == 1 .AND. ucn_j(2,i,j,k) .le. zero) then  !idem
    up(2,i,j,k) = um(2,i,j,k)
    um(2,i,j,k) = zero
  end if  

  if(boundary_obstacle(i,j,k,4) == 1 .AND. ucn_j(2,i,j,k) .ge. zero) then  !idem
    um(2,i,j,k) = up(2,i,j,k)
    up(2,i,j,k) = zero
  end if
  
  !zet direction    
  ucon=one_half*ucn_j(3,i,j,k)
  up(3,i,j,k)=ucon+abs(ucon)
  um(3,i,j,k)=ucon-abs(ucon)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          WARNING          !         
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Se asume que el eje z esta contenido en solo un procesador
  
  if(mydown == mpi_proc_null .AND. k == kl + kgp) then 
    up(3,i,j,k) = zero
    um(3,i,j,k) = ucn_j(3,i,j,k)
  end if
  
  if(myup == mpi_proc_null .AND. k == ku -kgp) then
    up(3,i,j,k) = ucn_j(3,i,j,k)
    um(3,i,j,k) = zero
  end if
      
end do
end do
end do


  
!call outputD3_real(ucn_j(1,:,:,:), 'ucn1')
      
!contribucion csi
!----------------------------------------------
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp


do i=iRHAD_sta, iRHAD_end
do j=jRHAD_sta, jRHAD_end
do k=kRHAD_sta, kRHAD_end
  if(obstacle_lsm_ad(i,j,k) == 1) &
   rightH_AD(i,j,k) = dc * (&
                      up(1,i,j,k) * (phi_AD(i,j,k) - phi_AD(i - 1,j,k)) + &
                      um(1,i,j,k) * (phi_AD(i + 1,j,k) - phi_AD(i,j,k))  )
                      
end do 
end do
end do

!contribucion eta
!----------------------------------------------
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp


do i=iRHAD_sta, iRHAD_end
do j=jRHAD_sta, jRHAD_end
do k=kRHAD_sta, kRHAD_end
  if(obstacle_lsm_ad(i,j,k) == 1) &
   rightH_AD(i,j,k) = de * (&
                      up(2,i,j,k) * (phi_AD(i,j,k) - phi_AD(i,j - 1,k)) +    &
                      um(2,i,j,k) * (phi_AD(i,j + 1,k) - phi_AD(i,j,k))  ) + &
                      rightH_AD(i,j,k)
                      
end do 
end do
end do

!contribucion zet
!----------------------------------------------
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

if(mydown == mpi_proc_null) kRHAD_sta = kl + kgp + 1

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp

if(myup == mpi_proc_null) kRHAD_end = ku - kgp -1


do i=iRHAD_sta, iRHAD_end
do j=jRHAD_sta, jRHAD_end
do k=kRHAD_sta, kRHAD_end
  if(obstacle_lsm_ad(i,j,k) == 1) &
   rightH_AD(i,j,k) = dz * (&
                      up(3,i,j,k) * (phi_AD(i,j,k) - phi_AD(i,j,k - 1)) +    &
                      um(3,i,j,k) * (phi_AD(i,j,k + 1) - phi_AD(i,j,k))  ) + &
                      rightH_AD(i,j,k)
                      
end do 
end do
end do

!call index_outputD3(116,91,61,real(obstacle_lsm_ad),debug_values(7))
!if(myid == 2) print *, ' '
!call index_outputD3(116,91,61,up(3,:,:,:),debug_values(1))
!call index_outputD3(116,91,61,um(3,:,:,:),debug_values(2))
!call index_outputD3(116,91,61,real(boundary_obstacle(:,:,:,4)),debug_values(3))
!call index_outputD3(116,91,61,phi_AD(:,:,:),debug_values(4))
!call index_outputD3(116,91,60,phi_AD(:,:,:),debug_values(5))
!call index_outputD3(116,91,61,rightH_AD(:,:,:),debug_values(6))

!if(myid == 2) print *, 'byhand', debug_values(1)*(debug_values(4) -debug_values(5)), myid
!--------------------------------------------------------------------------------------------- 

 deallocate (up, um)

end subroutine calc_RH_AD_obstacle
