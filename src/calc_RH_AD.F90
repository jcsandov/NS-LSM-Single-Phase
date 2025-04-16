subroutine calc_RH_AD(phi_AD,rightH_AD,AdvectionNodes)

!calcula lado derecho usando esquema upwind
implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD
integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: AdvectionNodes 


!local

real (kind = rdf) :: ucon
real (kind = rdf), dimension(:,:,:,:), allocatable :: up
real (kind = rdf), dimension(:,:,:,:), allocatable :: um

!index
integer :: iRHAD_sta,iRHAD_end
integer :: jRHAD_sta,jRHAD_end
integer :: kRHAD_sta,kRHAD_end

integer :: iRHAD,jRHAD,kRHAD
logical, dimension(il:iu,jl:ju,kl:ku) :: is_boundary
integer, dimension(il:iu,jl:ju,kl:ku) :: is_boundary_int  !para conversion en debug. Borrar
!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


allocate (up(1:3,il:iu,jl:ju,kl:ku), &
            um(1:3,il:iu,jl:ju,kl:ku))

is_boundary = .TRUE.

!Nodos incluyendo el borde
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp


!Por cada procesador determino si es que esta en el borde (borde +3 los que no son ENO3)
if(orderAD == 1) then !para determinar si se usa orden 3 o orden 1 (si entra es orden 3)
  if (myback == mpi_proc_null)  iRHAD_sta = il + igp + 3
  if (myleft == mpi_proc_null)  jRHAD_sta = jl + jgp + 3
  if (mydown == mpi_proc_null)  kRHAD_sta = kl + kgp + 3

  if (myfront == mpi_proc_null) iRHAD_end = iu - igp - 3
  if (myright == mpi_proc_null) jRHAD_end = ju - jgp - 3
  if (myup    == mpi_proc_null) kRHAD_end = ku - kgp - 3

  do iRHAD=iRHAD_sta, iRHAD_end
  do jRHAD=jRHAD_sta, jRHAD_end
  do kRHAD=kRHAD_sta, kRHAD_end
    is_boundary(iRHAD,jRHAD,kRHAD) = .FALSE.
  end do
  end do
  end do
end if

!is_boundary_int = is_boundary
!call outputD3_real(real(is_boundary_int),'is_boundary')


do kRHAD=kl,ku
do jRHAD=jl,ju
do iRHAD=il,iu
   
   if ( AdvectionNodes(iRHAD,jRHAD,kRHAD) > 0 ) then

        !csi direction
        ucon=one_half*ucn_j(1,iRHAD,jRHAD,kRHAD)
        up(1,iRHAD,jRHAD,kRHAD)=ucon+abs(ucon)
        um(1,iRHAD,jRHAD,kRHAD)=ucon-abs(ucon)

        !eta direction
        ucon=one_half*ucn_j(2,iRHAD,jRHAD,kRHAD)
        up(2,iRHAD,jRHAD,kRHAD)=ucon+abs(ucon)
        um(2,iRHAD,jRHAD,kRHAD)=ucon-abs(ucon)

        !zet direction
        ucon=one_half*ucn_j(3,iRHAD,jRHAD,kRHAD)
        up(3,iRHAD,jRHAD,kRHAD)=ucon+abs(ucon)
        um(3,iRHAD,jRHAD,kRHAD)=ucon-abs(ucon)

   end if

end do
end do
end do

!contribucion csi
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp

!bordes
if(myback == mpi_proc_null) then
        do jRHAD=jRHAD_sta, jRHAD_end
        do kRHAD=kRHAD_sta, kRHAD_end

         if ( AdvectionNodes(iRHAD_sta,jRHAD,kRHAD) > 0 ) then

                rightH_AD(iRHAD_sta,jRHAD,kRHAD) =dc * ucn_j(1,iRHAD_sta,jRHAD,kRHAD)  * &
                (phi_AD(iRHAD_sta + 1,jRHAD, kRHAD) - phi_AD(iRHAD_sta,jRHAD,kRHAD) )
         
         end if

        end do
        end do
end if

if(myfront == mpi_proc_null) then
        do jRHAD=jRHAD_sta, jRHAD_end
        do kRHAD=kRHAD_sta, kRHAD_end

         if ( AdvectionNodes(iRHAD_end,jRHAD,kRHAD) > 0 ) then

                rightH_AD(iRHAD_end,jRHAD,kRHAD) =dc * ucn_j(1,iRHAD_end,jRHAD,kRHAD)  * &
                (phi_AD(iRHAD_end,jRHAD, kRHAD) - phi_AD(iRHAD_end-1,jRHAD,kRHAD) )

         end if
        end do
        end do
end if

!no bordes
if(myback == mpi_proc_null)  iRHAD_sta = il + igp + 1
if(myfront == mpi_proc_null) iRHAD_end = iu - igp - 1

do iRHAD=iRHAD_sta, iRHAD_end
do jRHAD=jRHAD_sta, jRHAD_end
do kRHAD=kRHAD_sta, kRHAD_end

   if ( AdvectionNodes(iRHAD,jRHAD,kRHAD) > 0 ) then
   
     if(is_boundary(iRHAD,jRHAD,kRHAD)) &
            rightH_AD(iRHAD,jRHAD,kRHAD) = dc * (&
            up(1,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD,jRHAD,kRHAD) - phi_AD(iRHAD - 1,jRHAD,kRHAD)) + &
            um(1,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD + 1,jRHAD,kRHAD) - phi_AD(iRHAD,jRHAD,kRHAD))  )
   
   end if

end do 
end do
end do

!--------------------------------------------------------------------------------------------- 
!contribucion eta
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp

!bordes
if(myleft== mpi_proc_null) then
        do iRHAD=iRHAD_sta, iRHAD_end
        do kRHAD=kRHAD_sta, kRHAD_end
         if ( AdvectionNodes(iRHAD,jRHAD_sta,kRHAD) > 0 ) then

                rightH_AD(iRHAD,jRHAD_sta,kRHAD) =de * ucn_j(2,iRHAD,jRHAD_sta,kRHAD)  * &
                (phi_AD(iRHAD,jRHAD_sta + 1, kRHAD) - phi_AD(iRHAD,jRHAD_sta,kRHAD) )  + &
                rightH_AD(iRHAD,jRHAD_sta,kRHAD)
         end if
        end do
        end do
end if

if(myright == mpi_proc_null) then
        do iRHAD=iRHAD_sta, iRHAD_end
        do kRHAD=kRHAD_sta, kRHAD_end
               if ( AdvectionNodes(iRHAD,jRHAD_end,kRHAD) > 0 ) then

                rightH_AD(iRHAD,jRHAD_end,kRHAD) =de * ucn_j(2,iRHAD,jRHAD_end,kRHAD)  * &
                (phi_AD(iRHAD,jRHAD_end, kRHAD) - phi_AD(iRHAD,jRHAD_end - 1,kRHAD) )  + &
                rightH_AD(iRHAD,jRHAD_end,kRHAD)

               end if
        end do
        end do
end if

!no bordes
if(myleft == mpi_proc_null)  jRHAD_sta = jl + jgp + 1
if(myright == mpi_proc_null) jRHAD_end = ju - jgp - 1

do iRHAD=iRHAD_sta, iRHAD_end
do jRHAD=jRHAD_sta, jRHAD_end
do kRHAD=kRHAD_sta, kRHAD_end
   if ( AdvectionNodes(iRHAD,jRHAD,kRHAD) > 0 ) then
   
     if(is_boundary(iRHAD,jRHAD,kRHAD)) &    
            rightH_AD(iRHAD,jRHAD,kRHAD) = de * (&
            up(2,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD,jRHAD,kRHAD) - phi_AD(iRHAD,jRHAD - 1,kRHAD)) + &
            um(2,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD,jRHAD + 1,kRHAD) - phi_AD(iRHAD,jRHAD,kRHAD))  ) +  &
            rightH_AD (iRHAD,jRHAD,kRHAD)
   end if
end do 
end do
end do
   
!------------------------------------------------------------------------------------------------
!contribucion zet
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp

!bordes
if(mydown == mpi_proc_null) then
        do jRHAD=jRHAD_sta, jRHAD_end
        do iRHAD=iRHAD_sta, iRHAD_end

           if ( AdvectionNodes(iRHAD,jRHAD,kRHAD_sta) > 0 ) then

                rightH_AD(iRHAD,jRHAD,kRHAD_sta) =dz * ucn_j(3,iRHAD,jRHAD,kRHAD_sta)  * &
                (phi_AD(iRHAD,jRHAD, kRHAD_sta + 1) - phi_AD(iRHAD,jRHAD,kRHAD_sta) )  + &
                rightH_AD(iRHAD,jRHAD,kRHAD_sta)

           end if
        end do
        end do
end if

if(myup == mpi_proc_null) then
        do jRHAD=jRHAD_sta, jRHAD_end
        do iRHAD=iRHAD_sta, iRHAD_end
            if ( AdvectionNodes(iRHAD,jRHAD,kRHAD_end) > 0 ) then

                rightH_AD(iRHAD,jRHAD,kRHAD_end) =dz * ucn_j(3,iRHAD,jRHAD,kRHAD_end)  * &
                (phi_AD(iRHAD,jRHAD, kRHAD_end) - phi_AD(iRHAD,jRHAD,kRHAD_end - 1) )  + &
                rightH_AD(iRHAD,jRHAD,kRHAD_end)
            
            end if
        end do
        end do
end if

!no bordes
if(mydown == mpi_proc_null)  kRHAD_sta = kl + kgp + 1
if(myup == mpi_proc_null)    kRHAD_end = ku - kgp - 1

do iRHAD=iRHAD_sta, iRHAD_end
do jRHAD=jRHAD_sta, jRHAD_end
do kRHAD=kRHAD_sta, kRHAD_end
   if ( AdvectionNodes(iRHAD,jRHAD,kRHAD) > 0 ) then
   
     if(is_boundary(iRHAD,jRHAD,kRHAD)) &
            rightH_AD(iRHAD,jRHAD,kRHAD) = dz * (&
            up(3,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD,jRHAD,kRHAD) - phi_AD(iRHAD,jRHAD,kRHAD - 1)) + &
            um(3,iRHAD,jRHAD,kRHAD) * (phi_AD(iRHAD,jRHAD,kRHAD + 1) - phi_AD(iRHAD,jRHAD,kRHAD))  ) + &
            rightH_AD(iRHAD,jRHAD,kRHAD)
   
   end if
end do 
end do
end do

deallocate (up, um)

end subroutine calc_RH_AD
