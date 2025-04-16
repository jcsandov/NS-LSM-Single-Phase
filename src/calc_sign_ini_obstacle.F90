subroutine calc_sign_ini_obstacle(sgndf_RN,signo_RN,normagrad_RN,epslsm,M_Sphi0)

!calcula la norma del gradiente (mas y menos) de una signed distance function
!justo en los bordes del obstaculo

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: sgndf_RN
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: signo_RN, normagrad_RN
real (kind = rdf) :: epslsm

!-------------------------------------------------------------------------------------
! Switches para probar modificaciones en la reinicializacion

integer, intent(in) :: M_Sphi0

!local

real (kind = rdf)  :: normagradfmas, normagradfmenos, normagradf
!integer :: i,j,k
real (kind = rdf) :: dxb,dxf,dyb,dyf,dzb,dzf
real (kind = rdf) :: derx,dery,derz
real (kind = rdf), dimension (il:iu,jl:ju,kl:ku) :: dcsimas, dcsimenos
real (kind = rdf), dimension (il:iu,jl:ju,kl:ku) :: detamas, detamenos
real (kind = rdf), dimension (il:iu,jl:ju,kl:ku) :: dzetmas, dzetmenos
real (kind = rdf) :: dxbmas,dxbmenos, dxfmas,dxfmenos
real (kind = rdf) :: dybmas,dybmenos, dyfmas,dyfmenos
real (kind = rdf) :: dzbmas,dzbmenos, dzfmas,dzfmenos
real (kind = rdf) :: aux


!index
integer :: i_sta,i_end
integer :: j_sta,j_end
integer :: k_sta,k_end


!------------------------------------------------------
!direccion csi
!------------------------------------------------------
!Nodos incluyendo el borde
i_sta = il + igp
j_sta = jl + jgp
k_sta = kl + kgp

i_end = iu - igp        
j_end = ju - jgp
k_end = ku - kgp

do i = i_sta, i_end
do j = j_sta, j_end
do k = k_sta, k_end

        if(boundary_obstacle(i,j,k,1) == 1) then
          dcsimenos(i,j,k) = dc * (sgndf_RN(i,j,k) - sgndf_RN(i-1,j,k))
          dcsimas(i,j,k) = dcsimenos(i,j,k)
        end if
        
        if(boundary_obstacle(i,j,k,3) == 1) then
          dcsimas(i,j,k) = dc * (sgndf_RN(i+1,j,k) - sgndf_RN(i,j,k))
          dcsimenos(i,j,k) = dcsimas(i,j,k)
        end if
end do
end do
end do

!------------------------------------------------------
!direccion eta
!------------------------------------------------------
!Nodos incluyendo el borde
i_sta = il + igp
j_sta = jl + jgp
k_sta = kl + kgp

i_end = iu - igp        
j_end = ju - jgp
k_end = ku - kgp

do i = i_sta, i_end
do j = j_sta, j_end
do k = k_sta, k_end

        if(boundary_obstacle(i,j,k,2) == 1) then
          detamenos(i,j,k) = de * (sgndf_RN(i,j,k) - sgndf_RN(i,j-1,k))
          detamas(i,j,k) = detamenos(i,j,k)
        end if
        
        if(boundary_obstacle(i,j,k,4) == 1) then
          detamas(i,j,k) = de * (sgndf_RN(i,j+1,k) - sgndf_RN(i,j,k))
          detamenos(i,j,k) = detamas(i,j,k)
        end if
end do
end do
end do

!------------------------------------------------------
!direccion zet
!------------------------------------------------------
!Nodos incluyendo el borde
i_sta = il + igp
j_sta = jl + jgp
k_sta = kl + kgp

i_end = iu - igp        
j_end = ju - jgp
k_end = ku - kgp

if(mydown == mpi_proc_null) k_sta = k_sta + 1
if(myup == mpi_proc_null) k_end = k_end -1

do i = i_sta, i_end
do j = j_sta, j_end
do k = k_sta, k_end
        
        if(boundary_obstacle(i,j,k,1) == 1 .OR. boundary_obstacle(i,j,k,3) == 1 .OR. &
           boundary_obstacle(i,j,k,2) == 1 .OR. boundary_obstacle(i,j,k,4) == 1) then
           
          dzetmas(i,j,k) = dz * (sgndf_RN(i,j,k+1) - sgndf_RN(i,j,k))
          dzetmenos(i,j,k) = dz * (sgndf_RN(i,j,k) - sgndf_RN(i,j,k-1))
        end if
end do
end do
end do

k_sta = kl + kgp
k_end = ku - kgp

!bordes
if(mydown == mpi_proc_null) then
      do i =i_sta, i_end
        do j =j_sta, j_end
        
          if(boundary_obstacle(i,j,k,1) == 1 .OR. boundary_obstacle(i,j,k,3) == 1 .OR. &
             boundary_obstacle(i,j,k,2) == 1 .OR. boundary_obstacle(i,j,k,4) == 1) then
             
                  dzetmas(i,j,k_sta) = dz * (sgndf_RN(i,j,k_sta+1) - sgndf_RN(i,j,k_sta))
                  dzetmenos(i,j,k_sta) = dzetmas(i,j,k_sta)
          end if
          
        end do
      end do
end if

if(myup == mpi_proc_null) then
      do i =i_sta, i_end
        do j =j_sta, j_end
        
          if(boundary_obstacle(i,j,k,1) == 1 .OR. boundary_obstacle(i,j,k,3) == 1 .OR. &
             boundary_obstacle(i,j,k,2) == 1 .OR. boundary_obstacle(i,j,k,4) == 1) then        
                dzetmenos(i,j,k_end) = dz * (sgndf_RN(i,j,k_end) - sgndf_RN(i,j,k_end-1))
                dzetmas(i,j,k_end) = dzetmenos(i,j,k_end)
          end if
                 
        end do
      end do
end if

!ahora calculo la funcion signo
!-----------------------------------------------------
i_sta = il + igp
j_sta = jl + jgp
k_sta = kl + kgp

i_end = iu - igp        
j_end = ju - jgp
k_end = ku - kgp

do i = i_sta, i_end
do j = j_sta, j_end
do k = k_sta, k_end
  if(boundary_obstacle(i,j,k,1) == 1 .OR. boundary_obstacle(i,j,k,3) == 1 .OR. &
     boundary_obstacle(i,j,k,2) == 1 .OR. boundary_obstacle(i,j,k,4) == 1) then 
        !grad sgndf
        dxb = csi(1,i,j,k) * dcsimenos(i,j,k) + &
              eta(1,i,j,k) * detamenos(i,j,k) + &
              zet(1,i,j,k) * dzetmenos(i,j,k)

        dxf = csi(1,i,j,k) * dcsimas(i,j,k) + &
              eta(1,i,j,k) * detamas(i,j,k) + &
              zet(1,i,j,k) * dzetmas(i,j,k)

        dyb = csi(2,i,j,k) * dcsimenos(i,j,k) + &
              eta(2,i,j,k) * detamenos(i,j,k) + &
              zet(2,i,j,k) * dzetmenos(i,j,k)

        dyf = csi(2,i,j,k) * dcsimas(i,j,k) + &
              eta(2,i,j,k) * detamas(i,j,k) + &
              zet(2,i,j,k) * dzetmas(i,j,k)

        dzb = csi(3,i,j,k) * dcsimenos(i,j,k) + &
              eta(3,i,j,k) * detamenos(i,j,k) + &
              zet(3,i,j,k) * dzetmenos(i,j,k)

        dzf = csi(3,i,j,k) * dcsimas(i,j,k) + &
              eta(3,i,j,k) * detamas(i,j,k) + &
              zet(3,i,j,k) * dzetmas(i,j,k)

        
        aux = dxb/two ; dxbmas = aux + abs(aux) ; dxbmenos = aux - abs(aux)
        aux = dxf/two ; dxfmas = aux + abs(aux) ; dxfmenos = aux - abs(aux)
        aux = dyb/two ; dybmas = aux + abs(aux) ; dybmenos = aux - abs(aux)
        aux = dyf/two ; dyfmas = aux + abs(aux) ; dyfmenos = aux - abs(aux)
        aux = dzb/two ; dzbmas = aux + abs(aux) ; dzbmenos = aux - abs(aux)
        aux = dzf/two ; dzfmas = aux + abs(aux) ; dzfmenos = aux - abs(aux)

        normagradfmas = sqrt( max(dxbmas**2,dxfmenos**2) + &
                              max(dybmas**2,dyfmenos**2) + &
                              max(dzbmas**2,dzfmenos**2) )

        normagradfmenos = sqrt( max(dxbmenos**2,dxfmas**2) + &
                                max(dybmenos**2,dyfmas**2) + &
                                max(dzbmenos**2,dzfmas**2) )

        if(sgndf_RN(i,j,k) > zero) then
                normagradf = normagradfmas
        else if (sgndf_RN(i,j,k) < zero ) then
                normagradf = normagradfmenos
        else
                normagradf = one
        end if

        ! Calculo de s(phi0), se puede usar el modelo de Kang (piecewise)
        ! o el modelo original continuo.

        if (M_Sphi0==1) then ! s(phi0) Kang
            
            if(sgndf_RN(i,j,k).ge.epslsm) then
                signo_RN(i,j,k) = one
            else if (sgndf_RN(i,j,k).le.(-one)*epslsm) then
                signo_RN(i,j,k) = -one
            else
                signo_RN(i,j,k) = sgndf_RN(i,j,k)/epslsm-one/pi*sin(pi*sgndf_RN(i,j,k)/epslsm)
            end if

        else ! version original de Carvajal
            signo_RN(i,j,k) = sgndf_RN(i,j,k)/sqrt((sgndf_RN(i,j,k))**2 + (normagradf *epslsm)**2)
        end if

        normagrad_RN(i,j,k) = normagradf

  end if
end do
end do
end do


end subroutine calc_sign_ini_obstacle
