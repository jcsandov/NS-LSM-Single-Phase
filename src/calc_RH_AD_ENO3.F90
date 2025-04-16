subroutine calc_RH_AD_ENO3(phi_AD,rightH_AD)
!Calcula lado derecho de ecuacion de adveccion

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD


!local

real (kind = rdf) :: ucon
real (kind = rdf), dimension(:,:,:,:), allocatable :: up
real (kind = rdf), dimension(:,:,:,:), allocatable :: um
real (kind = rdf), dimension(:,:,:,:,:), allocatable :: Hf

real (kind = rdf), dimension(:,:,:,:), allocatable :: RSP !Roe Speed
integer, dimension(:,:,:,:,:), allocatable :: kmin !valores de kmin
real (kind = rdf), dimension(:,:,:,:), allocatable :: flux 
real (kind = rdf) :: aENO2, bENO2, cENO2, aENO3, bENO3, cENO3
integer :: kmn1,kmn2,kmn3
real (kind = rdf) :: aux1,aux2,aux3 !valores auxiliares para algunos calculos y hacer mas legible el codigo

!index
integer, dimension(4) :: iRHAD_sta,iRHAD_end
integer, dimension(4) :: jRHAD_sta,jRHAD_end
integer, dimension(4) :: kRHAD_sta,kRHAD_end

integer :: iRHAD,jRHAD,kRHAD
integer :: i_idxa,j_idxa,k_idxa

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


allocate (up(1:3,il:iu,jl:ju,kl:ku), &
            um(1:3,il:iu,jl:ju,kl:ku))

allocate (Hf(1:3,il:iu,jl:ju,kl:ku,1:3)  , & !ultima coordenada refiere al orden del ENO
          RSP(1:3,il:iu,jl:ju,kl:ku)     , &    
          kmin(1:3,il:iu,jl:ju,kl:ku,1:3), &
          flux(1:3,il:iu,jl:ju,kl:ku)      ) 


!Nodos incluyendo el borde
iRHAD_sta(1:4) = il + igp
jRHAD_sta(1:4) = jl + jgp
kRHAD_sta(1:4) = kl + kgp

iRHAD_end(1:4) = iu - igp        
jRHAD_end(1:4) = ju - jgp
kRHAD_end(1:4) = ku - kgp
!Nodos sin incluir borde

if (myback == mpi_proc_null)  iRHAD_sta(2) = il + igp + 1
if (myleft == mpi_proc_null)  jRHAD_sta(2) = jl + jgp + 1
if (mydown == mpi_proc_null)  kRHAD_sta(2) = kl + kgp + 1

if (myback == mpi_proc_null)  iRHAD_sta(3) = il + igp + 2
if (myleft == mpi_proc_null)  jRHAD_sta(3) = jl + jgp + 2
if (mydown == mpi_proc_null)  kRHAD_sta(3) = kl + kgp + 2

if (myback == mpi_proc_null)  iRHAD_sta(4) = il + igp + 3
if (myleft == mpi_proc_null)  jRHAD_sta(4) = jl + jgp + 3
if (mydown == mpi_proc_null)  kRHAD_sta(4) = kl + kgp + 3

if (myfront == mpi_proc_null) iRHAD_end(2) = iu - igp - 1
if (myright == mpi_proc_null) jRHAD_end(2) = ju - jgp - 1
if (myup    == mpi_proc_null) kRHAD_end(2) = ku - kgp - 1

if (myfront == mpi_proc_null) iRHAD_end(3) = iu - igp - 2
if (myright == mpi_proc_null) jRHAD_end(3) = ju - jgp - 2
if (myup    == mpi_proc_null) kRHAD_end(3) = ku - kgp - 2

if (myfront == mpi_proc_null) iRHAD_end(4) = iu - igp - 3
if (myright == mpi_proc_null) jRHAD_end(4) = ju - jgp - 3
if (myup    == mpi_proc_null) kRHAD_end(4) = ku - kgp - 3


!calcular tabla de diferencias ([l-1/2,..., l+ k+ 1/2])
!i index
!Eno 1
do iRHAD = iRHAD_sta(1),iRHAD_end(1)
do jRHAD = jRHAD_sta(2),jRHAD_end(2)
do kRHAD = kRHAD_sta(2),kRHAD_end(2)
        Hf(1,iRHAD,jRHAD,kRHAD,1) = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
        
end do
end do
end do

call rhs_exchng3_3d (Hf(1,:,:,:,1))

!Eno 2
do iRHAD = iRHAD_sta(1),iRHAD_end(2)
do jRHAD = jRHAD_sta(2),jRHAD_end(2)
do kRHAD = kRHAD_sta(2),kRHAD_end(2)
        Hf(1,iRHAD,jRHAD,kRHAD,2) = one/two * dc * &
                                   (phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD) - &
                                    phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD))
end do
end do
end do

call rhs_exchng3_3d (Hf(1,:,:,:,2))

!Eno 3

do iRHAD = iRHAD_sta(1), iRHAD_end(3)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do kRHAD = kRHAD_sta(2), kRHAD_end(2)
        Hf(1,iRHAD,jRHAD,kRHAD,3) = one/three * dc * &
                                    (Hf(1,iRHAD+1,jRHAD,kRHAD,2) - Hf(1,iRHAD,jRHAD,kRHAD,2))
end do
end do
end do

call rhs_exchng3_3d (Hf(1,:,:,:,3))
!calcular Roe Speed (podria ser incluido en el anterior). Almaceno en i+1/2
         
do iRHAD = iRHAD_sta(1), iRHAD_end(2)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do kRHAD = kRHAD_sta(2), kRHAD_end(2)


        RSP(1,iRHAD,jRHAD,kRHAD) = (ucn_j(1,iRHAD+1,jRHAD,kRHAD)*phi_AD(iRHAD+1,jRHAD,kRHAD) - &
                                   ucn_j(1,iRHAD,jRHAD,kRHAD)*phi_AD(iRHAD,jRHAD,kRHAD)  ) / &
                                    (phi_AD(iRHAD+1,jRHAD,kRHAD)  -   phi_AD(iRHAD,jRHAD,kRHAD) )   

end do
end do
end do

!ghost point de ambos (creo que en estricto rigor basta solo los de HF)



call rhs_exchng3_3d (RSP(1,:,:,:))

!call outputD3_real(RSP(1,:,:,:),'rsp')
!call outputD3_real(Hf(1,:,:,:,2),'HF2')
!ahora calculo flux
!bordes

!----------------------------------------------------------
!i back
if (myback == mpi_proc_null)  then

        !ENO1
        i_idxa = iRHAD_sta(1)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(1,i_idxa,jRHAD,kRHAD) .ge. zero) then
                        kmn1 = i_idxa
                else
                        kmn1 = i_idxa + 1
                end if
                
                flux(1,i_idxa,jRHAD,kRHAD) = Hf(1,kmn1,jRHAD,kRHAD,1)                        
        end do
        end do

        !ENO2
        i_idxa = iRHAD_sta(2)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(1,i_idxa,jRHAD,kRHAD) .ge. zero) then
                        kmn1 = i_idxa
                else
                        kmn1 = i_idxa + 1
                end if
                
                aENO2 = Hf(1,kmn1,jRHAD,kRHAD,2)
                bENO2 = Hf(1,kmn1-1,jRHAD,kRHAD,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(1,i_idxa,jRHAD,kRHAD) = Hf(1,kmn1,jRHAD,kRHAD,1) + &
                                            cENO2* (dc**(-1)*real((i_idxa - kmn1),rdf) + &
                                                   dc**(-1)*real((i_idxa -kmn1 + 1),rdf))

        end do
        end do
end if



!i front
if (myfront == mpi_proc_null)  then

        !ENO1
        i_idxa = iRHAD_end(2)  !un nodo antes

        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(1,i_idxa,jRHAD,kRHAD) .ge. zero) then
                        kmn1 = i_idxa
                else
                        kmn1 = i_idxa + 1
                end if

                flux(1,i_idxa,jRHAD,kRHAD) = Hf(1,kmn1,jRHAD,kRHAD,1)                        
        end do
        end do

        !ENO2
        i_idxa = iRHAD_end(3)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(1,i_idxa,jRHAD,kRHAD) .ge. zero) then
                        kmn1 = i_idxa
                else
                        kmn1 = i_idxa + 1
                end if
                
                aENO2 = Hf(1,kmn1,jRHAD,kRHAD,2)
                bENO2 = Hf(1,kmn1-1,jRHAD,kRHAD,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(1,i_idxa,jRHAD,kRHAD) = Hf(1,kmn1,jRHAD,kRHAD,1) + &
                                            cENO2* (dc**(-1)*real((i_idxa - kmn1),rdf) + &
                                                   dc**(-1)*real((i_idxa -kmn1 + 1),rdf))

        end do
        end do
end if

!call outputD3_real(flux(1,:,:,:),'flux1')


!ENO3
do KRHAD = KRHAD_sta(2), KRHAD_end(2)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do iRHAD = iRHAD_sta(3), iRHAD_end(4)
!if(myid == 1 .and.  jRHAD == 5 .and. kRHAD == 51) print *, iRHAD, myid 


        if(RSP(1,iRHAD,jRHAD,kRHAD) .ge. zero) then
                kmn1 = iRHAD                                
        else
                kmn1 = iRHAD + 1
        end if
                
        aENO2 = Hf(1,kmn1,jRHAD,kRHAD,2)
        bENO2 = Hf(1,kmn1-1,jRHAD,kRHAD,2)

        if(abs(aENO2) .ge. abs(bENO2)) then
                cENO2 = bENO2
                kmn2 = kmn1 - 1
        else
                cENO2 = aENO2
                kmn2 = kmn1
        end if

        aENO3 = Hf(1,kmn2,jRHAD,kRHAD,3)
        bENO3 = Hf(1,kmn2-1,jRHAD,kRHAD,3)

        if(abs(aENO3) .ge. abs(bENO3)) then
                cENO3 = bENO3
                kmn3 = kmn2 - 1
        else
                cENO3 = aENO3
                kmn3 = kmn2
        end if

        aux1 = dc**(-1) * real(iRHAD - kmn2 +1,rdf)
        aux2 = dc**(-1) * real(iRHAD - kmn2,rdf)
        aux3 = dc**(-1) * real(iRHAD - kmn2 -1,rdf)        

        flux(1,iRHAD,jRHAD,kRHAD) = Hf(1,kmn1,jRHAD,kRHAD,1) + &
                                    cENO2* (dc**(-1)*real((iRHAD - kmn1),rdf) + &
                                           dc**(-1)*real((iRHAD -kmn1 + 1),rdf)) + &
                                    cENO3* ( (aux2 +aux1)*aux3 + aux1*aux2 )    


end do
end do
end do

call rhs_exchng3_3d (flux(1,:,:,:))

!-----------------------------------------------------------------------
!j index

!Eno 1
do iRHAD = iRHAD_sta(2),iRHAD_end(2)
do jRHAD = jRHAD_sta(1),jRHAD_end(1)
do kRHAD = kRHAD_sta(2),kRHAD_end(2)
        Hf(2,iRHAD,jRHAD,kRHAD,1) = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
        
end do
end do
end do

call rhs_exchng3_3d (Hf(2,:,:,:,1))

!Eno 2
do iRHAD = iRHAD_sta(2),iRHAD_end(2)
do jRHAD = jRHAD_sta(1),jRHAD_end(2)
do kRHAD = kRHAD_sta(2),kRHAD_end(2)
        Hf(2,iRHAD,jRHAD,kRHAD,2) = one/two * de * &
                                   (phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD) - &
                                    phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD))
end do
end do
end do

call rhs_exchng3_3d (Hf(2,:,:,:,2))

!Eno 3

do iRHAD = iRHAD_sta(2), iRHAD_end(2)
do jRHAD = jRHAD_sta(1), jRHAD_end(3)
do kRHAD = kRHAD_sta(2), kRHAD_end(2)
        Hf(2,iRHAD,jRHAD,kRHAD,3) = one/three * de * &
                                    (Hf(2,iRHAD,jRHAD+1,kRHAD,2) - Hf(2,iRHAD,jRHAD,kRHAD,2))
end do
end do
end do

call rhs_exchng3_3d (Hf(2,:,:,:,3))
!calcular Roe Speed (podria ser incluido en el anterior). Almaceno en i+1/2
         
do iRHAD = iRHAD_sta(2), iRHAD_end(2)
do jRHAD = jRHAD_sta(1), jRHAD_end(2)
do kRHAD = kRHAD_sta(2), kRHAD_end(2)


        RSP(2,iRHAD,jRHAD,kRHAD) = (ucn_j(2,iRHAD,jRHAD+1,kRHAD)*phi_AD(iRHAD,jRHAD+1,kRHAD) - &
                                   ucn_j(2,iRHAD,jRHAD,kRHAD)*phi_AD(iRHAD,jRHAD,kRHAD)  ) / &
                                    (phi_AD(iRHAD,jRHAD+1,kRHAD)  -   phi_AD(iRHAD,jRHAD,kRHAD) )   

end do
end do
end do


call rhs_exchng3_3d (RSP(2,:,:,:))

if (myleft == mpi_proc_null)  then

        !ENO1
        j_idxa = jRHAD_sta(1)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
                if(RSP(2,iRHAD,j_idxa,kRHAD) .ge. zero) then
                        kmn1 = j_idxa
                else
                        kmn1 = j_idxa + 1
                end if
                
                flux(2,iRHAD,j_idxa,kRHAD) = Hf(2,iRHAD,kmn1,kRHAD,1)                        
        end do
        end do

        !ENO2
        j_idxa = jRHAD_sta(2)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
                if(RSP(2,iRHAD,j_idxa,kRHAD) .ge. zero) then
                        kmn1 = j_idxa
                else
                        kmn1 = j_idxa + 1
                end if
                
                aENO2 = Hf(2,iRHAD,kmn1,kRHAD,2)
                bENO2 = Hf(2,iRHAD,kmn1-1,kRHAD,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(2,iRHAD,j_idxa,kRHAD) = Hf(2,iRHAD,kmn1,kRHAD,1) + &
                                            cENO2* (de**(-1)*real((j_idxa - kmn1),rdf) + &
                                                   de**(-1)*real((j_idxa -kmn1 + 1),rdf))

        end do
        end do
end if

!call outputD3_real(RSP(2,:,:,:),'RSP')

!j right
if (myright == mpi_proc_null)  then

        !ENO1
        j_idxa = jRHAD_end(2)  !un nodo antes

        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
                if(RSP(2,iRHAD,j_idxa,kRHAD) .ge. zero) then
                        kmn1 = j_idxa
                else
                        kmn1 = j_idxa + 1
                end if

                flux(2,iRHAD,j_idxa,kRHAD) = Hf(2,iRHAD,kmn1,kRHAD,1)                        
        end do
        end do

        !ENO2
        j_idxa = jRHAD_end(3)
        do kRHAD = KRHAD_sta(2), kRHAD_end(2)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
                if(RSP(2,iRHAD,j_idxa,kRHAD) .ge. zero) then
                        kmn1 = j_idxa
                else
                        kmn1 = j_idxa + 1
                end if
                
                aENO2 = Hf(2,iRHAD,kmn1,kRHAD,2)
                bENO2 = Hf(2,iRHAD,kmn1-1,kRHAD,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(2,iRHAD,j_idxa,kRHAD) = Hf(2,iRHAD,kmn1,kRHAD,1) + &
                                            cENO2* (de**(-1)*real((j_idxa - kmn1),rdf) + &
                                                   de**(-1)*real((j_idxa -kmn1 + 1),rdf))

        end do
        end do
end if


!ENO3
do KRHAD = KRHAD_sta(2), KRHAD_end(2)
do jRHAD = jRHAD_sta(3), jRHAD_end(4)
do iRHAD = iRHAD_sta(2), iRHAD_end(2)
!if(myid == 1 .and.  jRHAD == 5 .and. kRHAD == 51) print *, iRHAD, myid 


        if(RSP(2,iRHAD,jRHAD,kRHAD) .ge. zero) then
                kmn1 = jRHAD                                
        else
                kmn1 = jRHAD + 1
        end if
                
        aENO2 = Hf(2,iRHAD,kmn1,kRHAD,2)
        bENO2 = Hf(2,iRHAD,kmn1-1,kRHAD,2)

        if(abs(aENO2) .ge. abs(bENO2)) then
                cENO2 = bENO2
                kmn2 = kmn1 - 1
        else
                cENO2 = aENO2
                kmn2 = kmn1
        end if

        aENO3 = Hf(2,iRHAD,kmn2,kRHAD,3)
        bENO3 = Hf(2,iRHAD,kmn2-1,kRHAD,3)

        if(abs(aENO3) .ge. abs(bENO3)) then
                cENO3 = bENO3
                kmn3 = kmn2 - 1
        else
                cENO3 = aENO3
                kmn3 = kmn2
        end if

        aux1 = de**(-1) * real(jRHAD - kmn2 +1,rdf)
        aux2 = de**(-1) * real(jRHAD - kmn2,rdf)
        aux3 = de**(-1) * real(jRHAD - kmn2 -1,rdf)        

        flux(2,iRHAD,jRHAD,kRHAD) = Hf(2,iRHAD,kmn1,kRHAD,1) + &
                                    cENO2* (de**(-1)*real((jRHAD - kmn1),rdf) + &
                                           de**(-1)*real((jRHAD -kmn1 + 1),rdf)) + &
                                    cENO3* ( (aux2 +aux1)*aux3 + aux1*aux2 )    


end do
end do
end do

call rhs_exchng3_3d (flux(2,:,:,:))


!---------------------------------------------
!k index
!Eno 1
do iRHAD = iRHAD_sta(2),iRHAD_end(2)
do jRHAD = jRHAD_sta(2),jRHAD_end(2)
do kRHAD = kRHAD_sta(1),kRHAD_end(1)
        Hf(3,iRHAD,jRHAD,kRHAD,1) = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
        
end do
end do
end do

call rhs_exchng3_3d (Hf(3,:,:,:,1))

!Eno 2
do iRHAD = iRHAD_sta(2),iRHAD_end(2)
do jRHAD = jRHAD_sta(2),jRHAD_end(2)
do kRHAD = kRHAD_sta(1),kRHAD_end(2)
        Hf(3,iRHAD,jRHAD,kRHAD,2) = one/two * dz * &
                                   (phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1) - &
                                    phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD))
end do
end do
end do

call rhs_exchng3_3d (Hf(3,:,:,:,2))

!Eno 3

do iRHAD = iRHAD_sta(2), iRHAD_end(2)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do kRHAD = kRHAD_sta(1), kRHAD_end(3)
        Hf(3,iRHAD,jRHAD,kRHAD,3) = one/three * dz * &
                                    (Hf(3,iRHAD,jRHAD,kRHAD+1,2) - Hf(3,iRHAD,jRHAD,kRHAD,2))
end do
end do
end do

call rhs_exchng3_3d (Hf(3,:,:,:,3))
!calcular Roe Speed (podria ser incluido en el anterior). Almaceno en i+1/2
         
do iRHAD = iRHAD_sta(2), iRHAD_end(2)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do kRHAD = kRHAD_sta(1), kRHAD_end(2)


        RSP(3,iRHAD,jRHAD,kRHAD) = (ucn_j(3,iRHAD,jRHAD,kRHAD+1)*phi_AD(iRHAD,jRHAD,kRHAD+1) - &
                                   ucn_j(3,iRHAD,jRHAD,kRHAD)*phi_AD(iRHAD,jRHAD,kRHAD)  ) / &
                                    (phi_AD(iRHAD,jRHAD,kRHAD+1)  -   phi_AD(iRHAD,jRHAD,kRHAD) )   

end do
end do
end do

!ghost point de ambos (creo que en estricto rigor basta solo los de HF)



call rhs_exchng3_3d (RSP(3,:,:,:))


!call outputD3_real(Hf(1,:,:,:,2),'HF2')
!ahora calculo flux
!bordes

!----------------------------------------------------------
!k down
if (mydown == mpi_proc_null)  then

        !ENO1
        k_idxa = kRHAD_sta(1)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(3,iRHAD,jRHAD,k_idxa) .ge. zero) then
                        kmn1 = k_idxa
                else
                        kmn1 = k_idxa + 1
                end if
                
                flux(3,iRHAD,jRHAD,k_idxa) = Hf(3,iRHAD,jRHAD,kmn1,1)                        
        end do
        end do

        !ENO2
        k_idxa = kRHAD_sta(2)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(3,iRHAD,jRHAD,k_idxa) .ge. zero) then
                        kmn1 = k_idxa
                else
                        kmn1 = k_idxa + 1
                end if
                
                aENO2 = Hf(3,iRHAD,jRHAD,kmn1,2)
                bENO2 = Hf(3,iRHAD,jRHAD,kmn1-1,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(3,iRHAD,jRHAD,k_idxa) = Hf(3,iRHAD,jRHAD,kmn1,1) + &
                                            cENO2* (dz**(-1)*real((k_idxa - kmn1),rdf) + &
                                                   dz**(-1)*real((k_idxa -kmn1 + 1),rdf))

        end do
        end do
end if



!k up
if (myup == mpi_proc_null)  then

        !ENO1
        k_idxa = kRHAD_end(2)  !un nodo antes

        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(3,iRHAD,jRHAD,k_idxa) .ge. zero) then
                        kmn1 = k_idxa
                else
                        kmn1 = k_idxa + 1
                end if

                flux(3,iRHAD,jRHAD,k_idxa) = Hf(3,iRHAD,jRHAD,kmn1,1)                        
        end do
        end do

        !ENO2
        k_idxa = kRHAD_end(3)
        do iRHAD = iRHAD_sta(2), iRHAD_end(2)
        do jRHAD = jRHAD_sta(2), jRHAD_end(2)
                if(RSP(3,iRHAD,jRHAD,k_idxa) .ge. zero) then
                        kmn1 = k_idxa
                else
                        kmn1 = k_idxa + 1
                end if
                
                aENO2 = Hf(3,iRHAD,jRHAD,kmn1,2)
                bENO2 = Hf(3,iRHAD,jRHAD,kmn1-1,2)

                if(abs(aENO2) .ge. abs(bENO2)) then
                        cENO2 = bENO2
                        kmn2 = kmn1 - 1
                else
                        cENO2 = aENO2
                        kmn2 = kmn1
                end if

                flux(3,iRHAD,jRHAD,k_idxa) = Hf(3,iRHAD,jRHAD,kmn1,1) + &
                                            cENO2* (dz**(-1)*real((k_idxa - kmn1),rdf) + &
                                                   dz**(-1)*real((k_idxa -kmn1 + 1),rdf))

        end do
        end do
end if

!call outputD3_real(flux(3,:,:,:),'flux3')
!call index_outputD3(100,2,2,Hf(3,:,:,:,1),debug_values(1))
!call index_outputD3(100,2,1,Hf(3,:,:,:,1),debug_values(2))

!ENO3
do KRHAD = KRHAD_sta(3), KRHAD_end(4)
do jRHAD = jRHAD_sta(2), jRHAD_end(2)
do iRHAD = iRHAD_sta(2), iRHAD_end(2)
!if(myid == 1 .and.  jRHAD == 5 .and. kRHAD == 51) print *, iRHAD, myid 


        if(RSP(3,iRHAD,jRHAD,kRHAD) .ge. zero) then
                kmn1 = kRHAD                                
        else
                kmn1 = kRHAD + 1
        end if
                
        aENO2 = Hf(3,iRHAD,jRHAD,kmn1,2)
        bENO2 = Hf(3,iRHAD,jRHAD,kmn1-1,2)

        if(abs(aENO2) .ge. abs(bENO2)) then
                cENO2 = bENO2
                kmn2 = kmn1 - 1
        else
                cENO2 = aENO2
                kmn2 = kmn1
        end if

        aENO3 = Hf(3,iRHAD,jRHAD,kmn2,3)
        bENO3 = Hf(3,iRHAD,jRHAD,kmn2-1,3)

        if(abs(aENO3) .ge. abs(bENO3)) then
                cENO3 = bENO3
                kmn3 = kmn2 - 1
        else
                cENO3 = aENO3
                kmn3 = kmn2
        end if

        aux1 = dz**(-1) * real(kRHAD - kmn2 +1,rdf)
        aux2 = dz**(-1) * real(kRHAD - kmn2,rdf)
        aux3 = dz**(-1) * real(kRHAD - kmn2 -1,rdf)        

        flux(3,iRHAD,jRHAD,kRHAD) = Hf(3,iRHAD,jRHAD,kmn1,1) + &
                                    cENO2* (dz**(-1)*real((kRHAD - kmn1),rdf) + &
                                           dz**(-1)*real((kRHAD -kmn1 + 1),rdf)) + &
                                    cENO3* ( (aux2 +aux1)*aux3 + aux1*aux2 )    


end do
end do
end do

call rhs_exchng3_3d (flux(3,:,:,:))


!call outputD3_real(flux(3,:,:,:),'flux3zx')
!call outputD3_real(flux(1,:,:,:),'flux1zx')



!ahora calculo el valor de RH
!11 de Julio 2018. Modificado para que solo calcule
!en los nodos que el flujo es calculado con orden 3

do iRHAD=iRHAD_sta(4), iRHAD_end(4)
do jRHAD=jRHAD_sta(4), jRHAD_end(4)
do kRHAD=kRHAD_sta(4), kRHAD_end(4)
        rightH_AD(iRHAD,jRHAD,kRHAD) = dc*(flux(1,iRHAD,jRHAD,kRHAD) - flux(1,iRHAD-1,jRHAD,kRHAD) )  + &
                                       de*(flux(2,iRHAD,jRHAD,kRHAD) - flux(2,iRHAD,jRHAD-1,kRHAD) )  + &
                                       dz*(flux(3,iRHAD,jRHAD,kRHAD) - flux(3,iRHAD,jRHAD,kRHAD-1) )
end do
end do
end do



!call outputD3_real(rightH_AD(:,:,:),'RH')
!-------------------------------------------------------------------------------------------------------------------------

deallocate(up,um,Hf,RSP,kmin,flux)


end subroutine calc_RH_AD_ENO3
