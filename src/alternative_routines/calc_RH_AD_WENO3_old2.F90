subroutine calc_RH_AD_WENO3(phi_AD,rightH_AD)
! Calcula lado derecho de ecuacion de adveccion usando reconstruccion WENO de 3er orden
! La adveccion es realizada usando la ecuacion de LS no conservativa en coordenadas
! curvilineas planteada en Kang & Sotiropoulos, 2010
!--------------------------------------------------------------------------------------
! TO DO:
! 
! * Quizas reemplazar los 0 y 1 para detectar si estoy en un borde y ajustar el stencil
!   por algo mas eficiente
! 
! * Analizar si esta implementacion es adecuada para trabajar con multiples obstaculos
!
! * Analizar si esta implementacion es adecuada para trabajar con mallas traslapadas
!
! * Verificar si los valores reales usados estan en precision rdf


implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD


!index
integer  :: iRHAD_sta,iRHAD_end
integer  :: jRHAD_sta,jRHAD_end
integer  :: kRHAD_sta,kRHAD_end

integer :: iRHAD,jRHAD,kRHAD
integer :: i,j,k


! Variables WENO 3 

real (kind = rdf) :: epsWENO !epsilon positivo para mantener acotada la velocidad de Roe (Durran p251)
real (kind = rdf) :: phiL,phiC,phiR ! Valores de phi segun el stencil escogido
real (kind = rdf) :: alpha1, alpha2 ! coeficientes de ponderacion no lineales WENO
real (kind = rdf) :: phi_flux_plus_i, phi_flux_minus_i,phi_flux_plus_j, phi_flux_minus_j,phi_flux_plus_k, phi_flux_minus_k
real (kind = rdf) :: Lc,Le,Lz ! Operador espacial eq HJ en las distintas direcciones
real (kind = rdf) :: dc2,de2,dz2 ! ( 1/(2Δξ) , 1/(2Δη) , 1/(2Δζ) )

! Switches (1 o 0) para forzar stencils bounded en los bordes del dominio. El codigo los determina
! automaticamente evaluando la posicion del procesador en el dominio computacional

integer :: bd_stencil_mb = 0 ! bounded stencil myback
integer :: bd_stencil_ml = 0
integer :: bd_stencil_md = 0

integer :: bd_stencil_mf = 0
integer :: bd_stencil_mr = 0
integer :: bd_stencil_mu = 0

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

dc2 = one_half * dc
de2 = one_half * de
dz2 = one_half * dz

epsWENO = ten_eminus_six ! 1x10^(-6) in rdf precision

!Nodos incluyendo el borde
iRHAD_sta = il + igp
jRHAD_sta = jl + jgp
kRHAD_sta = kl + kgp

iRHAD_end = iu - igp        
jRHAD_end = ju - jgp
kRHAD_end = ku - kgp

!Nodos sin incluir borde (el borde se actualiza en bcond_lsm.F90)

if (myback == mpi_proc_null)  iRHAD_sta = il + igp + 1
if (myleft == mpi_proc_null)  jRHAD_sta = jl + jgp + 1
if (mydown == mpi_proc_null)  kRHAD_sta = kl + kgp + 1

if (myfront == mpi_proc_null) iRHAD_end = iu - igp - 1
if (myright == mpi_proc_null) jRHAD_end = ju - jgp - 1
if (myup    == mpi_proc_null) kRHAD_end = ku - kgp - 1

if (is_obstacle==1) then

    do iRHAD = iRHAD_sta,iRHAD_end
    do jRHAD = jRHAD_sta,jRHAD_end
    do kRHAD = kRHAD_sta,kRHAD_end


        ! Definicion de los booleans para stencils bounded en los bordes

        if(myback.eq.mpi_proc_null.and.iRHAD.eq.iRHAD_sta) then
            bd_stencil_mb = 1
        else
            bd_stencil_mb = 0
        end if


        if(myleft.eq.mpi_proc_null.and.jRHAD.eq.jRHAD_sta) then
            bd_stencil_ml = 1
        else
            bd_stencil_ml = 0
        end if

        if(mydown.eq.mpi_proc_null.and.kRHAD.eq.kRHAD_sta) then
            bd_stencil_md = 1
        else
            bd_stencil_md = 0
        end if


        if(myfront.eq.mpi_proc_null.and.iRHAD.eq.iRHAD_end) then
            bd_stencil_mf = 1
        else
            bd_stencil_mf = 0
        end if

        if(myright.eq.mpi_proc_null.and.jRHAD.eq.jRHAD_end) then
            bd_stencil_mr = 1
        else
            bd_stencil_mr = 0
        end if

        if(myup.eq.mpi_proc_null.and.kRHAD.eq.kRHAD_end) then
            bd_stencil_mu = 1
        else
            bd_stencil_mu = 0
        end if
    
        ! if de si estoy dentro (incluyendo) del obstaculo
        if(iRHAD.ge.li_obs_ia.and.iRHAD.le.li_obs_ib.and.&
           jRHAD.ge.li_obs_ja.and.jRHAD.le.li_obs_jb.and.&
           kRHAD.ge.li_obs_ka.and.kRHAD.le.li_obs_kb) then
            
            ! Si estoy dentro del obstaculo el RHS de la ecuacion de 
            ! adveccion del LS es cero y por ende no hay evolucion
            ! temporal de phi dentro del obstaculo
    
            rightH_AD(iRHAD,jRHAD,kRHAD)=0.0
    
    
        else ! si estoy fuera del obstaculo

            !------------------------------------------------------------------------
            ! DIRECCION CSI
            !------------------------------------------------------------------------
    
            !-----------------
            ! i+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
            ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,1).eq.1.or.&
                bd_stencil_mf.eq.1) then ! upwind
    
                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_{i+1/2,j,k}
            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
    
            !-----------------
            ! i-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
            ! si estoy en la zona 3 del obstaculo, fuerzo downwind en i.
            ! Si estoy en el borde de entrada, evito stencil upwinded
            
            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,3).ne.1.and.&
                bd_stencil_mb.ne.1) then ! upwind
    
                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
    
            end if
    
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_{i-1/2,j,k}
            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
    
            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
    
    
            !------------------------------------------------------------------------
            ! DIRECCION ETA
            !------------------------------------------------------------------------
    
            !-----------------
            ! j+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
            ! si estoy en la zona 2, fuerzo upwind en j.
            ! Si estoy en el borde jm, fuerzo stencil centrado
    
            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,2).eq.1.or.&
                bd_stencil_mr.eq.1) then !upwind
    
                phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
                phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
                phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
                phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
                phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
            !-----------------
            ! j-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
            ! si estoy en la zona 4, fuerzo downwind en j
            ! Si estoy en j1, evito stencil upwinded
    
            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,4).ne.1.and.&
                bd_stencil_ml.ne.1) then ! upwind
    
                phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
                phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
                phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
    
            else ! downwind
    
                phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
                phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
                phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
    
            end if
    
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
            Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
    
            !------------------------------------------------------------------------
            ! DIRECCION ZETA
            !------------------------------------------------------------------------
            !-----------------
            ! k+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
            ! Si estoy en km, fuerzo stencil centrado

            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                bd_stencil_mu.eq.1) then ! upwind
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
    
            else ! downwind
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                              ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
            !-----------------
            ! k-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
            ! si estoy en k1, evito stencil upwinded

            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                bd_stencil_md.ne.1) then ! upwind
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
    
            else ! downwind
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
    
            end if
    
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    
            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
    

            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
    
        end if
    
    end do
    end do
    end do


else ! procesador sin parte del obstaculo

    do iRHAD = iRHAD_sta,iRHAD_end
    do jRHAD = jRHAD_sta,jRHAD_end
    do kRHAD = kRHAD_sta,kRHAD_end
    
        ! Definicion de los booleans para stencils bounded en los bordes

        if(myback.eq.mpi_proc_null.and.iRHAD.eq.iRHAD_sta) then
            bd_stencil_mb = 1
        else
            bd_stencil_mb = 0
        end if


        if(myleft.eq.mpi_proc_null.and.jRHAD.eq.jRHAD_sta) then
            bd_stencil_ml = 1
        else
            bd_stencil_ml = 0
        end if

        if(mydown.eq.mpi_proc_null.and.kRHAD.eq.kRHAD_sta) then
            bd_stencil_md = 1
        else
            bd_stencil_md = 0
        end if


        if(myfront.eq.mpi_proc_null.and.iRHAD.eq.iRHAD_end) then
            bd_stencil_mf = 1
        else
            bd_stencil_mf = 0
        end if

        if(myright.eq.mpi_proc_null.and.jRHAD.eq.jRHAD_end) then
            bd_stencil_mr = 1
        else
            bd_stencil_mr = 0
        end if

        if(myup.eq.mpi_proc_null.and.kRHAD.eq.kRHAD_end) then
            bd_stencil_mu = 1
        else
            bd_stencil_mu = 0
        end if


        ! I only advect the RHS of the level set function if I'm in the water
        ! phase (check if this is valid or not btw. For the moment, it's just
        ! to avoid some errors)

        !if(rsign(iRHAD,jRHAD,kRHAD) > one_half) then

        !------------------------------------------------------------------------
        ! DIRECCION CSI
        !------------------------------------------------------------------------
    
            !-----------------
            ! i+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado

            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                bd_stencil_mf.eq.1) then ! upwind
    
                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
    
            end if
    
            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_{i+1/2,j,k}
            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )

            !-----------------
            ! i-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
            ! Si estoy en el borde de entrada, evito stencil upwinded

            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                bd_stencil_mb.ne.1) then
    
                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
    
            else
    
                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
    
            end if
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_{i-1/2,j,k}
            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )


            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
    
    
        !------------------------------------------------------------------------
        ! DIRECCION ETA
        !------------------------------------------------------------------------
    
            !-----------------
            ! j+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
            ! Si estoy en el borde jm, fuerzo stencil centrado


            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                bd_stencil_mr.eq.1) then
    
                phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
                phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
                phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
    
            else
    
                phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
                phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
                phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )

            !-----------------
            ! j-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
            ! Si estoy en j1, evito stencil upwinded

            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                bd_stencil_ml.ne.1) then
    
                phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
                phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
                phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
    
            else
    
                phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
                phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
                phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
    
            end if
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )


            Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
    
        !------------------------------------------------------------------------
        ! DIRECCION ZETA
        !------------------------------------------------------------------------
            !-----------------
            ! k+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
            ! Si estoy en km, fuerzo stencil centrado

            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
                bd_stencil_mu.eq.1) then
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
    
            else
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_{i,j,k+1/2}
            phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                              ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )

            !-----------------
            ! k-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
            ! Si estoy en k1, evito stencil upwinded

            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
                bd_stencil_md.ne.1) then
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
    
            else
    
                phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
                phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
    
            end if
    
            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
    
    
            phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
    

            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
    
    
            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz

        !end if
    end do
    end do
    end do

!    ! - - - - - - - - - - - - - - - - - - - - - - - - 
!    ! BOUNDARY CONDITIONS
!    ! - - - - - - - - - - - - - - - - - - - - - - - - 
!
!    ! Solving the advection equation at the boundaries using a mixed 2nd-order UPWIND/CENTERED
!    ! scheme. This is a partial implementation (JS, 15 Feb 2023)
!    
!    if (myback == mpi_proc_null) then
!        
!        iRHAD_sta = il + igp
!
!        if (myleft == mpi_proc_null)  jRHAD_sta = jl + jgp + 1
!        if (mydown == mpi_proc_null)  kRHAD_sta = kl + kgp + 1
!
!        if (myfront == mpi_proc_null) iRHAD_end = iu - igp - 1
!        if (myright == mpi_proc_null) jRHAD_end = ju - jgp - 1
!        if (myup    == mpi_proc_null) kRHAD_end = ku - kgp - 1
!
!
!        do kRHAD=kRHAD_sta,kRHAD_end
!        do jRHAD=jRHAD_sta,jRHAD_end
!           
!           iRHAD=iRHAD_sta
!
!           ! UPWIND BIASED
!           Lc = ucn_j(1,iRHAD,jRHAD,kRHAD) * ( dc2 * ( - three * phi_AD( iRHAD  , jRHAD , kRHAD )  &
!                                                       + four  * phi_AD( iRHAD+1, jRHAD , kRHAD )  &
!                                                       - one   * phi_AD( iRHAD+2, jRHAD , kRHAD )    ) )
!
!            !------------------------------------------------------------------------
!            ! DIRECCION ETA
!            !------------------------------------------------------------------------
!        
!                !-----------------
!                ! j+1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para j+1/2
!                ! Si estoy en el borde jm, fuerzo stencil centrado
!    
!    
!                if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                    bd_stencil_mr.eq.1) then
!        
!                    phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
!                    phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
!                    phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!        
!                else
!        
!                    phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
!                    phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!                    phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
!        
!                end if
!        
!                ! Calculo los coeficientes de ponderacion no lineal
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                  (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!                !-----------------
!                ! j-1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para j-1/2
!                ! Si estoy en j1, evito stencil upwinded
!    
!                if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                    bd_stencil_ml.ne.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
!                    phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!                    phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
!                    phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
!                    phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!        
!                end if
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                   ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!    
!                Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
!        
!            !------------------------------------------------------------------------
!            ! DIRECCION ZETA
!            !------------------------------------------------------------------------
!                !-----------------
!                ! k+1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para k+1/2
!                ! Si estoy en km, fuerzo stencil centrado
!    
!                if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                    bd_stencil_mu.eq.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!                end if
!        
!                ! Calculo los coeficientes de ponderacion no lineal
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!                ! phi_{i,j,k+1/2}
!                phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                  ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!                !-----------------
!                ! k-1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para k-1/2
!                ! Si estoy en k1, evito stencil upwinded
!    
!                if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                    bd_stencil_md.ne.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!        
!                end if
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                   ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!        
!    
!                Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
!        
!        
!                rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (myback == mpi_proc_null)
!
!    if (myfront == mpi_proc_null) then
!        
!        iRHAD_end = iu - igp
!        
!        if (myback == mpi_proc_null)  iRHAD_sta = il + igp + 1
!        if (myleft == mpi_proc_null)  jRHAD_sta = jl + jgp + 1
!        if (mydown == mpi_proc_null)  kRHAD_sta = kl + kgp + 1
!
!        if (myright == mpi_proc_null) jRHAD_end = ju - jgp - 1
!        if (myup    == mpi_proc_null) kRHAD_end = ku - kgp - 1
!
!
!        do kRHAD=kRHAD_sta,kRHAD_end
!        do jRHAD=jRHAD_sta,jRHAD_end
!           
!           iRHAD=iRHAD_end
!
!           ! UPWIND BIASED
!           Lc = ucn_j(1,iRHAD,jRHAD,kRHAD) * ( dc2 * ( + three * phi_AD( iRHAD  , jRHAD , kRHAD )  &
!                                                       - four  * phi_AD( iRHAD-1, jRHAD , kRHAD )  &
!                                                       + one   * phi_AD( iRHAD-2, jRHAD , kRHAD )    ) )
!
!            !------------------------------------------------------------------------
!            ! DIRECCION ETA
!            !------------------------------------------------------------------------
!        
!                !-----------------
!                ! j+1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para j+1/2
!                ! Si estoy en el borde jm, fuerzo stencil centrado
!    
!    
!                if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                    bd_stencil_mr.eq.1) then
!        
!                    phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
!                    phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
!                    phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!        
!                else
!        
!                    phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
!                    phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!                    phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
!        
!                end if
!        
!                ! Calculo los coeficientes de ponderacion no lineal
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                  (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!                !-----------------
!                ! j-1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para j-1/2
!                ! Si estoy en j1, evito stencil upwinded
!    
!                if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                    bd_stencil_ml.ne.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
!                    phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!                    phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
!                    phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
!                    phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!        
!                end if
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                   ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!    
!                Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
!        
!            !------------------------------------------------------------------------
!            ! DIRECCION ZETA
!            !------------------------------------------------------------------------
!                !-----------------
!                ! k+1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para k+1/2
!                ! Si estoy en km, fuerzo stencil centrado
!    
!                if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                    bd_stencil_mu.eq.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!                end if
!        
!                ! Calculo los coeficientes de ponderacion no lineal
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!                ! phi_{i,j,k+1/2}
!                phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                  ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!                !-----------------
!                ! k-1/2
!                !-----------------
!                ! Escojo el stencil dependiendo de la velocidad local para k-1/2
!                ! Si estoy en k1, evito stencil upwinded
!    
!                if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                    bd_stencil_md.ne.1) then
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!                else
!        
!                    phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                    phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                    phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!        
!                end if
!        
!                alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!                alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!                phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                                   ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!        
!    
!                Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
!        
!        
!                rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (myfront == mpi_proc_null)
!
!
!    if (myleft == mpi_proc_null) then
!        
!        jRHAD_sta = jl + jgp
!        
!        if (myback  == mpi_proc_null)  iRHAD_sta = il + igp + 1
!        if (mydown  == mpi_proc_null)  kRHAD_sta = kl + kgp + 1
!
!        if (myfront == mpi_proc_null)  iRHAD_end = iu - igp - 1
!        if (myright == mpi_proc_null)  jRHAD_end = ju - jgp - 1
!        if (myup    == mpi_proc_null)  kRHAD_end = ku - kgp - 1
!
!
!        do kRHAD=kRHAD_sta,kRHAD_end
!        do iRHAD=iRHAD_sta,iRHAD_end
!           
!            jRHAD=jRHAD_sta
!
!            ! UPWIND BIASED
!            Le = ucn_j(2,iRHAD,jRHAD,kRHAD) * ( de2 * ( - three * phi_AD( iRHAD , jRHAD   , kRHAD )  &
!                                                        + four  * phi_AD( iRHAD , jRHAD+1 , kRHAD )  &
!                                                        - one   * phi_AD( iRHAD , jRHAD+2 , kRHAD )    ) )
!
!            !------------------------------------------------------------------------
!            ! DIRECCION CSI
!            !------------------------------------------------------------------------
!    
!            !-----------------
!            ! i+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
!            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mf.eq.1) then ! upwind
!    
!                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
!    
!            else ! downwind
!    
!                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i+1/2,j,k}
!            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! i-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
!            ! Si estoy en el borde de entrada, evito stencil upwinded
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_mb.ne.1) then
!    
!                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i-1/2,j,k}
!            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
!        
!            !------------------------------------------------------------------------
!            ! DIRECCION ZETA
!            !------------------------------------------------------------------------
!            !-----------------
!            ! k+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
!            ! Si estoy en km, fuerzo stencil centrado
!    
!            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mu.eq.1) then
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!        
!            else
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!            end if
!        
!            ! Calculo los coeficientes de ponderacion no lineal
!        
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!            ! phi_{i,j,k+1/2}
!            phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                              ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!            !-----------------
!            ! k-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
!            ! Si estoy en k1, evito stencil upwinded
!    
!            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_md.ne.1) then
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!            else
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!        
!            end if
!        
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!            phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!        
!    
!            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
!        
!        
!            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (myleft == mpi_proc_null)
!
!    if (myright == mpi_proc_null) then
!        
!        jRHAD_end = ju - jgp
!        
!        if (myback  == mpi_proc_null) iRHAD_sta = il + igp + 1
!        if (myleft  == mpi_proc_null) jRHAD_sta = jl + jgp + 1
!        if (mydown  == mpi_proc_null) kRHAD_sta = kl + kgp + 1
!
!        if (myfront == mpi_proc_null) iRHAD_end = iu - igp - 1
!        if (myup    == mpi_proc_null) kRHAD_end = ku - kgp - 1
!
!
!        do kRHAD=kRHAD_sta,kRHAD_end
!        do iRHAD=iRHAD_sta,iRHAD_end
!           
!            jRHAD=jRHAD_end
!
!            ! UPWIND BIASED
!            Le = ucn_j(2,iRHAD,jRHAD,kRHAD) * ( de2 * ( + three * phi_AD( iRHAD , jRHAD   , kRHAD )  &
!                                                        - four  * phi_AD( iRHAD , jRHAD-1 , kRHAD )  &
!                                                        + one   * phi_AD( iRHAD , jRHAD-2 , kRHAD )    ) )
!
!
!            !------------------------------------------------------------------------
!            ! DIRECCION CSI
!            !------------------------------------------------------------------------
!    
!            !-----------------
!            ! i+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
!            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mf.eq.1) then ! upwind
!    
!                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
!    
!            else ! downwind
!    
!                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i+1/2,j,k}
!            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! i-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
!            ! Si estoy en el borde de entrada, evito stencil upwinded
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_mb.ne.1) then
!    
!                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i-1/2,j,k}
!            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
!        
!            !------------------------------------------------------------------------
!            ! DIRECCION ZETA
!            !------------------------------------------------------------------------
!            !-----------------
!            ! k+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
!            ! Si estoy en km, fuerzo stencil centrado
!    
!            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mu.eq.1) then
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!        
!            else
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD+2 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!            end if
!        
!            ! Calculo los coeficientes de ponderacion no lineal
!        
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!            ! phi_{i,j,k+1/2}
!            phi_flux_plus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                              ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!    
!            !-----------------
!            ! k-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
!            ! Si estoy en k1, evito stencil upwinded
!    
!            if (ucn_j(3,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_md.ne.1) then
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD-2 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD   )
!        
!            else
!        
!                phiL = phi_AD( iRHAD , jRHAD , kRHAD+1 )
!                phiC = phi_AD( iRHAD , jRHAD , kRHAD   )
!                phiR = phi_AD( iRHAD , jRHAD , kRHAD-1 )
!        
!            end if
!        
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!        
!        
!            phi_flux_minus_k = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!        
!    
!            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD)*(dz*(phi_flux_plus_k-phi_flux_minus_k))
!        
!        
!            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (myright == mpi_proc_null)
!
!    if (mydown == mpi_proc_null) then
!        
!        kRHAD_sta = kl + kgp
!        
!        if (myback  == mpi_proc_null)  iRHAD_sta = il + igp + 1
!        if (myleft  == mpi_proc_null)  jRHAD_sta = jl + jgp + 1
!
!        if (myfront == mpi_proc_null)  iRHAD_end = iu - igp - 1
!        if (myright == mpi_proc_null)  jRHAD_end = ju - jgp - 1
!        if (myup    == mpi_proc_null)  kRHAD_end = ku - kgp - 1
!
!
!        do jRHAD=jRHAD_sta,jRHAD_end
!        do iRHAD=iRHAD_sta,iRHAD_end
!           
!            kRHAD=kRHAD_sta
!
!            ! UPWIND BIASED
!            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD) * ( dz2 * ( - three * phi_AD( iRHAD , jRHAD , kRHAD   )  &
!                                                        + four  * phi_AD( iRHAD , jRHAD , kRHAD+1 )  &
!                                                        - one   * phi_AD( iRHAD , jRHAD , kRHAD+2 )    ) )
!
!        !------------------------------------------------------------------------
!        ! DIRECCION CSI
!        !------------------------------------------------------------------------
!    
!            !-----------------
!            ! i+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
!            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mf.eq.1) then ! upwind
!    
!                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
!    
!            else ! downwind
!    
!                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i+1/2,j,k}
!            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! i-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
!            ! Si estoy en el borde de entrada, evito stencil upwinded
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_mb.ne.1) then
!    
!                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i-1/2,j,k}
!            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
!    
!    
!        !------------------------------------------------------------------------
!        ! DIRECCION ETA
!        !------------------------------------------------------------------------
!    
!            !-----------------
!            ! j+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
!            ! Si estoy en el borde jm, fuerzo stencil centrado
!
!
!            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mr.eq.1) then
!    
!                phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
!                phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
!                phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
!                phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!                phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
!    
!            end if
!    
!            ! Calculo los coeficientes de ponderacion no lineal
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!    
!            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! j-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
!            ! Si estoy en j1, evito stencil upwinded
!
!            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_ml.ne.1) then
!    
!                phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
!                phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!                phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
!    
!            else
!    
!                phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
!                phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
!                phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!    
!            phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
!        
!        
!            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (mydown == mpi_proc_null)
!
!    if (myup == mpi_proc_null) then
!        
!        kRHAD_end = ku - kgp
!        
!        if (myback  == mpi_proc_null) iRHAD_sta = il + igp + 1
!        if (myleft  == mpi_proc_null) jRHAD_sta = jl + jgp + 1
!        if (mydown  == mpi_proc_null) kRHAD_sta = kl + kgp + 1
!
!        if (myfront == mpi_proc_null) iRHAD_end = iu - igp - 1
!        if (myright == mpi_proc_null) jRHAD_end = ju - jgp - 1
!
!
!        do jRHAD=jRHAD_sta,jRHAD_end
!        do iRHAD=iRHAD_sta,iRHAD_end
!           
!            kRHAD=kRHAD_end
!
!            ! UPWIND BIASED
!            Lz = ucn_j(3,iRHAD,jRHAD,kRHAD) * ( dz2 * ( + three * phi_AD( iRHAD , jRHAD , kRHAD   )  &
!                                                        - four  * phi_AD( iRHAD , jRHAD , kRHAD-1 )  &
!                                                        + one   * phi_AD( iRHAD , jRHAD , kRHAD-2 )    ) )
!
!
!        !------------------------------------------------------------------------
!        ! DIRECCION CSI
!        !------------------------------------------------------------------------
!    
!            !-----------------
!            ! i+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
!            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mf.eq.1) then ! upwind
!    
!                phiL = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD+1 , jRHAD , kRHAD)  
!    
!            else ! downwind
!    
!                phiL = phi_AD(iRHAD+2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = (two/three) * (one / (((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = (one/three) * (one / (((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i+1/2,j,k}
!            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*( phiC/two + phiR/two       ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! i-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
!            ! Si estoy en el borde de entrada, evito stencil upwinded
!
!            if (ucn_j(1,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_mb.ne.1) then
!    
!                phiL = phi_AD(iRHAD-2 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD   , jRHAD , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD+1 , jRHAD , kRHAD)
!                phiC = phi_AD(iRHAD   , jRHAD , kRHAD)
!                phiR = phi_AD(iRHAD-1 , jRHAD , kRHAD)
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!            ! phi_{i-1/2,j,k}
!            phi_flux_minus_i = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Lc = ucn_j(1,iRHAD,jRHAD,kRHAD)*(dc*(phi_flux_plus_i-phi_flux_minus_i))
!    
!    
!        !------------------------------------------------------------------------
!        ! DIRECCION ETA
!        !------------------------------------------------------------------------
!    
!            !-----------------
!            ! j+1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
!            ! Si estoy en el borde jm, fuerzo stencil centrado
!
!
!            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.or.&
!                bd_stencil_mr.eq.1) then
!    
!                phiL = phi_AD(iRHAD , jRHAD-1 , kRHAD)
!                phiC = phi_AD(iRHAD , jRHAD   , kRHAD)
!                phiR = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!    
!            else
!    
!                phiL = phi_AD(iRHAD , jRHAD+2 , kRHAD)
!                phiC = phi_AD(iRHAD , jRHAD+1 , kRHAD)
!                phiR = phi_AD(iRHAD , jRHAD   , kRHAD)
!    
!            end if
!    
!            ! Calculo los coeficientes de ponderacion no lineal
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!    
!            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                              (alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!            !-----------------
!            ! j-1/2
!            !-----------------
!            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
!            ! Si estoy en j1, evito stencil upwinded
!
!            if (ucn_j(2,iRHAD,jRHAD,kRHAD).gt.zero.and.&
!                bd_stencil_ml.ne.1) then
!    
!                phiL = phi_AD( iRHAD , jRHAD-2 , kRHAD )
!                phiC = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!                phiR = phi_AD( iRHAD , jRHAD   , kRHAD )
!    
!            else
!    
!                phiL = phi_AD( iRHAD , jRHAD+1 , kRHAD )
!                phiC = phi_AD( iRHAD , jRHAD   , kRHAD )
!                phiR = phi_AD( iRHAD , jRHAD-1 , kRHAD )
!    
!            end if
!    
!            alpha1 = two/three * (one/(((phiC-phiR)**2+epsWENO)**2))
!            alpha2 = one/three * (one/(((phiL-phiC)**2+epsWENO)**2))
!    
!    
!            phi_flux_minus_j = ( alpha1/(alpha1+alpha2))*( phiC/two + phiR      /two ) + &
!                               ( alpha2/(alpha1+alpha2))*(-phiL/two + phiC*three/two )
!
!
!            Le = ucn_j(2,iRHAD,jRHAD,kRHAD)*(de*(phi_flux_plus_j-phi_flux_minus_j))
!        
!        
!            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
!
!        end do
!        end do
!
!    end if ! (myup == mpi_proc_null)

end if ! if isobstacle==1

end subroutine calc_RH_AD_WENO3
