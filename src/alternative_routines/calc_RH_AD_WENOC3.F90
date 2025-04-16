subroutine calc_RH_AD_WENOC3(phi_AD,rightH_AD)
! Calcula lado derecho de ecuacion de adveccion usando reconstruccion WENO de 3er orden.
! La adveccion es realizada usando la ecuacion de LS en forma conservativa y la variable
! U^j/J*phi es reconstruida en la cara mediante un esquema WENO3 como el propuesto
! por Kang (2010), en donde el stencil usado se decide tras evaluar el signo de la 
! velocidad contravariante en la cara analizada.
! La velocidad contravariante en la cara se obtiene simplemente como el promedio
! de la velocidad entre los dos nodos adyacentes.
!
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
! * ver la implementacion de Kang para weno3 que la define como una funcion nomas
! weno3(double f0, double f1, double f2, double f3, double wavespeed) y por dentro
! arregla lo del stencil
!

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD


!index
integer  :: iRHAD_sta,iRHAD_end
integer  :: jRHAD_sta,jRHAD_end
integer  :: kRHAD_sta,kRHAD_end

integer :: iRHAD,jRHAD,kRHAD


! Variables WENO 3 

real (kind = rdf) :: epsWENO !epsilon positivo para mantener acotada la velocidad de Roe (Durran p251)
real (kind = rdf) :: phiL,phiC,phiR ! Valores de phi segun el stencil escogido
real (kind = rdf) :: alpha1, alpha2 ! coeficientes de ponderacion no lineales WENO
real (kind = rdf) :: phi_flux_plus_i, phi_flux_minus_i,phi_flux_plus_j, phi_flux_minus_j,phi_flux_plus_k, phi_flux_minus_k
real (kind = rdf) :: Lc,Le,Lz ! Operador espacial eq HJ en las distintas direcciones

real (kind = rdf) :: U1plus,U1minus,U2plus,U2minus,U3plus,U3minus ! Auxiliar para evaluar el signo de la velocidad en la cara



! Switches (1 o 0) para forzar stencils bounded en los bordes del dominio

integer :: bd_stencil_mb = 0 ! bounded stencil myback
integer :: bd_stencil_ml = 0
integer :: bd_stencil_md = 0

integer :: bd_stencil_mf = 0
integer :: bd_stencil_mr = 0
integer :: bd_stencil_mu = 0

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


epsWENO = 1E-6

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
    
            U1plus = one/two*(ucn_j(1,iRHAD+1,jRHAD,kRHAD)+ucn_j(1,iRHAD,jRHAD,kRHAD))

            if (U1plus.gt.zero.or.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,1).eq.1.or.&
                bd_stencil_mf.eq.1) then ! upwind
    
                phiL = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)  
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+2,jRHAD,kRHAD)*ucn_j(1,iRHAD+2,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
        
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_plus_i = (U1/J*phi)_{i+1/2,j,k}
            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
    
            !-----------------
            ! i-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
            ! si estoy en la zona 3 del obstaculo, fuerzo downwind en i.
            ! Si estoy en el borde de entrada, evito stencil upwinded
            
            U1minus = one/two*(ucn_j(1,iRHAD,jRHAD,kRHAD)+ucn_j(1,iRHAD-1,jRHAD,kRHAD))

            if (U1minus.gt.zero.and.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,3).ne.1.and.&
                bd_stencil_mb.ne.1) then ! upwind
    
                phiL = phi_AD(iRHAD-2,jRHAD,kRHAD)*ucn_j(1,iRHAD-2,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_minus_i = (U1/J*phi)_{i-1/2,j,k}
            phi_flux_minus_i = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            Lc = dc*(phi_flux_plus_i-phi_flux_minus_i)
    
            !------------------------------------------------------------------------
            ! DIRECCION ETA
            !------------------------------------------------------------------------
    
            !-----------------
            ! j+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
            ! si estoy en la zona 2, fuerzo upwind en j.
            ! Si estoy en el borde jm, fuerzo stencil centrado
    
            U2plus = one/two*(ucn_j(2,iRHAD,jRHAD+1,kRHAD)+ucn_j(2,iRHAD,jRHAD,kRHAD))

            if (U2plus.gt.zero.or.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,2).eq.1.or.&
                bd_stencil_mr.eq.1) then !upwind
    
                phiL = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD+2,kRHAD)*ucn_j(2,iRHAD,jRHAD+2,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_plus_j = (U2/J*phi)_{i,j+1/2,k}
            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            !-----------------
            ! j-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
            ! si estoy en la zona 4, fuerzo downwind en j
            ! Si estoy en j1, evito stencil upwinded
     
            U2minus = one/two*(ucn_j(2,iRHAD,jRHAD,kRHAD)+ucn_j(2,iRHAD,jRHAD-1,kRHAD))

            if (U2minus.gt.zero.and.&
                obstacle_lsm_ad(iRHAD,jRHAD,kRHAD,4).ne.1.and.&
                bd_stencil_ml.ne.1) then ! upwind
    
                phiL = phi_AD(iRHAD,jRHAD-2,kRHAD)*ucn_j(2,iRHAD,jRHAD-2,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_minus_j = (U2/J*phi)_{i,j-1/2,k}
            phi_flux_minus_j = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            Le = de*(phi_flux_plus_j-phi_flux_minus_j)
    
            !------------------------------------------------------------------------
            ! DIRECCION ZETA
            !------------------------------------------------------------------------
            !-----------------
            ! k+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
            ! Si estoy en km, fuerzo stencil centrado

            U3plus = one/two*(ucn_j(3,iRHAD,jRHAD,kRHAD+1)+ucn_j(3,iRHAD,jRHAD,kRHAD))

            if (U3plus.gt.zero.or.&
                bd_stencil_mu.eq.1) then ! upwind
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD+2)*ucn_j(3,iRHAD,jRHAD,kRHAD+2)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))

            ! phi_flux_plus_k = (U3/J*phi)_{i,j,k+1/2}
            phi_flux_plus_k = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            !-----------------
            ! k-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
            ! si estoy en k1, evito stencil upwinded

            U3minus = one/two*(ucn_j(3,iRHAD,jRHAD,kRHAD)+ucn_j(3,iRHAD,jRHAD,kRHAD-1))

            if (U3minus.gt.zero.and.&
                bd_stencil_md.ne.1) then ! upwind
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD-2)*ucn_j(3,iRHAD,jRHAD,kRHAD-2)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))

            ! phi_flux_minus_k = (U3/J*phi)_{i,j,k-1/2}
            phi_flux_minus_k = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            Lz = dz*(phi_flux_plus_k-phi_flux_minus_k)
    
    
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


    !------------------------------------------------------------------------
    ! DIRECCION CSI
    !------------------------------------------------------------------------
    
            !-----------------
            ! i+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i+1/2
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado

            U1plus = one/two*(ucn_j(1,iRHAD+1,jRHAD,kRHAD)+ucn_j(1,iRHAD,jRHAD,kRHAD))

            if (U1plus.gt.zero.or.&
                bd_stencil_mf.eq.1) then ! upwind
    
                phiL = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)  

            else ! downwind
    
                phiL = phi_AD(iRHAD+2,jRHAD,kRHAD)*ucn_j(1,iRHAD+2,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_plus_i = (U1/J*phi)_{i+1/2,j,k}
            phi_flux_plus_i = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
    
            !-----------------
            ! i-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para i-1/2
            ! Si estoy en el borde de entrada, evito stencil upwinded

            U1minus = one/two*(ucn_j(1,iRHAD,jRHAD,kRHAD)+ucn_j(1,iRHAD-1,jRHAD,kRHAD))

            if (U1minus.gt.zero.and.&
                bd_stencil_mb.ne.1) then
    
                phiL = phi_AD(iRHAD-2,jRHAD,kRHAD)*ucn_j(1,iRHAD-2,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD+1,jRHAD,kRHAD)*ucn_j(1,iRHAD+1,jRHAD,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(1,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD-1,jRHAD,kRHAD)*ucn_j(1,iRHAD-1,jRHAD,kRHAD)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_minus_i = (U1/J*phi)_{i-1/2,j,k}
            phi_flux_minus_i = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
    
    
            Lc = dc*(phi_flux_plus_i-phi_flux_minus_i)
    
    
    !------------------------------------------------------------------------
    ! DIRECCION ETA
    !------------------------------------------------------------------------
    
            !-----------------
            ! j+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j+1/2
            ! Si estoy en el borde jm, fuerzo stencil centrado

            U2plus = one/two*(ucn_j(2,iRHAD,jRHAD+1,kRHAD)+ucn_j(2,iRHAD,jRHAD,kRHAD))

            if (U2plus.gt.zero.or.&
                bd_stencil_mr.eq.1) then
    
                phiL = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD+2,kRHAD)*ucn_j(2,iRHAD,jRHAD+2,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_plus_j = (U2/J*phi)_{i,j+1/2,k}
            phi_flux_plus_j = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            !-----------------
            ! j-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para j-1/2
            ! Si estoy en j1, evito stencil upwinded

            U2minus = one/two*(ucn_j(2,iRHAD,jRHAD,kRHAD)+ucn_j(2,iRHAD,jRHAD-1,kRHAD))

            if (U2minus.gt.zero.and.&
                bd_stencil_ml.ne.1) then
    
                phiL = phi_AD(iRHAD,jRHAD-2,kRHAD)*ucn_j(2,iRHAD,jRHAD-2,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD+1,kRHAD)*ucn_j(2,iRHAD,jRHAD+1,kRHAD)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(2,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD-1,kRHAD)*ucn_j(2,iRHAD,jRHAD-1,kRHAD)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_minus_j = (U2/J*phi)_{i,j-1/2,k}
            phi_flux_minus_j = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            Le = de*(phi_flux_plus_j-phi_flux_minus_j)
    
    !------------------------------------------------------------------------
    ! DIRECCION ZETA
    !------------------------------------------------------------------------
            !-----------------
            ! k+1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k+1/2
            ! Si estoy en km, fuerzo stencil centrado

            U3plus = one/two*(ucn_j(3,iRHAD,jRHAD,kRHAD+1)+ucn_j(3,iRHAD,jRHAD,kRHAD))

            if (U3plus.gt.zero.or.&
                bd_stencil_mu.eq.1) then
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
    
            else ! downwind
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD+2)*ucn_j(3,iRHAD,jRHAD,kRHAD+2)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
    
            end if
    
            ! Calculo los coeficientes de ponderacion no lineal
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_plus_k = (U3/J*phi)_{i,j,k+1/2}
            phi_flux_plus_k = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            !-----------------
            ! k-1/2
            !-----------------
            ! Escojo el stencil dependiendo de la velocidad local para k-1/2
            ! Si estoy en k1, evito stencil upwinded

            U3minus = one/two*(ucn_j(3,iRHAD,jRHAD,kRHAD)+ucn_j(3,iRHAD,jRHAD,kRHAD-1))

            if (U3minus.gt.zero.and.&
                bd_stencil_md.ne.1) then
    
                phiL = phi_AD(iRHAD,jRHAD,kRHAD-2)*ucn_j(3,iRHAD,jRHAD,kRHAD-2)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
    
            else ! downwind
         
                phiL = phi_AD(iRHAD,jRHAD,kRHAD+1)*ucn_j(3,iRHAD,jRHAD,kRHAD+1)
                phiC = phi_AD(iRHAD,jRHAD,kRHAD)*ucn_j(3,iRHAD,jRHAD,kRHAD)
                phiR = phi_AD(iRHAD,jRHAD,kRHAD-1)*ucn_j(3,iRHAD,jRHAD,kRHAD-1)
    
            end if
    
    
            alpha1 = two/three * (1/(((phiC-phiR)**2+epsWENO)**2))
            alpha2 = one/three * (1/(((phiL-phiC)**2+epsWENO)**2))
    
            ! phi_flux_minus_k = (U3/J*phi)_{i,j,k-1/2}
            phi_flux_minus_k = (alpha1/(alpha1+alpha2))*(phiC/two + phiR/two) + &
                            (alpha2/(alpha1+alpha2))*(-phiL/two + three*phiC/two)
    
            Lz = dz*(phi_flux_plus_k-phi_flux_minus_k)
    
    
            rightH_AD(iRHAD,jRHAD,kRHAD) = Lc+Le+Lz
    
    end do
    end do
    end do

end if ! if isobstacle==1

end subroutine calc_RH_AD_WENOC3
