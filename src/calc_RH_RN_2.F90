subroutine calc_RH_RN_2(sgndf_RN,signo_RN,rightH_RN)

!----------------------------------------------------------------------------------------
! Calcula lado derecho de ecuacion de Reinicializacion usando un esquema espacial ENO
! de segundo orden. Para ello calcula el gradiente usando la transformacion a coordenadas
! curvilineas planteada en Kang & Sotiropoulos, 2010
!----------------------------------------------------------------------------------------
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

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: sgndf_RN
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: signo_RN
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_RN
real (kind = rdf) :: Mphi, signPhi


!index
integer  :: iRN_sta,iRN_end
integer  :: jRN_sta,jRN_end
integer  :: kRN_sta,kRN_end

integer :: iRN,jRN,kRN


! Variables ENO2

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) :: dphi_dcsi, dphi_deta, dphi_dzet
real (kind = rdf) :: dphi_dcsi_plus,dphi_dcsi_minus, dphi_deta_plus, dphi_deta_minus, dphi_dzet_plus, dphi_dzet_minus  
real (kind = rdf) :: phiLL,phiL,phiC,phiR,phiRR 
real (kind = rdf) :: dphi_dx,dphi_dy,dphi_dz,norma_grad_phi 

! Switches (1 o 0) para forzar stencils bounded en los bordes del dominio

integer :: bd_stencil_mb = 0 ! bounded stencil myback
integer :: bd_stencil_ml = 0
integer :: bd_stencil_md = 0

integer :: bd_stencil_mf = 0
integer :: bd_stencil_mr = 0
integer :: bd_stencil_mu = 0

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!Nodos incluyendo el borde
iRN_sta = il + igp
jRN_sta = jl + jgp
kRN_sta = kl + kgp

iRN_end = iu - igp        
jRN_end = ju - jgp
kRN_end = ku - kgp

!Nodos sin incluir borde (el borde se actualiza en bcond_lsm.F90)

if (myback == mpi_proc_null)  iRN_sta = il + igp + 1
if (myleft == mpi_proc_null)  jRN_sta = jl + jgp + 1
if (mydown == mpi_proc_null)  kRN_sta = kl + kgp + 1

if (myfront == mpi_proc_null) iRN_end = iu - igp - 1
if (myright == mpi_proc_null) jRN_end = ju - jgp - 1
if (myup    == mpi_proc_null) kRN_end = ku - kgp - 1

if (is_obstacle==1) then

    do iRN = iRN_sta,iRN_end
    do jRN = jRN_sta,jRN_end
    do kRN = kRN_sta,kRN_end


        ! Definicion de los booleans para stencils bounded en los bordes

        if(myback.eq.mpi_proc_null.and.iRN.eq.iRN_sta) then
            bd_stencil_mb = 1
        else
            bd_stencil_mb = 0
        end if


        if(myleft.eq.mpi_proc_null.and.jRN.eq.jRN_sta) then
            bd_stencil_ml = 1
        else
            bd_stencil_ml = 0
        end if

        if(mydown.eq.mpi_proc_null.and.kRN.eq.kRN_sta) then
            bd_stencil_md = 1
        else
            bd_stencil_md = 0
        end if


        if(myfront.eq.mpi_proc_null.and.iRN.eq.iRN_end) then
            bd_stencil_mf = 1
        else
            bd_stencil_mf = 0
        end if

        if(myright.eq.mpi_proc_null.and.jRN.eq.jRN_end) then
            bd_stencil_mr = 1
        else
            bd_stencil_mr = 0
        end if

        if(myup.eq.mpi_proc_null.and.kRN.eq.kRN_end) then
            bd_stencil_mu = 1
        else
            bd_stencil_mu = 0
        end if
    
        ! if de si estoy dentro (incluyendo) del obstaculo
        if(iRN.ge.li_obs_ia.and.iRN.le.li_obs_ib.and.&
           jRN.ge.li_obs_ja.and.jRN.le.li_obs_jb.and.&
           kRN.ge.li_obs_ka.and.kRN.le.li_obs_kb) then
            
            ! Si estoy dentro del obstaculo el RHS de la ecuacion de 
            ! es cero y por ende, el valor de phi, no se actualiza por
            ! reinicializacion dentro del obstaculo
    
            rightH_RN(iRN,jRN,kRN)=0.0
    
    
        else ! si estoy fuera del obstaculo

            !------------------------------------------------------------------------
            ! DIRECCION CSI
            !------------------------------------------------------------------------

            ! Escojo el stencil para calcular el gradiente 
            ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
            if (obstacle_lsm_ad(iRN,jRN,kRN,1).eq.1.or.&
                bd_stencil_mf.eq.1) then
    
                phiLL = sgndf_RN(iRN-2,jRN,kRN)
                phiL  = sgndf_RN(iRN-1,jRN,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN+1,jRN,kRN)
                phiRR = 0.0 ! outbounded

                dphi_dcsi_plus  = phiR-phiC-one_half*(phiR-2*phiC+phiL)
                dphi_dcsi_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            else if(obstacle_lsm_ad(iRN,jRN,kRN,3).eq.1.or.&
                bd_stencil_mb.eq.1) then

                phiLL = 0.0 ! outbounded
                phiL  = sgndf_RN(iRN-1,jRN,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN+1,jRN,kRN)
                phiRR = sgndf_RN(iRN+2,jRN,kRN)

                dphi_dcsi_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_dcsi_minus = phiC-phiL+one_half*(phiR-2*phiC+phiL)

            else ! no tengo restricciones de stencil
                
                phiLL = sgndf_RN(iRN-2,jRN,kRN)
                phiL  = sgndf_RN(iRN-1,jRN,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN+1,jRN,kRN)
                phiRR = sgndf_RN(iRN+2,jRN,kRN)

                dphi_dcsi_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_dcsi_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            end if


            ! Con dphi_dcsi_plus y dphi_dcsi_minus calculadas, puedo calcular dphi_csi en el nodo

            if(signo_RN(iRN,jRN,kRN)*(phiR-phiC).lt.zero.and.&
               (signo_RN(iRN,jRN,kRN)*(phiC-phiL)+signo_RN(iRN,jRN,kRN)*(phiR-phiC)).lt.zero) then

                dphi_dcsi(iRN,jRN,kRN) = dphi_dcsi_plus

            else if(signo_RN(iRN,jRN,kRN)*(phiC-phiL).gt.zero.and.&
                    (signo_RN(iRN,jRN,kRN)*(phiR-phiC)+signo_RN(iRN,jRN,kRN)*(phiC-phiL)).gt.zero) then

                dphi_dcsi(iRN,jRN,kRN) = dphi_dcsi_minus

            else

                dphi_dcsi(iRN,jRN,kRN) = one_half*(dphi_dcsi_plus+dphi_dcsi_minus)
    
            end if

            !------------------------------------------------------------------------
            ! DIRECCION ETA
            !------------------------------------------------------------------------
            ! Escojo el stencil para calcular el gradiente 
            ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
            if (obstacle_lsm_ad(iRN,jRN,kRN,2).eq.1.or.&
                bd_stencil_mr.eq.1) then
    
                phiLL = sgndf_RN(iRN,jRN-2,kRN)
                phiL  = sgndf_RN(iRN,jRN-1,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN+1,kRN)
                phiRR = 0.0 ! outbounded

                dphi_deta_plus  = phiR-phiC-one_half*(phiR-2*phiC+phiL)
                dphi_deta_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            else if(obstacle_lsm_ad(iRN,jRN,kRN,4).eq.1.or.&
                bd_stencil_ml.eq.1) then

                phiLL = 0.0 ! outbounded
                phiL  = sgndf_RN(iRN,jRN-1,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN+1,kRN)
                phiRR = sgndf_RN(iRN,jRN+2,kRN)

                dphi_deta_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_deta_minus = phiC-phiL+one_half*(phiR-2*phiC+phiL)

            else ! no tengo restricciones de stencil
                
                phiLL = sgndf_RN(iRN,jRN-2,kRN)
                phiL  = sgndf_RN(iRN,jRN-1,kRN)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN+1,kRN)
                phiRR = sgndf_RN(iRN,jRN+2,kRN)

                dphi_deta_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_deta_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            end if


            ! Con dphi_deta_plus y dphi_deta_minus calculadas, puedo calcular dphi_eta en el nodo

            if(signo_RN(iRN,jRN,kRN)*(phiR-phiC).lt.zero.and.&
               (signo_RN(iRN,jRN,kRN)*(phiC-phiL)+signo_RN(iRN,jRN,kRN)*(phiR-phiC)).lt.zero) then

                dphi_deta(iRN,jRN,kRN) = dphi_deta_plus

            else if(signo_RN(iRN,jRN,kRN)*(phiC-phiL).gt.zero.and.&
                    (signo_RN(iRN,jRN,kRN)*(phiR-phiC)+signo_RN(iRN,jRN,kRN)*(phiC-phiL)).gt.zero) then

                dphi_deta(iRN,jRN,kRN) = dphi_deta_minus

            else

                dphi_deta(iRN,jRN,kRN) = one_half*(dphi_deta_plus+dphi_deta_minus)
    
            end if

            !------------------------------------------------------------------------
            ! DIRECCION ZETA
            !------------------------------------------------------------------------
            ! Escojo el stencil para calcular el gradiente 
            ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
            ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
            if (bd_stencil_mu.eq.1) then
    
                phiLL = sgndf_RN(iRN,jRN,kRN-2)
                phiL  = sgndf_RN(iRN,jRN,kRN-1)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN,kRN+1)
                phiRR = 0.0 ! outbounded

                dphi_dzet_plus  = phiR-phiC-one_half*(phiR-2*phiC+phiL)
                dphi_dzet_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            else if(bd_stencil_md.eq.1) then

                phiLL = 0.0 ! outbounded
                phiL  = sgndf_RN(iRN,jRN,kRN-1)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN,kRN+1)
                phiRR = sgndf_RN(iRN,jRN,kRN+2)

                dphi_dzet_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_dzet_minus = phiC-phiL+one_half*(phiR-2*phiC+phiL)

            else ! no tengo restricciones de stencil
                
                phiLL = sgndf_RN(iRN,jRN,kRN-2)
                phiL  = sgndf_RN(iRN,jRN,kRN-1)
                phiC  = sgndf_RN(iRN,jRN,kRN)
                phiR  = sgndf_RN(iRN,jRN,kRN+1)
                phiRR = sgndf_RN(iRN,jRN,kRN+2)

                dphi_dzet_plus  = phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC)
                dphi_dzet_minus = phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL)
    
            end if


            ! Con dphi_dzet_plus y dphi_dzet_minus calculadas, puedo calcular dphi_eta en el nodo

            if(signo_RN(iRN,jRN,kRN)*(phiR-phiC).lt.zero.and.&
               (signo_RN(iRN,jRN,kRN)*(phiC-phiL)+signo_RN(iRN,jRN,kRN)*(phiR-phiC)).lt.zero) then

                dphi_dzet(iRN,jRN,kRN) = dphi_dzet_plus

            else if(signo_RN(iRN,jRN,kRN)*(phiC-phiL).gt.zero.and.&
                    (signo_RN(iRN,jRN,kRN)*(phiR-phiC)+signo_RN(iRN,jRN,kRN)*(phiC-phiL)).gt.zero) then

                dphi_dzet(iRN,jRN,kRN) = dphi_dzet_minus

            else

                dphi_dzet(iRN,jRN,kRN) = one_half*(dphi_dzet_plus+dphi_dzet_minus)
    
            end if


            ! Calculo del gradiente en coordenadas curvilineas
        
        
            dphi_dx = csi(1,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
                      eta(1,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
                      zet(1,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)
        
            dphi_dy = csi(2,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
                      eta(2,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
                      zet(2,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)
        
            dphi_dx = csi(3,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
                      eta(3,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
                      zet(3,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)
        
            norma_grad_phi = sqrt((dphi_dx)**2+(dphi_dy)**2+(dphi_dz)**2)
        
            ! RHS ecuacion de reinicializacion: s(phi0)*(|grad(phi)|-1) 
        
            rightH_RN(iRN,jRN,kRN) = signo_RN(iRN,jRN,kRN) * (norma_grad_phi - one)     
        
        end if ! Fuera del obstaculo

    end do ! iRN
    end do ! jRN
    end do ! kRN


else ! procesador sin parte del obstaculo

    do iRN = iRN_sta,iRN_end
    do jRN = jRN_sta,jRN_end
    do kRN = kRN_sta,kRN_end


        ! Definicion de los booleans para stencils bounded en los bordes

        if(myback.eq.mpi_proc_null.and.iRN.eq.iRN_sta) then
            bd_stencil_mb = 1
        else
            bd_stencil_mb = 0
        end if


        if(myleft.eq.mpi_proc_null.and.jRN.eq.jRN_sta) then
            bd_stencil_ml = 1
        else
            bd_stencil_ml = 0
        end if

        if(mydown.eq.mpi_proc_null.and.kRN.eq.kRN_sta) then
            bd_stencil_md = 1
        else
            bd_stencil_md = 0
        end if


        if(myfront.eq.mpi_proc_null.and.iRN.eq.iRN_end) then
            bd_stencil_mf = 1
        else
            bd_stencil_mf = 0
        end if

        if(myright.eq.mpi_proc_null.and.jRN.eq.jRN_end) then
            bd_stencil_mr = 1
        else
            bd_stencil_mr = 0
        end if

        if(myup.eq.mpi_proc_null.and.kRN.eq.kRN_end) then
            bd_stencil_mu = 1
        else
            bd_stencil_mu = 0
        end if
    
        !------------------------------------------------------------------------
        ! DIRECCION CSI
        !------------------------------------------------------------------------

        ! Escojo el stencil para calcular el gradiente 
        ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
        ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
        if (bd_stencil_mf.eq.1) then
    
!            phiLL = sgndf_RN( iRN-2 , jRN , kRN )
!            phiL  = sgndf_RN( iRN-1 , jRN , kRN )
!            phiC  = sgndf_RN( iRN   , jRN , kRN )
!            phiR  = sgndf_RN( iRN+1 , jRN , kRN )
!            phiRR = zero ! outbounded
!
!            dphi_dcsi_plus  = dc * ( phiR-phiC-one_half*(phiR-2*phiC+phiL) )
!            dphi_dcsi_minus = dc * ( phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL))
!    
!            print *, 'myfront boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_dcsi_plus  = ', dphi_dcsi_plus
!            print *, 'dphi_dcsi_minus = ', dphi_dcsi_minus
!            print *, ' '

        else if(bd_stencil_mb.eq.1) then

!            phiLL = zero ! outbounded
!            phiL  = sgndf_RN( iRN-1 , jRN , kRN )
!            phiC  = sgndf_RN( iRN   , jRN , kRN )
!            phiR  = sgndf_RN( iRN+1 , jRN , kRN )
!            phiRR = sgndf_RN( iRN+2 , jRN , kRN )
!
!            dphi_dcsi_plus  = dc * ( phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC) )
!            dphi_dcsi_minus = dc * ( phiC-phiL+one_half*(phiR-2*phiC+phiL) )
!
!            print *, 'myback boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_dcsi_plus  = ', dphi_dcsi_plus
!            print *, 'dphi_dcsi_minus = ', dphi_dcsi_minus
!            print *, ' '

        else ! no tengo restricciones de stencil
            
            phiLL = sgndf_RN( iRN-2 , jRN , kRN )
            phiL  = sgndf_RN( iRN-1 , jRN , kRN )
            phiC  = sgndf_RN( iRN   , jRN , kRN )
            phiR  = sgndf_RN( iRN+1 , jRN , kRN )
            phiRR = sgndf_RN( iRN+2 , jRN , kRN )

            ! ---------- 
            ! (∂Φ/∂ξ)+ 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiRR-two*phiR+phiC)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiRR-two*phiR+phiC
            end if


            dphi_dcsi_plus  = dc * ( phiR - phiC - one_half * Mphi )
            

            ! ---------- 
            ! (∂Φ/∂ξ)- 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiC-two*phiL+phiLL)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiC-two*phiL+phiLL
            end if

            dphi_dcsi_minus = dc * ( phiC - phiL + one_half * Mphi ) 
    
        end if

        if      ( sgndf_RN(iRN,jRN,kRN) < -eps_sims) then       
            signPhi = -one
        else if ( sgndf_RN(iRN,jRN,kRN) >  eps_sims) then
            signPhi =  one
        else
            signPhi = zero            
        end if

        ! Con dphi_dcsi_plus y dphi_dcsi_minus calculadas, puedo calcular dphi_csi en el nodo

        if      ( signPhi*(phiR-phiC) <  zero .and. &
                  signPhi*(phiC-phiL) < -signPhi*(phiR-phiC) ) then

            dphi_dcsi(iRN,jRN,kRN) = dphi_dcsi_plus

        else if ( signPhi*(phiC-phiL) >  zero .and. &
                  signPhi*(phiR-phiC) > -signPhi*(phiC-phiL) ) then

            dphi_dcsi(iRN,jRN,kRN) = dphi_dcsi_minus

        else

            dphi_dcsi(iRN,jRN,kRN) = one_half * ( dphi_dcsi_plus + dphi_dcsi_minus )
    
        end if

        !------------------------------------------------------------------------
        ! DIRECCION ETA
        !------------------------------------------------------------------------
        ! Escojo el stencil para calcular el gradiente 
        ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
        ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
        if (bd_stencil_mr.eq.1) then
    
!            phiLL = sgndf_RN( iRN, jRN-2  , kRN )
!            phiL  = sgndf_RN( iRN, jRN-1  , kRN )
!            phiC  = sgndf_RN( iRN, jRN    , kRN )
!            phiR  = sgndf_RN( iRN, jRN+1  , kRN )
!            phiRR = zero ! outbounded
!
!            dphi_deta_plus  = de * ( phiR-phiC-one_half*(phiR-2*phiC+phiL) )
!            dphi_deta_minus = de * ( phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL) )
!    
!            print *, 'myright boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_deta_plus  = ', dphi_deta_plus
!            print *, 'dphi_deta_minus = ', dphi_deta_minus
!            print *, ' '

        else if(bd_stencil_ml.eq.1) then

!            phiLL = zero ! outbounded
!            phiL  = sgndf_RN( iRN , jRN-1 , kRN )
!            phiC  = sgndf_RN( iRN , jRN   , kRN )
!            phiR  = sgndf_RN( iRN , jRN+1 , kRN )
!            phiRR = sgndf_RN( iRN , jRN+2 , kRN )
!
!            dphi_deta_plus  = de * ( phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC) )
!            dphi_deta_minus = de * ( phiC-phiL+one_half*(phiR-2*phiC+phiL) )
!
!            print *, 'myleft boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_deta_plus  = ', dphi_deta_plus
!            print *, 'dphi_deta_minus = ', dphi_deta_minus
!            print *, ' '

        else ! no tengo restricciones de stencil
            
            phiLL = sgndf_RN( iRN , jRN-2 , kRN )
            phiL  = sgndf_RN( iRN , jRN-1 , kRN )
            phiC  = sgndf_RN( iRN , jRN   , kRN )
            phiR  = sgndf_RN( iRN , jRN+1 , kRN )
            phiRR = sgndf_RN( iRN , jRN+2 , kRN )

            ! ---------- 
            ! (∂Φ/∂η)+ 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiRR-two*phiR+phiC)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiRR-two*phiR+phiC
            end if


            dphi_deta_plus  = de * ( phiR - phiC - one_half * Mphi )
            

            ! ---------- 
            ! (∂Φ/∂η)- 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiC-two*phiL+phiLL)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiC-two*phiL+phiLL
            end if

            dphi_deta_minus = de * ( phiC - phiL + one_half * Mphi )
    
        end if

        if      ( sgndf_RN(iRN,jRN,kRN) < -eps_sims) then       
            signPhi = -one
        else if ( sgndf_RN(iRN,jRN,kRN) >  eps_sims) then
            signPhi =  one
        else
            signPhi = zero            
        end if

        ! Con dphi_dcsi_plus y dphi_dcsi_minus calculadas, puedo calcular dphi_csi en el nodo

        if      ( signPhi*(phiR-phiC) <  zero .and. &
                  signPhi*(phiC-phiL) < -signPhi*(phiR-phiC) ) then

            dphi_deta(iRN,jRN,kRN) = dphi_deta_plus

        else if ( signPhi*(phiC-phiL) >  zero .and. &
                  signPhi*(phiR-phiC) > -signPhi*(phiC-phiL) ) then

            dphi_deta(iRN,jRN,kRN) = dphi_deta_minus

        else

            dphi_deta(iRN,jRN,kRN) = one_half * ( dphi_deta_plus + dphi_deta_minus )
    
        end if

        !------------------------------------------------------------------------
        ! DIRECCION ZETA
        !------------------------------------------------------------------------
        ! Escojo el stencil para calcular el gradiente 
        ! si estoy en la zona 1 del obstaculo, fuerzo upwind en i.
        ! Si estoy en el borde de salida, tambien fuerzo stencil centrado
    
        if (bd_stencil_mu.eq.1) then
    
!            phiLL = sgndf_RN( iRN , jRN , kRN-2 )
!            phiL  = sgndf_RN( iRN , jRN , kRN-1 )
!            phiC  = sgndf_RN( iRN , jRN , kRN   )
!            phiR  = sgndf_RN( iRN , jRN , kRN+1 )
!            phiRR = zero ! outbounded
!
!            dphi_dzet_plus  = dz * ( phiR-phiC-one_half*(phiR-2*phiC+phiL) )
!            dphi_dzet_minus = dz * ( phiC-phiL+one_half*min(phiR-2*phiC+phiL,phiC-2*phiL+phiLL) )
!    
!            print *, 'myup boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_dzet_plus  = ', dphi_dzet_plus
!            print *, 'dphi_dzet_minus = ', dphi_dzet_minus
!            print *, ' '

        else if(bd_stencil_md.eq.1) then

!            phiLL = zero ! outbounded
!            phiL  = sgndf_RN( iRN , jRN , kRN-1 )
!            phiC  = sgndf_RN( iRN , jRN , kRN   )
!            phiR  = sgndf_RN( iRN , jRN , kRN+1 )
!            phiRR = sgndf_RN( iRN , jRN , kRN+2 )
!
!            dphi_dzet_plus  = dz * ( phiR-phiC-one_half*min(phiR-2*phiC+phiL,phiRR-2*phiR+phiC) )
!            dphi_dzet_minus = dz * ( phiC-phiL+one_half*(phiR-2*phiC+phiL) )
!
!            print *, 'mydown boundary'
!            print *, 'i,j,k  = ',iRN, jRN, kRN
!            print *, 'dphi_dzet_plus  = ', dphi_dzet_plus
!            print *, 'dphi_dzet_minus = ', dphi_dzet_minus
!            print *, ' '

        else ! no tengo restricciones de stencil
            
            phiLL = sgndf_RN( iRN , jRN , kRN-2 )
            phiL  = sgndf_RN( iRN , jRN , kRN-1 )
            phiC  = sgndf_RN( iRN , jRN , kRN   )
            phiR  = sgndf_RN( iRN , jRN , kRN+1 )
            phiRR = sgndf_RN( iRN , jRN , kRN+2 )

            ! ---------- 
            ! (∂Φ/∂ζ)+ 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiRR-two*phiR+phiC)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiRR-two*phiR+phiC
            end if


            dphi_dzet_plus  = dz * ( phiR - phiC - one_half * Mphi )
            

            ! ---------- 
            ! (∂Φ/∂ζ)- 
            ! ---------- 

            if(abs(phiR-two*phiC+phiL) < abs(phiC-two*phiL+phiLL)) then
                Mphi = phiR-two*phiC+phiL
            else
                Mphi = phiC-two*phiL+phiLL
            end if

            dphi_dzet_minus = dz * ( phiC - phiL + one_half * Mphi )
    
        end if

        if      ( sgndf_RN(iRN,jRN,kRN) < -eps_sims) then       
            signPhi = -one
        else if ( sgndf_RN(iRN,jRN,kRN) >  eps_sims) then
            signPhi =  one
        else
            signPhi = zero            
        end if

        ! Con dphi_dcsi_plus y dphi_dcsi_minus calculadas, puedo calcular dphi_csi en el nodo

        if      ( signPhi*(phiR-phiC) <  zero .and. &
                  signPhi*(phiC-phiL) < -signPhi*(phiR-phiC) ) then

            dphi_dzet(iRN,jRN,kRN) = dphi_dzet_plus

        else if ( signPhi*(phiC-phiL) >  zero .and. &
                  signPhi*(phiR-phiC) > -signPhi*(phiC-phiL) ) then

            dphi_dzet(iRN,jRN,kRN) = dphi_dzet_minus

        else

            dphi_dzet(iRN,jRN,kRN) = one_half * ( dphi_dzet_plus + dphi_dzet_minus )
    
        end if


    ! Calculo del gradiente en coordenadas curvilineas


    dphi_dx = csi(1,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
              eta(1,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
              zet(1,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)

    dphi_dy = csi(2,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
              eta(2,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
              zet(2,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)

    dphi_dz = csi(3,iRN,jRN,kRN) * dphi_dcsi(iRN,jRN,kRN) + &
              eta(3,iRN,jRN,kRN) * dphi_deta(iRN,jRN,kRN) + &
              zet(3,iRN,jRN,kRN) * dphi_dzet(iRN,jRN,kRN)

    norma_grad_phi = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

!    if(norma_grad_phi > 0.99_rdf) then
!        print *, 'i,j,k : ', iRN, jRN, kRN
!        print *, 'normagrad_RN(iRN,jRN,kRN) = ', norma_grad_phi
!  
!    end if

    ! RHS ecuacion de reinicializacion: s(phi0)*(|grad(phi)|-1) 

    rightH_RN(iRN,jRN,kRN) = signo_RN(iRN,jRN,kRN) * (norma_grad_phi - one)     
        
    end do ! iRN
    end do ! jRN
    end do ! kRN

end if ! if isobstacle==1

end subroutine calc_RH_RN_2
