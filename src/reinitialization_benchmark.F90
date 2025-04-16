subroutine reinitialization_benchmark(  phizero, phicorregido,                  &
                                        M_Sphi0, M_TI, M_Grad, M_eps, M_lambda, &
                                        InterfaceNodesID                          )

!phi que esta en la rutina levelsetmethod es considerado el valor inicial
!reinicializiacion primer orden
!rutina que sirve cuando es una funcion que esta cerca de una signed distance step function

use global_LSM, only : heaviside, deltaf
use global_app
use global_mpi
use global_debug
use DataTypes
!use LinkedListModule

implicit none

real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in)    :: phizero
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: phicorregido

!variables para reinicializacion
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: sgndf, sgndf_n
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: rh_RN
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: sign_ini, normagrad_ini
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: lambda_RN

integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: InterfaceNodesID 

integer :: iterReini, iterReiniMax
integer :: iRN,jRN,kRN  !index reinitialization
integer :: RK_RN        !runge kutta reinitialization index (RK-4)
real (kind = rdf),dimension(4) :: coef_RK_RN
real (kind = rdf),dimension(3) :: coef2,coef3
real (kind = rdf) :: timeStepRN
integer :: iRN_sta,iRN_end
integer :: jRN_sta,jRN_end
integer :: kRN_sta,kRN_end
character(len = 256) :: debugname
!-------------------------------------------------------------------------------------
! Switches para probar modificaciones en la reinicializacion

integer, intent(in) :: M_Sphi0,M_TI,M_Grad,M_eps,M_lambda
real (kind = rdf),dimension(3,3) :: coef_RKTVD_AD

coef_RKTVD_AD(1,1) = 1.0_rdf
coef_RKTVD_AD(1,2) = 3.0_rdf/4.0_rdf
coef_RKTVD_AD(1,3) = 1.0_rdf/3.0_rdf
coef_RKTVD_AD(2,1) = 0.0_rdf
coef_RKTVD_AD(2,2) = 1.0_rdf/4.0_rdf
coef_RKTVD_AD(2,3) = 2.0_rdf/3.0_rdf
coef_RKTVD_AD(3,1) = 1.0_rdf
coef_RKTVD_AD(3,2) = 1.0_rdf/4.0_rdf
coef_RKTVD_AD(3,3) = 2.0_rdf/3.0_rdf

!-------------------------------------------------------------------------------------
!Nodos incluyendo el borde
iRN_sta = il + igp
jRN_sta = jl + jgp
kRN_sta = kl + kgp

iRN_end = iu - igp
jRN_end = ju - jgp
kRN_end = ku - kgp

iterReiniMax =  RNitermax    !iteraciones máximas
timeStepRN   =  RNtimestep   !time step de la reinicializacion (mesh_size/10)

sgndf_n = phizero     !d0 (valorinicial)
sgndf   = sgndf_n     !inicializacion de el que va variando en el tiempo

!gp deben estar actualizados para la entrada

!calcular  s(d0)

if(M_Grad==1) then ! Calculo del gradiente de phi0 usando ENO2 de Kang 2010  
  call calc_sign_ini_2(phizero,sign_ini,normagrad_ini,epslsm,M_Sphi0)  !calcula s(d0) y |normad0| !requiere gp actualizado para phi
else
  call calc_sign_ini(phizero,sign_ini,normagrad_ini,epslsm,M_Sphi0)  !calcula s(d0) y |normad0| !requiere gp actualizado para phi
  if(is_obstacle == 1) call calc_sign_ini_obstacle(phizero,sign_ini,normagrad_ini,epslsm,M_Sphi0)
end if

!if(MOD(iteraciontiempo,phi_outputiter) == 0 ) then
        if(myid ==0) print *, 'iteracion tiempo = ', iteraciontiempo
        write(debugname, fmt ='(a,i6.6)') 'normgrad_ini',iteraciontiempo
        call outputD1_real(normagrad_ini,debugname) 
!end if


! Intercambio de gp con los otros procesadores

call rhs_exchng3_3d (sign_ini)
call rhs_exchng3_3d (normagrad_ini)

do iterReini=1, iterReiniMax
  
  if(M_TI==1) then ! RK3
  
    do RK = 1,3
      
      rh_RN = zero
  
      if(M_Grad==1) then ! Gradiente ENO2 de Kang 2010
  
        call calc_RH_RN_2(sgndf,sign_ini,rh_RN)
        if(is_obstacle == 1) call bcond_lsm_obstacle(sgndf) ! extrapolacion upwind de phi al obstaculo
  
      else
  
        call calc_RH_RN (sgndf,sign_ini,rh_RN)
        if(is_obstacle == 1) call calc_RH_RN_obstacle(sgndf,sign_ini,rh_RN) 
        if(is_obstacle == 1) call bcond_lsm_obstacle(sgndf) ! extrapolacion upwind de phi al obstaculo
  
      end if
  
      do iRN  = iRN_sta,iRN_end
        do jRN  = jRN_sta,jRN_end
          do kRN  = kRN_sta,kRN_end
            
            if(act_obstacle_ad(i,j,k) == 1) then
              sgndf(iRN,jRN,kRN) = coef_RKTVD_AD(1,RK)*sgndf_n(iRN,jRN,kRN) + &
                                   coef_RKTVD_AD(2,RK)*sgndf(iRN,jRN,kRN) - &
                                   coef_RKTVD_AD(3,RK)*timeStepRN*rh_RN(iRN,jRN,kRN)
            end if
            
          end do
        end do
      end do
    
      call rhs_exchng3_3d (sgndf)
  
    end do
  
  else ! Euler-Backward
  
    if(M_Grad==1) then ! Gradiente ENO2 de Kang 2010
  
      call calc_RH_RN_2(sgndf,sign_ini,rh_RN)
      if(is_obstacle == 1) call bcond_lsm_obstacle(sgndf) ! extrapolacion upwind de phi al obstaculo
  
    else
  
      call calc_RH_RN(sgndf,sign_ini,rh_RN)
      if(is_obstacle == 1) call calc_RH_RN_obstacle(sgndf,sign_ini,rh_RN) 
      if(is_obstacle == 1) call bcond_lsm_obstacle(sgndf) ! extrapolacion upwind de phi al obstaculo
  
    end if
  
    do iRN  = iRN_sta,iRN_end
    do jRN  = jRN_sta,jRN_end
    do kRN  = kRN_sta,kRN_end
        if(act_obstacle_rn(iRN,jRN,kRN) == 1) then
           sgndf(iRN,jRN,kRN) =  sgndf_n(iRN,jRN,kRN)  - timeStepRN * rh_RN(iRN,jRN,kRN)   
        end if
    end do
    end do
    end do
  
    call rhs_exchng3_3d (sgndf)
  
  end if

!---------------------------------------
! correcion lambda de SUSSMANN
!---------------------------------------

  call calc_lambda_RN(phizero,normagrad_ini,sgndf,lambda_RN,timeStepRN) !todos los gp deben tar en orden
 

  do iRN  = iRN_sta,iRN_end
  do jRN  = jRN_sta,jRN_end
  do kRN  = kRN_sta,kRN_end

      if(act_obstacle_rn(iRN,jRN,kRN) == 1) then
        sgndf(iRN,jRN,kRN) =  sgndf(iRN,jRN,kRN)  + &
                              timeStepRN * lambda_RN(iRN,jRN,kRN) * &
                              deltaf(phizero(iRN,jRN,kRN),epslsm) * normagrad_ini(iRN,jRN,kRN)
      end if

  end do
  end do
  end do

  call rhs_exchng3_3d (sgndf)        
          
  sgndf_n = sgndf
        
end do ! IterReini

! we correct the phi values only at the nodes where geometric 
! reinitialisation was not applied 

if ( hybrid_reinitialisation ) then

  where( InterfaceNodesID == 0 )
    phicorregido = sgndf 
  end where

!  do iRN  = iRN_sta,iRN_end
!  do jRN  = jRN_sta,jRN_end
!  do kRN  = kRN_sta,kRN_end
!    
!    if (InterfaceNodesID(iRN, jRN, kRN) == 0) then
!      phicorregido(iRN, jRN, kRN) = sgndf(iRN, jRN, kRN)
!    
!    else
!      print *, 'i,j,k = ',iRN,jRN,kRN
!      print *, 'phizero = ', phizero(iRN,jRN,kRN)
!      print *, ' '
!    end if


!  end do
!  end do
!  end do

else
  ! If we're using Sussman method only, then the signed distance function
  ! obtained from the PDE method is used to update the whole phi array
  
  phicorregido = sgndf   

end if

!if ( TotalVolumeComputation .and. .not.(hybrid_reinitialisation)) then
!  call TotalWaterVolumeSussman( phicorregido )
!end if

end subroutine reinitialization_benchmark
