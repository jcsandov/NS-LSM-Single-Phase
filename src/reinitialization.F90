subroutine reinitialization(phizero,phicorregido)

!phi que esta en la rutina levelsetmethod es considerado el valor inicial
!reinicializiacion primer orden
!rutina que sirve cuando es una funcion que esta cerca de una signed distance step function

use global_LSM, only : heaviside, deltaf
implicit none

real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(out) :: phicorregido
!variables para reinicializacion
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: sgndf, sgndf_n
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: rh_RN
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: sign_ini, normagrad_ini
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: lambda_RN

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
!Nodos incluyendo el borde
iRN_sta = il + igp
jRN_sta = jl + jgp
kRN_sta = kl + kgp

iRN_end = iu - igp
jRN_end = ju - jgp
kRN_end = ku - kgp

iterReiniMax =RNitermax !iteraciones máximas
timeStepRN =  RNtimestep   !time step de la reinicializacion (mesh_size/10)

sgndf_n = phizero       !d0 (valorinicial)
sgndf = sgndf_n     !inicializacion de el que va variando en el tiempo

!gp deben estar actualizados para la entrada
!calcular  s(d0)

call calc_sign_ini(phizero,sign_ini,normagrad_ini,epslsm)  !calcula s(d0) y |normad0| !requiere gp actualizado para phi
if(is_obstacle == 1) call calc_sign_ini_obstacle(phizero,sign_ini,normagrad_ini,epslsm)
call rhs_exchng3_3d (sign_ini);call rhs_exchng3_3d (normagrad_ini)

!Debugging------------------------
!  call outputD3_real(sign_ini,'sign_ini_new') !check
!--------------------------------

do iterReini=1, iterReiniMax

  
  call calc_RH_RN (sgndf,sign_ini,rh_RN) !requiere gp de sgndf actualizados
  if(is_obstacle == 1) call calc_RH_RN_obstacle(sgndf,sign_ini,rh_RN) !calculo de RH corregido en obstaculo

!Debug----
!  call outputD3_real(rh_RN, 'rh_RN_o') !check
!------
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

!---------------------------------------
! correcion lambda de SUSSMANN
!---------------------------------------

  call calc_lambda_RN(phizero,normagrad_ini,sgndf,lambda_RN,timeStepRN) !todos los gp deben tar en orden
 

  do iRN  = iRN_sta,iRN_end
  do jRN  = jRN_sta,jRN_end
  do kRN  = kRN_sta,kRN_end

      if(act_obstacle_rn(iRN,jRN,kRN) == 1) then
        sgndf(iRN,jRN,kRN) =  sgndf(iRN,jRN,kRN)  + &
                              timeStepRN * lambda_RN(iRN,jRN,kRN)* &
                              deltaf(phizero(iRN,jRN,kRN),epslsm) * normagrad_ini(iRN,jRN,kRN)
      end if

  end do
  end do
  end do

  call rhs_exchng3_3d (sgndf)        
          
  sgndf_n = sgndf

        
end do

phicorregido = sgndf 

end subroutine reinitialization
