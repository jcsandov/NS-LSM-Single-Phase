subroutine reinitialisation_test2 (  PhiZero, PhiReinitialised, InterfaceNodesID )

!phi que esta en la rutina levelsetmethod es considerado el valor inicial
!reinicializiacion primer orden
!rutina que sirve cuando es una funcion que esta cerca de una signed distance step function

use global_LSM, only : heaviside, deltaf, GetLambdaMassCorrection   , &
                       epslsm, epsReinitialisation, phi_outputiter  , &
                       OrderReinitialisationBoundaries, ENOBCReinitialisation
use global_app
use global_mpi
use global_debug
use AdvectionMethods

implicit none

real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in)    :: PhiZero
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: PhiReinitialised

!variables para reinicializacion
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: SignPhi0 , dn, dn1
!real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: rh_RN
real (kind = rdf) :: RHS_RN
real (kind = rdf) :: MaxNormGrad_d, MinNormGrad_d
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: NormGrad_d , NormGrad_d_convergence
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: lambda_RN

real(kind = rdf) :: LambdaMassCorrection, aux

integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: InterfaceNodesID 

integer :: iterReini, iterReiniMax
integer :: iRN,jRN,kRN  !index reinitialization
integer :: RK        !runge kutta reinitialization index (RK-4)
real (kind = rdf) :: timeStepRN
integer :: iRN_sta,iRN_end
integer :: jRN_sta,jRN_end
integer :: kRN_sta,kRN_end
character(len = 256) :: debugname
real (kind = rdf),dimension(3,3) :: RKTVDcoeff

RKTVDcoeff(1,1) = 1.0_rdf
RKTVDcoeff(1,2) = 3.0_rdf/4.0_rdf
RKTVDcoeff(1,3) = 1.0_rdf/3.0_rdf
RKTVDcoeff(2,1) = 0.0_rdf
RKTVDcoeff(2,2) = 1.0_rdf/4.0_rdf
RKTVDcoeff(2,3) = 2.0_rdf/3.0_rdf
RKTVDcoeff(3,1) = 1.0_rdf
RKTVDcoeff(3,2) = 1.0_rdf/4.0_rdf
RKTVDcoeff(3,3) = 2.0_rdf/3.0_rdf

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

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! Steps
!
! I have to solve
! ∂d / ∂τ + s(ϕ0) * ( |∇d| - 1 ) = λ * δ(d) * |∇d|
!
!
! 1) I receive ϕ0 from the advection step
! 2) I initialise d^(n+1) = ϕ0 and d^n = ϕ0  
! 3) The pseudo-time loop starts
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


! I get s(ϕ0) for the whole domain (including ghost points) using the piecewise 
! function described in Kang & Sotiropoulos (AWR, 2012)

do k = kl, ku
do j = jl, ju
do i = il, iu
   ! SignPhi0 = s(ϕ0)
   SignPhi0( i,j,k ) = GetSignPhi0( phizero( i,j,k ) , epslsm )
end do
end do
end do

NormGrad_d_convergence = zero
NormGrad_d             = zero
dn1                    = phizero

call calc_NormGrad_d ( dn1 , phizero , NormGrad_d )

do iterReini = 1, iterReiniMax ! pseudo-time iteration

   MaxNormGrad_d = maxval( NormGrad_d )
   MinNormGrad_d = minval( NormGrad_d )

   ! if | max( |∇d| )-1 | < ϵ and | min( |∇d| )-1 | < ϵ , it means
   ! that ϕ was correctly redistanced over the whole domain

   if (       abs( MaxNormGrad_d - one ) < epsReinitialisation  &
        .and. abs( MinNormGrad_d - one ) < epsReinitialisation    ) exit

   ! Otherwise, I keep iterating ...

   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   !
   ! RUNGE-KUTTA STEPS
   !
   ! First, I solve
   !  
   !  ∂d / ∂τ =  - s(ϕ0) * (|∇d| - 1) + λ * δ(d*^(n+1)) * |∇(d*^(n+1))|
   !
   ! RK iterations:
   !
   ! d^(1)    =       d^(n)               -       Δτ * L(d^(n)) 
   ! d^(2)    = 3/4 * d^(n) + 1/4 * d^(1) - 1/4 * Δτ * L(d^(1)) 
   ! d*^(n+1) = 1/3 * d^(n) + 2/3 * d^(1) - 2/3 * Δτ * L(d^(2)) 
   !
   ! Where the right hand side is
   ! L(d) = s(ϕ0) * ( |∇d| - 1 ) - λ * δ(d*^(n+1)) * |∇(d*^(n+1))|
   !
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

   dn = dn1

   ! I compute the norm of the gradient of the last RK loop to identify
   ! which nodes have to be updated in the next loop

   call calc_NormGrad_d ( dn , phizero , NormGrad_d_convergence )

   do RK = 1,3

      ! I compute the norm of the gradient of the updated ϕ distribution dn1
      ! at the current RK step. Based on the last update of the RK loop, the 
      ! convergence criterion | max( |∇d| )-1 | < ϵ and | min( |∇d| )-1 | < ϵ, 
      ! is evaluated

      call calc_NormGrad_d ( dn1 , phizero , NormGrad_d )

      do k = kRN_sta, kRN_end
      do j = jRN_sta, jRN_end
      do i = iRN_sta, iRN_end
         
         ! If the local norm of the gradient already converged in the previous
         ! RK loop, then I don't modify ϕ

         if ( NormGrad_d_convergence(i,j,k) - one > epsReinitialisation ) then 

            RHS_RN =  SignPhi0(i,j,k) * ( NormGrad_d(i,j,k) - one )

            ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
            !
            ! SUSSMAN & FATEMI CORRECTION (JCP, 1999)
            !
            ! To compute λ, I use the method proposed by Kang & Sotiropoulos 
            ! (AWR, 2012) 
            ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

            LambdaMassCorrection =  GetLambdaMassCorrection ( dn1        (i,j,k) ,  &
                                                              SignPhi0   (i,j,k) ,  &
                                                              NormGrad_d (i,j,k)       )

            aux = deltaf( dn1(i,j,k) , epslsm ) * NormGrad_d(i,j,k)

            RHS_RN = RHS_RN - LambdaMassCorrection * aux

            dn1(i,j,k) = RKTVDcoeff(1,RK) * dn  (i,j,k)         + &
                         RKTVDcoeff(2,RK) * dn1 (i,j,k)         - &
                         RKTVDcoeff(3,RK) * timeStepRN * RHS_RN 

         end if

      end do
      end do
      end do
   
   end do ! RK loop
end do ! Reinitialisation convergence loop

! we correct the values of ϕ only at the nodes where geometric 
! reinitialisation was not applied 

if ( hybrid_reinitialisation ) then

   where( InterfaceNodesID == 0 )
      PhiReinitialised = dn1 
   end where

   !do k = kRN_sta, kRN_end
   !do j = jRN_sta, jRN_end
   !do i = iRN_sta, iRN_end
   !   if( InterfaceNodesID(i,j,k) == 0 ) PhiReinitialised(i,j,k) = dn1(i,j,k)
   !end do
   !end do
   !end do

else

   ! If we're using Sussman method only, then the signed distance function
   ! obtained from the PDE method is used to update the whole phi array
   
   PhiReinitialised = dn1   

end if

!if ( TotalVolumeComputation .and. .not.(hybrid_reinitialisation) ) then
!   
!   call TotalWaterVolumeSussman( PhiReinitialised )
!
!end if

end subroutine reinitialisation_test2



