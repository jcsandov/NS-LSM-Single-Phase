subroutine calc_RH_AD_WENO3_test2( phi_AD , rightH_AD , AdvectionNodes )
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


use AdvectionMethods
implicit none

real (kind = rdf) , dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) , dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD
integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: AdvectionNodes 


!index
integer  :: iRHAD_sta,iRHAD_end
integer  :: jRHAD_sta,jRHAD_end
integer  :: kRHAD_sta,kRHAD_end

integer :: i,j,k

! Variables WENO 3 

real (kind = rdf) :: epsWENO !epsilon positivo para mantener acotada la velocidad de Roe (Durran p251)
real (kind = rdf) :: phiL,phiC,phiR ! Valores de phi segun el stencil escogido
real (kind = rdf) :: alpha1, alpha2 ! coeficientes de ponderacion no lineales WENO
real (kind = rdf) :: phi_flux_plus_i, phi_flux_minus_i,phi_flux_plus_j, phi_flux_minus_j,phi_flux_plus_k, phi_flux_minus_k
real (kind = rdf) :: Lc,Le,Lz ! Operador espacial eq HJ en las distintas direcciones
real (kind = rdf) :: dc2,de2,dz2 ! ( 1/(2Δξ) , 1/(2Δη) , 1/(2Δζ) )
real (kind = rdf) :: dc6,de6,dz6 ! ( 1/(6Δξ) , 1/(6Δη) , 1/(6Δζ) )

real (kind = rdf) :: WaveSpeed, phi0, phi1, phi2, phi3

logical           :: SecondOrderDerivativeBC , ThirdOrderDerivativeBC 
real (kind = rdf) :: dphi_dcsi_bc , dphi_deta_bc , dphi_dzet_bc 

! Switches (1 o 0) para forzar stencils bounded en los bordes del dominio. El codigo los determina
! automaticamente evaluando la posicion del procesador en el dominio computacional

integer :: bd_stencil_mb = 0 ! bounded stencil myback
integer :: bd_stencil_ml = 0
integer :: bd_stencil_md = 0

integer :: bd_stencil_mf = 0
integer :: bd_stencil_mr = 0
integer :: bd_stencil_mu = 0

integer :: iBiasDirection, jBiasDirection, kBiasDirection

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!dc2 = one_half * dc
!de2 = one_half * de
!dz2 = one_half * dz

!dc6 = one/six * dc
!de6 = one/six * de
!dz6 = one/six * dz

!SecondOrderDerivativeBC = .true.
!ThirdOrderDerivativeBC  = .false.

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

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! INTERIOR DOMAIN LOOP
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

print *, 'max phi   = ', maxval(phi_AD)
print *, 'min phi   = ', minval(phi_AD)
print *, ' '

do i = iRHAD_sta , iRHAD_end
do j = jRHAD_sta , jRHAD_end
do k = kRHAD_sta , kRHAD_end

if ( AdvectionNodes(i,j,k) > 0  ) then


   ! print *, '------------------------------'
   ! print *, 'i,j,k = ' , i,j,k
   ! print *, ' '

   ! Definicion de los booleans para stencils bounded en los bordes

   bd_stencil_mb = 0
   bd_stencil_ml = 0
   bd_stencil_md = 0

   bd_stencil_mf = 0
   bd_stencil_mr = 0
   bd_stencil_mu = 0
   
   if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
   if( myleft == mpi_proc_null  .and. j == jRHAD_sta ) bd_stencil_ml = 1
   if( mydown == mpi_proc_null  .and. k == kRHAD_sta ) bd_stencil_md = 1

   if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
   if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1
   if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

   !------------------------------------------------------------------------
   ! ξ - DIRECTION
   !------------------------------------------------------------------------

   !-----------------
   ! i+1/2
   !-----------------

   WaveSpeed = ucn_j(1,i,j,k)

   phi0 = phi_AD ( i-1 , j , k )
   phi1 = phi_AD ( i   , j , k )
   phi2 = phi_AD ( i+1 , j , k )

   ! Outbounded stencil
   if (bd_stencil_mf == 1) then 
       
      ! I impose positive WaveSpeed to force upwind stencil in 
      ! GetWENO3Reconstruction 

      WaveSpeed = one  ! dummy value
      phi3      = zero ! dummy value
   
   ! Not-outbounded stencil        
   else 
       phi3 = phi_AD ( i+2 , j , k )
   end if        

!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' i + 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if


   ! Φ_{i+1/2}
   phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   !-----------------
   ! i-1/2
   !-----------------

   WaveSpeed = ucn_j(1,i,j,k)

   ! Outbounded stencil
   if ( bd_stencil_mb == 1 ) then 
       
       ! I impose negative WaveSpeed to force downwind stencil in 
       ! GetWENO3Reconstruction 

       WaveSpeed = - one  ! dummy value
       phi0      =   zero ! dummy value
   
   ! Not-outbounded stencil
   else 
       phi0 = phi_AD ( i-2 , j , k )
   end if        

   phi1 = phi_AD ( i-1 , j , k )
   phi2 = phi_AD ( i   , j , k )
   phi3 = phi_AD ( i+1 , j , k )

!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' i - 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if


   ! Φ_{i-1/2}
   phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   ! Reset WaveSpeed to its actual value
   WaveSpeed = ucn_j(1,i,j,k)

   ! dphi_dcsi using WENO3
   !dphi_dcsi(i,j,k) = dc * ( phi_flux_plus_i - phi_flux_minus_i ) 

   ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
   Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

   !------------------------------------------------------------------------
   ! η - DIRECTION
   !------------------------------------------------------------------------
   
   !-----------------
   ! j+1/2
   !-----------------

   WaveSpeed = ucn_j(2,i,j,k)

   phi0 = phi_AD ( i, j-1 , k )
   phi1 = phi_AD ( i, j   , k )
   phi2 = phi_AD ( i, j+1 , k )

   ! Outbounded stencil
   if (bd_stencil_mr == 1) then 
       
       ! I impose positive WaveSpeed to force upwind stencil in 
       ! GetWENO3Reconstruction 

       WaveSpeed = one  ! dummy value
       phi3      = zero ! dummy value
   
   ! Not-outbounded stencil
   else 
       phi3 = phi_AD ( i , j+2 , k )
   end if        


!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' j + 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if

   ! Φ_{j+1/2}
   phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   !-----------------
   ! j-1/2
   !-----------------

   WaveSpeed = ucn_j(2,i,j,k)

   ! Outbounded stencil
   if (bd_stencil_ml == 1) then 
       
       ! I impose negative WaveSpeed to force downwind stencil in 
       ! GetWENO3Reconstruction 

       WaveSpeed = - one  ! dummy value
       phi0      =   zero ! dummy value
   
   ! Not-outbounded stencil
   else
       phi0 = phi_AD ( i , j-2 , k )
   end if        

   phi1 = phi_AD ( i , j-1 , k )
   phi2 = phi_AD ( i , j   , k )
   phi3 = phi_AD ( i , j+1 , k )


!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' j - 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if

   ! Φ_{j-1/2}
   phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   ! Reset WaveSpeed to its actual value
   WaveSpeed = ucn_j(2,i,j,k)

   ! dphi_deta using WENO3
   !dphi_deta (i,j,k) = de * ( phi_flux_plus_j - phi_flux_minus_j ) 

   ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
   Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

   !------------------------------------------------------------------------
   ! ζ - DIRECTION
   !------------------------------------------------------------------------

   !-----------------
   ! k+1/2
   !-----------------

   WaveSpeed = ucn_j(3,i,j,k)

   phi0 = phi_AD ( i, j , k-1 )
   phi1 = phi_AD ( i, j , k   )
   phi2 = phi_AD ( i, j , k+1 )

   ! Outbounded stencil
   if (bd_stencil_mu == 1) then 
       
       ! I impose positive WaveSpeed to force upwind stencil in 
       ! GetWENO3Reconstruction 

       WaveSpeed = one  ! dummy value
       phi3      = zero ! dummy value
   
   ! Not-outbounded stencil
   else 
       phi3 = phi_AD ( i , j , k+2 )
   end if        

!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' k + 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if

   ! Φ_{k+1/2}
   phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   !-----------------
   ! k-1/2
   !-----------------

   WaveSpeed = ucn_j(3,i,j,k)

   ! Outbounded stencil
   if (bd_stencil_md == 1) then 
       
       ! I impose negative WaveSpeed to force downwind stencil in 
       ! GetWENO3Reconstruction 

       WaveSpeed = - one  ! dummy value
       phi0      =   zero ! dummy value
   
   ! Not-outbounded stencil
   else 
       phi0 = phi_AD ( i , j , k-2 )
   end if        

   phi1 = phi_AD ( i , j , k-1 )
   phi2 = phi_AD ( i , j , k   )
   phi3 = phi_AD ( i , j , k+1 )

!   if ( abs(phi0)>0.15 .or. abs(phi1)>0.15 .or. abs(phi2)>0.15 .or. abs(phi3)>0.15 ) then
!
!      print *, ' '
!      print *, ' k - 1/2 face'
!      print *, 'i,j,k = ', i,j,k
!      print *, 'phi0 = ', phi0
!      print *, 'phi1 = ', phi1
!      print *, 'phi2 = ', phi2
!      print *, 'phi3 = ', phi3
!   
!      stop
!   
!   end if


   ! Φ_{k-1/2}
   phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

   ! Reset WaveSpeed to its actual value
   WaveSpeed = ucn_j(3,i,j,k)

   ! dphi_dzet using WENO3
   !dphi_dzet (i,j,k) = dz * ( phi_flux_plus_k - phi_flux_minus_k )

   ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
   Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )
          
   !-----------------------------------------------------
   ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
   !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
   !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ     
   !-----------------------------------------------------

   rightH_AD(i,j,k) = Lc + Le + Lz

end if
end do
end do
end do

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!  ######  ####### #     # #     # ######     #    ######  #     # 
!  #     # #     # #     # # #   # #     #  #   #  #     #   # #   
!  ######  #     # #     # #  #  # #     # #     # ######     #    
!  #     # #     # #     # #   # # #     # ####### #   #      #    
!  ######  #######  #####  #     # ######  #     # #     #    #    
!                                                                  
!  #######    #     #####  #######  #####  
!  #        #   #  #       #       #       
!  #####   #     # #       #####    #####  
!  #       ####### #       #             # 
!  #       #     #  #####  #######  #####  
!                                        
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!-------------------
! i = 1 
!-------------------


if ( myback == mpi_proc_null)  then
   
   i = il + igp

   iBiasDirection = 1

   do k = kRHAD_sta , kRHAD_end
   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_ml = 0
      bd_stencil_md = 0
   
      bd_stencil_mr = 0
      bd_stencil_mu = 0
      
      if( myleft == mpi_proc_null  .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( mydown == mpi_proc_null  .and. k == kRHAD_sta ) bd_stencil_md = 1
   
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1


      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
      
      !-----------------
      ! j+1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      phi0 = phi_AD ( i, j-1 , k )
      phi1 = phi_AD ( i, j   , k )
      phi2 = phi_AD ( i, j+1 , k )

      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        

      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! j-1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_ml == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else
          phi0 = phi_AD ( i , j-2 , k )
      end if        

      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )

      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)

      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )
             
      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

       ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
       ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )
      
      WaveSpeed = ucn_j(1,i,j,k)
      
      Lc = WaveSpeed * dphi_dcsi_bc

      !-----------------------------------------------------
      ! RHS = U^1/J *  ∂Φ / ∂ξ_BC                     +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if
   end do
   end do

end if

!-------------------
! i = im
!-------------------

if ( myfront == mpi_proc_null )  then
   
   i = iu - igp

   iBiasDirection = -1

   do k = kRHAD_sta , kRHAD_end
   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_ml = 0
      bd_stencil_md = 0
   
      bd_stencil_mr = 0
      bd_stencil_mu = 0
      
      if( myleft == mpi_proc_null  .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( mydown == mpi_proc_null  .and. k == kRHAD_sta ) bd_stencil_md = 1
   
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
      
      !-----------------
      ! j+1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      phi0 = phi_AD ( i, j-1 , k )
      phi1 = phi_AD ( i, j   , k )
      phi2 = phi_AD ( i, j+1 , k )

      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        

      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! j-1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_ml == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else
          phi0 = phi_AD ( i , j-2 , k )
      end if        

      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )

      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)

      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )
             
      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
      
      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !-----------------------------------------------------
      ! RHS = U^1/J *  ∂Φ / ∂ξ_BC                     +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if
   end do
   end do

end if

!-------------------
! j = 1
!-------------------

if ( myleft == mpi_proc_null )  then
   
   j = jl + jgp

   jBiasDirection = 1

   do k = kRHAD_sta , kRHAD_end
   do i = iRHAD_sta , iRHAD_end
   
   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_mb = 0
      bd_stencil_md = 0
   
      bd_stencil_mf = 0
      bd_stencil_mu = 0
      
      if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( mydown == mpi_proc_null  .and. k == kRHAD_sta ) bd_stencil_md = 1
   
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )
             
      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do
   end do

end if

!-------------------
! j = jm
!-------------------

if ( myright == mpi_proc_null )  then
   
   j = ju - jgp

   jBiasDirection = -1

   do k = kRHAD_sta , kRHAD_end
   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_mb = 0
      bd_stencil_md = 0
   
      bd_stencil_mf = 0
      bd_stencil_mu = 0
      
      if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( mydown == mpi_proc_null  .and. k == kRHAD_sta ) bd_stencil_md = 1
   
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )
             
      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do
   end do

end if

!-------------------
! k = 1
!-------------------

if ( mydown == mpi_proc_null )  then
   
   k = kl + kgp

   kBiasDirection = 1

   do j = jRHAD_sta , jRHAD_end
   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_mb = 0
      bd_stencil_ml = 0
   
      bd_stencil_mf = 0
      bd_stencil_mr = 0
      
      if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myleft == mpi_proc_null  .and. j == jRHAD_sta ) bd_stencil_ml = 1
   
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
      
      !-----------------
      ! j+1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      phi0 = phi_AD ( i, j-1 , k )
      phi1 = phi_AD ( i, j   , k )
      phi2 = phi_AD ( i, j+1 , k )

      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        

      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! j-1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_ml == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else
          phi0 = phi_AD ( i , j-2 , k )
      end if        

      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )

      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)

      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )
             
      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)

      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do
   end do

end if

!-------------------
! k = km
!-------------------

if ( myup == mpi_proc_null )  then
   
   k = ku - kgp

   kBiasDirection = -1

   do j = jRHAD_sta , jRHAD_end
   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      ! Definicion de los booleans para stencils bounded en los bordes
   
      bd_stencil_mb = 0
      bd_stencil_ml = 0
   
      bd_stencil_mf = 0
      bd_stencil_mr = 0
      
      if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myleft == mpi_proc_null  .and. j == jRHAD_sta ) bd_stencil_ml = 1
   
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)

      ! ∂Φ/∂ξ using WENO3
      ! dphi_dcsi(i,j,k) = dc * ( phi_flux_plus_i - phi_flux_minus_i )
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
      
      !-----------------
      ! j+1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      phi0 = phi_AD ( i, j-1 , k )
      phi1 = phi_AD ( i, j   , k )
      phi2 = phi_AD ( i, j+1 , k )

      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        

      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! j-1/2
      !-----------------

      WaveSpeed = ucn_j(2,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_ml == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else
          phi0 = phi_AD ( i , j-2 , k )
      end if        

      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )

      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)

      ! ∂Φ/Φη using WENO3
      !dphi_deta(i,j,k) = de * ( phi_flux_plus_j - phi_flux_minus_j )

      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )
             
      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )


      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do
   end do

end if



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!    
! ####### ######   #####  #######  #####  
! #       #     # #       #       #      
! #####   #     # #  #### #####    #####  
! #       #     # #     # #             # 
! ####### ######   #####  #######  #####  
!                                         
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!            i1,jm,km ---E11--- im,jm,km    
!                /|(8)           /|(7)                            
!        E12----/ |             /-|------E10                                   
!            i1,j1,km---E9---im,j1,km                   
!              |(5)           |(6)|                                   
!              | E8       E3  |   E7                                
!             E5  |       |   E6  |                                
!              | i1,jm,k1-----|-im,jm,k1           
!          E4--|-/(4)         |  /(3)                 
!              |/             | / -------E2
! ζ,k      i1,j1,k1---E1----im,j1,k1
! ^   η,j     (1)            (2)   
! |  7             
! | /
! |/
! +-----------> ξ,i

! Edge 1  : 1 - 2 --> i free ; j = 1  , k = 1
! Edge 1  : 2 - 3 --> j free ; i = im , k = 1
! Edge 3  : 3 - 4 --> i free ; j = jm , k = 1
! Edge 4  : 4 - 1 --> j free ; i = 1  , k = 1
!
! Edge 5  : 1 - 5 --> k free ; i = 1  , j = 1
! Edge 6  : 2 - 6 --> k free ; i = im , j = 1
! Edge 7  : 3 - 7 --> k free ; i = im , j = jm
! Edge 8  : 4 - 8 --> k free ; i = 1  , j = jm
!
! Edge 9  : 5 - 6 --> i free ; j = 1  , k = km
! Edge 10 : 6 - 7 --> j free ; i = im , k = km
! Edge 11 : 7 - 8 --> i free ; j = jm , k = km
! Edge 12 : 8 - 5 --> j free ; i = 1  , k = km
!
! i free : E1 , E3 , E9  , E11
! j free : E2 , E4 , E10 , E12
! k free : E5 , E6 , E7  , E8
!
! NOTE: This implementations DOESN'T consider
!       parallelisation. 
!
! - - - - - - - - - - - - - - - - - - - - - - - - 


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! i - FREE EDGES
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 1
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myleft == mpi_proc_null .and. mydown == mpi_proc_null)  then

   j = jl + jgp
   k = kl + kgp

   jBiasDirection = 1
   kBiasDirection = 1

   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then


      bd_stencil_mb = 0
      bd_stencil_mf = 0
      
      if( myback  == mpi_proc_null .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )
      
      WaveSpeed = ucn_j(2,i,j,k)

      Le = WaveSpeed * dphi_deta_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)

      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 3
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myright == mpi_proc_null .and. mydown == mpi_proc_null)  then

   j = ju - jgp
   k = kl + kgp

   jBiasDirection = -1
   kBiasDirection =  1

   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_mb = 0
      bd_stencil_mf = 0
      
      if( myback  == mpi_proc_null .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)

      Le = WaveSpeed * dphi_deta_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 9
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myleft == mpi_proc_null .and. myup == mpi_proc_null)  then

   j = jl + jgp
   k = ku - kgp

   jBiasDirection =  1
   kBiasDirection = -1

   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_mb = 0
      bd_stencil_mf = 0
      
      if( myback  == mpi_proc_null .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)

      Le = WaveSpeed * dphi_deta_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 11
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myright == mpi_proc_null .and. myup == mpi_proc_null)  then

   j = ju - jgp
   k = ku - kgp

   jBiasDirection = -1
   kBiasDirection = -1

   do i = iRHAD_sta , iRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_mb = 0
      bd_stencil_mf = 0
      
      if( myback  == mpi_proc_null .and. i == iRHAD_sta ) bd_stencil_mb = 1
      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! i+1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      phi0 = phi_AD ( i-1 , j , k )
      phi1 = phi_AD ( i   , j , k )
      phi2 = phi_AD ( i+1 , j , k )
   
      ! Outbounded stencil
      if (bd_stencil_mf == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i+2 , j , k )
      end if        
   
      ! Φ_{i+1/2}
      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! i-1/2
      !-----------------
   
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_mb == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i-2 , j , k )
      end if        
   
      phi1 = phi_AD ( i-1 , j , k )
      phi2 = phi_AD ( i   , j , k )
      phi3 = phi_AD ( i+1 , j , k )
   
      ! Φ_{i-1/2}
      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(1,i,j,k)
   
      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)

      Le = WaveSpeed * dphi_deta_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! j - FREE EDGES
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 2
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. mydown == mpi_proc_null)  then

   i = iu - igp
   k = kl + kgp

   iBiasDirection = -1
   kBiasDirection =  1

   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_ml = 0
      bd_stencil_mr = 0
      
      if( myleft  == mpi_proc_null .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! j+1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      phi0 = phi_AD ( i , j-1 , k )
      phi1 = phi_AD ( i , j   , k )
      phi2 = phi_AD ( i , j+1 , k )
   
      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        
   
      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! j-1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_ml == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j-2 , k )
      end if        
   
      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )
   
      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 4
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. mydown == mpi_proc_null)  then

   i = il + jgp
   k = kl + kgp

   iBiasDirection =  1
   kBiasDirection =  1

   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_ml = 0
      bd_stencil_mr = 0
      
      if( myleft  == mpi_proc_null .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! j+1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      phi0 = phi_AD ( i , j-1 , k )
      phi1 = phi_AD ( i , j   , k )
      phi2 = phi_AD ( i , j+1 , k )
   
      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        
   
      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! j-1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_ml == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j-2 , k )
      end if        
   
      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )
   
      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )
      
      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 10
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myup == mpi_proc_null)  then

   i = iu - jgp
   k = ku - kgp

   iBiasDirection = -1
   kBiasDirection = -1

   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_ml = 0
      bd_stencil_mr = 0
      
      if( myleft  == mpi_proc_null .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! j+1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      phi0 = phi_AD ( i , j-1 , k )
      phi1 = phi_AD ( i , j   , k )
      phi2 = phi_AD ( i , j+1 , k )
   
      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        
   
      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! j-1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_ml == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j-2 , k )
      end if        
   
      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )
   
      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 12
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myup == mpi_proc_null)  then

   i = il + igp
   k = ku - kgp

   iBiasDirection =  1
   kBiasDirection = -1

   do j = jRHAD_sta , jRHAD_end

   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_ml = 0
      bd_stencil_mr = 0
      
      if( myleft  == mpi_proc_null .and. j == jRHAD_sta ) bd_stencil_ml = 1
      if( myright == mpi_proc_null .and. j == jRHAD_end ) bd_stencil_mr = 1

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
   
      !-----------------
      ! j+1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      phi0 = phi_AD ( i , j-1 , k )
      phi1 = phi_AD ( i , j   , k )
      phi2 = phi_AD ( i , j+1 , k )
   
      ! Outbounded stencil
      if (bd_stencil_mr == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil        
      else 
          phi3 = phi_AD ( i , j+2 , k )
      end if        
   
      ! Φ_{j+1/2}
      phi_flux_plus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      !-----------------
      ! j-1/2
      !-----------------
   
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Outbounded stencil
      if ( bd_stencil_ml == 1 ) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 
   
          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j-2 , k )
      end if        
   
      phi1 = phi_AD ( i , j-1 , k )
      phi2 = phi_AD ( i , j   , k )
      phi3 = phi_AD ( i , j+1 , k )
   
      ! Φ_{j-1/2}
      phi_flux_minus_j = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
   
      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(2,i,j,k)
   
      ! Lη = U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη
      Le = WaveSpeed * ( de * ( phi_flux_plus_j - phi_flux_minus_j ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                            phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            kBiasDirection                                )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ_BC
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! k - FREE EDGES
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 5
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myleft == mpi_proc_null)  then

   i = il + igp
   j = jl + jgp

   iBiasDirection =  1
   jBiasDirection =  1

   do k = kRHAD_sta , kRHAD_end
   
   if ( AdvectionNodes(i,j,k) > 0  ) then

      bd_stencil_md = 0
      bd_stencil_mu = 0
      
      if( mydown  == mpi_proc_null .and. k == kRHAD_sta ) bd_stencil_md = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )
      
      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      
 
   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 6
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null)  then

   i = iu - igp
   j = jl + jgp

   iBiasDirection = -1
   jBiasDirection =  1

   do k = kRHAD_sta , kRHAD_end
   
   if ( AdvectionNodes(i,j,k) > 0 ) then

      bd_stencil_md = 0
      bd_stencil_mu = 0
      
      if( mydown  == mpi_proc_null .and. k == kRHAD_sta ) bd_stencil_md = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )


      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                       +
      !       U^2/J * ∂Φ / ∂η_BC                       +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end if

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 7
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null)  then

   i = iu - igp
   j = ju - jgp

   iBiasDirection = -1
   jBiasDirection = -1

   do k = kRHAD_sta , kRHAD_end

   if ( AdvectionNodes(i,j,k) > 0 ) then

      bd_stencil_md = 0
      bd_stencil_mu = 0
      
      if( mydown  == mpi_proc_null .and. k == kRHAD_sta ) bd_stencil_md = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------
      
      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                         +
      !       U^2/J * ∂Φ / ∂η_BC                         +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      
 
   end if

   end do

end if

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 8
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myright == mpi_proc_null)  then

   i = il + igp
   j = ju - jgp

   iBiasDirection =  1
   jBiasDirection = -1

   do k = kRHAD_sta , kRHAD_end

   if ( AdvectionNodes(i,j,k) > 0) then

      bd_stencil_md = 0
      bd_stencil_mu = 0
      
      if( mydown  == mpi_proc_null .and. k == kRHAD_sta ) bd_stencil_md = 1
      if( myup    == mpi_proc_null .and. k == kRHAD_end ) bd_stencil_mu = 1

      !------------------------------------------------------------------------
      ! ζ - DIRECTION
      !------------------------------------------------------------------------

      !-----------------
      ! k+1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      phi0 = phi_AD ( i, j , k-1 )
      phi1 = phi_AD ( i, j , k   )
      phi2 = phi_AD ( i, j , k+1 )

      ! Outbounded stencil
      if (bd_stencil_mu == 1) then 
          
          ! I impose positive WaveSpeed to force upwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = one  ! dummy value
          phi3      = zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi3 = phi_AD ( i , j , k+2 )
      end if        

      ! Φ_{k+1/2}
      phi_flux_plus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      !-----------------
      ! k-1/2
      !-----------------

      WaveSpeed = ucn_j(3,i,j,k)

      ! Outbounded stencil
      if (bd_stencil_md == 1) then 
          
          ! I impose negative WaveSpeed to force downwind stencil in 
          ! GetWENO3Reconstruction 

          WaveSpeed = - one  ! dummy value
          phi0      =   zero ! dummy value
      
      ! Not-outbounded stencil
      else 
          phi0 = phi_AD ( i , j , k-2 )
      end if        

      phi1 = phi_AD ( i , j , k-1 )
      phi2 = phi_AD ( i , j , k   )
      phi3 = phi_AD ( i , j , k+1 )

      ! Φ_{k-1/2}
      phi_flux_minus_k = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )

      ! Reset WaveSpeed to its actual value
      WaveSpeed = ucn_j(3,i,j,k)

      ! Lζ = U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      Lz = WaveSpeed * ( dz * ( phi_flux_plus_k - phi_flux_minus_k ) )

      !------------------------------------------------------------------------
      ! ξ - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                            phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                            OrderLSAdvectionBoundaries               , &
                                            iBiasDirection                               )

      WaveSpeed = ucn_j(1,i,j,k)

      Lc = WaveSpeed * dphi_dcsi_bc

      !------------------------------------------------------------------------
      ! η - DIRECTION
      !------------------------------------------------------------------------

      ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
      ! for the BiasedDerivative function

      dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                            phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                            OrderLSAdvectionBoundaries                , &
                                            jBiasDirection                                )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta_bc

      !-----------------------------------------------------
      ! RHS = U^1/J * ∂Φ / ∂ξ_BC                      +
      !       U^2/J * ∂Φ / ∂η_BC                      +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      
   
   end if
   
   end do

end if


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! #     # ####### ######  ####### ###  #####  #######  #####  
! #     # #       #     #    #     #  #       #       #       
! #     # #####   ######     #     #  #       #####    #####  
!  #   #  #       #   #      #     #  #       #             # 
!    #    ####### #    #     #    ###  #####  #######  #####  
!                                                            
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


   !    i,j+1,k+1 ----i+1,j+1,k+1    
   !     /|(8)           /|(7)                            
   !    / |             / |                                  
   ! i,j,k+1-------i+1,j,k+1                   
   !   |(5)           |(6)|                                   
   !   |  |           |   |                                
   !   |  |           |   |                                
   !   |  i,j+1,k-----|-i+1,j+1,k           
   !   | /(4)         |  /(3)                 
   !   |/             | /
   ! i,j,k-----------i+1,j,k
   !  (1)               (2)



! VERTEX 1

if ( myback == mpi_proc_null .and. myleft == mpi_proc_null &
                             .and. mydown == mpi_proc_null    )  then

   i = il + igp
   j = jl + jgp
   k = kl + kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection =  1
   jBiasDirection =  1
   kBiasDirection =  1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if


! VERTEX 2

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null &
                              .and. mydown == mpi_proc_null    )  then

   i = iu - igp
   j = jl + jgp
   k = kl + kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection = -1
   jBiasDirection =  1
   kBiasDirection =  1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

! VERTEX 3

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null &
                              .and. mydown  == mpi_proc_null    )  then

   i = iu - igp
   j = ju - jgp
   k = kl + kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection = -1
   jBiasDirection = -1
   kBiasDirection =  1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

! VERTEX 4

if ( myback == mpi_proc_null  .and. myright == mpi_proc_null &
                              .and. mydown  == mpi_proc_null    )  then

   i = il + igp
   j = ju - jgp
   k = kl + kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection =  1
   jBiasDirection = -1
   kBiasDirection =  1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if


! VERTEX 5

if ( myback == mpi_proc_null .and. myleft == mpi_proc_null &
                             .and. myup   == mpi_proc_null    )  then

   i = il + igp
   j = jl + jgp
   k = ku - kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection =  1
   jBiasDirection =  1
   kBiasDirection = -1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

! VERTEX 6

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null &
                              .and. myup   == mpi_proc_null    )  then

   i = iu - igp
   j = jl + jgp
   k = ku - kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection = -1
   jBiasDirection =  1
   kBiasDirection = -1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

! VERTEX 7

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null &
                              .and. myup    == mpi_proc_null    )  then

   i = iu - igp
   j = ju - jgp
   k = ku - kgp

   if ( AdvectionNodes(i,j,k) > 0  ) then

   iBiasDirection = -1
   jBiasDirection = -1
   kBiasDirection = -1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )

   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

! VERTEX 8

if ( myback == mpi_proc_null  .and. myright == mpi_proc_null &
                              .and. myup    == mpi_proc_null    )  then

   i = il + igp
   j = ju - jgp
   k = ku - kgp

   if ( AdvectionNodes(i,j,k) > 0 ) then

   iBiasDirection =  1
   jBiasDirection = -1
   kBiasDirection = -1

   ! if ( OrderLSAdvectionBoundaries == 2 ), Φ_{i,j,k ± 3} is a dummy variable 
   ! for the BiasedDerivative function

   dphi_dcsi_bc = dc * BiasedDerivative( phi_AD( i + iBiasDirection * 0 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 1 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 2 , j , k ) , &
                                         phi_AD( i + iBiasDirection * 3 , j , k ) , &
                                         OrderLSAdvectionBoundaries               , &
                                         iBiasDirection                                )

   dphi_deta_bc = de * BiasedDerivative( phi_AD( i , j + jBiasDirection * 0  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 1  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 2  , k ) , &
                                         phi_AD( i , j + jBiasDirection * 3  , k ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         jBiasDirection                                )

   dphi_dzet_bc = dz * BiasedDerivative( phi_AD( i , j  , k + kBiasDirection * 0 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 1 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 2 ) , &
                                         phi_AD( i , j  , k + kBiasDirection * 3 ) , &
                                         OrderLSAdvectionBoundaries                , &
                                         kBiasDirection                                )
   Lc = ucn_j(1,i,j,k) * dphi_dcsi_bc
   Le = ucn_j(2,i,j,k) * dphi_deta_bc
   Lz = ucn_j(3,i,j,k) * dphi_dzet_bc

   rightH_AD(i,j,k) = Lc + Le + Lz

   end if

end if

!print *, 'max rh    = ', maxval(rightH_AD)
!print *, 'min rh    = ', minval(rightH_AD)
!print *, ' '
!print *, 'max U^j/J = ', maxval(ucn_j)
!print *, 'min U^j/J = ', minval(ucn_j)
!print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - '


end subroutine calc_RH_AD_WENO3_test2
