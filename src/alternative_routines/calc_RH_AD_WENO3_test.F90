subroutine calc_RH_AD_WENO3_test(phi_AD,rightH_AD)
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

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: phi_AD
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: rightH_AD


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

real (kind = rdf) :: WaveSpeed, phi0, phi1, phi2, phi3

! Switches (1 o 0) para forzar stencils bounded en los bordes del dominio. El codigo los determina
! automaticamente evaluando la posicion del procesador en el dominio computacional

integer :: bd_stencil_mb = 0 ! bounded stencil myback
integer :: bd_stencil_ml = 0
integer :: bd_stencil_md = 0

integer :: bd_stencil_mf = 0
integer :: bd_stencil_mr = 0
integer :: bd_stencil_mu = 0

real (kind = rdf), dimension(3) :: csivecB , etavecB , zetvecB ! Metrics at the boundary
real (kind = rdf)               :: dphi_dcsi , dphi_deta , dphi_dzet ! gradient components
integer                         :: boundary ! 1: i , 2: j, 3: k

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

do i = iRHAD_sta , iRHAD_end
do j = jRHAD_sta , jRHAD_end
do k = kRHAD_sta , kRHAD_end

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
   ! DIRECCION CSI
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
   ! DIRECCION ETA
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
   ! DIRECCION ZETA
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
          
   !-----------------------------------------------------
   ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
   !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
   !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ     
   !-----------------------------------------------------

   rightH_AD(i,j,k) = Lc + Le + Lz

end do
end do
end do


!---------------------------------------------------------------------------------------------------
!    BOUNDARIES
!---------------------------------------------------------------------------------------------------

!-------------------
! i = 1 
!-------------------

if ( myback == mpi_proc_null)  then
   
   i = il + igp

   do k = kRHAD_sta , kRHAD_end
   do j = jRHAD_sta , jRHAD_end

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
      ! DIRECCION ETA
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
      ! DIRECCION ZETA
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
      ! DIRECCION CSI
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at i boundary

      boundary  = 1
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      ! dphi_dcsi to be updated in BCSignedDistanceFunction
      
      dphi_dcsi = dc2 * (   - one   * phi_AD( i+2 , j , k ) &
                          & + four  * phi_AD( i+1 , j , k ) &
                          & - three * phi_AD( i   , j , k )   )

      dphi_deta = de * ( phi_flux_plus_j - phi_flux_minus_j )
      dphi_dzet = dz * ( phi_flux_plus_k - phi_flux_minus_k )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(1,i,j,k)
      
      Lc = WaveSpeed * dphi_dcsi

      !-----------------------------------------------------
      ! RHS = U^1/J *  ∂Φ / ∂ξ                        +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !
      ! where ∂Φ / ∂ξ is obtained as a function of ∂Φ / ∂η 
      ! and ∂Φ / ∂ζ and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

!-------------------
! i = im
!-------------------

if ( myfront == mpi_proc_null )  then
   
   i = iu + igp

   do k = kRHAD_sta , kRHAD_end
   do j = jRHAD_sta , jRHAD_end

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
      ! DIRECCION ETA
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
      ! DIRECCION ZETA
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
      ! DIRECCION CSI
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at i boundary

      boundary  = 2
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      ! dphi_dcsi to be updated in BCSignedDistanceFunction
      
      dphi_dcsi = dc2 * (     one   * phi_AD( i-2 , j , k ) &
                          & - four  * phi_AD( i-1 , j , k ) &
                          & + three * phi_AD( i   , j , k )   )

      dphi_deta = de * ( phi_flux_plus_j - phi_flux_minus_j )
      dphi_dzet = dz * ( phi_flux_plus_k - phi_flux_minus_k )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(1,i,j,k)
      
      Lc = WaveSpeed * dphi_dcsi

      !-----------------------------------------------------
      ! RHS = U^1/J *  ∂Φ / ∂ξ                        +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !
      ! where ∂Φ / ∂ξ is obtained as a function of ∂Φ / ∂η 
      ! and ∂Φ / ∂ζ and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

!-------------------
! j = 1
!-------------------

if ( myleft == mpi_proc_null )  then
   
   j = jl + jgp

   do k = kRHAD_sta , kRHAD_end
   do i = iRHAD_sta , iRHAD_end

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
      ! DIRECCION CSI
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
      ! DIRECCION ZETA
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
      ! DIRECCION ETA
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at j boundary

      boundary  = 3
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      dphi_dcsi = dc * ( phi_flux_plus_i - phi_flux_minus_i )

      ! dphi_deta to be updated in BCSignedDistanceFunction
      
      dphi_deta = de2 * (   - one   * phi_AD( i , j+2 , k ) &
                          & + four  * phi_AD( i , j+1 , k ) &
                          & - three * phi_AD( i , j   , k )   )

      dphi_dzet = dz * ( phi_flux_plus_k - phi_flux_minus_k )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η                         +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !
      ! where ∂Φ / ∂η is obtained as a function of ∂Φ / ∂ξ 
      ! and ∂Φ / ∂ζ and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

!-------------------
! j = jm
!-------------------

if ( myright == mpi_proc_null )  then
   
   j = ju + jgp

   do k = kRHAD_sta , kRHAD_end
   do i = iRHAD_sta , iRHAD_end

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
      ! DIRECCION CSI
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
      ! DIRECCION ZETA
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
      ! DIRECCION ETA
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at j boundary

      boundary  = 4
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      dphi_dcsi = dc * ( phi_flux_plus_i - phi_flux_minus_i )

      ! dphi_deta to be updated in BCSignedDistanceFunction
      
      dphi_deta = de2 * (     one   * phi_AD( i , j-2 , k ) &
                          & - four  * phi_AD( i , j-1 , k ) &
                          & + three * phi_AD( i , j   , k )   )

      dphi_dzet = dz * ( phi_flux_plus_k - phi_flux_minus_k )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(2,i,j,k)
      
      Le = WaveSpeed * dphi_deta

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ∂Φ / ∂η                         +
      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
      !
      ! where ∂Φ / ∂η is obtained as a function of ∂Φ / ∂ξ 
      ! and ∂Φ / ∂ζ and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

!-------------------
! k = 1
!-------------------

if ( mydown == mpi_proc_null )  then
   
   k = kl + kgp

   do j = jRHAD_sta , jRHAD_end
   do i = iRHAD_sta , iRHAD_end

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
      ! DIRECCION CSI
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
      ! DIRECCION ETA
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
      ! DIRECCION ZET
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at k boundary

      boundary  = 5
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      dphi_dcsi = dc * ( phi_flux_plus_i - phi_flux_minus_i )
      dphi_deta = de * ( phi_flux_plus_j - phi_flux_minus_j )

      ! dphi_dzet to be updated in BCSignedDistanceFunction
      
      dphi_dzet = dz2 * (   - one   * phi_AD( i , j , k+2 ) &
                          & + four  * phi_AD( i , j , k+1 ) &
                          & - three * phi_AD( i , j , k   )   )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ
      !
      ! where ∂Φ / ∂ζ is obtained as a function of ∂Φ / ∂ξ 
      ! and ∂Φ / ∂η and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

!-------------------
! k = km
!-------------------

if ( myup == mpi_proc_null )  then
   
   k = ku + kgp

   do j = jRHAD_sta , jRHAD_end
   do i = iRHAD_sta , iRHAD_end

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
      ! DIRECCION CSI
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
      ! DIRECCION ETA
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
      ! DIRECCION ZET
      !------------------------------------------------------------------------

      ! We impose the condition ǁ ∇Φ ǁ = 1

      ! metrics at k boundary

      boundary  = 6
      
      csivecB   = csi(:,i,j,k)
      etavecB   = eta(:,i,j,k)
      zetvecB   = zet(:,i,j,k)
      
      dphi_dcsi = dc * ( phi_flux_plus_i - phi_flux_minus_i )
      dphi_deta = de * ( phi_flux_plus_j - phi_flux_minus_j )

      ! dphi_dzet to be updated in BCSignedDistanceFunction
      
      dphi_dzet = dz2 * (     one   * phi_AD( i , j , k-2 ) &
                          & - four  * phi_AD( i , j , k-1 ) &
                          & + three * phi_AD( i , j , k   )   )

      call BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                      dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                      boundary )

      WaveSpeed = ucn_j(3,i,j,k)
      
      Lz = WaveSpeed * dphi_dzet

      !-----------------------------------------------------
      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
      !       U^2/J * ( Φ_{j+1/2} - Φ_{j-1/2} ) / Δη  +
      !       U^3/J * ∂Φ / ∂ζ
      !
      ! where ∂Φ / ∂ζ is obtained as a function of ∂Φ / ∂ξ 
      ! and ∂Φ / ∂η and the metrics from the condition 
      ! ǁ ∇Φ ǁ = 1     
      !-----------------------------------------------------

      rightH_AD(i,j,k) = Lc + Le + Lz      

   end do
   end do

end if

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGES
! - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
!         i1,jm,km ---E11--- im,jm,km    
!             /|(8)           /|(7)                            
!     E12----/ |             /-|------E10                                   
!         i1,j1,km---E9---im,j1,km                   
!           |(5)           |(6)|                                   
!           | E8       E3  |   E7                                
!          E5  |       |   E6  |                                
!           | i1,jm,k1-----|-im,jm,k1           
!       E4--|-/(4)         |  /(3)                 
!           |/             | / -------E2
! ζ     i1,j1,k1---E1----im,j1,k1
! ^   η    (1)            (2)   
! |  7             
! | /
! |/
! +-----------> ξ

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
! NOTE: This implementations DOESN'T consider
!       parallelisation. 
!
! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 1
! - - - - - - - - - - - - - - - - - - - - - - - - 

!if ( myleft == mpi_proc_null)  then
!
!   j = jl + jgp
!   k = kl + kgp
!
!   do i = iRHAD_sta , iRHAD_end
!
!      bd_stencil_mb = 0
!      bd_stencil_mf = 0
!      
!      if( myback == mpi_proc_null  .and. i == iRHAD_sta ) bd_stencil_mb = 1
!      if( myfront == mpi_proc_null .and. i == iRHAD_end ) bd_stencil_mf = 1
!
!      !------------------------------------------------------------------------
!      ! DIRECCION CSI
!      !------------------------------------------------------------------------
!   
!      !-----------------
!      ! i+1/2
!      !-----------------
!   
!      WaveSpeed = ucn_j(1,i,j,k)
!   
!      phi0 = phi_AD ( i-1 , j , k )
!      phi1 = phi_AD ( i   , j , k )
!      phi2 = phi_AD ( i+1 , j , k )
!   
!      ! Outbounded stencil
!      if (bd_stencil_mf == 1) then 
!          
!          ! I impose positive WaveSpeed to force upwind stencil in 
!          ! GetWENO3Reconstruction 
!   
!          WaveSpeed = one  ! dummy value
!          phi3      = zero ! dummy value
!      
!      ! Not-outbounded stencil        
!      else 
!          phi3 = phi_AD ( i+2 , j , k )
!      end if        
!   
!      ! Φ_{i+1/2}
!      phi_flux_plus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
!   
!      !-----------------
!      ! i-1/2
!      !-----------------
!   
!      WaveSpeed = ucn_j(1,i,j,k)
!   
!      ! Outbounded stencil
!      if ( bd_stencil_mb == 1 ) then 
!          
!          ! I impose negative WaveSpeed to force downwind stencil in 
!          ! GetWENO3Reconstruction 
!   
!          WaveSpeed = - one  ! dummy value
!          phi0      =   zero ! dummy value
!      
!      ! Not-outbounded stencil
!      else 
!          phi0 = phi_AD ( i-2 , j , k )
!      end if        
!   
!      phi1 = phi_AD ( i-1 , j , k )
!      phi2 = phi_AD ( i   , j , k )
!      phi3 = phi_AD ( i+1 , j , k )
!   
!      ! Φ_{i-1/2}
!      phi_flux_minus_i = GetWENO3Reconstruction ( phi0, phi1, phi2, phi3, WaveSpeed )
!   
!      ! Reset WaveSpeed to its actual value
!      WaveSpeed = ucn_j(1,i,j,k)
!   
!      ! Lξ = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ
!      Lc = WaveSpeed * ( dc * ( phi_flux_plus_i - phi_flux_minus_i ) )
!
!      !------------------------------------------------------------------------
!      ! DIRECCION ETA
!      !------------------------------------------------------------------------
!
!
!      !------------------------------------------------------------------------
!      ! DIRECCION ZETA
!      !------------------------------------------------------------------------
!
!
!
!
!      !-----------------------------------------------------
!      ! RHS = U^1/J * ( Φ_{i+1/2} - Φ_{i-1/2} ) / Δξ  +
!      !       U^2/J * ∂Φ / ∂η                         +
!      !       U^3/J * ( Φ_{k+1/2} - Φ_{k-1/2} ) / Δζ
!      !
!      ! where ∂Φ / ∂η is obtained as a function of ∂Φ / ∂ξ 
!      ! and ∂Φ / ∂ζ and the metrics from the condition 
!      ! ǁ ∇Φ ǁ = 1     
!      !-----------------------------------------------------
!
!      rightH_AD(i,j,k) = Lc + Le + Lz      
!
!
!   end do
!
!end if



end subroutine calc_RH_AD_WENO3_test
