subroutine calc_RH_AD_WENO3( phi_AD , rightH_AD , AdvectionNodes )
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
integer  :: i_mysta , i_myend
integer  :: j_mysta , j_myend
integer  :: k_mysta , k_myend

integer  :: ista , iend
integer  :: jsta , jend
integer  :: ksta , kend

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

logical :: BlankingFlag

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

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!Nodos incluyendo el borde
i_mysta = il + igp
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp        
j_myend = ju - jgp
k_myend = ku - kgp

!Nodos sin incluir borde (el borde se actualiza en bcond_lsm.F90)

if (myback == mpi_proc_null)  i_mysta = il + igp + 1
if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

if (myfront == mpi_proc_null) i_myend = iu - igp - 1
if (myright == mpi_proc_null) j_myend = ju - jgp - 1
if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

! Physical boundaries of the processor
   
ista = il ; jsta = jl ; ksta = kl 
iend = iu ; jend = ju ; kend = ku 

if ( myback  == mpi_proc_null )  ista = il + igp 
if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

if ( myfront == mpi_proc_null )  iend = iu - igp
if ( myright == mpi_proc_null )  jend = ju - jgp
if ( myup    == mpi_proc_null )  kend = ku - kgp


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! INTERIOR DOMAIN LOOP
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!print *, 'max phi   = ', maxval(phi_AD)
!print *, 'min phi   = ', minval(phi_AD)
!print *, ' '

do i = i_mysta , i_myend
do j = j_mysta , j_myend
do k = k_mysta , k_myend

   if ( AdvectionNodes(i,j,k) == 0  ) cycle

   ! Definicion de los booleans para stencils bounded en los bordes

   bd_stencil_mb = 0
   bd_stencil_ml = 0
   bd_stencil_md = 0

   bd_stencil_mf = 0
   bd_stencil_mr = 0
   bd_stencil_mu = 0
   
   if( myback == mpi_proc_null  .and. i == i_mysta ) bd_stencil_mb = 1
   if( myleft == mpi_proc_null  .and. j == j_mysta ) bd_stencil_ml = 1
   if( mydown == mpi_proc_null  .and. k == k_mysta ) bd_stencil_md = 1

   if( myfront == mpi_proc_null .and. i == i_myend ) bd_stencil_mf = 1
   if( myright == mpi_proc_null .and. j == j_myend ) bd_stencil_mr = 1
   if( myup    == mpi_proc_null .and. k == k_myend ) bd_stencil_mu = 1

   ! Blanking stencil biasing
   BlankingFlag = .false.

   if ( nblk/=0 ) then

      do nb = 1,nblk

         if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. & 
              j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
              k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
   
              BlankingFlag = .true.
   
         end if
   
         ! Blanking stencil biasing
         if( i==li_blk_ia(1,nb)-1 .and. &
             j>=li_blk_ja(1,nb)   .and. j<=li_blk_jb(1,nb) ) bd_stencil_mf = 1
   
         if( i==li_blk_ib(1,nb)+1 .and. &
             j>=li_blk_ja(1,nb)   .and. j<=li_blk_jb(1,nb) ) bd_stencil_mb = 1
   
         if( j==li_blk_ja(1,nb)-1 .and. &
             i>=li_blk_ia(1,nb)   .and. i<=li_blk_ib(1,nb) ) bd_stencil_mr = 1
   
         if( j==li_blk_jb(1,nb)+1 .and. &
             i>=li_blk_ia(1,nb)   .and. i<=li_blk_ib(1,nb) ) bd_stencil_ml = 1
   
      end do

   end if

   ! Cycle if I'm within the Blanking region
   if ( BlankingFlag ) cycle

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

end do
end do
end do



end subroutine calc_RH_AD_WENO3
