subroutine calc_h_gradient_ENO2( il,iu          ,&
                                 jl,ju          ,&
                                 kl,ku          ,&
                                 igp, jgp, kgp  ,&
                                 dc, de, dz     ,&
                                 csi            ,&
                                 eta            ,&
                                 zet            ,&
                                 h              ,&
                                 h_gradient      &
                                )

! Calculo de |grad(h0)| mediante ENO2, el cual es usado posteriormente para la correccion
! de Sussman. Tambien se calcula funcion S(h0).

use global_mpi
use global_app
use global_lsm, only : OrderReinitialisationBoundaries , ENOBCReinitialisation

use AdvectionMethods

implicit none


! Input variables

integer, intent(in) :: il,iu,jl,ju,kl,ku ! external nodes
integer, intent(in) :: igp, jgp, kgp ! # of ghostpoint at local processor 
real (kind = rdf), intent(in) :: dc,de,dz ! dcsi, deta, dzet
real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi , eta, zet ! metrics
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: h

! Output variable
real (kind = rdf), dimension(1:2,il:iu,jl:ju,kl:ku) , intent(inout) :: h_gradient


!local

!index
integer  :: i_mysta,i_myend
integer  :: j_mysta,j_myend
integer  :: k_mysta,k_myend

integer  :: ista , iend
integer  :: jsta , jend
integer  :: ksta , kend

integer :: i , j , k 

! Variables ENO2

real (kind = rdf) :: dh_dcsi, dh_deta, dh_dzet
real (kind = rdf) :: h0 , hLL , hL , hC , hR , hRR 
real (kind = rdf) :: dh_dx , dh_dy , dh_dz
real (kind = rdf) :: dc2,de2,dz2 ! ( 1/(2Δξ) , 1/(2Δη) , 1/(2Δζ) )

! dummy variables

real (kind = rdf) :: dummy

! Switches (logicals) para forzar stencils bounded en los bordes del dominio

logical :: BackOutbounded   , LeftOutbounded   , DownOutbounded 
logical :: FrontOutbounded  , RightOutbounded  , UpOutbounded   

integer :: iBiasDirection, jBiasDirection, kBiasDirection

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


dc2 = one_half * dc
de2 = one_half * de
dz2 = one_half * dz
    
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! INTERIOR DOMAIN LOOP
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

do k = k_mysta , k_myend
do j = j_mysta , j_myend
do i = i_mysta , i_myend

   BackOutbounded  = .false.
   LeftOutbounded  = .false.
   DownOutbounded  = .false.

   FrontOutbounded = .false.
   RightOutbounded = .false.
   UpOutbounded    = .false.

   if( myback  == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
   if( myleft  == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
   if( mydown  == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.

   if( myfront == mpi_proc_null  .and. i == i_myend ) FrontOutbounded = .true.
   if( myright == mpi_proc_null  .and. j == j_myend ) RightOutbounded = .true.
   if( myup    == mpi_proc_null  .and. k == k_myend ) UpOutbounded    = .true.

   if ( nblk/=0 ) then

      do nb = 1,nblk

         if ( i > li_blk_ia(1,nb) .and. i < li_blk_ib(1,nb) .and. & 
              j > li_blk_ja(1,nb) .and. j < li_blk_jb(1,nb) .and. &
              k > li_blk_ka(1,nb) .and. k < li_blk_kb(1,nb) ) then
   
              h_gradient(:,i,j,k) = zero
              cycle
   
         end if
   
         ! Blanking stencil biasing
         if( i==li_blk_ia(1,nb)-1 .and. &
             j>=li_blk_ja(1,nb)   .and. j<=li_blk_jb(1,nb) ) FrontOutbounded  = .true.
   
         if( i==li_blk_ib(1,nb)+1 .and. &
             j>=li_blk_ja(1,nb)   .and. j<=li_blk_jb(1,nb) ) BackOutbounded   = .true.
   
         if( j==li_blk_ja(1,nb)-1 .and. &
             i>=li_blk_ia(1,nb)   .and. i<=li_blk_ib(1,nb) ) RightOutbounded  = .true.
   
         if( j==li_blk_jb(1,nb)+1 .and. &
             i>=li_blk_ia(1,nb)   .and. i<=li_blk_ib(1,nb) ) LeftOutbounded   = .true.

      end do

   end if

   ! local ϕ0 value
   ! h0 = hzero(i,j,k)

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   hLL = zero ! dummy
   hL  = h( i-1 , j , k ) 
   hC  = h( i   , j , k )
   hR  = h( i+1 , j , k )
   hRR = zero ! dummy

   if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
   if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )

   dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                           BackOutbounded, FrontOutbounded )


   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   hLL = zero ! dummy
   hL  = h( i , j-1 , k )
   hC  = h( i , j   , k )
   hR  = h( i , j+1 , k )
   hRR = zero ! dummy

   if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
   if (.not. RightOutbounded ) hRR = h( i , j+2 , k )

   dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                           LeftOutbounded, RightOutbounded )


   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   hLL = zero ! dummy
   hL  = h( i , j , k-1 )
   hC  = h( i , j , k   )
   hR  = h( i , j , k+1 )
   hRR = zero ! dummy

   if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
   if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )

   dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                           DownOutbounded, UpOutbounded ) 


   ! Blanking boundary condition
   if (nblk /= 0) then

      do nb = 1,nblk

         if( i==li_blk_ia(1,nb) .and. &
             j>=li_blk_ja(1,nb) .and. j<=li_blk_jb(1,nb) ) dh_dcsi = zero
   
         if( i==li_blk_ib(1,nb) .and. &
             j>=li_blk_ja(1,nb) .and. j<=li_blk_jb(1,nb) ) dh_dcsi = zero
   
         if( j==li_blk_ja(1,nb) .and. &
             i>=li_blk_ia(1,nb) .and. i<=li_blk_ib(1,nb) ) dh_deta = zero
   
         if( j==li_blk_jb(1,nb) .and. &
             i>=li_blk_ia(1,nb) .and. i<=li_blk_ib(1,nb) ) dh_deta = zero

      end do
   
   end if


   ! Calculo del gradiente en coordenadas curvilineas

   dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
           eta( 1 , i , j , k ) * dh_deta + &
           zet( 1 , i , j , k ) * dh_dzet

   dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
           eta( 2 , i , j , k ) * dh_deta + &
           zet( 2 , i , j , k ) * dh_dzet

  ! ∇h_ijk

  h_gradient(1,i,j,k) = dh_dx
  h_gradient(2,i,j,k) = dh_dy

end do 
end do
end do

! Boundaries, edges and corners

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
   
   !i = il + igp
   i = ista

   iBiasDirection = 1

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      DownOutbounded  = .false.
   
      RightOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myleft == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
   
      if( myright == mpi_proc_null .and. j == j_myend ) RightOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == k_myend ) UpOutbounded    = .true.
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries       , &
                                       iBiasDirection                        , &
                                       ENOBCReinitialisation                     )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                              LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                              DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
      
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
   end do
   end do

end if

!-------------------
! i = im 
!-------------------


if ( myfront == mpi_proc_null)  then
   
   !i = iu - igp

   i = iend

   iBiasDirection = -1

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      DownOutbounded  = .false.
   
      RightOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myleft == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
   
      if( myright == mpi_proc_null .and. j == j_myend ) RightOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == k_myend ) UpOutbounded    = .true.
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                              DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do
   end do

end if


!-------------------
! j = 1 
!-------------------

if (myleft == mpi_proc_null) then

   !j = jl + jgp
   j = jsta

   jBiasDirection = 1
   
   do k = k_mysta , k_myend
   do i = i_mysta , i_myend
   
      BackOutbounded  = .false.
      DownOutbounded  = .false.
   
      FrontOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myback == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == i_myend ) FrontOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == k_myend ) UpOutbounded    = .true.
   
      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
   end do 
   end do

end if

!-------------------
! j = jm 
!-------------------

if (myright == mpi_proc_null) then

   !j = ju - jgp
    j = jend

   jBiasDirection = -1

   do k = k_mysta , k_myend
   do i = i_mysta , i_myend
   
      BackOutbounded  = .false.
      DownOutbounded  = .false.
   
      FrontOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myback == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == i_myend ) FrontOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == k_myend ) UpOutbounded    = .true.
   
      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                         h( i , j + jBiasDirection * 1 , k ) , &
                                         h( i , j + jBiasDirection * 2 , k ) , &
                                         h( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries     , &
                                         jBiasDirection                      , &
                                         ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
   end do 
   end do

end if

!-------------------
! k = 1 
!-------------------

if (mydown == mpi_proc_null) then

   !k = kl + kgp
   k = ksta

   kBiasDirection = 1

   do j = j_mysta , j_myend
   do i = i_mysta , i_myend
   
      BackOutbounded  = .false.
      LeftOutbounded  = .false.
   
      FrontOutbounded = .false.
      RightOutbounded = .false.
   
      if( myback == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myleft == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == i_myend ) FrontOutbounded = .true.
      if( myright == mpi_proc_null .and. j == j_myend ) RightOutbounded = .true.
   
      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
   end do 
   end do

end if

!-------------------
! k = km 
!-------------------

if (myup == mpi_proc_null) then

   !k = ku - kgp
   k = kend

   kBiasDirection = -1

   do j = j_mysta , j_myend
   do i = i_mysta , i_myend
   
      BackOutbounded  = .false.
      LeftOutbounded  = .false.
   
      FrontOutbounded = .false.
      RightOutbounded = .false.
   
      if( myback == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myleft == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == i_myend ) FrontOutbounded = .true.
      if( myright == mpi_proc_null .and. j == j_myend ) RightOutbounded = .true.
   
      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
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


!**************************
!
! i - free edges
!
!**************************

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 1
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myleft == mpi_proc_null .and. mydown == mpi_proc_null)  then

   !j = jl + jgp
   !k = kl + kgp

   j = jsta
   k = ksta

   jBiasDirection = 1
   kBiasDirection = 1

   do i = i_mysta , i_myend

      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == i_myend ) FrontOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz
   
   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 3
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myright == mpi_proc_null .and. mydown == mpi_proc_null)  then

   !j = ju - jgp
   !k = kl + kgp

   j = jend
   k = ksta

   jBiasDirection = -1
   kBiasDirection =  1

   do i = i_mysta , i_myend
   
      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == i_myend ) FrontOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 9
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myleft == mpi_proc_null .and. myup == mpi_proc_null)  then

   !j = jl + jgp
   !k = ku - kgp

   j = jsta
   k = kend

   jBiasDirection =  1
   kBiasDirection = -1

   do i = i_mysta , i_myend
 
      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == i_myend ) FrontOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 11
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myright == mpi_proc_null .and. myup == mpi_proc_null)  then

   !j = ju - jgp
   !k = ku - kgp

   j = jend
   k = kend

   jBiasDirection = -1
   kBiasDirection = -1

   do i = i_mysta , i_myend

      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == i_mysta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == i_myend ) FrontOutbounded = .true.
   
      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i-1 , j , k )
      hC  = h( i   , j , k )
      hR  = h( i+1 , j , k )
      hRR = zero ! dummy
   
      if (.not. BackOutbounded  ) hLL = h( i-2 , j , k )
      if (.not. FrontOutbounded ) hRR = h( i+2 , j , k )
   
      dh_dcsi = dc * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                              BackOutbounded, FrontOutbounded )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if

!**************************
!
! j - free edges
!
!**************************

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 2
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. mydown == mpi_proc_null)  then

   !i = iu - igp
   !k = kl + kgp

   i = iend
   k = ksta

   iBiasDirection = -1
   kBiasDirection =  1

   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == j_myend ) RightOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if

! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 4
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. mydown == mpi_proc_null)  then

   !i = il + jgp
   !k = kl + kgp

   i = ista
   k = ksta

   iBiasDirection =  1
   kBiasDirection =  1

   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == j_myend ) RightOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dcsi = dc2 * (   - one   * h( i+2 , j , k ) &
                        & + four  * h( i+1 , j , k ) &
                        & - three * h( i   , j , k )   )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------

      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 10
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myup == mpi_proc_null)  then

   !i = iu - jgp
   !k = ku - kgp

   i = iend
   k = kend

   iBiasDirection = -1
   kBiasDirection = -1

   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == j_myend ) RightOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 12
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myup == mpi_proc_null)  then

   !i = il + igp
   !k = ku - kgp

   i = ista
   k = kend

   iBiasDirection =  1
   kBiasDirection = -1

   do j = j_mysta , j_myend

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == j_mysta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == j_myend ) RightOutbounded = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j-1 , k )
      hC  = h( i , j   , k )
      hR  = h( i , j+1 , k )
      hRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) hLL = h( i , j-2 , k )
      if (.not. RightOutbounded ) hRR = h( i , j+2 , k )
   
      dh_deta = de * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            LeftOutbounded, RightOutbounded )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dzet = dz * BiasedDerivative( h( i , j , k + kBiasDirection * 0 ) , &
                                       h( i , j , k + kBiasDirection * 1 ) , &
                                       h( i , j , k + kBiasDirection * 2 ) , &
                                       h( i , j , k + kBiasDirection * 3 ) , &
                                       OrderReinitialisationBoundaries     , &
                                       kBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if

!**************************
!
! k-free edges
!
!**************************


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 5
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myleft == mpi_proc_null)  then

   !i = il + igp
   !j = jl + jgp

   i = ista
   j = jsta

   iBiasDirection =  1
   jBiasDirection =  1

   do k = k_mysta , k_myend

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == k_myend ) UpOutbounded    = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 6
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null)  then

   !i = iu - igp
   !j = jl + jgp

   i = iend
   j = jsta

   iBiasDirection = -1
   jBiasDirection =  1

   do k = k_mysta , k_myend

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == k_myend ) UpOutbounded    = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 7
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null)  then

   !i = iu - igp
   !j = ju - jgp

   i = iend
   j = jend

   iBiasDirection = -1
   jBiasDirection = -1

   do k = k_mysta , k_myend

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == k_myend ) UpOutbounded    = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if


! - - - - - - - - - - - - - - - - - - - - - - - - 
! EDGE 8
! - - - - - - - - - - - - - - - - - - - - - - - - 

if ( myback == mpi_proc_null .and. myright == mpi_proc_null)  then

   !i = il + igp
   !j = ju - jgp

   i = ista
   j = jend

   iBiasDirection =  1
   jBiasDirection = -1

   do k = k_mysta , k_myend

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == k_mysta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == k_myend ) UpOutbounded    = .true.

      ! local ϕ0 value
      !h0 = hzero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dh_dcsi = dc * BiasedDerivative( h( i + iBiasDirection * 0 , j , k ) , &
                                       h( i + iBiasDirection * 1 , j , k ) , &
                                       h( i + iBiasDirection * 2 , j , k ) , &
                                       h( i + iBiasDirection * 3 , j , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       iBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dh_deta = de * BiasedDerivative( h( i , j + jBiasDirection * 0 , k ) , &
                                       h( i , j + jBiasDirection * 1 , k ) , &
                                       h( i , j + jBiasDirection * 2 , k ) , &
                                       h( i , j + jBiasDirection * 3 , k ) , &
                                       OrderReinitialisationBoundaries     , &
                                       jBiasDirection                      , &
                                       ENOBCReinitialisation                   )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      hLL = zero ! dummy
      hL  = h( i , j , k-1 )
      hC  = h( i , j , k   )
      hR  = h( i , j , k+1 )
      hRR = zero ! dummy
   
      if (.not. DownOutbounded  ) hLL = h( i , j , k-2 )
      if (.not. UpOutbounded    ) hRR = h( i , j , k+2 )
   
      dh_dzet = dz * GetENO2Reconstruction( hC, hLL, hL, hC, hR, hRR , &
                                            DownOutbounded, UpOutbounded ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dh_dx = csi( 1 , i , j , k ) * dh_dcsi + &
              eta( 1 , i , j , k ) * dh_deta + &
              zet( 1 , i , j , k ) * dh_dzet
   
      dh_dy = csi( 2 , i , j , k ) * dh_dcsi + &
              eta( 2 , i , j , k ) * dh_deta + &
              zet( 2 , i , j , k ) * dh_dzet
   
      !dh_dz = csi( 3 , i , j , k ) * dh_dcsi + &
      !        eta( 3 , i , j , k ) * dh_deta + &
      !        zet( 3 , i , j , k ) * dh_dzet
   
      ! ∇h_ijk
      h_gradient(1,i,j,k) = dh_dx
      h_gradient(2,i,j,k) = dh_dy
      !h_gradient(3,i,j,k) = dh_dz

   end do

end if

! Update the ghost nodes
call rhs_exchng3_4d( h_gradient )

contains

include 'rhs_exchng3_4d.F90'


end subroutine calc_h_gradient_ENO2
