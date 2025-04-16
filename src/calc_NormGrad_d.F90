subroutine calc_NormGrad_d( phi , phizero , NormGrad_RN )

! Calculo de |grad(phi0)| mediante ENO2, el cual es usado posteriormente para la correccion
! de Sussman. Tambien se calcula funcion S(phi0).

use AdvectionMethods

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in)  :: phi
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in)  :: phizero
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out) :: NormGrad_RN

!local

!index
integer  :: iRN_sta,iRN_end
integer  :: jRN_sta,jRN_end
integer  :: kRN_sta,kRN_end

integer :: i , j , k 

! Variables ENO2

real (kind = rdf) :: dphi_dcsi, dphi_deta, dphi_dzet
real (kind = rdf) :: phi0 , phiLL , phiL , phiC , phiR , phiRR 
real (kind = rdf) :: dphi_dx , dphi_dy , dphi_dz
real (kind = rdf) :: dc2,de2,dz2 ! ( 1/(2Δξ) , 1/(2Δη) , 1/(2Δζ) )

! Switches (logicals) para forzar stencils bounded en los bordes del dominio

logical :: BackOutbounded   , LeftOutbounded   , DownOutbounded 
logical :: FrontOutbounded  , RightOutbounded  , UpOutbounded   

integer :: iBiasDirection, jBiasDirection, kBiasDirection

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

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

dc2 = one_half * dc
de2 = one_half * de
dz2 = one_half * dz

NormGrad_RN = zero
    
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
! INTERIOR DOMAIN LOOP
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

do k = kRN_sta , kRN_end
do j = jRN_sta , jRN_end
do i = iRN_sta , iRN_end

   BackOutbounded  = .false.
   LeftOutbounded  = .false.
   DownOutbounded  = .false.

   FrontOutbounded = .false.
   RightOutbounded = .false.
   UpOutbounded    = .false.

   if( myback  == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
   if( myleft  == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
   if( mydown  == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.

   if( myfront == mpi_proc_null  .and. i == iRN_end ) FrontOutbounded = .true.
   if( myright == mpi_proc_null  .and. j == jRN_end ) RightOutbounded = .true.
   if( myup    == mpi_proc_null  .and. k == kRN_end ) UpOutbounded    = .true.

   ! local ϕ0 value
   phi0 = phizero(i,j,k)

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   phiLL = zero ! dummy
   phiL  = phi( i-1 , j , k )
   phiC  = phi( i   , j , k )
   phiR  = phi( i+1 , j , k )
   phiRR = zero ! dummy

   if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
   if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )

   dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                           BackOutbounded, FrontOutbounded          )


   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   phiLL = zero ! dummy
   phiL  = phi( i , j-1 , k )
   phiC  = phi( i , j   , k )
   phiR  = phi( i , j+1 , k )
   phiRR = zero ! dummy

   if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
   if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )

   dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                           LeftOutbounded, RightOutbounded          )


   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   phiLL = zero ! dummy
   phiL  = phi( i , j , k-1 )
   phiC  = phi( i , j , k   )
   phiR  = phi( i , j , k+1 )
   phiRR = zero ! dummy

   if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
   if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )

   dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                           DownOutbounded, UpOutbounded             ) 

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |∇(Φ)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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
   
   i = il + igp

   iBiasDirection = 1

   do k = kRN_sta , kRN_end
   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      DownOutbounded  = .false.
   
      RightOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myleft == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
   
      if( myright == mpi_proc_null .and. j == jRN_end ) RightOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == kRN_end ) UpOutbounded    = .true.
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
   end do
   end do

end if

!-------------------
! i = im 
!-------------------


if ( myfront == mpi_proc_null)  then
   
   i = iu - igp

   iBiasDirection = -1

   do k = kRN_sta , kRN_end
   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      DownOutbounded  = .false.
   
      RightOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myleft == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
   
      if( myright == mpi_proc_null .and. j == jRN_end ) RightOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == kRN_end ) UpOutbounded    = .true.
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

   end do
   end do

end if


!-------------------
! j = 1 
!-------------------

if (myleft == mpi_proc_null) then

   j = jl + jgp

   jBiasDirection = 1
   
   do k = kRN_sta , kRN_end
   do i = iRN_sta , iRN_end
   
      BackOutbounded  = .false.
      DownOutbounded  = .false.
   
      FrontOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myback == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == iRN_end ) FrontOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == kRN_end ) UpOutbounded    = .true.
   
      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
   end do 
   end do

end if

!-------------------
! j = jm 
!-------------------

if (myright == mpi_proc_null) then

   j = ju - jgp
   
   jBiasDirection = -1

   do k = kRN_sta , kRN_end
   do i = iRN_sta , iRN_end
   
      BackOutbounded  = .false.
      DownOutbounded  = .false.
   
      FrontOutbounded = .false.
      UpOutbounded    = .false.
   
      if( myback == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( mydown == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == iRN_end ) FrontOutbounded = .true.
      if( myup    == mpi_proc_null .and. k == kRN_end ) UpOutbounded    = .true.
   
      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
   end do 
   end do

end if

!-------------------
! k = 1 
!-------------------

if (mydown == mpi_proc_null) then

   k = kl + kgp

   kBiasDirection = 1

   do j = jRN_sta , jRN_end
   do i = iRN_sta , iRN_end
   
      BackOutbounded  = .false.
      LeftOutbounded  = .false.
   
      FrontOutbounded = .false.
      RightOutbounded = .false.
   
      if( myback == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myleft == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == iRN_end ) FrontOutbounded = .true.
      if( myright == mpi_proc_null .and. j == jRN_end ) RightOutbounded = .true.
   
      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
   end do 
   end do

end if

!-------------------
! k = km 
!-------------------

if (myup == mpi_proc_null) then

   k = ku - kgp

   kBiasDirection = -1

   do j = jRN_sta , jRN_end
   do i = iRN_sta , iRN_end
   
      BackOutbounded  = .false.
      LeftOutbounded  = .false.
   
      FrontOutbounded = .false.
      RightOutbounded = .false.
   
      if( myback == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myleft == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
   
      if( myfront == mpi_proc_null .and. i == iRN_end ) FrontOutbounded = .true.
      if( myright == mpi_proc_null .and. j == jRN_end ) RightOutbounded = .true.
   
      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
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

   j = jl + jgp
   k = kl + kgp

   jBiasDirection = 1
   kBiasDirection = 1

   do i = iRN_sta , iRN_end

      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == iRN_end ) FrontOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   
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

   do i = iRN_sta , iRN_end
   
      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == iRN_end ) FrontOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do i = iRN_sta , iRN_end
 
      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == iRN_end ) FrontOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do i = iRN_sta , iRN_end

      BackOutbounded  = .false.
      FrontOutbounded = .false.
   
      if( myback  == mpi_proc_null  .and. i == iRN_sta ) BackOutbounded  = .true.
      if( myfront == mpi_proc_null  .and. i == iRN_end ) FrontOutbounded = .true.
   
      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i-1 , j , k )
      phiC  = phi( i   , j , k )
      phiR  = phi( i+1 , j , k )
      phiRR = zero ! dummy
   
      if (.not. BackOutbounded  ) phiLL = phi( i-2 , j , k )
      if (.not. FrontOutbounded ) phiRR = phi( i+2 , j , k )
   
      dphi_dcsi = dc * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              BackOutbounded, FrontOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   i = iu - igp
   k = kl + kgp

   iBiasDirection = -1
   kBiasDirection =  1

   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == jRN_end ) RightOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == jRN_end ) RightOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc2 * (   - one   * phi( i+2 , j , k ) &
                          & + four  * phi( i+1 , j , k ) &
                          & - three * phi( i   , j , k )   )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == jRN_end ) RightOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do j = jRN_sta , jRN_end

      LeftOutbounded  = .false.
      RightOutbounded = .false.
   
      if( myleft  == mpi_proc_null  .and. j == jRN_sta ) LeftOutbounded  = .true.
      if( myright == mpi_proc_null  .and. j == jRN_end ) RightOutbounded = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)

      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j-1 , k )
      phiC  = phi( i , j   , k )
      phiR  = phi( i , j+1 , k )
      phiRR = zero ! dummy
   
      if (.not. LeftOutbounded  ) phiLL = phi( i , j-2 , k )
      if (.not. RightOutbounded ) phiRR = phi( i , j+2 , k )
   
      dphi_deta = de * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              LeftOutbounded, RightOutbounded          )
   
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         phi( i , j , k + kBiasDirection * 3 ) , &
                                         OrderReinitialisationBoundaries          , &
                                         kBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   i = il + igp
   j = jl + jgp

   iBiasDirection =  1
   jBiasDirection =  1

   do k = kRN_sta , kRN_end

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == kRN_end ) UpOutbounded    = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do k = kRN_sta , kRN_end

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == kRN_end ) UpOutbounded    = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------

      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do k = kRN_sta , kRN_end

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == kRN_end ) UpOutbounded    = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   do k = kRN_sta , kRN_end

      DownOutbounded  = .false.
      UpOutbounded    = .false.
   
      if( mydown  == mpi_proc_null  .and. k == kRN_sta ) DownOutbounded  = .true.
      if( myup    == mpi_proc_null  .and. k == kRN_end ) UpOutbounded    = .true.

      ! local ϕ0 value
      phi0 = phizero(i,j,k)
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------

      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         phi( i + iBiasDirection * 3 , j , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         iBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         phi( i , j + jBiasDirection * 3 , k ) , &
                                         OrderReinitialisationBoundaries          , &
                                         jBiasDirection                           , &
                                         ENOBCReinitialisation                        )

      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      phiLL = zero ! dummy
      phiL  = phi( i , j , k-1 )
      phiC  = phi( i , j , k   )
      phiR  = phi( i , j , k+1 )
      phiRR = zero ! dummy
   
      if (.not. DownOutbounded  ) phiLL = phi( i , j , k-2 )
      if (.not. UpOutbounded    ) phiRR = phi( i , j , k+2 )
   
      dphi_dzet = dz * GetENO2Reconstruction( phi0, phiLL, phiL, phiC, phiR, phiRR , &
                                              DownOutbounded, UpOutbounded             ) 
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      ! |grad(phi)| 
   
      NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

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

   iBiasDirection =  1
   jBiasDirection =  1
   kBiasDirection =  1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )


   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

end if


! VERTEX 2

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null &
                              .and. mydown == mpi_proc_null    )  then

   
   i = iu - igp
   j = jl + jgp
   k = kl + kgp

   iBiasDirection = -1
   jBiasDirection =  1
   kBiasDirection =  1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

end if

! VERTEX 3

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null &
                              .and. mydown  == mpi_proc_null    )  then

   i = iu - igp
   j = ju - jgp
   k = kl + kgp

   iBiasDirection = -1
   jBiasDirection = -1
   kBiasDirection =  1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )


end if

! VERTEX 4

if ( myback == mpi_proc_null  .and. myright == mpi_proc_null &
                              .and. mydown  == mpi_proc_null    )  then

   i = il + igp
   j = ju - jgp
   k = kl + kgp

   iBiasDirection =  1
   jBiasDirection = -1
   kBiasDirection =  1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

end if


! VERTEX 5

if ( myback == mpi_proc_null .and. myleft == mpi_proc_null &
                             .and. myup   == mpi_proc_null    )  then

   i = il + igp
   j = jl + jgp
   k = ku - kgp

   iBiasDirection =  1
   jBiasDirection =  1
   kBiasDirection = -1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )


   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )
   

end if

! VERTEX 6

if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null &
                              .and. myup   == mpi_proc_null    )  then
   
   i = iu - igp
   j = jl + jgp
   k = ku - kgp

   iBiasDirection = -1
   jBiasDirection =  1
   kBiasDirection = -1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )


end if

! VERTEX 7

if ( myfront == mpi_proc_null .and. myright == mpi_proc_null &
                              .and. myup    == mpi_proc_null    )  then

   i = iu - igp
   j = ju - jgp
   k = ku - kgp
   
   iBiasDirection = -1
   jBiasDirection = -1
   kBiasDirection = -1

   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )


end if

! VERTEX 8

if ( myback == mpi_proc_null  .and. myright == mpi_proc_null &
                              .and. myup    == mpi_proc_null    )  then

   i = il + igp
   j = ju - jgp
   k = ku - kgp

   iBiasDirection =  1
   jBiasDirection = -1
   kBiasDirection = -1
   
   !--------------------------------------------------------------------
   ! ξ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                      phi( i + iBiasDirection * 1 , j , k ) , &
                                      phi( i + iBiasDirection * 2 , j , k ) , &
                                      phi( i + iBiasDirection * 3 , j , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      iBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! η - DIRECTION
   !--------------------------------------------------------------------

   dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                      phi( i , j + jBiasDirection * 1 , k ) , &
                                      phi( i , j + jBiasDirection * 2 , k ) , &
                                      phi( i , j + jBiasDirection * 3 , k ) , &
                                      OrderReinitialisationBoundaries          , &
                                      jBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   !--------------------------------------------------------------------
   ! ζ - DIRECTION
   !--------------------------------------------------------------------

   dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                      phi( i , j , k + kBiasDirection * 1 ) , &
                                      phi( i , j , k + kBiasDirection * 2 ) , &
                                      phi( i , j , k + kBiasDirection * 3 ) , &
                                      OrderReinitialisationBoundaries          , &
                                      kBiasDirection                           , &
                                      ENOBCReinitialisation                        )

   ! Calculo del gradiente en coordenadas curvilineas

   dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
             eta( 1 , i , j , k ) * dphi_deta + &
             zet( 1 , i , j , k ) * dphi_dzet

   dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
             eta( 2 , i , j , k ) * dphi_deta + &
             zet( 2 , i , j , k ) * dphi_dzet

   dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
             eta( 3 , i , j , k ) * dphi_deta + &
             zet( 3 , i , j , k ) * dphi_dzet

   ! |grad(phi)| 

   NormGrad_RN(i,j,k) = sqrt( (dphi_dx)**2 + (dphi_dy)**2 + (dphi_dz)**2 )

end if


end subroutine calc_NormGrad_d
