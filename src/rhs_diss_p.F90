
subroutine rhs_diss_p
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate dissipation for continutiy equation from fourth-order
  ! derivative of the pressure in all three directions

  ! input
  !     pdiss_coef (from input file) (global_param)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     q(1,ijk)
  !     dtau(ijk)

  ! output
  !     diss(1,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global_app
  implicit none

  ! local variables
  real (kind = rdf) :: tmp1, tmp2, tmp3

  real (kind = rdf) :: d2p_dcsi2 , d2p_deta2 , d2p_dzet2


  real (kind = rdf) :: dcsq, desq, dzsq
  real (kind = rdf) :: g11, g22, g33
  !real (kind = rdf) :: rsignc ! 1 if the nodes for a central scheme are
                              ! all within the water phase. 0 otherwise
  
  real (kind = rdf) :: thickness , damping

  integer :: ista , iend , jsta , jend , ksta , kend
  integer :: i,j,k
  integer :: im,jm,km,ip,jp,kp

  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 

  if ( myback  == mpi_proc_null )  ista = il + igp 
  if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
  if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

  if ( myfront == mpi_proc_null )  iend = iu - igp
  if ( myright == mpi_proc_null )  jend = ju - jgp
  if ( myup    == mpi_proc_null )  kend = ku - kgp

  fv   = zero
  diss = zero

  dcsq = dc * dc
  desq = de * de
  dzsq = dz * dz

  ! inner second derivatives and scaling

  ! csi-derivatives
  ! this scheme requires two ghost points on each side
  !
!   do k = k_mysta - 1, k_myend + 1
!   do j = j_mysta - 1, j_myend + 1
!   do i = i_mysta - 1, i_myend + 1

  ! use ista, iend to keep from having
  ! a array-bound access error when process
  ! owns the boundary

  ! ! I apply a damping function to the pressure directly
  !  do k = ksta, kend
  !  do j = jsta, jend
  !  do i = ista, iend
   
  !    ! I apply a Heavised-Sigmoid type function to make the dissipation zero
  !    ! at the free-surface. I use as the tickness, twice the minimum between
  !    ! dx, dy, dz
     
  !    im = max( i-1 , ista )
  !    jm = max( j-1 , jsta )
  !    km = max( k-1 , ksta )

  !    ip = min( i+1 , iend )
  !    jp = min( j+1 , jend )
  !    kp = min( k+1 , kend )

  !    thickness = min( abs ( ( z( i  , j  , kp ) -  z( i  , j  , km ) ) ) , &
  !                     abs ( ( y( i  , jp , k  ) -  y( i  , jm , k  ) ) ) , &
  !                     abs ( ( x( ip , j  , k  ) -  x( im , j  , k  ) ) )   &
  !                   ) 
    
  !    damping   =  heaviside( phi(i,j,k) , thickness ) 
       
  !    q(1,i,j,k) = damping * q(1,i,j,k)
        
  !  end do
  !  end do
  !  end do


   ista = i_mysta - 1
   iend = i_myend + 1

   if ( myback  == mpi_proc_null ) ista = i_mysta
   if ( myfront == mpi_proc_null ) iend = i_myend

   do k = k_mysta-1 , k_myend+1
   do j = j_mysta-1 , j_myend+1
   do i = ista      , iend
  
      g11   = csi(1,i,j,k) * csi(1,i,j,k) + &
              csi(2,i,j,k) * csi(2,i,j,k) + &
              csi(3,i,j,k) * csi(3,i,j,k)

      ! tmp1 = g11 * dtau / J
      tmp1 = g11 * dtau(i,j,k) / aj(i,j,k) 

      ! ∂2p / ∂ξ2
      d2p_dcsi2 = dcsq * (q(1,i+1,j,k) - two * q(1,i,j,k) + q(1,i-1,j,k))
      
      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      ! I'm testi what happens if I use the total pressue into the 
      ! dissipation term
      !
      ! d2p_dcsi2 = dcsq * (         ( q(1,i+1,j,k) + h(i+1,j,k) ) - &
      !                        two * ( q(1,i,j,k)   + h(i,j,k)   ) + &
      !                              ( q(1,i-1,j,k) + h(i-1,j,k) ) )
      !
      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

      fv(1,i,j,k) = tmp1 * d2p_dcsi2
   
   end do
   end do
   end do

   jsta = j_mysta - 1
   jend = j_myend + 1
   
   if ( myleft  == mpi_proc_null ) jsta = j_mysta
   if ( myright == mpi_proc_null ) jend = j_myend

   do k = k_mysta-1 , k_myend+1
   do j = jsta      , jend
   do i = i_mysta-1 , i_myend+1

      g22   = eta(1,i,j,k) * eta(1,i,j,k) + &
              eta(2,i,j,k) * eta(2,i,j,k) + &
              eta(3,i,j,k) * eta(3,i,j,k)

      tmp2      = g22 * dtau(i,j,k) / aj(i,j,k) 

      ! ∂2p / ∂η2
      d2p_deta2 = desq * (q(1,i,j+1,k) - two * q(1,i,j,k) + q(1,i,j-1,k)) 
      
      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      ! I'm testi what happens if I use the total pressue into the 
      ! dissipation term
      !
      ! d2p_deta2 = dcsq * (         ( q(1,i,j+1,k) + h(i,j+1,k) ) - &
      !                        two * ( q(1,i,j,k)   + h(i,j,k)   ) + &
      !                              ( q(1,i,j-1,k) + h(i,j-1,k) ) )

      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

      fv(2,i,j,k) = tmp2 * d2p_deta2
  
   end do
   end do
   end do


   ksta = k_mysta - 1
   kend = k_myend + 1
  
   if ( mydown == mpi_proc_null ) ksta = k_mysta
   if ( myup   == mpi_proc_null ) kend = k_myend

   do k = ksta      , kend
   do j = j_mysta-1 , j_myend+1
   do i = i_mysta-1 , i_myend+1
   
      g33   = zet(1,i,j,k) * zet(1,i,j,k) + &
              zet(2,i,j,k) * zet(2,i,j,k) + &
              zet(3,i,j,k) * zet(3,i,j,k)
 
      tmp3      = g33 * dtau(i,j,k) / aj(i,j,k) 

      ! ∂2p / ∂ζ2
      !d2p_dzet2 = dzsq * (q(1,i,j,k+1) - two * q(1,i,j,k) + q(1,i,j,k-1))

      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      ! I'm testi what happens if I use the total pressue into the 
      ! dissipation term
      !
      ! d2p_dzet2 = dcsq * (         ( q(1,i,j,k+1) + h(i,j,k+1) ) - &
      !                        two * ( q(1,i,j,k)   + h(i,j,k)   ) + &
      !                              ( q(1,i,j,k-1) + h(i,j,k-1) ) )

      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

      fv(3,i,j,k) = tmp3 * d2p_dzet2
  
   end do
   end do
   end do

   ! boundary values
   ! 

   ! Put valid axu for nodes on domain (not process) boundary
   ! 
   if ( myback  == mpi_proc_null ) fv(1,i_mysta-1,:,:) = fv(1,i_mysta,:,:)
   if ( myfront == mpi_proc_null ) fv(1,i_myend+1,:,:) = fv(1,i_myend,:,:)
   if ( myleft  == mpi_proc_null ) fv(2,:,j_mysta-1,:) = fv(2,:,j_mysta,:)
   if ( myright == mpi_proc_null ) fv(2,:,j_myend+1,:) = fv(2,:,j_myend,:)
   if ( mydown  == mpi_proc_null ) fv(3,:,:,k_mysta-1) = fv(3,:,:,k_mysta)
   if ( myup    == mpi_proc_null ) fv(3,:,:,k_myend+1) = fv(3,:,:,k_myend)

   !call outputD3_real(fv(1,:,:,:),'disscsi_bmod')
   !if(myid == 0) print *, 'output_disscsi'

   ! boundary treatment of blanking nodes
   if (nblk /= 0) then
      
      do nb = 1, nblk

        if (blktype(1,nb,myzone) == 0) then
      
          i = li_blk_ia(n,nb)
     
          do k = li_blk_ka(n,nb) , li_blk_kb(n,nb)
          do j = li_blk_ja(n,nb) , li_blk_jb(n,nb)
            fv(1,i,j,k) = fv(1, max(i-1,ista) ,j,k)
          end do
          end do
      
        end if

        if (blktype(2,nb,myzone) == 0) then
      
          i = li_blk_ib(n,nb)
     
          do k = li_blk_ka(n,nb) , li_blk_kb(n,nb)
          do j = li_blk_ja(n,nb) , li_blk_jb(n,nb)
            fv(1,i,j,k) = fv(1, min(i+1,iend) ,j,k)
          end do
          end do
     
        end if

        if (blktype(3,nb,myzone) == 0) then
          
          j = li_blk_ja(n,nb)
          
          do k = li_blk_ka(n,nb) , li_blk_kb(n,nb)
          do i = li_blk_ia(n,nb) , li_blk_ib(n,nb)
            fv(2,i,j,k) = fv(2,i, max(j-1,jsta) ,k)
          end do
          end do
        
        end if

        if (blktype(4,nb,myzone) == 0) then
          
          j = li_blk_jb(n,nb)
          
          do k = li_blk_ka(n,nb) , li_blk_kb(n,nb)
          do i = li_blk_ia(n,nb) , li_blk_ib(n,nb)
            fv(2,i,j,k) = fv(2,i, min(j+1,jend) ,k)
          end do
          end do
        
        end if

        if (blktype(5,nb,myzone) == 0) then
        
          k = li_blk_ka(n,nb)
          
          do j = li_blk_ja(n,nb) , li_blk_jb(n,nb)
          do i = li_blk_ia(n,nb) , li_blk_ib(n,nb)
            fv(3,i,j,k) = fv(3,i,j, max(k-1,ksta) )
          end do
          end do
        
        end if

        if (blktype(6,nb,myzone) == 0) then
          
          k = li_blk_kb(n,nb)
          
          do j = li_blk_ja(n,nb) , li_blk_jb(n,nb)
          do i = li_blk_ia(n,nb) , li_blk_ib(n,nb)
            fv(3,i,j,k) = fv(3,i,j,min(k+1,kend))
          end do
          end do
        
        end if
    
      end do
    
    end if

   ! I copy the dissipation at the first air layer
   ! TO DO: maybe there's a more efficient way to do this by including 
   ! this computations within some of the other loops

!   do k = k_mysta, k_myend
!   do j = j_mysta, j_myend
!   do i = i_mysta, i_myend
!
!
!      fv(1,i,j,k) =                          rsign(i  ,j,k) * fv(1,i  ,j,k) &
!                    + (one - rsign(i,j,k)) * rsign(i-1,j,k) * fv(1,i-1,j,k) &
!                    + (one - rsign(i,j,k)) * rsign(i+1,j,k) * fv(1,i+1,j,k) 
!
!
!      fv(2,i,j,k) =                          rsign(i,j  ,k) * fv(2,i,j  ,k) &
!                    + (one - rsign(i,j,k)) * rsign(i,j-1,k) * fv(2,i,j-1,k) &
!                    + (one - rsign(i,j,k)) * rsign(i,j+1,k) * fv(2,i,j+1,k) 
!
!
!      fv(3,i,j,k) =                          rsign(i,j,k  ) * fv(3,i,j,k  ) &
!                    + (one - rsign(i,j,k)) * rsign(i,j,k-1) * fv(3,i,j,k-1) &
!                    + (one - rsign(i,j,k)) * rsign(i,j,k+1) * fv(3,i,j,k+1) 
!
!   end do
!   end do
!   end do

  ! Clearing out dissipation fluxes
  !fv(1,:,:,:) = fv(1,:,:,:) * rsign(:,:,:) 
  !fv(2,:,:,:) = fv(2,:,:,:) * rsign(:,:,:) 
  !fv(3,:,:,:) = fv(3,:,:,:) * rsign(:,:,:) 

  !! Update ghost nodes
  !call rhs_exchng3_4d( fv )
  !! I extrapolate the dissipative flux to the first air layer:
  !call rhs_ghost_fluid_nodes_extp_4d ( fv )
  !! Update ghost nodes
  !call rhs_exchng3_4d( fv )

  ! outer second derivative, interior nodes only
  !

   ! This loop is to apply a heaviside damping near the free-surface
   do k = k_mysta, k_myend
   do j = j_mysta, j_myend
   do i = i_mysta, i_myend
   
     ! I apply a Heavised-Sigmoid type function to make the dissipation zero
     ! at the free-surface. I use as the tickness, twice the minimum between
     ! dx, dy, dz
     
      !thickness = min( abs ( ( z( i   , j   , k+1 ) -  z( i   , j   , k-1 ) ) ) , &
      !                 abs ( ( y( i   , j+1 , k   ) -  y( i   , j-1 , k   ) ) ) , &
      !                 abs ( ( x( i+1 , j   , k   ) -  x( i-1 , j   , k   ) ) )   &
      !              ) 
    
      !damping   =  heaviside( phi(i,j,k) , thickness ) 
       
      !diss(1,i,j,k) =  pdiss_coef * damping * ( &
      !    ( fv(1,i+1,j,k) - two * fv(1,i,j,k) + fv(1,i-1,j,k) ) + &
      !    ( fv(2,i,j+1,k) - two * fv(2,i,j,k) + fv(2,i,j-1,k) ) + &
      !    ( fv(3,i,j,k+1) - two * fv(3,i,j,k) + fv(3,i,j,k-1) ) )
        
      diss(i,j,k) =  pdiss_coef *  ( &
                     ( fv(1,i+1,j,k) - two * fv(1,i,j,k) + fv(1,i-1,j,k) ) + &
                     ( fv(2,i,j+1,k) - two * fv(2,i,j,k) + fv(2,i,j-1,k) ) + &
                     ( fv(3,i,j,k+1) - two * fv(3,i,j,k) + fv(3,i,j,k-1) ) )

   end do
   end do
   end do

  ! Clearing out dissipation 
  !diss(1,:,:,:) = diss(1,:,:,:) * rsign(:,:,:) 

  ! Update ghost nodes
  call rhs_exchng3_3d( diss )


  !! extrapolate the dissipation to the first air layer
  !call rhs_ghost_fluid_nodes_extp_3d ( diss(1,:,:,:) )
  !! Update ghost nodes
  !call rhs_exchng3_3d( diss(1,:,:,:) )


!  ! THIS IS THE MODIFIED VERSION WITH THE DISSIPATION ESTIMATED USING
!  ! DIFFERENT SCHEMES. I COMMENTED IT TO TRY THE ORIGINAL ONE FOR THE
!  ! SLOSHING TANK BENCHMARK
!  
!  ista = i_mysta - 1
!  iend = i_myend + 1
!  
!  if (myback == mpi_proc_null)  ista = i_mysta
!  if (myfront == mpi_proc_null) iend = i_myend
!
!  do k = k_mysta-1, k_myend+1
!  do j = j_mysta-1, j_myend+1
!  do i = ista, iend
!
!   tmp1 = zero
!   tmp2 = zero
!
!   if( rsign(i,j,k) > one_half ) then
!
!      g1   = csi(1,i,j,k) * csi(1,i,j,k) + &
!             csi(2,i,j,k) * csi(2,i,j,k) + &
!             csi(3,i,j,k) * csi(3,i,j,k)
!
!      tmp1 = g1 * dtau(i,j,k) / aj(i,j,k) 
!
!      ! I compute tmp2 using water phase nodes only
!
!      !if i-1,i,i+1 = water --> rsignc = 1. Otherwise, rsignc = 0      
!      rsignc = rsign(i-1,j,k) * rsign(i,j,k) * rsign(i+1,j,k)
!
!      if ( rsignc > one_half ) then
!
!         tmp2 = dcsq * ( q(1,i+1,j,k) - two * q(1,i,j,k) + q(1,i-1,j,k) )
!   
!      else ! --> rsign(i-1,j,k) = 0 or rsign(i+1,j,k) = 0, or both
!
!         if ( i /= ista ) then
!
!            if ( rsign(i-1,j,k) > one_half .and. rsign(i-2,j,k) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!
!               tmp2 = dcsq * (q(1,i-2,j,k) - two * q(1,i-1,j,k) + q(1,i,j,k) )
!
!            end if
!
!         end if
!
!         if ( i /= iend ) then
!
!            if ( rsign(i+1,j,k) > one_half .and. rsign(i+2,j,k) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!               
!               tmp2 = dcsq * (q(1,i+2,j,k) - two * q(1,i+1,j,k) + q(1,i,j,k))
!
!            end if
!
!         end if
!
!      end if
!
!   end if
!
!   fv(1,i,j,k) = tmp1 * tmp2
!  
!  end do
!  end do
!  end do
!
!  jsta = j_mysta - 1
!  jend = j_myend + 1
!  if (myleft == mpi_proc_null)  jsta = j_mysta
!  if (myright == mpi_proc_null) jend = j_myend
!
!
!  ! eta-derivatives
!!   do k = k_mysta - 1, k_myend + 1
!!   do j = j_mysta - 1, j_myend + 1
!!   do i = i_mysta - 1, i_myend + 1
!
!  do k = k_mysta-1, k_myend+1
!  do j = jsta, jend
!  do i = i_mysta-1, i_myend+1
!
!   tmp1 = zero
!   tmp2 = zero
!
!   if( rsign(i,j,k) > one_half ) then
!
!      g2   = eta(1,i,j,k) * eta(1,i,j,k) + &
!             eta(2,i,j,k) * eta(2,i,j,k) + &
!             eta(3,i,j,k) * eta(3,i,j,k)
!
!      tmp1 = g2 * dtau(i,j,k) / aj(i,j,k)
!
!      ! I compute tmp2 using water phase nodes only
!
!      !if j-1,j,j+1 = water --> rsignc = 1. Otherwise, rsignc = 0      
!      rsignc = rsign(i,j-1,k) * rsign(i,j,k) * rsign(i,j+1,k)
!
!      if ( rsignc > one_half ) then
!
!         tmp2 = desq * (q(1,i,j+1,k) - two * q(1,i,j,k) + q(1,i,j-1,k))
!   
!      else ! --> rsign(i,j-1,k) = 0 or rsign(i,j+1,k) = 0, or both
!
!         if ( j /= jsta ) then
!
!            if ( rsign(i,j-1,k) > one_half .and. rsign(i,j-2,k) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!
!               tmp2 = desq * (q(1,i,j-2,k) - two * q(1,i,j-1,k) + q(1,i,j,k))
!
!            end if
!
!         end if
!
!         if ( j /= jend ) then
!
!            if ( rsign(i,j+1,k) > one_half .and. rsign(i,j+2,k) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!               
!               tmp2 = desq * (q(1,i,j+2,k) - two * q(1,i,j+1,k) + q(1,i,j,k))
!
!            end if
!
!         end if
!
!      end if
!
!   end if 
!
!   fv(2,i,j,k) = tmp1 * tmp2
!  
!  end do
!  end do
!  end do
!
!  ! zeta-derivatives
!!   do k = k_mysta - 1, k_myend + 1
!!   do j = j_mysta - 1, j_myend + 1
!!   do i = i_mysta - 1, i_myend + 1
!
!  ksta = k_mysta - 1
!  kend = k_myend + 1
!
!  if (mydown == mpi_proc_null)  ksta = k_mysta
!  if (myup   == mpi_proc_null)  kend = k_myend
!
!  do k = ksta      , kend
!  do j = j_mysta-1 , j_myend+1
!  do i = i_mysta-1 , i_myend+1
!
!   tmp1 = zero
!   tmp2 = zero
!
!   if( rsign(i,j,k) > one_half ) then
!
!      g3   = zet(1,i,j,k) * zet(1,i,j,k) + &
!             zet(2,i,j,k) * zet(2,i,j,k) + &
!             zet(3,i,j,k) * zet(3,i,j,k)
!
!      tmp1 = g3 * dtau(i,j,k) / aj(i,j,k) 
!
!      ! I compute tmp2 using water phase nodes only
!
!      !if k-1,k,k+1 = water --> rsignc = 1. Otherwise, rsignc = 0      
!      rsignc = rsign(i,j,k-1) * rsign(i,j,k) * rsign(i,j,k+1)
!
!      if ( rsignc > one_half ) then
!
!         tmp2 = dzsq * (q(1,i,j,k+1) - two * q(1,i,j,k) + q(1,i,j,k-1))
!   
!      else ! --> rsign(i,j,k-1) = 0 or rsign(i,j,k+1) = 0, or both
!
!         if ( k /= ksta ) then
!
!            if ( rsign(i,j,k-1) > one_half .and. rsign(i,j,k-2) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!
!               tmp2 = dzsq * (q(1,i,j,k-2) - two * q(1,i,j,k-1) + q(1,i,j,k))
!
!            end if
!
!         end if
!
!         if ( k /= kend ) then
!
!            if ( rsign(i,j,k+1) > one_half .and. rsign(i,j,k+2) > one_half ) then
!
!               ! First order accuracy one-sided second order derivative
!               ! https://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
!               
!               tmp2 = dzsq * (q(1,i,j,k+2) - two * q(1,i,j,k+1) + q(1,i,j,k))
!
!            end if
!
!         end if
!
!      end if
!
!   end if
!
!   fv(3,i,j,k) = tmp1 * tmp2
!
!  end do
!  end do
!  end do
!
!  ! boundary values
!  ! 
!
!  ! Put valid axu for nodes on domain (not process) boundary
!  ! 
!  if (myback  == mpi_proc_null)  fv(1,i_mysta-1,:,:) = fv(1,i_mysta,:,:)
!  if (myfront == mpi_proc_null)  fv(1,i_myend+1,:,:) = fv(1,i_myend,:,:)
!  if (myleft  == mpi_proc_null)  fv(2,:,j_mysta-1,:) = fv(2,:,j_mysta,:)
!  if (myright == mpi_proc_null)  fv(2,:,j_myend+1,:) = fv(2,:,j_myend,:)
!  if (mydown  == mpi_proc_null)  fv(3,:,:,k_mysta-1) = fv(3,:,:,k_mysta)
!  if (myup    == mpi_proc_null)  fv(3,:,:,k_myend+1) = fv(3,:,:,k_myend)
!
!  !call outputD3_real(fv(1,:,:,:),'disscsi_bmod')
!  !if(myid == 0) print *, 'output_disscsi'
!
!  ! boundary treatment of blanking nodes
!  if (nblk /= 0) then
!  do nb = 1, nblk
!
!     if (blktype(1,nb,myzone) == 0) then
!        i = li_blk_ia(n,nb)
!     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
!     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
!        fv(1,i,j,k) = fv(1,i-1,j,k)
!     end do
!     end do
!     end if
!
!     if (blktype(2,nb,myzone) == 0) then
!        i = li_blk_ib(n,nb)
!     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
!     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
!        fv(1,i,j,k) = fv(1,i+1,j,k)
!     end do
!     end do
!     end if
!
!     if (blktype(3,nb,myzone) == 0) then
!        j = li_blk_ja(n,nb)
!     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
!     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
!        fv(2,i,j,k) = fv(2,i,j-1,k)
!     end do
!     end do
!     end if
!
!     if (blktype(4,nb,myzone) == 0) then
!        j = li_blk_jb(n,nb)
!     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
!     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
!        fv(2,i,j,k) = fv(2,i,j+1,k)
!     end do
!     end do
!     end if
!
!     if (blktype(5,nb,myzone) == 0) then
!        k = li_blk_ka(n,nb)
!     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
!     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
!        fv(3,i,j,k) = fv(3,i,j,k-1)
!     end do
!     end do
!     end if
!
!     if (blktype(6,nb,myzone) == 0) then
!        k = li_blk_kb(n,nb)
!     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
!     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
!        fv(3,i,j,k) = fv(3,i,j,k+1)
!     end do
!     end do
!     end if
!  end do
!  end if
!
!  ! outer second derivative, interior nodes only
!  !
!
!  do k = k_mysta, k_myend
!  do j = j_mysta, j_myend
!  do i = i_mysta, i_myend
!   
!     diss(1,i,j,k) =     rsign(i,j,k) * dissLSM( one )                       * & 
!                     ( ( fv(1,i+1,j,k) - two * fv(1,i,j,k) + fv(1,i-1,j,k) ) + &
!                       ( fv(2,i,j+1,k) - two * fv(2,i,j,k) + fv(2,i,j-1,k) ) + &
!                       ( fv(3,i,j,k+1) - two * fv(3,i,j,k) + fv(3,i,j,k-1) )     )
!  
!
!  end do
!  end do
!  end do

end subroutine rhs_diss_p



