subroutine p_correction_post_advection( )

  use InterpolationMethods

  implicit none

  ! local variables

  integer :: i,j,k
  integer :: ii, jj, kk

  real(kind=rdf) :: rIntp(3)

  real(kind=rdf) :: PhiIntp, pIntp
  real(kind=rdf) :: pAux, LocalWeight, TotalWeight
  real(kind=rdf) , parameter :: extp_threshold = 1.0E-05
  real(kind=rdf) , parameter :: pmax = two

  integer :: i_mysta , j_mysta , k_mysta , &
             i_myend , j_myend , k_myend 
  
  integer :: ista, iend, jsta, jend, ksta, kend


  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp
  
  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp
           
  ! processes on the domain boundaries
  
  if ( myback  == mpi_proc_null )  i_mysta = il + igp + 1
  if ( myleft  == mpi_proc_null )  j_mysta = jl + jgp + 1
  if ( mydown  == mpi_proc_null )  k_mysta = kl + kgp + 1
  
  if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1
  if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1
  if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1
  
  ! Physical boundaries
  
  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 
  
  if ( myback  == mpi_proc_null )  ista = il + igp 
  if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
  if ( mydown  == mpi_proc_null )  ksta = kl + kgp 
  
  if ( myfront == mpi_proc_null )  iend = iu - igp
  if ( myright == mpi_proc_null )  jend = ju - jgp
  if ( myup    == mpi_proc_null )  kend = ku - kgp

  !---------------------------------------------------------------------------------

  do k = k_mysta - 1 , k_myend + 1 !ksta , kend !k_mysta - 1 , k_myend + 1
  do j = j_mysta - 1 , j_myend + 1 !jsta , jend !j_mysta - 1 , j_myend + 1
  do i = i_mysta - 1 , i_myend + 1 !ista , iend !i_mysta - 1 , i_myend + 1
          
    if ( rsign(i,j,k) > one_half .and. phi(i,j,k) * phi_n(i,j,k) < zero ) then !eps_sims ) then

      ! If the node is closer to the free-surface than the machine-precision, then
      ! I'm on the free surface and the pressure is zero  
      if ( abs( phi(i,j,k) ) < eps_sims ) then

        q(1,i,j,k)   = zero
        rsign(i,j,k) = -one

        cycle

      end if

      ! variables initialisation
      
      pAux        = zero
      TotalWeight = zero
    
      do kk = max( k-2 , ksta ) , min( k+2 , kend )
      do jj = max( j-2 , jsta ) , min( j+2 , jend )
      do ii = max( i-2 , ista ) , min( i+2 , iend )
  
        ! If the node is too close to the fs, I don't use
        ! if for interpolation. Also, if the gradient is too large, maybe
        ! I'm taking nodes away from the geometric reinitialisation range,
        ! so I discard those nodes too.

        !if ( ( ii==i .and. jj==j .and. kk == k )                         .or. &
        !       phi(ii,jj,kk)                       < extp_threshold      .or. &
        !       q(1,ii,jj,kk)                       < extp_threshold      .or. &
        !       phi(ii,jj,kk) * q(1,ii,jj,kk)       < eps_sims            .or. &
        !       norm2( phi_gradient(1:3,ii,jj,kk) ) > two                 .or. &
        !       norm2( phi_gradient(1:3,ii,jj,kk) ) < pt1                        ) cycle
  
        if ( ( ii==i .and. jj==j .and. kk == k )                         .or. &
               phi(ii,jj,kk)                       < extp_threshold      .or. &
               phi(ii,jj,kk) * phi_n(ii,jj,kk)     < eps_sims            .or. &
               norm2( phi_gradient(1:3,ii,jj,kk) ) > two                 .or. &
               norm2( phi_gradient(1:3,ii,jj,kk) ) < pt1                        ) cycle

        PhiIntp = phi(ii,jj,kk)
        pIntp   = q(1,ii,jj,kk)
  
        rIntp = (/ x(ii,jj,kk) - x(i,j,k) , &
                   y(ii,jj,kk) - y(i,j,k) , &
                   z(ii,jj,kk) - z(i,j,k)    /)

        !if ( norm2( phi_gradient(1:3,ii,jj,kk) )> two .or. norm2( phi_gradient(1:3,ii,jj,kk)) < 0.1 ) then
        !  print *,myid,ii,jj,kk,phi(ii,jj,kk),phi_gradient(1:3,ii,jj,kk)
        !end if

        LocalWeight = abs( dot_product ( rIntp  / norm2( rIntp ) , &
                           phi_gradient(1:3,ii,jj,kk) / norm2( phi_gradient(1:3,ii,jj,kk) ) ) )!**two
  
        ! If the locally extrapolated pressure is too big, then it's not considered
        if ( ( phi(i,j,k) / PhiIntp ) * pIntp > pmax ) cycle

        pAux        = pAux + LocalWeight * ( phi(i,j,k) / PhiIntp ) * pIntp 
        TotalWeight = TotalWeight + LocalWeight
  
      end do
      end do
      end do
  
      if ( abs( TotalWeight ) > eps_sims ) then
        
        q(1,i,j,k) = pAux / TotalWeight
        
      else

        q(1,i,j,k) = two * q(1,i,j,k-1) - q(1,i,j,k-2)
        
      end if
      

    end if ! ( rsign(i,j,k) > one_half .and. phi(i,j,k) * q(1,i,j,k) < zero )

  end do
  end do
  end do

end subroutine p_correction_post_advection