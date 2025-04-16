subroutine pressure_extrapolation2( )

  ! local variables

  integer :: i,j,k
  integer :: ii, jj, kk
  integer :: isearch, jsearch, ksearch
  integer :: iic, jic, kic
  integer :: is, js, ks
  integer :: ioc, joc, koc
  integer :: iaux, jaux, kaux

  real(kind=rdf) :: rAir2Water(3)

  real(kind=rdf) :: PhiAir, PhiWater, pWater
  real(kind=rdf) :: xIntpNode, yIntpNode, zIntpNode
  real(kind=rdf) :: dWater, dAir, dTotal
  real(kind=rdf) :: pAux, LocalWeight, TotalWeight
  real(kind=rdf) :: dx   
  real(kind=rdf) :: extp_threshold

  !integer :: i_mysta, &
  !           j_mysta, &
  !           k_mysta, &
  !           i_myend, &
  !           j_myend, &
  !           k_myend

  integer :: nWaterNodesExtp

  ! Interior nodes including domain boundaries
  !i_mysta = il + igp
  !j_mysta = jl + jgp
  !k_mysta = kl + kgp

  !i_myend = iu - igp
  !j_myend = ju - jgp
  !k_myend = ku - kgp

  extp_threshold  = 1.0E-08

  do k = ksta , kend !k_mysta - 1 , k_myend + 1
  do j = jsta , jend !j_mysta - 1 , j_myend + 1
  do i = ista , iend !i_mysta - 1 , i_myend + 1
            
    ! ----------------------------------------------------------------------------
    ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
    ! ----------------------------------------------------------------------------
    ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
    ! phase (rsign = 1), then this if is true. This condition identifies air nodes
    ! next to the free-surface
    ! ----------------------------------------------------------------------------

    ! If the node is closer to the free-surface than the machine-precision, then
    ! I'm on the free surface and the pressure is zero  
    if ( abs( phi(i,j,k) ) < eps_sims ) then

      q(1,i,j,k) = zero
      cycle

    end if

    ! local lenght scale
    !dx = aj(iaux,jaux,kaux)**( -one_third )
    dx = aj(i,j,k)**( -one_third )

    ! If I'm in the water phase or too far from the fs, I just cycle
    if ( rsign(i,j,k) > one_half .or. abs( phi(i,j,k) ) > radius_lsqm * dx ) cycle ! water-phase

    ! If I'm too deep into the water phase or too far from the fs, I just cycle
    !if ( phi(i,j,k) > extp_threshold .or. abs( phi(i,j,k) ) > radius_lsqm * dx ) cycle ! water-phase
               
    nWaterNodesExtp = 0

    ! Check if there are water or near water nodes in the neighbourhood 

    searchloop2 : &
    do kk = max( k-1 , ksta ) , min( k+1 , kend )
    do jj = max( j-1 , jsta ) , min( j+1 , jend )
    do ii = max( i-1 , ista ) , min( i+1 , iend )
      
      !if ( phi (ii,jj,kk) > -extp_threshold ) then

      if ( rsign (ii,jj,kk) > one_half ) then
            
        nWaterNodesExtp = 1
            
        exit searchloop2
            
      end if

    end do
    end do
    end do searchloop2

    if ( nWaterNodesExtp == 0 ) cycle

    ! variables initialisation
    
    pAux        = zero
    TotalWeight = zero

    PhiAir = abs( phi(i,j,k) )

    do kk = max( k-2 , ksta ) , min( k+2 , kend )
    do jj = max( j-2 , jsta ) , min( j+2 , jend )
    do ii = max( i-2 , ista ) , min( i+2 , iend )

      ! If the node is too close to the fs, I don't use
      ! if for interpolation
      if ( phi(ii,jj,kk) < extp_threshold .or. q(1,ii,jj,kk) < zero) cycle

      PhiWater = abs ( phi(ii,jj,kk) )
      pWater   = q(1,ii,jj,kk)

      rAir2Water = (/ x(ii,jj,kk) - x(i,j,k) , &
                      y(ii,jj,kk) - y(i,j,k) , &
                      z(ii,jj,kk) - z(i,j,k)    /)

      dTotal = norm2( rAir2Water )


      dWater = dTotal / ( one + PhiAir / PhiWater )
      dAir   = dTotal - dWater

      !print *, myid,ista,iend,jsta,jend,ksta,kend,ii,jj,kk, phi_gradient(1:3,ii,jj,kk)

      LocalWeight = (  dot_product ( rAir2Water  / norm2( rAir2Water ) , &
                      phi_gradient(1:3,ii,jj,kk) / norm2( phi_gradient(1:3,ii,jj,kk) ) ) )!**two

      pAux        = pAux - LocalWeight * ( dAir / dWater ) * pWater 
      TotalWeight = TotalWeight + LocalWeight

    end do
    end do
    end do

    if ( TotalWeight > eps_sims ) q(1,i,j,k) = pAux / TotalWeight

  end do
  end do
  end do


end subroutine pressure_extrapolation2