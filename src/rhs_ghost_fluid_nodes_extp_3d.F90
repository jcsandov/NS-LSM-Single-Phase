subroutine rhs_ghost_fluid_nodes_extp_3d( var , InteriorNodesOnly)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Description:
  !
  ! extrapolate a variable var to the first ghost-nodes layer using the
  ! geometric procedure used for pressure. 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use AdvectionMethods

  implicit none

  ! variable to extrapolate
  ! 
  real (kind = rdf), dimension (il:iu,jl:ju,kl:ku), intent(inout) :: var
  logical, intent(in), optional :: InteriorNodesOnly

  ! local variables
  real (kind = rdf) :: tmp
  real (kind = rdf) :: var_aux
  real (kind = rdf) :: search_radius
  real (kind = rdf) :: phase
  real (kind = rdf) :: weight , TotalWeight

  real(kind=rdf) :: rAir2Water(3) , phi_grad_avg(3)

  real(kind=rdf) :: alpha1 , alpha2
  real(kind=rdf) :: num , denom 
  real(kind=rdf), parameter :: eps_denom = 1E-6 

  ! neighbourhood phase counter
  integer :: nWaterNodesExtp

  ! stencil indexes  
  integer :: i1,j1,k1,i2,j2,k2 

  ! Interior nodes
  integer :: i_mysta, j_mysta, k_mysta
  integer :: i_myend, j_myend, k_myend

  ! physical boundaries indexes
  integer :: ista , iend , jsta , jend , ksta , kend ! 

  integer :: i,j,k,ii,jj,kk

  ! Loop boundaries
  integer :: ilbound, iubound, jlbound, jubound, klbound, kubound

  ! Interior nodes
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


  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 

  if ( myback  == mpi_proc_null )  ista = il + igp 
  if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
  if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

  if ( myfront == mpi_proc_null )  iend = iu - igp
  if ( myright == mpi_proc_null )  jend = ju - jgp
  if ( myup    == mpi_proc_null )  kend = ku - kgp

  ! Loop boundaries
  ilbound = i_mysta - 1 ; iubound = i_myend + 1 
  jlbound = j_mysta - 1 ; jubound = j_myend + 1 
  klbound = k_mysta - 1 ; kubound = k_myend + 1 

  ! Search restriction to interior nodes only
  if ( present ( InteriorNodesOnly ) ) then

    if ( InteriorNodesOnly ) then

      ! Loop boundaries
      ilbound = i_mysta ; iubound = i_myend 
      jlbound = j_mysta ; jubound = j_myend 
      klbound = k_mysta ; kubound = k_myend 
  
      if ( myback  == mpi_proc_null )  ista = i_mysta 
      if ( myleft  == mpi_proc_null )  jsta = j_mysta 
      if ( mydown  == mpi_proc_null )  ksta = k_mysta 
    
      if ( myfront == mpi_proc_null )  iend = i_myend
      if ( myright == mpi_proc_null )  jend = j_myend
      if ( myup    == mpi_proc_null )  kend = k_myend

    end if

  end if

  do k = klbound , kubound
  do j = jlbound , jubound
  do i = ilbound , iubound

    ! search radius for extrapolation
    !search_radius = three * aj(i,j,k)**( -one_third )
    search_radius = ten * aj(i,j,k)**( -one_third )

    ! If I'm in the water phase or too far from the fs, I just cycle
    if ( rsign(i,j,k) > one_half .or. abs(phi(i,j,k)) > search_radius ) cycle 

    ! Check if there are water or near water nodes in the neighbourhood 
    nWaterNodesExtp = 0

    searchloop2 : &
    do kk = max( k-1 , ksta ) , min( k+1 , kend )
    do jj = max( j-1 , jsta ) , min( j+1 , jend )
    do ii = max( i-1 , ista ) , min( i+1 , iend )
      
      if ( rsign (ii,jj,kk) > one_half ) then
            
        nWaterNodesExtp = 1
            
        exit searchloop2
            
      end if

    end do
    end do
    end do searchloop2

    if ( nWaterNodesExtp == 0 ) cycle

    ! Total weight and auxiliary vars initialisation
    TotalWeight = zero
    var_aux     = zero

    do kk = -1,1
    do jj = -1,1
    do ii = -1,1

      ! Local weight initialisation
      weight = zero
      tmp    = zero

      ! Extended stencil
      i1 = i +   ii
      i2 = i + 2*ii

      j1 = j +   jj
      j2 = j + 2*jj

      k1 = k +   kk
      k2 = k + 2*kk

      ! skipping central node and avoiding outbounding
      if ( all( [ii,jj,kk] == 0 )   .or. &
           i2 < ista .or. i2 > iend .or. &
           j2 < jsta .or. j2 > jend .or. &
           k2 < ksta .or. k2 > kend ) cycle

      ! checking the phase of the stencil
      phase = rsign(i1,j1,k1) * rsign(i2,j2,k2)

      ! skipping non-water or mixed stencils
      if ( phase < one_half ) cycle

      ! distance vector between the average position of the stencil and
      ! the air node to be extrapolated
      rAir2Water = (/ ( x(i1,j1,k1) + x(i2,j2,k2) ) / two - x(i,j,k) , &
                      ( y(i1,j1,k1) + y(i2,j2,k2) ) / two - y(i,j,k) , &
                      ( z(i1,j1,k1) + z(i2,j2,k2) ) / two - z(i,j,k)    /)
  
      phi_grad_avg = ( phi_gradient(:,i1,j1,k1) + phi_gradient(:,i2,j2,k2) ) / two

      ! weight = Δr * ∇ϕ (normalised)
      num   = abs( dot_product ( rAir2Water   / norm2( rAir2Water   ) , &
                                 phi_grad_avg / norm2( phi_grad_avg )     ) )**(three/two)

      denom = norm2( rAir2Water )**(three/two) + eps_denom

      weight = num / denom

      ! tmp = 2v1-v2
      !tmp   = two * var(i1,j1,k1) - var(i2,j2,k2)
  
      call GetLinearExtrapolationCoefficients( x(i1,j1,k1)      ,  y(i1,j1,k1) ,  z(i1,j1,k1) , &
                                               x(i2,j2,k2)      ,  y(i2,j2,k2) ,  z(i2,j2,k2) , &
                                               x(i,j,k)         ,  y(i,j,k)    ,  z(i,j,k)    , & 
                                               alpha1 , alpha2                                  &
                                             )

      tmp   = alpha1 * var(i1,j1,k1) + alpha2 * var(i2,j2,k2)
      
      ! if phase = 0, var_aux doesn't change
      ! if phase = 1, var_aux = var_aux + weight * tmp
      var_aux     =  var_aux     + weight * tmp 
      TotalWeight =  TotalWeight + weight

    end do
    end do
    end do

    ! DEBUG
    if ( TotalWeight > eps_sims ) var(i,j,k) = var_aux / TotalWeight
    
    !if ( TotalWeight > eps_sims ) then
    !  var(i,j,k) = two * var(i,j,k-1) - var(i,j,k-2)
    !end if

  end do
  end do
  end do
  
end subroutine rhs_ghost_fluid_nodes_extp_3d


