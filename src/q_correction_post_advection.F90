subroutine q_correction_post_advection( )

  use InterpolationMethods
  use AdvectionMethods

  implicit none

  ! local variables

  integer :: i,j,k
  integer :: ii, jj, kk
  integer :: i1, j1, k1
  integer :: i2, j2, k2

  real(kind=rdf) :: rIntp(3)

  real (kind = rdf), dimension(1:5) :: tmp
  real (kind = rdf), dimension(1:5) :: var_aux

  real(kind=rdf) :: PhiIntp, pIntp
  real(kind=rdf) :: weight, TotalWeight
  real(kind=rdf) :: alpha1 , alpha2
  real(kind=rdf) :: phase
  real(kind=rdf) :: num , denom 

  real(kind=rdf) :: rAir2Water(3) , phi_grad_avg(3)
  real(kind=rdf) :: uvec1(3), uvec2(3) , umaxvec(3)
  real(kind=rdf) :: exsign 

  logical :: BlankingFlag

  real(kind=rdf) , parameter :: extp_threshold = 1.0E-05
  real(kind=rdf) , parameter :: pmax = two

  integer :: i_mysta , j_mysta , k_mysta , &
             i_myend , j_myend , k_myend 
  
  ! indexes of the maximum velocity
  integer :: iumax, jumax, kumax
  real(kind=rdf) :: max_u_norm 

  integer :: ista, iend, jsta, jend, ksta, kend

  real(kind=rdf), parameter :: eps_denom = 1E-08

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

      ! Total weight and auxiliary vars initialisation
      TotalWeight = zero
      var_aux     = zero
      max_u_norm  = zero
      umaxvec     = zero

      do kk = -1,1
      do jj = -1,1
      do ii = -1,1

        !---------------------------------------------------------------------------------
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

        !---------------------------------------------------------------------------------
        ! skipping central node and avoiding outbounding
        if ( all( [ii,jj,kk] == 0 )   .or. &
             i2 < ista .or. i2 > iend .or. &
             j2 < jsta .or. j2 > jend .or. &
             k2 < ksta .or. k2 > kend ) cycle

        !---------------------------------------------------------------------------------
        ! Skipping Blanking area
        BlankingFlag = .false.

        if ( nblke /= 0 ) then
          
          do nb = 1,nblke

            if ( i2 > le_blk_ia(1,nb) .and. i2 < le_blk_ib(1,nb) .and. &  
                 j2 > le_blk_ja(1,nb) .and. j2 < le_blk_jb(1,nb) .and. &
                 k2 > le_blk_ka(1,nb) .and. k2 < le_blk_kb(1,nb) ) then

              BlankingFlag = .true.
        
            end if

          end do  
        
        end if

        if ( BlankingFlag ) cycle

        !---------------------------------------------------------------------------------
        ! checking the phase of the stencil
        phase = rsign(i1,j1,k1) * rsign(i2,j2,k2)

        ! skipping non-water stencils
        if ( phase < one_half ) cycle

        ! checking that the node wasnt a node that needs correction too
        phase = ( sign( one , phi(i1,j1,k1) * phi_n(i1,j1,k1) ) + one ) * &
                ( sign( one , phi(i2,j2,k2) * phi_n(i2,j2,k2) ) + one )      

        ! skipping non-water stencils
        if ( phase < one_half ) cycle

        !---------------------------------------------------------------------------------
        ! distance vector between the average position of the stencil and
        ! the air node to be extrapolated
        rAir2Water = (/ ( x(i1,j1,k1) + x(i2,j2,k2) ) / two - x(i,j,k) , &
                        ( y(i1,j1,k1) + y(i2,j2,k2) ) / two - y(i,j,k) , &
                        ( z(i1,j1,k1) + z(i2,j2,k2) ) / two - z(i,j,k)    /)
  
        phi_grad_avg = ( phi_gradient(:,i1,j1,k1) + &
                         phi_gradient(:,i2,j2,k2) ) / two

        ! I discard phi gradients where it's not well defined
        if ( norm2( phi_grad_avg ) < eps_denom .or. norm2( phi_grad_avg ) > five ) cycle

        !---------------------------------------------------------------------------------
        ! weight = Δr * ∇ϕ (normalised)
        num   = abs( dot_product ( rAir2Water   / norm2( rAir2Water   ) , &
                                   phi_grad_avg / norm2( phi_grad_avg )     ) )**(three/two)

        denom = one !norm2( rAir2Water )**(three/two) + eps_denom

        weight = num / denom
               
        ! tmp = 2v1-v2
        ! tmp   = two * var(:,i1,j1,k1) - var(:,i2,j2,k2)

        call GetLinearExtrapolationCoefficients( x(i1,j1,k1)      ,  y(i1,j1,k1) ,  z(i1,j1,k1) , &
                                                 x(i2,j2,k2)      ,  y(i2,j2,k2) ,  z(i2,j2,k2) , &
                                                 x(i,j,k)         ,  y(i,j,k)    ,  z(i,j,k)    , & 
                                                 alpha1 , alpha2                                  &
                                               )

        !---------------------------------------------------------------------------------
        ! I keep the track of the maximum input velocity to generate
        ! the extrapolation to limit it if necessary
        uvec1 = (/q(2,i1,j1,k1) , q(3,i1,j1,k1) , q(4,i1,j1,k1)/)
        uvec2 = (/q(2,i2,j2,k2) , q(3,i2,j2,k2) , q(4,i2,j2,k2)/)

        ! i1,j1,k1
        ! |uvec1| > |umax| --> exsign = 1 . exsign = 0 otherwise 
        exsign     = ( sign( one , norm2( uvec1 ) - max_u_norm ) + one ) / two
        max_u_norm = max_u_norm + exsign * ( norm2( uvec1 ) - max_u_norm )
        umaxvec    = umaxvec    + exsign * ( uvec1 - umaxvec ) 

        ! i2,j2,k2
        ! |uvec2| > |umax| --> exsign = 1 . exsign = 0 otherwise 
        exsign     = ( sign( one , norm2( uvec2 ) - max_u_norm ) + one ) / two
        max_u_norm = max_u_norm + exsign * ( norm2( uvec2 ) - max_u_norm )
        umaxvec    = umaxvec    + exsign * ( uvec2 - umaxvec ) 

        ! pressure and velocities
        !tmp(2:4)   = alpha1 * q(2:4,i1,j1,k1) + alpha2 * q(2:4,i2,j2,k2)
        tmp(1:4)   = alpha1 * q(1:4,i1,j1,k1) + alpha2 * q(1:4,i2,j2,k2)
        ! xnut
        tmp(5)     = alpha1 * xnut(i1,j1,k1)  + alpha2 * xnut(i2,j2,k2)
      
        var_aux     =  var_aux     + weight * tmp 
        TotalWeight =  TotalWeight + weight

      end do
      end do
      end do

      ! DEBUG
      if ( TotalWeight > eps_sims ) then 

        !q(2:4,i,j,k) = var_aux(2:4) / TotalWeight
        q(1:4,i,j,k) = var_aux(1:4) / TotalWeight
        xnut(i,j,k)  = var_aux(5)   / TotalWeight
  
        ! limit the velocity if necessary
        exsign = ( sign( one , norm2( q(2:4,i,j,k) ) - max_u_norm ) + one ) / two
        !
        q(2:4,i,j,k) = q(2:4,i,j,k) + exsign * ( umaxvec - q(2:4,i,j,k) )
        
      end if
    
    end if ! ( rsign(i,j,k) > one_half .and. phi(i,j,k) * q(1,i,j,k) < zero )

  end do
  end do
  end do

end subroutine q_correction_post_advection