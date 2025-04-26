subroutine get_weight_LSM_GFM( xs , ys , zs , nvec , dmax , ii , jj , kk , weight )

  implicit none

  real(kind = rdf) , intent(in) :: xs,ys,zs
  real(kind=rdf),dimension(1:3) , intent(in) :: nvec
  integer , intent(in) :: ii,jj,kk
  real(kind = rdf) , intent(in) :: dmax
  real(kind = rdf) :: weight

  !local variables
  real(kind=rdf),dimension(1:3) :: delta_r
  real(kind=rdf) :: num , denom
  real(kind=rdf), parameter :: eps_denom = 1.0E-6

  if ( rsign( phi(ii,jj,kk) ) < one_half ) return

  !---------------------------------------------------
  !   
  !          | r̂ • n̂ |
  !   w  = -------------- 
  !            d̄ + ε   
  !
  ! where:  n̂ : normal vector at the free surface
  !         r̂ = r / ‖ r ‖
  !         d̄ = d / dmax
  !         ε : parameter to avoid division by zero  
  !--------------------------------------------------

  delta_r(1) =  xs - x(ii,jj,kk) 
  delta_r(2) =  ys - y(ii,jj,kk) 
  delta_r(3) =  zs - z(ii,jj,kk) 

  num   = abs( dot_product( delta_r / norm2(delta_r)  , nvec ) )**two!**(three/two)
  denom = one !( norm2(delta_r) / dmax ) + eps_denom

  weight = num / denom

end subroutine get_weight_LSM_GFM