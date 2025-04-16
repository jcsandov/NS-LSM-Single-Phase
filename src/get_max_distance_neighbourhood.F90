subroutine get_max_distance_neighbourhood(xs,ys,zs,im,ip,jm,jp,km,kp,dmax)

  !------------------------------------------------------------------------
  ! Description:
  ! -------------
  ! get_max_distance_neighbourhood returns the largest distance
  ! among the water nodes in the region defined by the nodes
  ! im,ip,jm,jp,km,kp
  !------------------------------------------------------------------------

  implicit none
  ! input/output variables
  real(kind = rdf) , intent(in) :: xs,ys,zs
  integer , intent(in) :: im,ip,jm,jp,km,kp
  real(kind = rdf) :: dmax

  !local variables
  real(kind=rdf),dimension(1:3) :: delta_r
  integer :: ii,jj,kk
  real(kind=rdf) :: aux_sign

  dmax = zero

  do kk = km , kp
  do jj = jm , jp
  do ii = im , ip

    aux_sign = one

    ! If I check an air node, the Δr vector is set to zero
    if( rsign(ii,jj,kk) < one_half ) aux_sign = zero

    delta_r(1) = aux_sign * ( xs - x(ii,jj,kk) )
    delta_r(2) = aux_sign * ( ys - y(ii,jj,kk) )
    delta_r(3) = aux_sign * ( zs - z(ii,jj,kk) )

    dmax = max( dmax , norm2(delta_r) )

  end do
  end do
  end do

end subroutine get_max_distance_neighbourhood
