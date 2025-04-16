
subroutine narrowband()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none
  real (kind = rdf) :: tol
  integer :: iNB,jNB,kNB

  tol = 1.0E-5   
  nband = 0           
  do kNB = kl, ku
  do jNB = jl, ju
  do iNB = il, iu
        if(abs(phi_n(iNB,jNB,kNB))<narrowcoef*epslsm + tol) nband(iNB,jNB,kNB) = 1
  end do
  end do
  end do
  
  if(isnarrow == 0) nband = 1  
end subroutine narrowband



