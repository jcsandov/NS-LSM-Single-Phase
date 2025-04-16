
subroutine rhs_contra_j
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate the contravariant velocity (U) divided by
  ! the Jacobian (J) at all nodes. U/J appears (slightly) more
  ! often than U so it makes sense to use U/J instead of U

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  real (kind = rdf) :: tmp

  integer :: ista , iend , jsta , jend , ksta , kend
  integer :: i,j,k

  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 

  if ( myback  == mpi_proc_null )  ista = il + igp 
  if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
  if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

  if ( myfront == mpi_proc_null )  iend = iu - igp
  if ( myright == mpi_proc_null )  jend = ju - jgp
  if ( myup    == mpi_proc_null )  kend = ku - kgp

  do k = ksta , kend
  do j = jsta , jend
  do i = ista , iend

      tmp = one / aj(i,j,k)

      ucn_j(1,i,j,k) = tmp * ( csi(1,i,j,k) * q(2,i,j,k) + &
                               csi(2,i,j,k) * q(3,i,j,k) + &
                               csi(3,i,j,k) * q(4,i,j,k)      )

      ucn_j(2,i,j,k) = tmp * ( eta(1,i,j,k) * q(2,i,j,k) + &
                               eta(2,i,j,k) * q(3,i,j,k) + &
                               eta(3,i,j,k) * q(4,i,j,k)      )
      
      ucn_j(3,i,j,k) = tmp * ( zet(1,i,j,k) * q(2,i,j,k) + &
                               zet(2,i,j,k) * q(3,i,j,k) + &
                               zet(3,i,j,k) * q(4,i,j,k)     )
  
  end do
  end do
  end do
  
end subroutine rhs_contra_j



