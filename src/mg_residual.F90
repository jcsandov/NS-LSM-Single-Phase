
subroutine mg_residual (n, residual)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! This routine computes the L2-norm of the residual
  ! on the finest grid at the end of each pseudo time step

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  !real (kind = rdf) :: erp, eru, erv, erw, erx
  real (kind = rdf) :: tnu

  real (kind = rdf), dimension (me) :: resd
  real (kind = rdf), dimension (me) :: resd_tmp
  real (kind = rdf), dimension (me) :: residual
  integer :: i, j, k, l
  integer :: m
  integer :: n

  integer :: k_mysta
  integer :: j_mysta
  integer :: i_mysta

  integer :: k_myend
  integer :: j_myend
  integer :: i_myend

  integer, dimension (0:nproc-1) :: recvcounts

  tnu = real((img(n,myzone)-2)*(jmg(n,myzone)-2)*(kmg(n,myzone)-2), kind = rdf)

  k_mysta = li_ka(n)
  j_mysta = li_ja(n)
  i_mysta = li_ia(n)

  k_myend = li_kb(n)
  j_myend = li_jb(n)
  I_myend = li_ib(n)

  ! use interior nodes only
  ! 
  if (myback  == mpi_proc_null) i_mysta = i_mysta + 1
  if (myleft  == mpi_proc_null) j_mysta = j_mysta + 1
  if (mydown  == mpi_proc_null) k_mysta = k_mysta + 1

  if (myfront == mpi_proc_null) i_myend = i_myend - 1
  if (myright == mpi_proc_null) j_myend = j_myend - 1
  if (myup    == mpi_proc_null) k_myend = k_myend - 1

  ! calculate L2 norm
  ! 
  resd_tmp(:) = zero
  resd(:)     = zero
  residual(:) = zero

  do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           l = le_idx(i,j,k,n)
           !do m = 1, me
           !   resd_tmp(m) = resd_tmp(m) + abs(q(m,l) - qold_mg(m,l))
           !end do

           do m = 1, 4
              resd_tmp(m) = resd_tmp(m) + abs( q(m,l)  - qold_mg(m,l))
           end do

           resd_tmp(5)    = resd_tmp(5) + abs( xnut(l) - qold_mg(5,l) )

        end do
     end do
  end do

  recvcounts(:) = me

  call mpi_reduce (resd_tmp, resd, me, MPI_REAL_TYPE, mpi_sum,&
       & root, mpi_comm_world, ierr)

  !   call mpi_reduce_scatter (resd_tmp, resd, recvcounts, MPI_REAL_TYPE, mpi_sum,&
  !        & mpi_comm_world, ierr)

  ! calculate average residual
  !
  if ( myid == root ) then

!      resd(:) = sqrt(resd(:) / tnu)
      resd(:) = resd(:) / tnu

      do m = 1,me

         if ( resd(m) > eps_sims ) then
            residual(m) = log10(resd(m)) ! zero !log10(resd(:))
         else
            residual(m) = -9999.9
         end if

      end do
!      erp = log10(resd(1))
!      eru = log10(resd(2))
!      erv = log10(resd(3))
!      erw = log10(resd(4))

!      if (me == 5) erx = log10(resd(5))

  else

!      erp = zero
!      eru = zero
!      erv = zero
!      erw = zero

     residual(:) = zero

  end if

end subroutine mg_residual

