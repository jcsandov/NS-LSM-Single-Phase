
subroutine conservation(il, iu, jl, ju, kl, ku, & !  1
			       igp, jgp, kgp,          & !  2
                         xc, yc, zc,        &
                         xe, ye, ze,        &
                         xz, yz, zz,        &
                         csi,               & !  8
                         eta,               & !  9
                         zet,               & ! 10
			             aj,                &
			             phi,		        &
                         rsign,             &
                         q ,                &
                         x,y,z)


  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! gathers grid and re-calculates interpolation coeff and distance
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global_param
  use global_app
  use global_mpi
  use global_lsm, only :  rhoLSM, k_surface
  !use global_debug

  implicit none
  !-------debug-----------
  character (len=256) :: filename
  integer :: myunit

  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: var_aux
  !-------debug-----------


  integer :: il,jl,kl,iu,ju,ku

  integer :: igp,jgp,kgp

!!  integer, dimension(nzones) :: im,jm,km

  integer :: i,  j,  k, l, nv

  integer :: np

  integer :: i_mysta, &
             j_mysta, &
             k_mysta, &
             i_myend, &
             j_myend, &
             k_myend

  !Level set method
  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi, rsign !phi_n

  ! local arrays
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: q 
  !real (kind = rdf), dimension(    il:iu,jl:ju,kl:ku) :: rsign 
  real (kind = rdf), dimension(    il:iu,jl:ju,kl:ku) :: aj,ajn,ajnm1
  real (kind = rdf), dimension(    il:iu,jl:ju,kl:ku) :: xc, yc, zc, xe, ye, ze, xz, yz, zz
  real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) :: csi, eta, zet

  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), optional :: x,y,z


  ! allocatable buffers
  real (kind = rdf), dimension (:), allocatable :: sbuf1, sbuf2, rbuf


  ! buffer construction
  integer :: nz,num_vars,sbuf_size1,sbuf_size2
  integer :: rbuf_size,imax,jmax,kmax
  integer :: istatus(mpi_status_size)

  real (kind = rdf), dimension(:,:,:,:), allocatable :: ucn_j
  real (kind = rdf):: volchange, vol, voln, volnm1, dumflux, dumvol
  real (kind = rdf):: fluxin, fluxout, delta_e

  ! For LSM scalation

  real (kind = rdf), dimension (:), allocatable :: sbuf3, sbuf4
  integer :: sbuf_size3, sbuf_size4
  real (kind = rdf) :: fluxin_w, fluxout_w, fluxin_a, fluxout_a 
  real (kind = rdf) :: delta_e_w, delta_e_a

  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp


  ! allocate variables 
  allocate (ucn_j(1:3,il:iu,jl:ju,kl:ku))

  ucn_j = zero

  ! calculate contravariant velocities (that is, U/J)
  ! I commented this call to rhs_contra_j for debbuging LSM
  ! uncomment afterwards.
  ! call rhs_contra_j ()

!Debug------------------------------------------
  !call outputD3_real(ucn_j(1,:,:,:),"ucn_j")
!-----------------------------------------------


num_vars = 1  ! ucn_j(1,i,:,:) 

! local size of buffer:

!Water

sbuf_size1 = (j_myend-j_mysta+1) *&
		(k_surface-k_mysta+1)

sbuf_size2 = (j_myend-j_mysta+1) *&
		(k_surface-k_mysta+1)

!Air
sbuf_size3 = (j_myend-j_mysta+1) *&
		(k_myend-k_surface)

sbuf_size4 = (j_myend-j_mysta+1) *&
		(k_myend-k_surface)



allocate ( sbuf1(sbuf_size1), sbuf2(sbuf_size2)  )

allocate ( sbuf3(sbuf_size3), sbuf4(sbuf_size4) ) 

sbuf1=zero
sbuf2=zero
sbuf3=zero
sbuf4=zero

! pack send bufs with current contravariant velocity in each processor
! we'll broadcast only i_mysta and i_myend planes

!Water
l = 0
do k = k_mysta, k_surface; do j = j_mysta, j_myend 
    l=l+1; sbuf1(l) = ucn_j(1,i_mysta,j,k)
end do; end do 

l = 0
do k = k_mysta, k_surface; do j = j_mysta, j_myend 
    l=l+1; sbuf2(l) = ucn_j(1,i_myend,j,k)
end do; end do 


!Air
l = 0
do k = k_surface + 1, k_myend; do j = j_mysta, j_myend 
    l=l+1; sbuf3(l) = ucn_j(1,i_mysta,j,k)
end do; end do 

l = 0
do k = k_surface + 1, k_myend; do j = j_mysta, j_myend 
    l=l+1; sbuf4(l) = ucn_j(1,i_myend,j,k)
end do; end do 


! every processor has 2 "i" planes
! only boundary processors send the corresponding domains

if (myback.ne.mpi_proc_null)   sbuf1 = zero
if (myfront.ne.mpi_proc_null)  sbuf2 = zero

if (myback.ne.mpi_proc_null)   sbuf3 = zero
if (myfront.ne.mpi_proc_null)  sbuf4 = zero



dumflux=zero
dumflux=sum(sbuf1)
call mpi_allreduce(dumflux, fluxin_w, 1, MPI_REAL_TYPE, &
                   mpi_sum, mpi_comm_world, ierr)

dumflux=zero
dumflux=sum(sbuf2)
call mpi_allreduce(dumflux, fluxout_w, 1, MPI_REAL_TYPE, &
                   mpi_sum, mpi_comm_world, ierr)

dumflux=zero
dumflux=sum(sbuf3)
call mpi_allreduce(dumflux, fluxin_a, 1, MPI_REAL_TYPE, &
                   mpi_sum, mpi_comm_world, ierr)

dumflux=zero
dumflux=sum(sbuf4)
call mpi_allreduce(dumflux, fluxout_a, 1, MPI_REAL_TYPE, &
                   mpi_sum, mpi_comm_world, ierr)


delta_e_w= (fluxin_w-fluxout_w)	
delta_e_a= (fluxin_a-fluxout_a)    

! now delta_e contains the total error, which is proportionally distributed

  !Water
        if (myfront == mpi_proc_null) then
                do k = k_mysta, k_surface; do j = j_mysta, j_myend
                        ucn_j(1,i_myend,j,k)=ucn_j(1,i_myend,j,k)+ delta_e_w*ucn_j(1,i_myend,j,k)/fluxout_w

                        ! transform back:
                        i = i_myend

                        q(2,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*xc(i,j,k) +&
                                ucn_j(2,i,j,k)*xe(i,j,k)  + ucn_j(3,i,j,k)*xz(i,j,k) )
                        q(3,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*yc(i,j,k) +&
                                ucn_j(2,i,j,k)*ye(i,j,k)  + ucn_j(3,i,j,k)*yz(i,j,k) )
                        q(4,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*zc(i,j,k) +&
                                ucn_j(2,i,j,k)*ze(i,j,k)  + ucn_j(3,i,j,k)*zz(i,j,k) )

                end do; end do
        end if


  !Air
        if (myfront == mpi_proc_null) then
                do k = k_surface + 1, k_myend; do j = j_mysta, j_myend
                        ucn_j(1,i_myend,j,k)=ucn_j(1,i_myend,j,k)+ delta_e_a*ucn_j(1,i_myend,j,k)/fluxout_a

                        ! transform back:
                        i = i_myend

                        q(2,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*xc(i,j,k) +&
                                ucn_j(2,i,j,k)*xe(i,j,k)  + ucn_j(3,i,j,k)*xz(i,j,k) )
                        q(3,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*yc(i,j,k) +&
                                ucn_j(2,i,j,k)*ye(i,j,k)  + ucn_j(3,i,j,k)*yz(i,j,k) )
                        q(4,i,j,k)=aj(i,j,k)*(ucn_j(1,i,j,k)*zc(i,j,k) +&
                                ucn_j(2,i,j,k)*ze(i,j,k)  + ucn_j(3,i,j,k)*zz(i,j,k) )

                end do; end do
        end if

     call rhs_exchng3_4d (q)      ! update ghost points: q

deallocate ( sbuf1, sbuf2, ucn_j )


contains

  include 'rhs_contra_j.F90'
  include 'rhs_exchng3_4d.F90'


end subroutine conservation



