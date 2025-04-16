
subroutine rhs_daf_dtau
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculates dtau(ijk) at each node (local, pseudo time stepping)
  ! eigenvalues of Jacobian matrices

  ! input
  !     ren
  !     cfl1
  !     vnn1
  !     csi(3,ijk)
  !     aj(ijk)
  !     ucn(3,ijk)

  ! output
  !     dtau(ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  real (kind = rdf) :: c1, c2, c3
  real (kind = rdf) :: e1, e2, e3
  real (kind = rdf) :: g11, g22, g33
  real (kind = rdf) :: dti1, dti2, dti3
  real (kind = rdf) :: dtv1, dtv2, dtv3
  real (kind = rdf) :: dcsq, desq, dzsq

  real (kind = rdf) :: cfl, vnn
  real (kind = rdf) :: dtemp
  real (kind = rdf) :: rei, ret

  real (kind = rdf) :: min_aux , min_dtau_global

  integer :: ista , iend
  integer :: jsta , jend
  integer :: ksta , kend

  integer :: i,j,k

  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 

  if (myback  == mpi_proc_null)  ista = il + igp 
  if (myleft  == mpi_proc_null)  jsta = jl + jgp 
  if (mydown  == mpi_proc_null)  ksta = kl + kgp 

  if (myfront == mpi_proc_null)  iend = iu - igp
  if (myright == mpi_proc_null)  jend = ju - jgp
  if (myup    == mpi_proc_null)  kend = ku - kgp


  cfl=cfl1(myzone)
  vnn=vnn1(myzone)
  rei = one / ren ! 1/Re

  ! we have made the inverse of the grid spacing
  ! universal; therefore, dsp = dc^(-1) in old code
  ! this caused a bug in the translations

  ! Arnone et al. (1995)
  dtemp = delti / (one_pt_five * two ** (3 - 1))

  dcsq=dc*dc
  desq=de*de
  dzsq=dz*dz

  ! Note, that each direction is treated independently

  !do k=kl,ku
  !do j=jl,ju
  !do i=il,iu

  min_aux = ten * ten ;

  do k = ksta , kend
  do j = jsta , jend
  do i = ista , iend
         
      !if( rsign(i,j,k) > one_half ) then
      
         g11 = csi(1,i,j,k) * csi(1,i,j,k) + &
               csi(2,i,j,k) * csi(2,i,j,k) + &
               csi(3,i,j,k) * csi(3,i,j,k)

         g22 = eta(1,i,j,k) * eta(1,i,j,k) + &
               eta(2,i,j,k) * eta(2,i,j,k) + &
               eta(3,i,j,k) * eta(3,i,j,k)

         g33 = zet(1,i,j,k) * zet(1,i,j,k) + &
               zet(2,i,j,k) * zet(2,i,j,k) + &
               zet(3,i,j,k) * zet(3,i,j,k)
      
         !
         !
         !

         !dti1=cfl/sqrt(g11)/dc
         !dti2=cfl/sqrt(g22)/de
         !dti3=cfl/sqrt(g33)/dz

         if ( dynamic_dtau ) then
            
            ! Compute ∆τ using the local velocities
            c1 = abs( ucn_j(1,i,j,k) * aj(i,j,k) ) + sqrt( (ucn_j(1,i,j,k)*aj(i,j,k))**two + beta * g11)            
            c2 = abs( ucn_j(2,i,j,k) * aj(i,j,k) ) + sqrt( (ucn_j(2,i,j,k)*aj(i,j,k))**two + beta * g22)            
            c3 = abs( ucn_j(3,i,j,k) * aj(i,j,k) ) + sqrt( (ucn_j(3,i,j,k)*aj(i,j,k))**two + beta * g33)            
         
         else

            ! Compute ∆τ using just the local grid size
            c1 = sqrt( g11 )
            c2 = sqrt( g22 )
            c3 = sqrt( g33 )
         
         end if  

         dti1 = cfl / c1 / dc
         dti2 = cfl / c2 / de
         dti3 = cfl / c3 / dz

!         ret = rei*muLSM(phi_n(i,j,k))/rhoLSM(phi_n(i,j,k))+xnut(i,j,k)
         
         ! I put this back as the original formulation
         ret = rei + abs( xnut(i,j,k) )      
         
         dtv1 = vnn / ret / dcsq / g11
         dtv2 = vnn / ret / desq / g22
         dtv3 = vnn / ret / dzsq / g33
      
         dtau(i,j,k) = min( dti1 , dti2 , dti3 , dtv1 , dtv2 , dtv3 )
         !min_aux = min( min_aux , dti1 , dti2 , dti3 , dtv1 , dtv2 , dtv3 )
         !if (unsteady) dtau(i,j,k) = min(dtau(i,j,k), dtemp)

      !end if

  end do
  end do
  end do

  !call mpi_allreduce( min_aux, min_dtau_global , &
  !                    1 , MPI_REAL_TYPE , mpi_min , mpi_comm_world, ierr )

  ! I use the same dtau for the whole domain
  !dtau = min_dtau_global

  ! blanking area
  !
  if (nblk /= 0) then
  do nb = 1, nblk
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        dtau(i,j,k) = zero
     end do
     end do
     end do
  end do
  end if

!call outputD3_real(dtau(:,:,:),'dtau_dmod')
!if(myid == 0) print *, 'output_dtau0'
end subroutine rhs_daf_dtau


