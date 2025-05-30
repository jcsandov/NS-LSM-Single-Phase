!============================================================
!============================================================
!
program main
  use global
  use global_param
  use global_mpi
  use global_lsm
  use global_debug, only : dirRead_debug

  use phi_inlet 


  implicit none

  real (kind = rdf) :: t1, t2

  integer :: itt
  integer :: l, m
  integer :: is_nan
  
  character (len = 256) :: dirOUT

  call init_mpi ()

#ifdef DEBUG
  print*, 'mpi initialized ', myid
#endif

  call directoryRead() !solo lee donde escribire las soluciones
  call dirRead_debug() !lo mismo que antes pero para las variables de debug (y superficie libre)

  call init ()

  call init_LSM()

  call init_obstacle()


!  if ( rfgt ) call rfg_inlet_init ()

  ! init time variables and solution monitoring
  ! 
  total_time = zero
  nt1        = nt1 + 1
  time0      = time0 + delti
  itt        = 1 + itm(ns) / icn

  !NaN Check
  is_nan = 0

   ! physical time loop
   !
   do ntime = nt1, nt2         

    if ( myid == root ) then

        print *, ' '
        print *, '============================================================'
        print *, ' TIME STEP = ', ntime
        print *, '============================================================'
        print *, ' '

    end if

      time = time0 + (ntime - nt1) * delti

      call init_output_conver ()

      !if ( phi_inlet_prescribed ) call phi_mpi_inlet_forcing ()

      call cpu_time (t1)         ! monitor cpu time

      call mg_driver (ntime)

      call mpi_barrier (mpi_comm_world, ierr)

      call cpu_time (t2)         ! monitor cpu time
      
      ! output required
      ! 
      call output_cpu_performance (total_time) ! report cpu time
      call output_probe ()
      call output_restart_lsm ()  
      call output_restart ()
      call output_real_residual ()

      !nan check
      if( is_nan == 1 ) exit

   end do 
  
  if(myid == root) print *, 'succesful'
  call final_check ()
  call finish_mpi ()            ! shut down mpi cleanly

contains

  ! --

  subroutine output_cpu_performance ( total_time )

    real (kind = rdf)                :: time_step_time
    real (kind = rdf), intent(inout) :: total_time

    time_step_time = t2 - t1
    total_time = total_time + time_step_time

    ! later write to a file instead of STDOUT
    !
    if (myid == root) then
      !print *, 'CPU time = ', time_step_time
       open (unit = 1100, file = 'cpu_performance', position = 'append')
       write(unit = 1100, fmt = '(2(1x,i5.5),2(1X,G15.7))') &
            ntime, itc, time_step_time, total_time
       close(unit = 1100)
    end if

  end subroutine output_cpu_performance
  
  ! --

  subroutine init_output_conver ()

    ! header for convergence file
    !
    if (myid == root) then
       open  (8,   file = 'conver',      position = 'append')
       write (8, *) '-------------------------------------------------------------'
       write (8, *) ' Time = ', time,'  Time step No. = ', ntime
       write (8, *) '-------------------------------------------------------------'
       close (unit = 8)
    end if
  end subroutine init_output_conver

  ! --

  subroutine final_check ()

    use checksum

    implicit none 

    real (kind = rdf), dimension (:,:), allocatable :: qt

    integer :: n
    integer :: k
    integer :: j
    integer :: i
    integer :: la
    integer :: lb

    character (len = 80) :: annotation

    allocate(qt(1:4, 1:li_idx_b(1)))
    
    n = 1
    do k = li_ka(n), li_kb(n)
    do j = li_ja(n), li_jb(n)
    do i = li_ia(n), li_ib(n)
       la = li_idx(i,j,k,n)
       lb = le_idx(i,j,k,n)
       qt(1:4,la) = q(1:4,lb)
    end do
    end do
    end do
    
    write (annotation, fmt = '(a12,i5.5,a9)') &
         'iteration = ', ntime, 'result = '

    call checksum_2d_par (trim(annotation), qt)

    deallocate (qt)

  end subroutine final_check
  
  ! --

  subroutine output_real_residual ()

    implicit none 

    real (kind = rdf), dimension(1:me) :: residual


    ! residual calculations
    !
    call real_residual (1, residual)

    if (myid == root) then
       open  (unit = 9, file = 'conver_real', position = 'append', form = 'formatted')
       write (unit = 9, fmt = '(i5.5,1x,5(g13.6,1x))') ntime, residual
       close (unit = 9)
    end if



   if(myid == root) then
	if(residual(1) /= residual(1)) is_nan = 1 !NaN is not even equal itself
   end if
 
   call mpi_bcast(is_nan,1 ,MPI_INTEGER, root, MPI_COMM_WORLD, ierr)


   

  end subroutine output_real_residual

  ! --

  subroutine real_residual (n, residual)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This routine computes the L2-norm of the residual
    ! on the finest grid at the end of each pseudo time step

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    implicit none

    real (kind = rdf), dimension(1:me) :: residual 

!    real (kind = rdf) :: erp, eru, erv, erw, erx, erw
    real (kind = rdf) :: tnu

    real (kind = rdf), dimension (1:me) :: resd
    real (kind = rdf), dimension (1:me) :: resd_tmp

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

    tnu = real((img(n,myzone) - 2) * (jmg(n,myzone) - 2) * (kmg(n,myzone) - 2), kind = rdf)

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
    resd_tmp = zero
    resd     = zero

    do k = k_mysta, k_myend
    do j = j_mysta, j_myend
    do i = i_mysta, i_myend
       l = le_idx(i,j,k,n)
       do m = 1, me
          resd_tmp(m) = resd_tmp(m) + (q(m,l) - qn(m,l))**2
       end do
    end do
    end do
    end do

    recvcounts(:) = me

    call mpi_reduce (resd_tmp, resd, me, MPI_REAL_TYPE, mpi_sum, &
                     root, mpi_comm_world, ierr)
    !call mpi_reduce (resd_tmp, resd, me, MPI_REAL_TYPE, mpi_sum, &
    !                 root, comm3d, ierr)

    ! call mpi_reduce_scatter (resd_tmp, resd, recvcounts, MPI_REAL_TYPE, mpi_sum,&
    !      & mpi_comm_world, ierr)

    ! calculate average residual
    !
    if ( myid == root ) then

       resd(:) = sqrt(resd(:) / tnu)

       !residual(:) = log10(resd) !zero !log10(resd)

      do m = 1,me

         if ( abs(resd(m)) > eps_sims ) then
            residual(m) = log10(resd(m)) ! zero !log10(resd(:))
         else
            residual(m) = -9999.9
         end if

      end do

!        erp = log10(resd(1))
!        eru = log10(resd(2))
!        erv = log10(resd(3))
!        erw = log10(resd(4))

!        if (me == 5) erx = log10(resd(5))

    else

       residual(:) = zero

    end if

  end subroutine real_residual

  ! --
  
  subroutine output_restart ()

    !real (kind = rdf), dimension(li_ia(1):li_ib(1), &
    !                             li_ja(1):li_jb(1), &
    !                             li_ka(1):li_kb(1) ) :: vartmp
    real (kind = rdf), dimension(:,:,:), allocatable  :: vartmp

    integer :: myunit
    integer :: m
    integer :: la, lb
    integer :: offset

    character (len = 255) :: filename

    ! Steps
    ! =====
    !
    ! * each process opens file
    ! * strip ghost cells from each variable to output
    ! * write entire variable to file one on top of the other
    !   using f90 array notation
    ! 

    allocate(vartmp(li_ia(1):li_ib(1),li_ja(1):li_jb(1),li_ka(1):li_kb(1)))

    la = le_idx_a(1)
    lb = le_idx_b(1)

    if ( (ntime / checkpoint) * checkpoint == ntime ) then
       
       ! check via output file
       !
       offset = 60
       myunit = myid + offset

       write (filename, fmt = '(2a,i3.3)') &
                trim(dirOUT),'solufiles/solu.', myunit - offset
       !write (filename, fmt = '(a,i2.2)') 'solu.', myunit - offset

       open (unit = myunit, file = trim(filename), form = 'unformatted')
       
       do m = 1, me
          call strip_ghost_pts (q(m,la:lb), vartmp)
          write(unit = myunit) vartmp
       end do
       
       if (turbulence) then

          call strip_ghost_pts (xnut(la:lb), vartmp)
          write(unit = myunit) vartmp

       end if 

       do m = 1, me
          call strip_ghost_pts (qn(m,la:lb), vartmp)
          write(unit = myunit) vartmp
       end do

        if (nlinc) then
           do m = 1, 6
              call strip_ghost_pts (uij(m,la:lb), vartmp)
              write(unit = myunit) vartmp
           end do
        end if
       
       close (unit = myunit)

       if ( myid == root ) then
          open (unit = 200, file = 'filestat')
          write(unit = 200, fmt = *) time, ntime
          close(unit = 200)
       end if

    end if

    deallocate (vartmp)
    
  end subroutine output_restart
  
  ! --
  subroutine output_restart_lsm ()


    real (kind = rdf), dimension(:,:,:), allocatable  :: vartmp

    integer :: myunit
    integer :: m
    integer :: la, lb
    integer :: offset

    character (len = 255) :: filename

    ! Steps
    ! =====
    !
    ! * each process opens file
    ! * strip ghost cells from each variable to output
    ! * write entire variable to file one on top of the other
    !   using f90 array notation
    ! 

    allocate(vartmp(li_ia(1):li_ib(1),li_ja(1):li_jb(1),li_ka(1):li_kb(1)))

    la = le_idx_a(1)
    lb = le_idx_b(1)

    if ( (ntime / checkpoint) * checkpoint == ntime ) then
       
       ! check via output file
       !
       offset = 60
       myunit = myid + offset

       write (filename, fmt = '(2a,i3.3)') &
                trim(dirOUT),'solufiles/phi.', myunit - offset
       !write (filename, fmt = '(a,i2.2)') 'solu.', myunit - offset

       open (unit = myunit, file = trim(filename), form = 'unformatted')
       
          call strip_ghost_pts (phi(la:lb), vartmp)
          write(unit = myunit) vartmp

          call strip_ghost_pts (phi_n(la:lb), vartmp)
          write(unit = myunit) vartmp
       
       close (unit = myunit)

    end if

    deallocate (vartmp)
    
  end subroutine output_restart_lsm
  
  ! --

  subroutine strip_ghost_pts ( varin, varout )

    real (kind = rdf), dimension(le_idx_a(1):le_idx_b(1)), intent(in) :: varin
    real (kind = rdf), dimension(:,:,:), intent(out) :: varout

    integer :: &
         i, &
         j, &
         k, &
         l

    do k = li_ka(1), li_kb(1)
    do j = li_ja(1), li_jb(1)
    do i = li_ia(1), li_ib(1)
       l = le_idx(i,j,k,1)
       varout(i,j,k) = varin(l)
    end do
    end do
    end do

  end subroutine strip_ghost_pts
  
  ! --
  subroutine strip_ghost_kplane ( kfl, varin, varout_pl )
  ! this routine strongly depends on applied case

    integer, intent(in) :: kfl
    real (kind = rdf), dimension(:), intent(in) :: varin
    real (kind = rdf), dimension(:,:), intent(out) :: varout_pl
    integer :: i, j, k, l

       k = kfl - gi_ka(1)+1
    do j = li_ja(1), li_jb(1)
    do i = li_ia(1), li_ib(1)
       l = le_idx(i,j,k,1)
       varout_pl(i,j) = varin(l)
    end do
    end do

  end subroutine strip_ghost_kplane
  
  ! --
  subroutine strip_ghost_jplane ( jfl, varin, varout_pl )
  ! this routine strongly depends on applied case

    integer, intent(in) :: jfl
    real (kind = rdf), dimension(:), intent(in) :: varin
    real (kind = rdf), dimension(:,:), intent(out) :: varout_pl
    integer :: i, j, k, l

       j = jfl - gi_ja(1)+1
    do k = li_ka(1), li_kb(1)
    do i = li_ia(1), li_ib(1)
       l = le_idx(i,j,k,1)
       varout_pl(i,k) = varin(l)
    end do
    end do

  end subroutine strip_ghost_jplane

  ! --
  subroutine strip_ghost_iplane ( ifl, varin, varout_pl )
  ! this routine strongly depends on applied case

    integer, intent(in) :: ifl
    real (kind = rdf), dimension(:), intent(in) :: varin
    real (kind = rdf), dimension(:,:), intent(out) :: varout_pl
    integer :: i, j, k, l

       i = ifl - gi_ia(1)+1
    do k = li_ka(1), li_kb(1)
    do j = li_ja(1), li_jb(1)
       l = le_idx(i,j,k,1)
       varout_pl(j,k) = varin(l)
    end do
    end do

  end subroutine strip_ghost_iplane


  ! --

  subroutine output_probe ()
    !
    !============================================================
    ! Probe
    ! =====
    !
    ! at every time step output three planes
    !
    ! 1) Station D1 from Kim expt. (needs interpolation)
    ! 2) Station D2 from Kim expt. (needs interpolation)
    ! 3) \theta = 45\degree plane with cross-flow velocity
    !    components parallel in this plane and axial component
    !    perpendicular
    !
    ! at every icnw time steps output the whole solution for
    ! computing time averages later
    ! 
    real (kind = rdf), dimension(:,:,:), allocatable :: vartmp
    real (kind = rdf), dimension(:,:),   allocatable :: vartmp_kpl
    real (kind = rdf), dimension(:,:),   allocatable :: vartmp_jpl
    real (kind = rdf), dimension(:,:),   allocatable :: vartmp_ipl
    real (kind = rdf), dimension(:,:),   allocatable :: nutmp
    real (kind = rdf), dimension(:,:,:), allocatable :: utmp

    real (kind = rdf), parameter :: th45 = pi / four

    integer :: myunit
    integer :: m
    integer :: la, lb

    integer :: offset           ! use offset for unit numbers
                                ! because unit number = 0 may
                                ! cause a problem
    integer :: l
    integer :: i
    integer :: j
    integer :: k

    integer :: npl

    character (len = 255) :: filename

    ! Assumptions
    ! ==========
    !
    ! each plane is located on the same processor therefore
    ! we don't need to search to find the processor which
    ! owns a particular plane but assume that every
    ! processor owns a part of the plane that needs to be
    ! output
    !
    ! Notes
    ! =====
    !
    ! we don't need to strip the ghost points off because
    ! we do that using the tmp variables which are only
    ! defined on the interior of each processor's nodes
    !
    ! Steps
    ! =====
    !
    ! * for each plane compute the quantities on that plane
    ! * each processor writes a file containing its part
    ! * if right time step, each processor writes its part
    !   of the entire domain

    allocate ( vartmp(li_ia(1):li_ib(1),li_ja(1):li_jb(1),li_ka(1):li_kb(1)) )
    allocate ( vartmp_kpl(li_ia(1):li_ib(1),li_ja(1):li_jb(1)) )
    allocate ( vartmp_jpl(li_ia(1):li_ib(1),li_ka(1):li_kb(1)) )
    allocate ( vartmp_ipl(li_ja(1):li_jb(1),li_ka(1):li_kb(1)) )

    if ( nlinc )      allocate (utmp(1:6,li_ja(1):li_jb(1),li_ka(1):li_kb(1)) )
    if ( turbulence ) allocate (   nutmp(li_ja(1):li_jb(1),li_ka(1):li_kb(1)) )

    ! ----------
    ! monitoring
    do i = 1, monitor_num_points
    if (gi_ka(1) <= monitor_point_ijk(i,3) .and. &
                    monitor_point_ijk(i,3) <= gi_kb(1) .and. &
        gi_ja(1) <= monitor_point_ijk(i,2) .and. &
                    monitor_point_ijk(i,2) <= gi_jb(1) .and. &
        gi_ia(1) <= monitor_point_ijk(i,1) .and. &
                    monitor_point_ijk(i,1) <= gi_ib(1) .and. &
          myzone == monitor_point_ijk(i,4) ) then

        l = gi_2_le_idx(monitor_point_ijk(i,1), &
                        monitor_point_ijk(i,2), &
                        monitor_point_ijk(i,3),1)

       open (unit = 88, file = 'history', position = 'append')
       write(unit = 88, fmt = '(i5.5,1x,6(g12.6,1x))') ntime, q(1:5,l), xnut(l)
       close(unit = 88)
    end if
    end do

    ! ----------------------------------------------------------
    ! store selected planes, selected column or whole flow field
    la = le_idx_a(1)
    lb = le_idx_b(1)
    offset = 70
    myunit = myid + offset


    ! - - - - - - - - -
    ! store selected planes at every time step (CHANGED FOR NOW)
    ! only for Dargahi's flow
    !la=le_idx_a(1)
    !lb=le_idx_b(1)
    !write(filename,fmt='(a,i6.6,a1,i2.2)') &
    !      '/scratch2/cristian/coarse/data/slice/sl.',ntime,'.', myunit - offset


    !open (unit=myunit, file=trim(filename), form='unformatted')

       ! print k-plane
    !   npl=2
    !   if (gi_ka(1) <= npl .and. npl <= gi_kb(1)) then
    !      do m=1,me
    !      call strip_ghost_kplane (npl, q(m,la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !      end do
    !      call strip_ghost_kplane (npl,xnut(la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !   end if

    !   npl=14
    !   if (gi_ka(1) <= npl .and. npl <= gi_kb(1)) then
    !      do m=1,me
    !      call strip_ghost_kplane (npl, q(m,la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !      end do
    !      call strip_ghost_kplane (npl,xnut(la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !   end if

    !   npl=68
    !   if (gi_ka(1) <= npl .and. npl <= gi_kb(1)) then
    !      do m=1,me
    !      call strip_ghost_kplane (npl, q(m,la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !      end do
    !      call strip_ghost_kplane (npl,xnut(la:lb), vartmp_kpl)
    !      write(unit = myunit) vartmp_kpl
    !   end if

       ! depend on nzone
    !   if (myzone == 1) then
    !      !j-plane
    !      npl=(jmg(1,myzone)+1)/2
    !      if (gi_ja(1) <= npl .and. npl <= gi_jb(1)) then
    !         do m=1,me
    !         call strip_ghost_jplane (npl, q(m,la:lb), vartmp_jpl)
    !         write(unit = myunit) vartmp_jpl
    !         end do
    !         call strip_ghost_jplane (npl, xnut(la:lb), vartmp_jpl)
    !         write(unit = myunit) vartmp_jpl
    !      end if
    !   else if (myzone == 2) then
    !      !i-planes
    !      npl=(img(1,myzone))/2
    !      if (gi_ia(1) <= npl .and. npl <= gi_ib(1)) then
    !         do m=1,me
    !         call strip_ghost_iplane (npl, q(m,la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !         end do
    !         call strip_ghost_iplane (npl, xnut(la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !      end if

    !   else if (myzone == 3) then

    !      npl=(img(1,myzone)+1)/2 ! 180 degree
    !      if (gi_ia(1) <= npl .and. npl <= gi_ib(1)) then
    !         do m=1,me
    !         call strip_ghost_iplane (npl, q(m,la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !         end do
    !         call strip_ghost_iplane (npl,xnut(la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !      end if

	!   else if (myzone == 4) then
          !i-planes
    !      npl=(img(1,myzone)+1)/2
    !      if (gi_ia(1) <= npl .and. npl <= gi_ib(1)) then
    !         do m=1,me
    !         call strip_ghost_iplane (npl, q(m,la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !         end do
    !         call strip_ghost_iplane (npl, xnut(la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !      end if

	!   else if (myzone == 5) then
          !i-planes
    !      npl=(img(1,myzone)+1)/2
    !      if (gi_ia(1) <= npl .and. npl <= gi_ib(1)) then
    !         do m=1,me
    !         call strip_ghost_iplane (npl, q(m,la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !         end do
    !         call strip_ghost_iplane (npl, xnut(la:lb), vartmp_ipl)
    !         write(unit = myunit) vartmp_ipl
    !      end if

    !   end if

    !close (unit = myunit)

    ! - - - - - - - - -
    ! store whole flow field (that is usolufiles) at regular intervals
    ! 
    if ((ntime/icnw)*icnw == ntime .or. ntime == current_NSdebug_tstep) then

       la=le_idx_a(1)
       lb=le_idx_b(1)
       !write (filename,fmt='(a,i6.6,a1,i2.2)') &
       !      '/scratch2/cristian/coarse/data/usolufiles/usolu.',ntime,'.', myunit - offset

       write (filename,fmt='(2a,i6.6,a1,i2.2)') &
             trim(dirOUT),'usolufiles/usolu.',ntime,'.', myunit - offset

       open  (unit = myunit, file = trim(filename), form = 'unformatted')

       do m=1,me
          call strip_ghost_pts (q(m,la:lb), vartmp)
          write(unit = myunit) vartmp
       end do
       
       if (turbulence) then
          call strip_ghost_pts (xnut(la:lb), vartmp)
          write(unit = myunit) vartmp
       end if 

       close (unit = myunit)

    end if

    deallocate (vartmp)
    if ( allocated(utmp) ) deallocate (utmp)
    if ( allocated(nutmp)) deallocate (nutmp)


    if (ntime == current_NSdebug_tstep) then
        ! we then update the value of current_LSdebug_tstep
        current_NSdebug_counter = current_NSdebug_counter + 1

        ! check if I already saved all the debug solutions
        if (current_NSdebug_counter <= n_save_NSdebug_tsteps) then
            current_NSdebug_tstep   = save_NSdebug_tsteps(current_NSdebug_counter)
        end if
    end if

  end subroutine output_probe


subroutine directoryRead()

  open(unit = 41, file = "directory", form = "formatted")
    read(unit = 41, fmt = "(a)") dirOUT
  close(unit = 41)

end subroutine directoryRead
  
end program main





