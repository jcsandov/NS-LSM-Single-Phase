
! Exchange ghost points for rank-3 array
! 
! limited to 1d data decomposition
! 
subroutine rhs_exchng3_3d (var)

  implicit none

  integer :: buf_kmx           ! total size of buf
  integer :: buf_jmx           ! total size of buf
  integer :: buf_imx           ! total size of buf

  ! variable to exchange
  ! 
  real (kind = rdf), dimension (il:,jl:,kl:), intent(inout) :: var

  ! send and receive buffers
  ! 
  real (kind = rdf), dimension (:), allocatable :: bufu
  real (kind = rdf), dimension (:), allocatable :: bufd
  real (kind = rdf), dimension (:), allocatable :: ubuf
  real (kind = rdf), dimension (:), allocatable :: dbuf

  real (kind = rdf), dimension (:), allocatable :: bufr
  real (kind = rdf), dimension (:), allocatable :: bufl
  real (kind = rdf), dimension (:), allocatable :: rbuf
  real (kind = rdf), dimension (:), allocatable :: lbuf

  real (kind = rdf), dimension (:), allocatable :: buff
  real (kind = rdf), dimension (:), allocatable :: bufb
  real (kind = rdf), dimension (:), allocatable :: fbuf
  real (kind = rdf), dimension (:), allocatable :: bbuf

  integer :: l

  ! set size (based on ghost layers)
  ! 
  buf_kmx = ( (iu - il + 1) * (ju - jl + 1) * kgp )
  buf_jmx = ( (iu - il + 1) * (ku - kl + 1) * jgp )
  buf_imx = ( (ju - jl + 1) * (ku - kl + 1) * igp )

  allocate(bufu(buf_kmx), bufd(buf_kmx), &
           ubuf(buf_kmx), dbuf(buf_kmx) )
  allocate(bufr(buf_jmx), bufl(buf_jmx), &
           rbuf(buf_jmx), lbuf(buf_jmx) )
  allocate(buff(buf_imx), bufb(buf_imx), &
           fbuf(buf_imx), bbuf(buf_imx) )

  ! PACK THE DATA GOING TO PROCESSES BACK AND FRONT

  ! pack data going to process back
  !
  if (myback /= mpi_proc_null) then
    
    l = 0
    
    do i = i_mysta , i_mysta + igp - 1
    do k = kl      , ku
    do j = jl      , ju
      
      l = l + 1
      bufb(l) = var(i,j,k)

    end do
    end do
    end do

  end if

  ! pack data going to process front
  ! 
  if (myfront /= mpi_proc_null) then
    
    l = 0
    
    do i = i_myend - igp + 1 , i_myend
    do k = kl                , ku
    do j = jl                , ju
    
      l = l + 1
      buff(l) = var(i,j,k)
    
    end do
    end do
    end do
  
  end if


  ! sync before communication
  ! 
  call mpi_barrier(mpi_comm_world, ierr)


  ! send to myback & receive from myfront
  ! 
  call mpi_sendrecv (bufb, buf_imx, MPI_REAL_TYPE, myback , 299, &
                     fbuf, buf_imx, MPI_REAL_TYPE, myfront, 299, &
                     mpi_comm_world, status, ierr)

  ! send to myfront and receive from myback
  ! 
  call mpi_sendrecv (buff, buf_imx, MPI_REAL_TYPE, myfront, 300, &
                     bbuf, buf_imx, MPI_REAL_TYPE, myback , 300, &
                     mpi_comm_world, status, ierr)


  ! unpack arrray from myback
  ! 
  if (myback /= mpi_proc_null) then
    
    l = 0
    
    do i = il , il + igp - 1
    do k = kl , ku
    do j = jl , ju
      
      l = l + 1
      var(i,j,k) = bbuf(l)
    
    end do
    end do
    end do
  
  end if

  ! unpack array from myfront
  ! 
  if (myfront /= mpi_proc_null) then
    
    l = 0
    
    do i = iu - igp + 1 , iu
    do k = kl           , ku
    do j = jl           , ju
    
      l = l + 1
      var(i,j,k) = fbuf(l)
    
    end do
    end do
    end do
  
  end if

  deallocate(buff, bufb, fbuf, bbuf)

  !-----------------------------------------------------------------------------

  ! PACK THE DATA GOING TO PROCESSES LEFT AND RIGHT

  ! pack data going to process left
  !
  if (myleft /= mpi_proc_null) then
    
    l = 0
    
    do j = j_mysta , j_mysta + jgp - 1
    do k = kl      , ku
    do i = il      , iu
    
      l = l + 1
      bufl(l) = var(i,j,k)
    
    end do
    end do
    end do
  
  end if

  ! pack data going to process right
  ! 
  if (myright /= mpi_proc_null) then
    
    l = 0
    
    do j = j_myend - jgp + 1 , j_myend
    do k = kl                , ku
    do i = il                , iu
    
        l = l + 1
        bufr(l) = var(i,j,k)
    
    end do
    end do
    end do
  
  end if

  ! sync before communication
  ! 
  call mpi_barrier(mpi_comm_world, ierr)


  ! send to myleft & receive from myright
  ! 
  call mpi_sendrecv (bufl, buf_jmx, MPI_REAL_TYPE, myleft , 199, &
                     rbuf, buf_jmx, MPI_REAL_TYPE, myright, 199, &
                     mpi_comm_world, status, ierr)

  ! send to myright and receive from myleft
  ! 
  call mpi_sendrecv (bufr, buf_jmx, MPI_REAL_TYPE, myright, 200, &
                     lbuf, buf_jmx, MPI_REAL_TYPE, myleft , 200, &
                     mpi_comm_world, status, ierr)


  ! unpack arrray from myleft
  ! 
  if (myleft /= mpi_proc_null) then
    
    l = 0
    
    do j = jl , jl + jgp - 1
    do k = kl , ku
    do i = il , iu

        l = l + 1
        var(i,j,k) = lbuf(l)
    
    end do
    end do
    end do

  end if

  ! unpack array from myright
  ! 
  if (myright /= mpi_proc_null) then
    
    l = 0
    
    do j = ju - jgp + 1 , ju
    do k = kl           , ku
    do i = il           , iu
      
      l = l + 1
      var(i,j,k) = rbuf(l)
    
    end do
    end do
    end do
  
  end if

  deallocate(bufr, bufl, rbuf, lbuf)


  !-----------------------------------------------------------------------------

  ! PACK THE DATA GOING TO PROCESSES DOWN AND UP

  ! pack data going to process below
  !
  if (mydown /= mpi_proc_null) then
    
    l = 0
    
    do k = k_mysta , k_mysta + kgp - 1
    do j = jl      , ju
    do i = il      , iu
    
      l = l + 1
      bufd(l) = var(i,j,k)
    
    end do
    end do
    end do
  
  end if

  ! pack data going to process above
  ! 
  if (myup /= mpi_proc_null) then
    
    l = 0
    
    do k = k_myend - kgp + 1 , k_myend
    do j = jl                , ju
    do i = il                , iu
    
      l = l + 1
      bufu(l) = var(i,j,k)
    
    end do
    end do
    end do
  
  end if

  ! sync before communication
  ! 
  call mpi_barrier(mpi_comm_world, ierr)

  ! send to mydown & receive from myup
  ! 
  call mpi_sendrecv (bufd, buf_kmx, MPI_REAL_TYPE, mydown, 99, &
                     ubuf, buf_kmx, MPI_REAL_TYPE, myup  , 99, &
                     mpi_comm_world, status, ierr)

  ! send to myup and receive from mydown
  ! 
  call mpi_sendrecv (bufu, buf_kmx, MPI_REAL_TYPE, myup  , 100, &
                     dbuf, buf_kmx, MPI_REAL_TYPE, mydown, 100, &
                     mpi_comm_world, status, ierr)


  ! *******************
  ! potential bug
  ! 
  !    don't use ka & kb or similar notation below;
  !    because ka & kb, etc. vary depending on
  !    whether we are on an edge process or a inner process!!!
  !
  ! *******************

  ! unpack arrray from mydown
  ! 
  if (mydown /= mpi_proc_null) then
    
    l = 0
    
    do k = kl , kl + kgp - 1
    do j = jl , ju
    do i = il , iu
    
      l = l + 1
      var(i,j,k) = dbuf(l)
  
    end do
    end do
    end do
  
  end if

  ! unpack array from myup
  ! 
  if (myup /= mpi_proc_null) then
    
    l = 0
    
    do k = ku - kgp + 1, ku
    do j = jl, ju
    do i = il, iu
    
      l = l + 1
      var(i,j,k) = ubuf(l)
    
    end do
    end do
    end do
  
  end if

  deallocate(bufu, bufd, ubuf, dbuf)


end subroutine rhs_exchng3_3d
