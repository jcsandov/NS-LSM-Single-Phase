        
!==================================================================================
! make_phi_restart: 
! -----------------
!
! compilation: 
! --------------
!
! gfortran -o mk_phi.restart make_phi.restart.f90 str_int.o f2kcli.o
!
!
! execution: 
! ----------
! 
! ./mk_phi.restart nproc 
!
! where nproc is the amount of processors used for the simulation
! 
!==================================================================================
! 

! ifort -o postproc postproc.ussl.f90 f2kcli.o str_int.o
program postprocess

   use f2kcli
   use str_int

  implicit none 

  integer, parameter :: rdf = selected_real_kind(p=7)
  

  character (len = 256) :: line
  character (len = 256) :: exe
  character (len =  40) :: cmd
  character (len = 256) :: filename

  integer, parameter :: nzone = 1 
  integer, parameter :: me = 5
  integer, parameter :: me_save =5
  integer :: icnw

  integer, dimension(nzone) :: im, jm, km
  integer :: nz
  integer :: i
  integer :: j
  integer :: k
  integer :: m
  integer :: ntime

  integer :: narg
  integer :: iarg

  integer :: ntimeStart
  integer :: ntimeEnd
  integer :: ntimeSkip
  integer :: nproc  
  integer :: np

  integer :: myunit
  integer :: offset

  integer :: dummy_int
  real (kind = rdf) :: dummy_real

  integer, dimension (:), allocatable :: &
       ks, &
       ke, &
       js, &
       je, &
       is, &
       ie

  integer :: npa, npb

  integer, dimension (:), allocatable :: myzone, nproc_nz

  integer, dimension (:,:), allocatable :: dims, coords

  character (len = 256) :: debugname


  integer :: num_debug, i_deb
  character (len = 256) :: debug_name

  character (len = 256) :: dirname

  ! get command-line input
  ! 
   narg = command_argument_count ()

   call get_command (line)
   call get_command_argument (0, exe)

  if (narg /= 4) then
	 print *,'nproc,ntimeStart,ntimeEnd,ntimeSkip'
     !print *, 'Need number of processors for ', trim(exe)
     !print *, 'command syntax: ', trim(exe), 'nproc'
     stop
  end if




   call get_command_argument (1, cmd); call decode(nproc, cmd)
   call get_command_argument (2, cmd); call decode(ntimeStart, cmd)
   call get_command_argument (3, cmd); call decode(ntimeEnd, cmd)
   call get_command_argument (4, cmd); call decode(ntimeSkip, cmd)

  open (unit = 21, file = '../ind3dmg.dat', form = 'formatted')
  read (unit = 21, fmt = *) 
  do nz = 1, nzone
  read (unit = 21, fmt = *) im(nz), jm(nz), km(nz)
  end do
  do nz = 1, nzone
  read (unit = 21, fmt = *)
  end do
  !read (unit = 21, fmt = *) 
  !read (unit = 21, fmt = *)
  read (unit = 21, fmt = *) 
  read (unit = 21, fmt = *) 
  read (unit = 21, fmt = *) 
  read (unit = 21, fmt = *) dummy_real, dummy_real, dummy_int, icnw
  close(unit = 21)
  print*, 'icnw=',icnw

  open (unit = 41, file = "../directory", form = "formatted")
	read(unit = 41, fmt = "(a)") dirname
  close (unit = 41)


  ! compute 3d decomposition
  ! =============================
  ! 
  allocate (ks(0:nproc-1), &
            ke(0:nproc-1), &
            js(0:nproc-1), &
            je(0:nproc-1), &
            is(0:nproc-1), &
            ie(0:nproc-1)  )

  allocate (  dims(1:3,0:nproc-1), &
            coords(1:3,0:nproc-1) )

  allocate (nproc_nz(nzone))
  allocate (myzone(0:nproc-1))

  open (2, file = '../nproc.dat')
  read (2,*) (nproc_nz(nz), nz = 1, nzone)
  close(2)

  ! establish the myzone to which this processor belongs
  !
  do np = 0, nproc - 1
     npa = 0
     npb = -1
     do nz = 1, nzone
        npb = npb + nproc_nz(nz)
        if ( np >= npa .and. np <= npb ) then
           myzone(np) = nz
        end if
        npa = npb + 1
     end do
  end do

  ! read dims
  ! read dims
  do np = 0, nproc - 1

     offset = 10
     myunit = np + offset

     write (filename, fmt = '(a,i3.3)') '../decofiles/dims_coords', myunit
     open (unit = myunit, file = trim(filename))
     read (unit = myunit, fmt = *) dims(1:3,np)
     read (unit = myunit, fmt = *) coords(1:3,np)
     close(unit = myunit)

  end do



  do np = 0, nproc - 1
     if (dims(1,np) == 1) then
        is(np) = 1
        ie(np) = im(myzone(np))
     else
        call mpe_decomp1d (im(myzone(np)),dims(1,np),coords(1,np),is(np),ie(np))
     end if
     if (dims(2,np) == 1) then
        js(np) = 1
        je(np) = jm(myzone(np))
     else
        call mpe_decomp1d (jm(myzone(np)),dims(2,np),coords(2,np),js(np),je(np))
     end if
     if (dims(3,np) == 1) then
        ks(np) = 1
        ke(np) = km(myzone(np))
     else
        call mpe_decomp1d (km(myzone(np)),dims(3,np),coords(3,np),ks(np),ke(np))
     end if
  end do








  do ntime = ntimeStart, ntimeEnd, ntimeSkip

	!aqui proceso phi
     if ((ntime / ntimeSkip) * ntimeSkip == ntime) then

 	         write(debugname, fmt ='(a,i6.6)') 'phi',ntime
                 call process_debug_files (debugname)
                !write(debugname, fmt ='(a,i6.6)') 'nband',ntime
                !call process_debug_files (debugname)
     
        print *, 'processing phi',ntime
	print *, " "
     end if

 end do
  deallocate (is, ie, js, je, ks, ke, myzone)

contains

! --

subroutine mpe_decomp1d (n, numprocs, myid, s, e)

  ! This file contains a routine for producing a decomposition of a
  ! 1-d array when given a number of processors.  It may be used in
  ! "direct" product decomposition.  The values returned assume
  ! a "global" domain in [1:n]
  
  ! From the book:
  ! 
  ! "Using MPI", Gropp, Lusk, & Skjellum
  ! 2nd edition, MIT Press, 1999
  
  ! global variables
  !
  integer, intent(in) :: n
  integer, intent(in) :: numprocs
  integer, intent(in) :: myid
  integer, intent(out) :: s
  integer, intent(out) :: e

  ! local variables
  ! 
  integer :: nlocal
  integer :: deficit

  nlocal = n / numprocs
  s = myid * nlocal + 1
  deficit = mod(n,numprocs)
  s = s + min(myid,deficit)
  if (myid < deficit) then
     nlocal = nlocal + 1
  end if
    
  e = s + nlocal - 1
  if (e > n .or. myid == numprocs-1) e = n

end subroutine mpe_decomp1d

! --


! --

subroutine process_debug_files (filenameD)

  ! process usolu.ntime.np restart files
  ! ====================================
  ! 

  integer :: ts

  real (kind = rdf), dimension(:,:,:,:), allocatable :: varD

  real (kind = rdf), dimension(:,:,:), allocatable :: varDtmp
  integer :: imx, jmx, kmx

  character (len = *) :: filenameD
  character (len = 256) :: filenameOUT

  imx = maxval(im)
  jmx = maxval(jm)
  kmx = maxval(km)

  ! process solu.np restart files
  ! =============================
  ! 
  allocate (varD(1:imx,1:jmx,1:kmx,nzone))
  do np = 0, nproc - 1

     allocate (varDtmp(is(np):ie(np),js(np):je(np),ks(np):ke(np))  )

     offset = 70
     myunit = np + offset

     write (filenameOUT,fmt='(2a,a,a,i2.2)') &
                 trim(dirname),'debugfiles/', &
                 trim(filenameD),'.',myunit -offset


     print *,'process file ',trim(filenameD)
     open (unit = myunit, file = trim(filenameOUT), form = 'unformatted')
     
     read (unit = myunit) varDtmp(:,:,:)

     close (unit = myunit)

     ! store this part in global vectors
     ! ---------------------------------
     ! 
     varD(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = varDtmp(:,:,:)
     deallocate (varDtmp)

  end do

  ! changes required in mk_tec etc.
  ! if writing of files stay as below
  ! 
  


  write(filenameOUT, fmt = '(2a,a)') &
                trim(dirname),'surface/',trim(filenameD)

  open  (unit = 22, file = trim(filenameOUT), form = 'unformatted')
  do nz = 1, nzone

        write (unit = 22) varD(1:im(nz),1:jm(nz),1:km(nz),nz)
  end do
  close (unit = 22)


  deallocate (varD )

end subroutine process_debug_files

! --



end program postprocess
