
!============================================================
!
!
! postprocess : post-process output of mpi program by
! -----------   combining into one file the output written 
!               separately by each processor
!               
! author      : J Paik <joongcheol.paik@ce.gatech.edu>
! ------
!
! date        : 2004 September
! ----
! 
! version     : 3.0
! -------
! 
!============================================================
! 
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
  
  character (len = 256) :: dirname

  ! get command-line input
  ! 
  narg = command_argument_count ()

  call get_command (line)
  call get_command_argument (0, exe)

  if (narg /= 1) then
     print *, 'Need number of processors for ', trim(exe)
     print *, 'command syntax: ', trim(exe), 'nproc'
     stop
  end if

  call get_command_argument (1, cmd); call decode(nproc, cmd)

! nproc = 25
print *,'start', nproc
open (unit = 21, file = '../ind3dmg.dat', form = 'formatted')
  read (unit = 21, fmt = *)
  do nz = 1, nzone
  read (unit = 21, fmt = *) im(nz), jm(nz), km(nz)
  end do
  do nz = 1, nzone
  read (unit = 21, fmt = *)
  end do
  read (unit = 21, fmt = *)
  read (unit = 21, fmt = *)
  read (unit = 21, fmt = *)
 ! read (unit = 21, fmt = *)
  read (unit = 21, fmt = *) dummy_real, dummy_real, dummy_int, icnw
  close(unit = 21)

print *,'ind3dmg read',icnw

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
print *,'nproc.dat read'
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
  do np = 0, nproc - 1

     offset = 10 !70
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


  ! postprocessing solu
  call process_solu_files ()

  deallocate (is, ie, js, je, ks, ke)

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

subroutine process_solu_files ()

  real (kind = rdf), dimension(:,:,:,:,:), allocatable :: q
  real (kind = rdf), dimension(:,:,:,:,:), allocatable :: qn
  real (kind = rdf), dimension(:,:,:,:),   allocatable :: nu
  

  real (kind = rdf), dimension(:,:,:,:), allocatable :: qtmp
  real (kind = rdf), dimension(:,:,:,:), allocatable :: qntmp
  real (kind = rdf), dimension(:,:,:),   allocatable :: nutmp


  integer :: imx, jmx, kmx
  
  imx = maxval(im)
  jmx = maxval(jm)
  kmx = maxval(km)

  ! process solu.np restart files
  ! =============================
  ! 
  allocate (q(1:me,1:imx,1:jmx,1:kmx,nzone), &
           qn(1:me,1:imx,1:jmx,1:kmx,nzone), &
                nu(1:imx,1:jmx,1:kmx,nzone)  )

  

  do np = 0, nproc - 1

     allocate (qtmp(1:me,is(np):ie(np),js(np):je(np),ks(np):ke(np)), &
              qntmp(1:me,is(np):ie(np),js(np):je(np),ks(np):ke(np)), &
                   nutmp(is(np):ie(np),js(np):je(np),ks(np):ke(np))  )

  
     offset = 60!
     myunit = np + offset



     !read solu
     !-------------------------------------------------------------------------
     !write (filename, fmt = '(a,i3.3)') '/home/matias/output/solufiles/solu.', myunit - offset

     
     write (filename,fmt='(2a,i3.3)') &
     trim(dirname),'solufiles/solu.', myunit - offset
     open (unit = myunit, file = trim(filename), form = 'unformatted')
          
       do m = 1, me
          read (unit = myunit) qtmp(m,:,:,:)
       end do

       read (unit = myunit) nutmp(:,:,:)

       do m = 1, me
          read (unit = myunit) qntmp(m,:,:,:)
       end do
     
     close (unit = myunit)
     !---------------------------------------------------------------------------

     ! store this part in global vectors
     ! ---------------------------------
     ! 
       q(:,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = qtmp
      qn(:,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = qntmp
        nu(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = nutmp

     deallocate(qtmp, qntmp, nutmp)


  end do

  ! write so they can be read by cduct
  ! 

 
  write(filename, fmt = '(2a)') trim(dirname),'solufiles/solu.restart'
  open  (unit = 22, file = trim(filename), form = 'unformatted')
  do nz = 1, nzone
     do m = 1, me
        write (unit = 22) q(m,1:im(nz),1:jm(nz),1:km(nz),nz)
     end do
     write (unit = 22) nu(1:im(nz),1:jm(nz),1:km(nz),nz)
     do m = 1, me
        write (unit = 22) qn(m,1:im(nz),1:jm(nz),1:km(nz),nz)
     end do
  end do
  close (unit = 22)

  deallocate (q, qn, nu)
  

end subroutine process_solu_files

! --



end program postprocess
