program get_tseries

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! This routine read usolufiles, phi files and all the debug files declared in debugnames.txt
   ! and extract tseries of the points given within tseries_control.txt
   !  
   !
   ! COMPILATION
   ! ---------------
   ! The compilation of this files is as follows 
   !
   ! gfortran -o mktseries get_tseries.f90 f2kcli.o str_int.o
   ! 
   ! * (f2kcli.f90 and str_int.f90 must be in the same folder where compilation is made)
   !
   ! INPUTS
   ! ---------
   ! 
   ! The input files for this routines are (considering this routine is located in a folder called
   ! postproc, which is located in the main folder of the simulations):
   !
   ! 1. ind3dmg.dat (from ../ind3dmg.dat) --> main information about grid dimensions
   ! 2. grid (from ../grid) --> grid coordinates
   ! 3. directory (from ../directory) --> to obtain the output folder path
   ! 4. nproc.dat (from ../nproc.dat) --> to obtain the processors distribution
   ! 5. dim_coords files (from ../decofiles/dim_coords) --> to obtain the subdomain of every processor
   ! 6. tseries_control.txt (in the same folder of this routine) --> to obtain the total amount of 
   !                                                                 variables (and their names) and 
   !                                                                 total amount of points 
   !
   ! RUNNING
   ! ---------
   ! 
   ! ./mktseries nproc ntimeStart ntimeEnd ntimeSkip
   !
   ! For example, if we have a simulation that was run using 32 procs and I want to generate
   ! output debug files between time steps 100 and 200 each 10 time steps, it has to be run as
   ! follows
   ! 
   ! ./mktseries 32 100 200 10
   !  
   ! The output files are called debugvars000100.plt, debugvars000110.plt, etc
   ! 
   ! and they contains the following variables:
   ! X, Y, Z, P, U, V, W, PHI, DV1, DV2, .... DVN
   !
   ! where DV1, DV2, ... are the names of the variables declared in debugnames.txt
   !
   ! Jorge Sandoval, UoE/PUC. Edinburgh, February 13th, 2025.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
   integer :: counter
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
  
   real (kind = rdf) , dimension(:,:,:,:), allocatable :: x,y,z,wd

   ! debugging variables

   integer :: numdebugvars ! total number of debugged variables
   character (len = 256), dimension(30) :: debugnames
   integer :: ideb

   ! out folder path
   character (len = 256) :: dirname

   ! flags to export variables
   logical :: pressureflag, uflag, vflag, wflag, nutflag, phiflag

   ! indexes for the time-series extraction
   integer, allocatable :: points(:,:)
   integer :: num_points

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
      print *, 'nz = ', nz, ', im = ', im(nz), ', jm = ', jm(nz), ', km =  ', km(nz)
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


   allocate(  x(1:maxval(im),1:maxval(jm),1:maxval(km),1:nzone) ,  &
              y(1:maxval(im),1:maxval(jm),1:maxval(km),1:nzone) ,  &
              z(1:maxval(im),1:maxval(jm),1:maxval(km),1:nzone) ,  &
             wd(1:maxval(im),1:maxval(jm),1:maxval(km),1:nzone))

   ! Grid reading
   open(1,file="../grid",form="unformatted")
   
   do nz = 1, nzone
      read(1) (((x(i,j,k,nz),  i = 1,im(nz)), j = 1,jm(nz)), k = 1,km(nz))
      read(1) (((y(i,j,k,nz),  i = 1,im(nz)), j = 1,jm(nz)), k = 1,km(nz))
      read(1) (((z(i,j,k,nz),  i = 1,im(nz)), j = 1,jm(nz)), k = 1,km(nz))
      read(1) (((wd(i,j,k,nz), i = 1,im(nz)), j = 1,jm(nz)), k = 1,km(nz))
   end do
   
   close(1)

   ! Here we read the path to the output folder

   open (unit = 41, file = "../directory", form = "formatted")
      read(unit = 41, fmt = "(a)") dirname
   close (unit = 41)

   ! =============================
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

   ! I read all the variables to be exported
   call read_input_file( )

   ! Loop exporting tec files
   do ntime = ntimeStart, ntimeEnd, ntimeSkip
      print *, 'ts = ', ntime
      call process_comb_files_debug (ntime)
   end do

  deallocate (x, y, z, wd, is, ie, js, je, ks, ke, myzone)

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

subroutine process_comb_files_debug_old (ts)
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This is the subroutine that reads every file for the tecplot master file construction
   ! it receives the time iteration (ts) and read, gather and write the processed files
   !
   ! steps: 
   ! 1. usolufiles reading from output/usolufiles
   ! 2. phi files reading from output/debugfiles
   ! 3. debug files reading from output/debugfiles
   ! 4. gathering of all the data in the tecplotfiles
   ! 5. writing of .plt files onto output/debug
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   integer :: ts

   real (kind = rdf), dimension(:,:,:,:,:) ,  allocatable :: q
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: nu
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: phi

   !---------------------------------------------------------------------
   ! Debug Variables

   real (kind = rdf), dimension(:,:,:,:,:), allocatable :: vardebug_master
   character (len = 512) :: tecvarnames

   !---------------------------------------------------------------------
   !  TEMPORAL ARRAYS
   !---------------------------------------------------------------------

   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: qtmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: nutmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: phitmp

   !---------------------------------------------------------------------
   ! Temporal arrays for debug variables

   real (kind = rdf), dimension(:,:,:,:), allocatable :: vardebug_master_tmp

   !---------------------------------------------------------------------

   integer :: imx, jmx, kmx
   integer :: idummy

   !-------------------------------------------------
   ! TECPLOT VARIABLES
   !- - - - - - - - - - - - - - - - - - - - - - - - - 
   ! character (len = 256) :: filename
   ! variables to enable writing of TecPlot binary (*.plt) files
   
   integer (kind = 4)           :: TecIni, TecDat, TecZne
   integer (kind = 4)           :: TecEnd
   integer (kind = 4)           :: VIsDouble = 0
   integer (kind = 4)           :: Debug = 1
   
   integer (kind = 4)           :: III
   
   character (len = 1)          :: nullchr = char(0)
   !-----------------------------------------------------


   imx = maxval(im)
   jmx = maxval(jm)
   kmx = maxval(km)


   allocate (q(1:me,1:imx,1:jmx,1:kmx,nzone)  ,  &
             nu(1:imx,1:jmx,1:kmx,nzone)      ,  &
             phi(1:imx,1:jmx,1:kmx,nzone))

   allocate (vardebug_master(1:numdebugvars,1:imx,1:jmx,1:kmx,nzone))
   

   q                    =  0.0_rdf
   nu                   =  0.0_rdf
   phi                  =  0.0_rdf
   vardebug_master      =  0.0_rdf
   
   !------------------------------------------------------------------------------------------------
   ! usolufiles reading
   do np = 0, nproc - 1

      allocate (qtmp(1:me,is(np):ie(np),js(np):je(np),ks(np):ke(np)),  &
                   nutmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))

      qtmp   =  0.0_rdf
      nutmp  =  0.0_rdf

      offset = 70
      myunit = np + offset

      write (filename,fmt='(2a,i6.6,a1,i2.2)') &
            trim(dirname),'usolufiles/usolu.',ts,'.', myunit - offset


      print *,'processing ', filename
      open (unit = myunit, file = trim(filename), form = 'unformatted')
     
      do m = 1, me
         read (unit = myunit) qtmp(m,:,:,:)
      end do
      read (unit = myunit) nutmp(:,:,:)

      close (unit = myunit)

      ! store this part in global vectors
      ! ---------------------------------

      q(:,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = qtmp(:,:,:,:)
      nu(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  nutmp(:,:,:)
     
      deallocate (qtmp, nutmp)

   end do

   !------------------------------------------------------------------------------------------------
   ! phi files reading

   do np = 0, nproc - 1

      allocate (phitmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))

      phitmp = 0.0_rdf

      offset = 60
      myunit = np + offset

      write (filename,fmt='(2a,i6.6,a1,i2.2)') &
            trim(dirname),'debugfiles/phi',ts,'.', myunit - offset


      print *,'processing ', filename
      open (unit = myunit, file = trim(filename), form = 'unformatted')
      read (unit = myunit) phitmp(:,:,:)
      close (unit = myunit)

      ! ---------------------------------
      ! store this part in global vectors
      ! --------------------------------- 
      
      ! incorporating all the data from the processor np to the global array

      phi(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  phitmp(:,:,:)

      deallocate (phitmp)

   end do

   ! --------------------------------------------------------------------------------------------
   ! Reading every debug variable
  
   do np = 0, nproc - 1

      allocate (vardebug_master_tmp(1:numdebugvars,is(np):ie(np),js(np):je(np),ks(np):ke(np)))
   
      vardebug_master_tmp  =  0.0_rdf

      offset = 60
      myunit = np + offset


      do ideb = 1, numdebugvars

         write (filename,fmt='(3a,i6.6,a1,i2.2)') &
               trim(dirname),'debugfiles/',trim(debugnames(ideb)),ts,'.', myunit - offset

         print *,'processing ', filename
         open (unit = myunit, file = trim(filename), form = 'unformatted')     
         read (unit = myunit) vardebug_master_tmp(ideb,:,:,:)
         close (unit = myunit)

      end do

      ! ---------------------------------
      ! store this part in global vectors
      ! ---------------------------------
      ! 

      do ideb = 1, numdebugvars
         vardebug_master(ideb,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  vardebug_master_tmp(ideb,:,:,:)
      end do

      deallocate (vardebug_master_tmp)

   end do

   ! -------------------------------------------------------------------------------------------- 
   ! Writing Tecplot File
   
   write(filename, fmt = '(2a,i6.6)') &
        trim(dirname),'debug/debugvars',ts  

   ! tecvarnames is the string to generate the variable set for tecplot file
   tecvarnames = ''

   ! concatenation loop: we format and concatenate all the variables name from debugnames.txt
   do ideb = 1, numdebugvars-1
      write(tecvarnames,fmt = '(3a)') &
            trim(tecvarnames),trim(debugnames(ideb)),','
   end do

   ! the last name doesn't need a final comma
   write(tecvarnames,fmt = '(2a)') &
         trim(tecvarnames),trim(debugnames(numdebugvars))

   ! TO DO: in the future, we could generate in the debugnames.txt file an instance where
   ! we can choose which variables from q or phi vectors we want 

   I = TecIni('debugfile'//NULLCHR,    &   ! title of file
     'X, Y, Z, P, U, V, W, PHI,'// trim(tecvarnames) //NULLCHR,   &   ! list of variables
     trim(filename)//'.plt'//NULLCHR, &   ! output file name
     '.'//NULLCHR,                 &
     Debug,                        &
     VIsDouble)
   
   do nz = 1, nzone
   
      I = TecZne('Zone'//NULLCHR,    &     
      im(nz),                  &
      jm(nz),                  &
      km(nz),                  &
      'BLOCK'//NULLCHR,       &
      NULLCHR//NULLCHR)
      
      ! total number of points
      III = km(nz) * jm(nz) * im(nz)
      
      I   = TecDat(III,x(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,y(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,z(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(1,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(2,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(3,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(4,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,phi(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      
      ! writing the debug variables
      
      do ideb = 1, numdebugvars
         I   = TecDat(III,vardebug_master(ideb,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      end do      

   end do !nz

   I   = TecEnd()

   deallocate (q, nu, phi, vardebug_master)

end subroutine process_comb_files_debug_old

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

subroutine process_comb_files_debug (ts)
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This is the subroutine that reads every file for the tecplot master file construction
   ! it receives the time iteration (ts) and read, gather and write the processed files
   !
   ! steps: 
   ! 1. usolufiles reading from output/usolufiles
   ! 2. phi files reading from output/debugfiles
   ! 3. debug files reading from output/debugfiles
   ! 4. gathering of all the data in the tecplotfiles
   ! 5. writing of .plt files onto output/debug
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   integer :: ts

   real (kind = rdf), dimension(:,:,:,:,:) ,  allocatable :: q
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: nu
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: phi

   !---------------------------------------------------------------------
   ! Debug Variables

   real (kind = rdf), dimension(:,:,:,:,:), allocatable :: vardebug_master
   character (len = 256) :: tseries_filename
   character (len = 256) :: tmp_string
   logical :: FileExist
   integer :: iout , jout , kout
   real (kind = rdf) :: var

   !---------------------------------------------------------------------
   !  TEMPORAL ARRAYS
   !---------------------------------------------------------------------

   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: qtmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: nutmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: phitmp

   !---------------------------------------------------------------------
   ! Temporal arrays for debug variables

   real (kind = rdf), dimension(:,:,:,:), allocatable :: vardebug_master_tmp

   !---------------------------------------------------------------------

   integer :: imx, jmx, kmx
   integer :: idummy
   integer :: ipoints

   !-------------------------------------------------
   ! TECPLOT VARIABLES
   !- - - - - - - - - - - - - - - - - - - - - - - - - 
   ! character (len = 256) :: filename
   ! variables to enable writing of TecPlot binary (*.plt) files
   
   integer (kind = 4)           :: TecIni, TecDat, TecZne
   integer (kind = 4)           :: TecEnd
   integer (kind = 4)           :: VIsDouble = 0
   integer (kind = 4)           :: Debug = 1
   
   integer (kind = 4)           :: III
   
   character (len = 1)          :: nullchr = char(0)
   !-----------------------------------------------------

   imx = maxval(im)
   jmx = maxval(jm)
   kmx = maxval(km)

   ! if I don't need any information from usolufiles, I don't even
   ! allocate memory for them

   if ( pressureflag .or. uflag .or. vflag .or. wflag .or. nutflag ) then
      allocate (q(1:me,1:imx,1:jmx,1:kmx,nzone)  ,  &
                nu(1:imx,1:jmx,1:kmx,nzone)            )
      
      ! variables initialisation
      q                    =  0.0_rdf
      nu                   =  0.0_rdf

   end if

   ! the same for phi
   if ( phiflag ) then
      allocate ( phi(1:imx,1:jmx,1:kmx,nzone) )
                 phi = 0.0_rdf
   end if

   if ( numdebugvars > 0 ) then
      allocate (vardebug_master(1:numdebugvars,1:imx,1:jmx,1:kmx,nzone))
                vardebug_master =  0.0_rdf
   end if
      
   !------------------------------------------------------------------------------------------------
   ! usolufiles reading

   if ( pressureflag .or. uflag .or. vflag .or. wflag .or. nutflag ) then
      do np = 0, nproc - 1
   
         allocate (qtmp(1:me,is(np):ie(np),js(np):je(np),ks(np):ke(np)),  &
                      nutmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))
   
         qtmp   =  0.0_rdf
         nutmp  =  0.0_rdf
   
         offset = 70
         myunit = np + offset
   
         write (filename,fmt='(2a,i6.6,a1,i2.2)') &
               trim(dirname),'usolufiles/usolu.',ts,'.', myunit - offset
   
   
         print *,'processing ', filename
         open (unit = myunit, file = trim(filename), form = 'unformatted')
        
         do m = 1, me
            read (unit = myunit) qtmp(m,:,:,:)
         end do
         read (unit = myunit) nutmp(:,:,:)
   
         close (unit = myunit)
   
         ! store this part in global vectors
         ! ---------------------------------
   
         q(:,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = qtmp(:,:,:,:)
         nu(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  nutmp(:,:,:)
        
         deallocate (qtmp, nutmp)
   
      end do
   end if

   !------------------------------------------------------------------------------------------------
   ! phi files reading

   if ( phiflag ) then
      do np = 0, nproc - 1
   
         allocate (phitmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))
   
         phitmp = 0.0_rdf
   
         offset = 60
         myunit = np + offset
   
         write (filename,fmt='(2a,i6.6,a1,i2.2)') &
               trim(dirname),'debugfiles/phi',ts,'.', myunit - offset
   
   
         print *,'processing ', filename
         open (unit = myunit, file = trim(filename), form = 'unformatted')
         read (unit = myunit) phitmp(:,:,:)
         close (unit = myunit)
   
         ! ---------------------------------
         ! store this part in global vectors
         ! --------------------------------- 
         
         ! incorporating all the data from the processor np to the global array
   
         phi(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  phitmp(:,:,:)
   
         deallocate (phitmp)
   
      end do
   end if

   ! --------------------------------------------------------------------------------------------
   ! Reading every debug variable

   if ( numdebugvars > 0 ) then     
      do np = 0, nproc - 1
   
         allocate (vardebug_master_tmp(1:numdebugvars,is(np):ie(np),js(np):je(np),ks(np):ke(np)))
      
         vardebug_master_tmp  =  0.0_rdf
   
         offset = 60
         myunit = np + offset
   
   
         do ideb = 1, numdebugvars
   
            write (filename,fmt='(3a,i6.6,a1,i2.2)') &
                  trim(dirname),'debugfiles/',trim(debugnames(ideb)),ts,'.', myunit - offset
   
            print *,'processing ', filename
            open (unit = myunit, file = trim(filename), form = 'unformatted')     
            read (unit = myunit) vardebug_master_tmp(ideb,:,:,:)
            close (unit = myunit)
   
         end do
   
         ! ---------------------------------
         ! store this part in global vectors
         ! ---------------------------------
         ! 
   
         do ideb = 1, numdebugvars
            vardebug_master(ideb,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  vardebug_master_tmp(ideb,:,:,:)
         end do
   
         deallocate (vardebug_master_tmp)
   
      end do
   end if

   ! -------------------------------------------------------------------------------------------- 
   ! Writing time series 

   ! chenge if there more than 1 zone
   nz = 1

   if ( pressureflag ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('p', iout, jout , kout , tmp_string )

         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = q(1,iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)

      end do

   end if

   if ( uflag        ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('u', iout, jout , kout , tmp_string )

         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = q(2,iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)



      end do


   end if

   if ( vflag        ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('u', iout, jout , kout , tmp_string )
         
         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = q(3,iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)


      end do

   end if

   if ( wflag        ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('w', iout, jout , kout , tmp_string )

         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = q(4,iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)


      end do

   
   end if

   if ( nutflag      ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('nut', iout, jout , kout , tmp_string )

         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = nu(iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)


      end do

   
   end if

   if ( phiflag      ) then

      do ipoints = 1 , num_points

         iout = points(ipoints,1)
         jout = points(ipoints,2)
         kout = points(ipoints,3)

         call format_string('phi', iout, jout , kout , tmp_string )

         tseries_filename = trim(tmp_string)

         print *, 'processing ', tseries_filename

         var = phi(iout,jout,kout,nz)

         ! Writting the file        
         inquire(file = tseries_filename, exist = FileExist)
                  
         if (FileExist) then
            open( 12 , file     = tseries_filename   , status = "old"      , &
                       position = "append"           , action = "write"    )
         else
            
            open( 12 , file     = tseries_filename   , status = "new"      , & 
                       action = "write"                                    )
         end if
                  
         write(12, *) var
         
         close(12)


      end do


   end if

   ! concatenation loop: we format and concatenate all the variables name from debugnames.txt
   if( numdebugvars > 0) then
      do ideb = 1, numdebugvars

         do ipoints = 1 , num_points


            iout = points(ipoints,1)
            jout = points(ipoints,2)
            kout = points(ipoints,3)
            
            call format_string( trim( debugnames(ideb) ), iout, jout , kout , tmp_string )

            tseries_filename = trim(tmp_string)

            print *, 'processing ', tseries_filename

            var = vardebug_master(ideb,iout,jout,kout,nz)

            ! Writting the file        
            inquire(file = tseries_filename, exist = FileExist)
                  
            if (FileExist) then
               open( 12 , file     = tseries_filename   , status = "old"      , &
                          position = "append"           , action = "write"    )
            else
            
               open( 12 , file     = tseries_filename   , status = "new"      , & 
                          action = "write"                                    )
            end if
                  
            write(12, *) var
         
            close(12)

         end do
      end do
   end if


   if ( allocated (q)               ) deallocate ( q )
   if ( allocated (nu)              ) deallocate ( nu )
   if ( allocated (phi)             ) deallocate ( phi )
   if ( allocated (vardebug_master) ) deallocate ( vardebug_master )

end subroutine process_comb_files_debug

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

subroutine process_comb_files (ts)
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This is the subroutine that reads every file for the tecplot master file construction
   ! it receives the time iteration (ts) and read, gather and write the processed files
   !
   ! steps: 
   ! 1. usolufiles reading from output/usolufiles
   ! 2. phi files reading from output/debugfiles
   ! 3. gathering of all the data in the tecplotfiles
   ! 4. writing of .plt files onto output/debug
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   integer :: ts

   real (kind = rdf), dimension(:,:,:,:,:) ,  allocatable :: q
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: nu
   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: phi

   !---------------------------------------------------------------------
   !  TEMPORAL ARRAYS
   !---------------------------------------------------------------------

   real (kind = rdf), dimension(:,:,:,:)   ,  allocatable :: qtmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: nutmp
   real (kind = rdf), dimension(:,:,:)     ,  allocatable :: phitmp

   integer :: imx, jmx, kmx
   integer :: idummy

   !-------------------------------------------------
   ! TECPLOT VARIABLES
   !- - - - - - - - - - - - - - - - - - - - - - - - - 
   ! character (len = 256) :: filename
   ! variables to enable writing of TecPlot binary (*.plt) files
   
   integer (kind = 4)           :: TecIni, TecDat, TecZne
   integer (kind = 4)           :: TecEnd
   integer (kind = 4)           :: VIsDouble = 0
   integer (kind = 4)           :: Debug = 1
   
   integer (kind = 4)           :: III
   
   character (len = 1)          :: nullchr = char(0)
   !-----------------------------------------------------


   imx = maxval(im)
   jmx = maxval(jm)
   kmx = maxval(km)


   allocate (q(1:me,1:imx,1:jmx,1:kmx,nzone)  ,  &
             nu(1:imx,1:jmx,1:kmx,nzone)      ,  &
             phi(1:imx,1:jmx,1:kmx,nzone))

   q                    =  0.0_rdf
   nu                   =  0.0_rdf
   phi                  =  0.0_rdf

   !------------------------------------------------------------------------------------------------
   ! usolufiles reading
   do np = 0, nproc - 1

      allocate (qtmp(1:me,is(np):ie(np),js(np):je(np),ks(np):ke(np)),  &
                   nutmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))

      qtmp   =  0.0_rdf
      nutmp  =  0.0_rdf

      offset = 70
      myunit = np + offset

      write (filename,fmt='(2a,i6.6,a1,i2.2)') &
            trim(dirname),'usolufiles/usolu.',ts,'.', myunit - offset


      print *,'processing ', filename
      open (unit = myunit, file = trim(filename), form = 'unformatted')
     
      do m = 1, me
         read (unit = myunit) qtmp(m,:,:,:)
      end do
      read (unit = myunit) nutmp(:,:,:)

      close (unit = myunit)

      ! store this part in global vectors
      ! ---------------------------------

      q(:,is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) = qtmp(:,:,:,:)
      nu(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  nutmp(:,:,:)
     
      deallocate (qtmp, nutmp)

   end do

   !------------------------------------------------------------------------------------------------
   ! phi files reading

   do np = 0, nproc - 1

      allocate (phitmp(is(np):ie(np),js(np):je(np),ks(np):ke(np)))

      phitmp = 0.0_rdf

      offset = 60
      myunit = np + offset

      write (filename,fmt='(2a,i6.6,a1,i2.2)') &
            trim(dirname),'debugfiles/phi',ts,'.', myunit - offset


      print *,'processing ', filename
      open (unit = myunit, file = trim(filename), form = 'unformatted')
      read (unit = myunit) phitmp(:,:,:)
      close (unit = myunit)

      ! ---------------------------------
      ! store this part in global vectors
      ! --------------------------------- 
      
      ! incorporating all the data from the processor np to the global array

      phi(is(np):ie(np),js(np):je(np),ks(np):ke(np),myzone(np)) =  phitmp(:,:,:)

      deallocate (phitmp)

   end do


   ! -------------------------------------------------------------------------------------------- 
   ! Writing Tecplot File
   
   write(filename, fmt = '(2a,i6.6)') &
        trim(dirname),'debug/debugvars',ts  

   I = TecIni('debugfile'//NULLCHR,    &   ! title of file
     'X, Y, Z, P, U, V, W, PHI'//NULLCHR,   &   ! list of variables
     trim(filename)//'.plt'//NULLCHR, &   ! output file name
     '.'//NULLCHR,                 &
     Debug,                        &
     VIsDouble)
   
   do nz = 1, nzone
   
      I = TecZne('Zone'//NULLCHR,    &     
      im(nz),                  &
      jm(nz),                  &
      km(nz),                  &
      'BLOCK'//NULLCHR,       &
      NULLCHR//NULLCHR)
      
      ! total number of points
      III = km(nz) * jm(nz) * im(nz)
      
      I   = TecDat(III,x(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,y(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,z(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(1,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(2,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(3,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,q(4,1:im(nz),1:jm(nz),1:km(nz),nz),0)
      I   = TecDat(III,phi(1:im(nz),1:jm(nz),1:km(nz),nz),0)
      
   end do !nz

   I   = TecEnd()

   deallocate (q, nu, phi)

end subroutine process_comb_files

function remove_comma(string)

   character(len=*), intent(in) :: string
   character(len=len_trim(string)) :: new_string
   character(len=len_trim(string)) :: remove_comma

   if (string(len_trim(string):len_trim(string)) == ',') then

     new_string = adjustl(string(1:len_trim(string)-1))

   else

     new_string = string

   end if

   remove_comma = new_string

end function remove_comma

subroutine read_input_file( )

   character(512) :: line, var_name, value_str
   integer :: pos, iunit, ierr, stat, counter

   debugnames = ''

   ! Open input file
   open(unit=10, file='tseries_control.txt', status='old', action='read', iostat=ierr)

   if (ierr /= 0) then
      print *, "Error opening input file!"
      stop
   end if

   ! Read variable values
   do while (.not. is_iostat_end(stat) )

      read(10, '(A)' , iostat = stat) line
      
      line      = trim(line)    
      pos       = index(line, "=")
      var_name  = trim(line(1:pos-1))
      value_str = trim(line(pos+1:))
        
      select case (trim(var_name))
         case ("pressure")
            read(value_str, *) pressureflag
         case ("u-velocity")
            read(value_str, *) uflag
         case ("v-velocity")
            read(value_str, *) vflag
         case ("w-velocity")
            read(value_str, *) wflag
         case ("nu-t")
            read(value_str, *) nutflag
         case ("phi")
            read(value_str, *) phiflag
         case ("Number of variables")
            read(value_str, *) numdebugvars
         case ("Number of points")
            read(value_str, *) num_points            
      end select

      if (line == "NAME OF THE VARIABLES TO BE EXPORTED") then 
         if( numdebugvars>0 ) then
            counter = 1
            do while (.not. is_iostat_end(stat) .and. counter <= numdebugvars)
               read(10, '(A)' , iostat = stat) line
               if(trim(line) == '' ) then
                  cycle
               else
                  line = trim(line)    
                  debugnames(counter) = line
                  counter = counter + 1
               end if
            end do
         end if
      end if

      if (line == "LIST OF POINTS TO EXTRACT TIME SERIES") then
         
         if (num_points > 0) then
         
            allocate(points(num_points,3))
         
            counter = 1
            do while (.not. is_iostat_end(stat) .and. counter <= num_points)
         
               read(10, *, iostat = stat) points(counter,1), points(counter,2), points(counter,3)
               
               if (stat /= 0) then
                  cycle  ! Skip to the next iteration if an error occurs
               end if
         
               counter = counter + 1
         
            end do
         
         end if
      
      end if

    end do

   ! Close input file
   close(unit=10)

   ! Print variable values
   print *, "pressure     = ", pressureflag
   print *, "u-velocity   = ", uflag
   print *, "v-velocity   = ", vflag
   print *, "w-velocity   = ", wflag
   print *, "nu-t         = ", nutflag
   print *, "phi          = ", phiflag
   print *, "numdebugvars = ", numdebugvars
   print *, ' '
   print *, 'debugnames: '
   print *, ' '

   do counter = 1, numdebugvars
      print *, '   ' , trim(debugnames(counter))
   end do


end subroutine read_input_file


SUBROUTINE format_string(base, i, j, k, formatted)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: base
    INTEGER, INTENT(IN)           :: i, j, k
    CHARACTER(LEN=256), INTENT(OUT) :: formatted
    CHARACTER(LEN=256) :: temp  ! Temporary variable to hold the full formatted string

    ! Write to temporary variable
    WRITE(temp, '(A,"_i",I4.4,"_j",I4.4,"_k",I4.4,".txt")') TRIM(base), i, j, k

    ! Trim the output and store it in formatted
    formatted = TRIM(temp)

END SUBROUTINE format_string


! --

end program get_tseries

