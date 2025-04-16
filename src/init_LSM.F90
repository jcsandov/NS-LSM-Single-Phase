subroutine init_LSM ()

  use global
  use global_param
  use global_app
  use global_mpi
  use global_osg
  use checksum

  use global_lsm
  use global_debug

implicit none

integer :: i,j,k,l,nz

!cosas debugeo integral
integer :: lla,llb
integer :: kka,kkb,jja,jjb,iia,iib
real (kind = rdf) :: px,py,pz, valortest
real (kind = rdf),dimension(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)) :: funcionprobar
real (kind = rdf) :: valorintegral 
integer :: tsdbg ! counter variable to loop over the debug tsteps


call init_mg_lsm()

call init_input_fileLSM() !inicializa valor de epslsm de acuerdo a espaciamiento de la malla

!print *, epslsm
!print *, thetainc
!print *, drho, dmu
!print *, FrLSM, WeLSM
!print *, myid
!print *, btype_lsm(:,:), myid
!print *, RNfreq
!print *, numiter_lsm
!print *, delti_lsm
!print *, isnarrow, narrowcoef
!print *, orderAD

call init_phizero_dat() !inicializa phizero
!call init_phistatic()
!call probarromu()


call outputD1_real(phi_zero,'phi_zero')
call outputD1_real(phi_static,'phi_static')
!call outputD3_real(funcionprobar,'dmu')


!------------------------------------------------------------------------------------------

contains

!rutina para debugear integral3d
subroutine crearFuncion()
        implicit none


        funcionprobar =zero



        do k=le_ka(1),le_kb(1)
        do j=le_ja(1),le_jb(1)
        do i=le_ia(1),le_ib(1)
                l = le_idx(i,j,k,1)
                call gfunctionint(x(l),y(l),z(l),funcionprobar(i,j,k))
        end do
        end do
        end do

        

end subroutine

subroutine probarromu()
        implicit none


        funcionprobar =zero


        do k=le_ka(1),le_kb(1)
        do j=le_ja(1),le_jb(1)
        do i=le_ia(1),le_ib(1)
                l = le_idx(i,j,k,1)
                !funcionprobar(i,j,k) = rhoLSM(phi_zero(l))
                funcionprobar(i,j,k) = muLSM(phi_zero(l))
        end do
        end do
        end do

        

end subroutine probarromu

subroutine init_mg_lsm()

  implicit none
  integer :: la,lb
  la = le_idx_a(1)
  lb = le_idx_b(ng)

  !allocate variables
  allocate( phi(la:lb)                , &
            phi_n(la:lb)              , &
            phi_gradient(1:3,la:lb)   , &
            h(la:lb)                  , &            
            hn(la:lb)                 , &            
            h_gradient(1:2,la:lb)     , &
            rsign(la:lb)                &
        )

  !inicializarlas en cero
  !        sgndf = zero ; sgndf_n = zero
  phi          = zero 
  phi_n        = zero
  rsign        = zero
  phi_gradient = zero
  h            = zero 
  hn           = zero 
  h_gradient   = zero
  

  allocate( btype_lsm ( 6 , nzone ) )
         
end subroutine init_mg_lsm

subroutine init_input_fileLSM()

  implicit none

  !Lo pondre como input, de esa manera no me olvidare de modificarlo...
   open(60, file = 'ind3dmg_LSM.dat')
  
      read(60,*) epslsm, numiter_lsm
      read(60,*) thetainc
      read(60,*) drho, dmu,diss_w,diss_a
      read(60,*) FrLSM, WeLSM
      read(60,*) IG,IGitermax, IGtimestep,IGstep2
      read(60,*) RNitermax, RNtimestep, RNfreq
      read(60,*) phi_outputiter
      read(60,*) k_surface
      read(60,*) ((btype_lsm(i,nz), i = 1, 6), nz = 1, nzone)
      read(60,*) bc_extrapolation_order
      read(60,*) isnarrow, narrowcoef
      read(60,*) orderAD
      read(60,*) nobstacles  
      read(60,*) sussman_correction_method
      read(60,*) call_levelsetmethod
      read(60,*) call_reinitialisation
      read(60,*) hybrid_reinitialisation
      read(60,*) sweep_lsqm
      read(60,*) radius_lsqm
      read(60,*) ConvergenceToleranceGeomReini
      read(60,*) TotalVolumeComputation
      read(60,*) GlobalMassCorrection
      read(60,*) n_save_LSdebug_tsteps
      if( n_save_LSdebug_tsteps > 0 ) then
        read(60,*) (save_LSdebug_tsteps(tsdbg), tsdbg = 1, n_save_LSdebug_tsteps)
      end if
      read(60,*) epsReinitialisation
      read(60,*) OrderLSAdvectionBoundaries
      read(60,*) OrderReinitialisationBoundaries
      read(60,*) ENOBCReinitialisation
      read(60,*) BigPhi
      read(60,*) limit_ghost_velocities
      read(60,*) zero_pressure_fs


   close(60)
   
  ! we initialise the current time debugger (the time itereation when a 
  ! debug solution is gonna be saved)

  if(n_save_LSdebug_tsteps > 0) then
    current_LSdebug_counter = 1
    current_LSdebug_tstep = save_LSdebug_tsteps(1)
  end if

  delti_lsm = delti / real( numiter_lsm , kind = rdf )
                
end subroutine init_input_fileLSM



!subroutine init_phizeroBD3D()
!        !Breaking dam 3d
!        implicit none
!
!        real (kind = rdf) :: Lzero,Azero
!        integer :: la, lb
!        integer :: ka,kb,ja,jb,ia,ib
!        
!
!
!        Lzero = 1.0_rdf !nivel referencia
!        Azero = 0.4_rdf !amplitud inicial
!        la = le_idx_a(1)
!        lb = le_idx_b(ng)
!
!        allocate(phi_zero(la:lb))
!
!        ka = li_ka(1)
!        kb = li_kb(1)
!        ja = li_ja(1)
!        jb = li_jb(1)
!        ia = li_ia(1)
!        ib = li_ib(1)
!
!        do k=ka,kb
!        do j=ja,jb
!        do i=ia,ib
!                l = le_idx(i,j,k,1)
!                phi_zero(l) = Azero/cosh(sqrt(3.0_rdf*Azero)*x(l)/2.0_rdf)**2 +  Lzero - z(l) 
!                !phi_zero(l) = Lzero - z(l)
!        end do
!        end do
!        end do
!
!        call mg_exchng3_1d (1, phi_zero(:))
!        
!        !call outputD1_real(phi_zero,'phi_zeroushijima_antes')   
!
!        la = le_idx_a(1)
!        lb = le_idx_b(1)
!        if(IG == 1) then
!                call IG_LSM(le_ia(1),le_ib(1)           ,&
!                            le_ja(1),le_jb(1)           ,&
!                            le_ka(1),le_kb(1)           ,&
!                            igp(1), jgp(1), kgp(1)      ,& 
!                            dc(1), de(1), dz(1)         ,& 
!                            csi(1:3,la:lb)              ,& 
!                            eta(1:3,la:lb)              ,& 
!                            zet(1:3,la:lb)              ,&
!                            aj(la:lb)                   ,&
!                            phi_zero(la:lb)             ,&
!                            epslsm) 
!        end if
!
!        !call outputD1_real(phi_zero,'phi_zeroushijima_despues')  
!
!        phi(la:lb) = phi_zero(la:lb)
!        phi_n(la:lb) = phi_zero(la:lb)
!
!        ! no gp actualizado, solo la parte interior
!
!
!
!
!end subroutine init_phizeroBD3D      

subroutine init_phistatic()
        !Yue 2D laminar Open channel flow
        implicit none

        real (kind = rdf) :: Lzero
        integer :: la, lb
        integer :: ka,kb,ja,jb,ia,ib
        



        Lzero = 1.0_rdf
        la = le_idx_a(1)
        lb = le_idx_b(ng)

        allocate(phi_static(la:lb))

        ka = li_ka(1)
        kb = li_kb(1)
        ja = li_ja(1)
        jb = li_jb(1)
        ia = li_ia(1)
        ib = li_ib(1)

        do k=ka,kb
        do j=ja,jb
        do i=ia,ib
                l = le_idx(i,j,k,1)
                phi_static(l) = Lzero - z(l)
        end do
        end do
        end do

        call mg_exchng3_1d (1, phi_static(:))
        
      
 



end subroutine init_phistatic

   

subroutine init_phizero_dat()
     
	!Lee el phi_inicial del archivo phi_ini
	
     implicit none

    integer ::  tmp_size
    integer :: sbuf_size

    integer, dimension(:,:,:,:,:), allocatable :: sbuf_idx

    real (kind = rdf), dimension (:), allocatable :: sbuf


    real (kind = rdf), dimension(:),   allocatable :: phitmp , phi_n_tmp

    integer :: num_vars, me_start
    integer :: ia, ib
    integer :: ja, jb
    integer :: ka, kb
    integer :: kk, jj, ii

    integer :: nv
    integer :: l1, l2,la,lb
    integer :: nodes_total

    integer :: myunit
    character (len=256) :: filename

    integer :: imax, jmax, kmax

    la = le_idx_a(1)
    lb = le_idx_b(ng)

    allocate(phi_zero(la:lb))
    allocate(phi_static(la:lb))
    
    phi_zero   = zero
    phi_static = zero

    num_vars = 2

    nodes_total = 0
    do nz = 1, nzone
       nodes_total = nodes_total + img(1,nz) * jmg(1,nz) * kmg(1,nz)
    end do

    sbuf_size = num_vars * nodes_total

    imax = maxval(img)
    jmax = maxval(jmg)
    kmax = maxval(kmg)

    allocate (sbuf(sbuf_size), &
              sbuf_idx(num_vars,imax,jmax,kmax,nzone))

    l = 0
    do nz = 1, nzone
    do k = 1, kmg(1,nz)
       do j = 1, jmg(1,nz)
          do i = 1, img(1,nz)
             do nv = 1, num_vars
                l = l + 1; sbuf_idx(nv,i,j,k,nz) = l
             end do
          end do
       end do
    end do
    end do

    if (myid == root) then


       ! get solution from 'phi_ini' file
       !


          open  (unit = 42, file = 'phi_ini', form = 'unformatted')
          do nz = 1, nzone

             ! temporary arrays on root
             ! 
             tmp_size = img(1,nz) * jmg(1,nz) * kmg(1,nz)
      
             allocate(phitmp(tmp_size))
             allocate(phi_n_tmp(tmp_size))             
                     

             read  (unit = 42)  (phitmp(l),l = 1, tmp_size)
             read  (unit = 42)  (phi_n_tmp(l),l = 1, tmp_size)
             
             ! pack send buf
             ! 
             l1 = 0
             if (nz == 1) l2 = 0
             do k = 1, kmg(1,nz)
                do j = 1, jmg(1,nz)
                   do i = 1, img(1,nz)
                         l1 = l1 + 1
                         l2 = l2 + 1; sbuf(l2) =  phitmp(l1)
                         l2 = l2 + 1; sbuf(l2) =  phi_n_tmp(l1)
                   end do
                end do
             end do

             deallocate (phitmp)
             deallocate (phi_n_tmp) 
          end do
          close (unit = 42)

       end if


    ! broadcast to all processes
    ! 
    call mpi_bcast (sbuf, sbuf_size, MPI_REAL_TYPE, root, mpi_comm_world, ierr)

    ! upack subuf including ghost celss into
    ! local copies of phi, phi_static
    !
    ia = le_ia(1)
    ja = le_ja(1)
    ka = le_ka(1)

    ib = le_ib(1)
    jb = le_jb(1)
    kb = le_kb(1)

    ! 3d decomp
    ! 
    if (myup    == mpi_proc_null) kb = li_kb(1)
    if (mydown  == mpi_proc_null) ka = li_ka(1)
    if (myright == mpi_proc_null) jb = li_jb(1)
    if (myleft  == mpi_proc_null) ja = li_ja(1)
    if (myfront == mpi_proc_null) ib = li_ib(1)
    if (myback  == mpi_proc_null) ia = li_ia(1)

    do k = ka, kb
       do j = ja, jb
          do i = ia, ib
             ii = i + gi_ia(1) - 1
             jj = j + gi_ja(1) - 1
             kk = k + gi_ka(1) - 1

             phi_zero(le_idx(i,j,k,1))   = sbuf(sbuf_idx(1,ii,jj,kk,myzone))
             phi_static(le_idx(i,j,k,1)) = sbuf(sbuf_idx(2,ii,jj,kk,myzone))
             
          end do
       end do
    end do

    deallocate(sbuf, sbuf_idx)
    

    call mg_exchng3_1d ( 1 , phi_zero   (:) )
    call mg_exchng3_1d ( 1 , phi_static (:) )
        

    la = le_idx_a(1)
    lb = le_idx_b(1)

    phi(la:lb)   = phi_zero(la:lb)

    ! I use phi_static just as an array for saving the phi_n variable when it's read 
    ! from the file
    phi_n(la:lb) = phi_static(la:lb)

    ! no gp actualizado, solo la parte interior

    ! Initialisation of h
    if ( hydraulic_mode )then

        call calc_h ( le_ia(1),le_ib(1)           , &
                      le_ja(1),le_jb(1)           , &
                      le_ka(1),le_kb(1)           , &
                      igp(1), jgp(1), kgp(1)      , &
                      z(la:lb)                    , &
                      phi(la:lb)                  , &
                      h(la:lb)                      &
                    )
    end if

end subroutine init_phizero_dat


include 'mg_exchng3_1d.F90'


end subroutine init_LSM
