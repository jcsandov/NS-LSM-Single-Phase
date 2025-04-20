
subroutine mg_driver(iteraciontiempo)
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! Multigrid driver (v-cycle only) using solver_rk (runge-kutta)
   ! as a smoother. The logic is intertwined between these two routines.
   ! And although, it might be nice I am not sure separate the logic
   ! so that the routines can be truely modular.
   !
   ! Improved version using solver_daf (daigonal approx. factorization)
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   use global
   use global_param
   use global_mpi
   use global_osg
   use global_app
   !use wf_mpi

   use global_lsm
   use global_debug
  
   !for debug
   use global_obstacle, only : act_obstacle_ad
  
   implicit none

   ! integers
    
   integer :: decide_recalc_rh=0

   ! counters & placeholders
   ! 
   integer :: n
   integer :: l, m
   integer :: i, j, k
   integer :: la
   integer :: lb
   integer :: iteraciontiempo
   logical :: converged
   character(len = 256) :: debugname

   ! defaults for v-cycle multigrid
  
   nstop = ng
   itc = 0
   itmax = itm(1)               ! remove connect to input file
   
   ! store 'n' and 'n-1' levels for new time step 
   
   la = le_idx_a(1)
   lb = le_idx_b(1)

   qnm1(:,la:lb) = qn(:,la:lb)
   qn  (:,la:lb) = q (:,la:lb)

   ! initialize q on coarse grids (mp should work)
   !
   do m = 1, me
      call mg_inject (   qn (m,:) )
      call mg_inject ( qnm1 (m,:) )
   end do

   !-------------------------------------------------------
   ! Level Set

   ! store 'n' and 'n-1' levels for new time step
    
   la = le_idx_a(1)
   lb = le_idx_b(1)

   !phi_n(la:lb) =  phi(la:lb)
   !hn(:) = h(:)

   ! Compute rsign flag-variable
   ! eps_sims: smallest representable number for a given
   !           precision   
   ! rsign = 1 ---> water - phase
   ! rsign = 0 ---> air   - phase
   
   !rsign = zero
   !where ( phi > eps_sims ) rsign = one

   ! Set rsign = 0 within the obstacle if needed
   !call rsign_blanking ( le_ia(1),le_ib(1)           , &
   !                      le_ja(1),le_jb(1)           , &
   !                      le_ka(1),le_kb(1)           , &
   !                      igp(1), jgp(1), kgp(1)      , & 
   !                      rsign(la:lb)                  &
   !                     )

   !call calc_phi_gradient( le_ia(1),le_ib(1)           ,&
   !                        le_ja(1),le_jb(1)           ,&
   !                        le_ka(1),le_kb(1)           ,&
   !                        igp(1), jgp(1), kgp(1)      ,&
   !                        dc(1), de(1), dz(1)         ,&
   !                        csi(1:3,la:lb)              ,&
   !                        eta(1:3,la:lb)              ,&
   !                        zet(1:3,la:lb)              ,&
   !                        phi(la:lb)                  ,&
   !                        phi_gradient(1:3,la:lb)      &
   !                      )

   call calc_phi_gradient_ENO2 ( le_ia(1),le_ib(1)           , &
                                 le_ja(1),le_jb(1)           , &
                                 le_ka(1),le_kb(1)           , &
                                 igp(1), jgp(1), kgp(1)      , &
                                 dc(1), de(1), dz(1)         , &
                                 csi(1:3,la:lb)              , &
                                 eta(1:3,la:lb)              , &
                                 zet(1:3,la:lb)              , &
                                 phi(la:lb)                  , &
                                 phi_gradient(1:3,la:lb)       &
                               )

   if ( hydraulic_mode ) then

      ! Make the pressure, total pressure for GFM using ϕn (hn) 
      ! to apply the post-advection correction (ϕn was the one used
      ! for the NS solver, so it's the one I need to use to apply the
      ! correction)
      !q(1,:) = q(1,:) + ( one / FrLSM**two ) * hn(:) 

      ! Get h(x,y,z) from the ϕ distribution and store it 
      call calc_h ( le_ia(1),le_ib(1)           , &
                    le_ja(1),le_jb(1)           , &
                    le_ka(1),le_kb(1)           , &
                    igp(1), jgp(1), kgp(1)      , &
                    z(la:lb)                    , &
                    phi(la:lb)                  , &
                    h(la:lb)                      &
                   )
      
      ! Get ∂h/∂x, and ∂h/∂y for the momentum equations
      call calc_h_gradient_ENO2 ( le_ia(1),le_ib(1)           , &
                                  le_ja(1),le_jb(1)           , &
                                  le_ka(1),le_kb(1)           , &
                                  igp(1), jgp(1), kgp(1)      , &
                                  dc(1), de(1), dz(1)         , &
                                  csi(1:3,la:lb)              , &
                                  eta(1:3,la:lb)              , &
                                  zet(1:3,la:lb)              , &
                                  h(la:lb)                    , &
                                  h_gradient(1:2,la:lb)         &
                                 )
   end if

   ! initialize q on coarse grids (mp should work)
   !
   call mg_inject ( phi_n (:) )
   call mg_inject ( phi   (:) )
   !call mg_inject ( rsign (:) )

   !gp
   do n=1,ng

      !call mg_exchng3_2d ( n , phi_gradient(1:3,:) ) 
      call mg_exchng3_1d ( n , phi_n(:) )
      call mg_exchng3_1d ( n , phi  (:) )
      !call mg_exchng3_1d ( n , rsign(:) )

   end do
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! 
   ! ########   ######  ######## ##     ## ########   #######          ######## #### ##     ## ########    
   ! ##     ## ##    ## ##       ##     ## ##     ## ##     ##            ##     ##  ###   ### ##          
   ! ##     ## ##       ##       ##     ## ##     ## ##     ##            ##     ##  #### #### ##          
   ! ########   ######  ######   ##     ## ##     ## ##     ## #######    ##     ##  ## ### ## ######      
   ! ##              ## ##       ##     ## ##     ## ##     ##            ##     ##  ##     ## ##          
   ! ##        ##    ## ##       ##     ## ##     ## ##     ##            ##     ##  ##     ## ##          
   ! ##         ######  ########  #######  ########   #######             ##    #### ##     ## ########    
   ! 
   ! #### ######## ######## ########     ###    ######## ####  #######  ##    ##                           
   !  ##     ##    ##       ##     ##   ## ##      ##     ##  ##     ## ###   ##                           
   !  ##     ##    ##       ##     ##  ##   ##     ##     ##  ##     ## ####  ##                           
   !  ##     ##    ######   ########  ##     ##    ##     ##  ##     ## ## ## ##                           
   !  ##     ##    ##       ##   ##   #########    ##     ##  ##     ## ##  ####                           
   !  ##     ##    ##       ##    ##  ##     ##    ##     ##  ##     ## ##   ###                           
   ! ####    ##    ######## ##     ## ##     ##    ##    ####  #######  ##    ##                           
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   inner_iteration : &
   do it = 1, itmax
  
      ! check if stop.now is present, and if so, it stops the simulation
      call stop_file_check()

      converged = .false.

      itc = itc + 1

      ! save fine grid solution for
      ! various residual calculations

      !do l = le_idx_a(1) , le_idx_b(1)
      !   qold_mg(1:4,l) = q(1:4,l)
      !   qold_mg(5,l)   = xnut(l)
      !end do

      qold_mg(1:4,:) = q(1:4,:)
      qold_mg(5,:)   = xnut(:)

      ! multigrid v-cycle
      do n = ns, nstop

         la = le_idx_a(n)
         lb = le_idx_b(n)

         if (n  > 1) then
            ! exchange ghost points for q
            call mg_exchng3_2d (n-1, rh)

            ! mg_inject; mg_calc_pk
            
            !call mg_nlevel (n)

            ! update ghost points for uij
            
            if (turbulence .and. nlinc) then
               call mg_exchng3_2d(n, uij)
            end if
           
            call mg_zero_pk_boundary (n, pk)
!           call mg_cksum ('nl-pk',n, pk(:,:))

         end if

         do itr = 1, iter(n)

            if (n /= ng) then
               if (itr /= iter(n)) then
                  decide_recalc_rh = 0
               else
                  decide_recalc_rh = 1
               end if
            end if

            ! exchange ghost points for q
            ! call mg_exchng3_2d (n, q)

            !if (turbulence) then
            !   call mg_exchng3_1d (n, xnut)
            !end if


            if ( des .and. n == ns ) then

               !call des_eddy (le_ia(n),le_ib(n),  & ! 1
               !  le_ja(n),le_jb(n),            & ! 1
               !  le_ka(n),le_kb(n),            & ! 1
               !  igp(n), jgp(n), kgp(n),       & ! 2
               !  dc(n), de(n), dz(n),          & ! 3
               !  dtev(la:lb),           & ! 4
               !  q(1:me,la:lb),         & ! 5
               !  qn(me,la:lb),          & ! 6
               !  qnm1(me,la:lb),        & ! 7
               !  csi(1:3,la:lb),        & ! 8
               !  eta(1:3,la:lb),        & ! 9
               !  zet(1:3,la:lb),        & ! 10
               !  aj(la:lb),             & ! 11
               !  wd(la:lb),             & ! 12
               !  xnut(la:lb),           & ! 13
               !  rsign(la:lb))            
                
               if( nzone > 1 ) call mg_bintp_tur

            end if  ! end of des_eddy

                                                                                                                                                           
            ! -------------------------------------------------------------------------------------------
            !   ####  #    #  ####   ####  #####         ###### #      #    # # #####      
            !  #    # #    # #    # #        #           #      #      #    # # #    #     
            !  #      ###### #    #  ####    #           #####  #      #    # # #    #     
            !  #  ### #    # #    #      #   #           #      #      #    # # #    #     
            !  #    # #    # #    # #    #   #           #      #      #    # # #    #     
            !   ####  #    #  ####   ####    #           #      ######  ####  # #####      
            !                                                                                                                                                           
            !
            !  ###### #    # ##### #####    ##   #####   ####  #        ##   ##### #  ####  #    #
            !  #       #  #    #   #    #  #  #  #    # #    # #       #  #    #   # #    # ##   #
            !  #####    ##     #   #    # #    # #    # #    # #      #    #   #   # #    # # #  #
            !  #        ##     #   #####  ###### #####  #    # #      ######   #   # #    # #  # #
            !  #       #  #    #   #   #  #    # #      #    # #      #    #   #   # #    # #   ##
            !  ###### #    #   #   #    # #    # #       ####  ###### #    #   #   #  ####  #    #
            ! --------------------------------------------------------------------------------------------


            ! ghost fluid method to extrapolate velocites and pressure to the air nodes
            ! next to the free-surface (based on the method described in Watanabe, Saruwatari 
            ! & Ingram, JCP, 2008)

            if ( call_solver_daf ) then
            !if ( mod( itc+1 , 2 ) == 0 .and. call_solver_daf ) then

                ! I set all the flow variables to zero within the air-phase               
                !q(1,:) = q(1,:) * rsign
                !q(2,:) = q(2,:) * rsign
                !q(3,:) = q(3,:) * rsign
                !q(4,:) = q(4,:) * rsign
                !q(5,:) = q(5,:) * rsign
               
               do l = la,lb
                  q(:,l) = rsign ( phi(l) ) * q(:,l)
               end do

               ! We extrapolate the velocity and the pressure field using ghost fluid method
               call ghost_fluid_extrapolation( le_ia(n),le_ib(n)           , &
                                               le_ja(n),le_jb(n)           , &
                                               le_ka(n),le_kb(n)           , &
                                               igp(n), jgp(n), kgp(n)      , & 
                                               dc(n), de(n), dz(n)         , & 
                                               q(1:4,la:lb)                , &
                                               xnut(la:lb)                 , &
                                               csi(1:3,la:lb)              , & 
                                               eta(1:3,la:lb)              , & 
                                               zet(1:3,la:lb)              , &
                                               aj(la:lb)                   , &
                                               phi(la:lb)                  , &
                                               phi_n(la:lb)                , &
                                               phi_gradient(1:3, la:lb)    , &
                                               x(la:lb),y(la:lb),z(la:lb)    &
                                             )

                  
               ! write(debugname, fmt ='(a,i6.6)') 'p_postgfm',itc
               ! call outputD1_real(q(1,:),debugname) 

               ! I clear up rsign and return it to 1 and 0 definition (in ghost_fluid_extrapolation
               ! I set the first layer of pressure extrapolated nodes to -1 for some operations )
               !rsign = zero
               !where ( phi > eps_sims ) rsign = one

               
            end if ! if ( call_solver_daf )

            ! processor ghost-points update
            call mg_exchng3_2d (n, q)

            ! Apply boundary conditions
            call bcond_fm ( le_ia(n)     , le_ib(n)          , &
                            le_ja(n)     , le_jb(n)          , &
                            le_ka(n)     , le_kb(n)          , &
                            igp(n)       , jgp(n)   , kgp(n) , &
                            dc(n)        , de(n)    , dz(n)  , &
                            q(1:4,la:lb)                     , &
                            eta(1:3,la:lb)                   , &
                            aj(la:lb)                        , &
                            x(la:lb),y(la:lb),z(la:lb)       , &
                            xnut(la:lb)                      , &
                            phi(la:lb)                         &
                           )


            ! LES dynamic Smagorinksy
            call les_dynamic_smagorinsky ( le_ia(n),le_ib(n)          , &
                                           le_ja(n),le_jb(n)          , &
                                           le_ka(n),le_kb(n)          , &
                                           igp(n), jgp(n), kgp(n)     , &
                                           dc(n), de(n), dz(n)        , &
                                           x(la:lb),y(la:lb),z(la:lb) , &
                                           q(1:5,la:lb)               , &
                                           csi(1:3,la:lb)             , &
                                           eta(1:3,la:lb)             , &
                                           zet(1:3,la:lb)             , &
                                           aj(la:lb)                  , &
                                           xnut(la:lb)                  & 
                                         ) 

            !xnut = zero

            ! ---------------------------------------------------------------------------------------
            !                                                                                                                                  
            ! #    #   ##   #    # # ###### #####         ####  #####  ####  #    # ######  ####   
            ! ##   #  #  #  #    # # #      #    #       #        #   #    # #   #  #      #       
            ! # #  # #    # #    # # #####  #    # #####  ####    #   #    # ####   #####   ####   
            ! #  # # ###### #    # # #      #####             #   #   #    # #  #   #           #  
            ! #   ## #    #  #  #  # #      #   #        #    #   #   #    # #   #  #      #    #  
            ! #    # #    #   ##   # ###### #    #        ####    #    ####  #    # ######  ####   
            !
            !
            !   ####   ####  #      #    # ###### #####  
            !  #      #    # #      #    # #      #    # 
            !   ####  #    # #      #    # #####  #    # 
            !       # #    # #      #    # #      #####  
            !  #    # #    # #       #  #  #      #   #  
            !   ####   ####  ######   ##   ###### #    #                                                                                                                                            
            !
            ! ---------------------------------------------------------------------------------------

            if ( call_solver_daf ) then
            
               ! we assume that velocity field is updated near the free-surface using
               ! ghost fluid method
            
               call solver_daf ( me, n, decide_recalc_rh, itr                          , & 
                                 le_ia(n),le_ib(n),le_ja(n),le_jb(n),le_ka(n),le_kb(n) , & 
                                 igp(n), jgp(n), kgp(n)                                , & 
                                 dc(n), de(n), dz(n)                                   , & 
                                 q(1:4,la:lb)                                          , & 
                                 qn(1:4,la:lb)                                         , & 
                                 qnm1(1:4,la:lb)                                       , & 
                                 csi(1:3,la:lb)                                        , & 
                                 eta(1:3,la:lb)                                        , & 
                                 zet(1:3,la:lb)                                        , & 
                                 aj(la:lb)                                             , & 
                                 xnut(la:lb)                                           , & 
                                 rh(1:4,la:lb)                                         , & 
                                 phi(la:lb)                                            , &
                                 phi_gradient(1:3, la:lb)                              , &
                                 h_gradient(1:2, la:lb)                                , &
                                 x(la:lb),y(la:lb),z(la:lb)                            , &
                                 iteraciontiempo                                         &
                               )   
            
            end if ! if ( call_solver_daf )

         end do ! itr

         !call mg_zero_pk_boundary (n, pk)

      end do ! do n = ns, nstop

      ! prolong from grid n+1 to grid n

      do n = nstop-1, ns, -1

         ! prolong from coarse grid to fine grid
         ! 
         call mg_prolong (n)

      end do


      ! finished V-cycle; apply boundary conditions
      ! 
      n  = 1                      
      la = le_idx_a(1)
      lb = le_idx_b(1)
     
      !call bcond_fm ( le_ia(n)     , le_ib(n)          , &
      !                le_ja(n)     , le_jb(n)          , &
      !                le_ka(n)     , le_kb(n)          , &
      !                igp(n)       , jgp(n)   , kgp(n) , &
      !                q(1:4,la:lb)                     , &
      !                rsign(la:lb)                       &
      !               )

      call bcond_fm ( le_ia(n)     , le_ib(n)          , &
                      le_ja(n)     , le_jb(n)          , &
                      le_ka(n)     , le_kb(n)          , &
                      igp(n)       , jgp(n)   , kgp(n) , &
                      dc(n)        , de(n)    , dz(n)  , &
                      q(1:4,la:lb)                     , &
                      eta(1:3,la:lb)                   , &
                      aj(la:lb)                        , &
                      x(la:lb),y(la:lb),z(la:lb)       , &
                      xnut(la:lb)                      , &
                      phi(la:lb)                         &
                     )

      call mg_exchng3_2d (n, q) 

      !call conservation(le_ia(n),le_ib(n),&
      !                  le_ja(n),le_jb(n),&
      !                  le_ka(n),le_kb(n),&
      !                  igp(n), jgp(n), kgp(n), &
      !                  xc(la:lb), yc(la:lb), zc(la:lb),        &
      !                  xe(la:lb), ye(la:lb), ze(la:lb),        &
      !                  xz(la:lb), yz(la:lb), zz(la:lb),        &
      !                  csi(1:3,la:lb),                         & ! 8
      !                  eta(1:3,la:lb),                         & ! 9
      !                  zet(1:3,la:lb),                         & ! 10
      !                  aj(la:lb) ,                             & ! 11
      !                  phi_n(la:lb) ,                          &
      !                  rsign(la:lb) ,                          &
      !                  q(1:4,la:lb) )


      if( nzone > 1 ) call mg_bintp_mom

      ! check for convergence during this time step
      ! 
      call mg_ck_conver ( converged )

      if ( converged ) then

         ! If converged == true, I call gfm a last time to go into the level-set
         ! method with latest velocity field extrapolated to the ghost nodes

         ! q(1,:) = q(1,:) * rsign
         ! q(2,:) = q(2,:) * rsign
         ! q(3,:) = q(3,:) * rsign
         ! q(4,:) = q(4,:) * rsign
         ! q(5,:) = q(5,:) * rsign

         do l = la,lb
            q(:,l) = rsign ( phi(l) ) * q(:,l)
         end do
                                             
         ! We extrapolate the velocity and the pressure field using ghost fluid method
         call ghost_fluid_extrapolation( le_ia(n),le_ib(n)           , &
                                         le_ja(n),le_jb(n)           , &
                                         le_ka(n),le_kb(n)           , &
                                         igp(n), jgp(n), kgp(n)      , & 
                                         dc(n), de(n), dz(n)         , & 
                                         q(1:4,la:lb)                , &
                                         xnut(la:lb)                 , &
                                         csi(1:3,la:lb)              , & 
                                         eta(1:3,la:lb)              , & 
                                         zet(1:3,la:lb)              , &
                                         aj(la:lb)                   , &
                                         phi(la:lb)                  , &
                                         phi_n(la:lb)                , &
                                         phi_gradient(1:3, la:lb)    , &
                                         x(la:lb),y(la:lb),z(la:lb)    &
                                       )

         ! Make sure rsign is properly defined before going to the level-set method
         !rsign = zero
         !where ( phi > eps_sims ) rsign = one

         ! I save xnut here to compute conver_real
         q(5,:) = xnut(:)

         exit inner_iteration

      end if


   end do inner_iteration


! ======================================================================================


   n = 1
   la = le_idx_a(1)
   lb = le_idx_b(1)
        
   ! -----------------------------------------------------------------  
   !                                                                                                          
   ! #      ###### #    # ###### #             ####  ###### #####    
   ! #      #      #    # #      #            #      #        #      
   ! #      #####  #    # #####  #      #####  ####  #####    #      
   ! #      #      #    # #      #                 # #        #      
   ! #      #       #  #  #      #            #    # #        #      
   ! ###### ######   ##   ###### ######        ####  ######   #      
   !
   ! #    # ###### ##### #    #  ####  #####  
   ! ##  ## #        #   #    # #    # #    # 
   ! # ## # #####    #   ###### #    # #    # 
   ! #    # #        #   #    # #    # #    # 
   ! #    # #        #   #    # #    # #    # 
   ! #    # ######   #   #    #  ####  #####  
   !
   ! -----------------------------------------------------------------  
                                                                                                          
   if( call_levelsetmethod .or. call_reinitialisation ) then

      phi_n = phi

      call levelsetmethod( le_ia(n) , le_ib(n)                    , &
                           le_ja(n) , le_jb(n)                    , &
                           le_ka(n) , le_kb(n)                    , &
                           igp(n), jgp(n), kgp(n)                 , & 
                           dc(n), de(n), dz(n)                    , & 
                           q(1:4,la:lb)                           , & 
                           csi(1:3,la:lb)                         , & 
                           eta(1:3,la:lb)                         , & 
                           zet(1:3,la:lb)                         , &
                           aj(la:lb)                              , &
                           phi(la:lb), phi_n(la:lb)               , &
                           x(la:lb),y(la:lb),z(la:lb)             , &
                           iteraciontiempo                          &
                         )
    
   end if

   if( MOD( iteraciontiempo , phi_outputiter ) == 0 ) then

      write(debugname, fmt ='(a,i6.6)') 'phi',iteraciontiempo
      call outputD1_real(phi,debugname) 

      ! saving h   
      write(debugname, fmt ='(a,i6.6)') 'h',iteraciontiempo
      call outputD1_real( h , debugname ) 

   end if

    ! check if we have to save debugging solutions

   if( iteraciontiempo == current_LSdebug_tstep ) then

      ! saving the debug solution   
      write(debugname, fmt ='(a,i6.6)') 'phi',iteraciontiempo
      call outputD1_real(phi,debugname) 
   
      ! we then update the value of current_LSdebug_tstep
      current_LSdebug_counter = current_LSdebug_counter + 1

      ! check if I already saved all the debug solutions
      if (current_LSdebug_counter <= n_save_LSdebug_tsteps) then
         current_LSdebug_tstep   = save_LSdebug_tsteps(current_LSdebug_counter)
      end if
   
   end if

contains

  include 'mg_exchng3_2d.F90'
  include 'mg_exchng3_1d.F90'
  include 'mg_cksum.F90'
  include 'mg_prolong.F90'
  include 'mg_residual.F90'
  include 'mg_bintp_tur.F90'
  include 'mg_bintp_mom.F90'

   subroutine stop_file_check()

      ! check if the stop file was created to stop the simulation
      ! Only rank 0 checks for the presence of the file
      if (myid == root) then
         inquire(file='stop.now', exist = stop_now )
      end if

      ! Broadcast the result to all processes
      call MPI_Bcast(stop_now, 1, MPI_LOGICAL, root , MPI_COMM_WORLD, ierr)

      ! If file 'stop.now' exists, terminate the simulation
      if ( stop_now ) then

         if ( myid == root ) print *, 'File stop.now found. Terminating simulation.'

         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

      end if

   end subroutine stop_file_check
   
  subroutine mg_ck_conver (converged)

    implicit none 

    !real (kind = rdf) :: erp
    real (kind = rdf) :: eru
    real (kind = rdf) :: erv
    real (kind = rdf) :: erw
    !real (kind = rdf) :: erx

    real (kind = rdf), save :: ermax
    real (kind = rdf), save :: eomax

    ! residual norm calcuations
    !
    real (kind = rdf), dimension(1:me), save :: init_resd
    real (kind = rdf), dimension(1:me) :: residual

    logical :: converged

    ! if (itc == 1 .or. ((itc / icn) * icn == itc)) then

    ! residual calculations
    !
    call mg_residual (1, residual)

    if (myid == root) then
       open  (unit = 8, file = 'conver', position = 'append', form = 'formatted')
       write (unit = 8, fmt = '(i7.7,1x,6(g13.6,1x))') itc, residual(:)
       close (unit = 8)
       !end if

       ! used mpi_reduce_scatter () in mg_resid; thus
       ! compute logical converged on all processes

       ! save initial residual
       !
       if (itc == 1) then

          init_resd(:) = residual(:)

!           init_resd(1) = erp
!           init_resd(2) = eru
!           init_resd(3) = erv
!           init_resd(4) = erw
!           init_resd(5) = erx
       end if

       ! check converged & go to next time step?
       !
       eomax = max(residual(2),residual(3),residual(4))

       eru = residual(2) - init_resd(2)
       erv = residual(3) - init_resd(3)
       erw = residual(4) - init_resd(4)

       ermax = max(eru, erv, erw)

       ! if (duct) ermax = max(eru, erv)
       !$ermax = erv

       converged = .false.
       if ((itc > it_min) .and. (ermax < er_min)) converged = .true.
       if ((itc > it_min) .and. (eomax < eo_min)) converged = .true.

    end if

    call mpi_bcast (converged, 1, mpi_logical, root, mpi_comm_world, ierr)
    !call mpi_bcast (converged, 1, mpi_logical, root, comm3d, ierr)

  end subroutine mg_ck_conver
  
  subroutine mg_zero_pk_boundary(n, var)

    integer , intent(in) :: n
    real (kind = rdf), dimension (:,:), intent(inout) :: var

    if (myback == mpi_proc_null)  then
       do k = le_ka(n), le_kb(n)
       do j = le_ja(n), le_jb(n)
          l = le_idx(li_ia(n),j,k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myfront == mpi_proc_null) then
       do k = le_ka(n), le_kb(n)
       do j = le_ja(n), le_jb(n)
          l = le_idx(li_ib(n),j,k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myleft == mpi_proc_null)  then
       do k = le_ka(n), le_kb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,li_ja(n),k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myright == mpi_proc_null) then
       do k = le_ka(n), le_kb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,li_jb(n),k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (mydown == mpi_proc_null)  then
       do j = le_ja(n), le_jb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,j,li_ka(n),n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myup == mpi_proc_null) then
       do j = le_ja(n), le_jb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,j,li_kb(n),n)
          var(:,l) = zero
       end do
       end do
    end if
    
  end subroutine mg_zero_pk_boundary


!  subroutine calc_rsign ( il, iu        , &
!                          jl, ju        , &
!                          kl, ku        , &
!                          igp, jgp, kgp , &
!                          phi, rsign      &
!                        )
!
!     integer, intent(in) :: il, iu, jl, ju, kl, ku
!     integer, intent(in) :: igp , jgp , kgp
!
!     real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in)    :: phi     
!     real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(inout) :: rsign 
!
!     integer :: ista, iend, jsta, jend, ksta, kend
!     integer :: i,j,k
!
!     ista = il ; jsta = jl ; ksta = kl 
!     iend = iu ; jend = ju ; kend = ku 
!
!     if (myback  == mpi_proc_null)  ista = il + igp 
!     if (myleft  == mpi_proc_null)  jsta = jl + jgp 
!     if (mydown  == mpi_proc_null)  ksta = kl + kgp 
!
!     if (myfront == mpi_proc_null)  iend = iu - igp
!     if (myright == mpi_proc_null)  jend = ju - jgp
!     if (myup    == mpi_proc_null)  kend = ku - kgp
!
!
!     do i = ista, iend
!     do j = jsta, jend
!     do k = ksta, kend
!
!         rsign(i,j,k) = ( sign( one , phi(i,j,k) - eps_sims ) + one ) / two
!
!     end do
!     end do
!     end do
!
!  end subroutine calc_rsign

!  subroutine rsign_blanking ( il, iu        , &
!                              jl, ju        , &
!                              kl, ku        , &
!                              igp, jgp, kgp , &
!                              rsign           &
!                            )
!
!
!
!     ! This routine set rsign to zero within the blanking regions
!     integer, intent(in) :: il, iu, jl, ju, kl, ku
!     integer, intent(in) :: igp , jgp , kgp
!
!     real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(inout) :: rsign 
!
!     integer :: ista, iend, jsta, jend, ksta, kend
!     integer :: i,j,k
!
!     ista = il ; jsta = jl ; ksta = kl 
!     iend = iu ; jend = ju ; kend = ku 
!
!     if (myback  == mpi_proc_null)  ista = il + igp 
!     if (myleft  == mpi_proc_null)  jsta = jl + jgp 
!     if (mydown  == mpi_proc_null)  ksta = kl + kgp 
!
!     if (myfront == mpi_proc_null)  iend = iu - igp
!     if (myright == mpi_proc_null)  jend = ju - jgp
!     if (myup    == mpi_proc_null)  kend = ku - kgp
!
!     if ( nblke == 0 ) return
!
!     do i = ista, iend
!     do j = jsta, jend
!     do k = ksta, kend
!
!         do nb = 1,nblke
!            
!            ! We allow rsign to be 1 on the blanking region, because
!            ! phi and velocity/pressure values on the obstacles are
!            ! valid entries
!            if ( i > le_blk_ia(1,nb) .and. i < le_blk_ib(1,nb) .and. & 
!                 j > le_blk_ja(1,nb) .and. j < le_blk_jb(1,nb) .and. &
!                 k > le_blk_ka(1,nb) .and. k < le_blk_kb(1,nb) ) then
!
!                 rsign(i,j,k) = zero
!      
!            end if
!
!         end do  
!
!     end do
!     end do
!     end do
!
!  end subroutine rsign_blanking

end subroutine mg_driver

!
! save for later use
!
!



!!$           if (rkm) then
!!$           call solver_rk (me, n, decide_recalc_rh, itr, & 
!!$                img(n), jmg(n), kmg(n), &
!!$                dc, de, dz,             &
!!$                dtau(ls(n):le(n)),      &
!!$                dtev(ls(n):le(n)),      &
!!$                q(1:me,ls(n):le(n)),     &
!!$                qn(1:me,ls(n):le(n)),    &
!!$                qnm1(1:me,ls(n):le(n)),  &
!!$                csi(1:3,ls(n):le(n)),   &
!!$                eta(1:3,ls(n):le(n)),   &
!!$                zet(1:3,ls(n):le(n)),   &
!!$                aj(ls(n):le(n)),        &
!!$                wd(ls(n):le(n)),        &
!!$                xnut(ls(n):le(n)),      &
!!$                pk(1:4,ls(n):le(n)),    &
!!$                rh(1:4,ls(n):le(n)) )
!!$           else if (daf) then
!            call solver_daf (me, n, decide_recalc_rh, itr, & 
!                 img(n), jmg(n), kmg(n), &
!                 dc, de, dz,             &
!                 dtau(ls(n):le(n)),      &
!                 dtev(ls(n):le(n)),      &
!                 q(1:me,ls(n):le(n)),     &
!                 qn(1:me,ls(n):le(n)),    &
!                 qnm1(1:me,ls(n):le(n)),  &
!                 csi(1:3,ls(n):le(n)),   &
!                 eta(1:3,ls(n):le(n)),   &
!                 zet(1:3,ls(n):le(n)),   &
!                 mai(1:4,1:4,ls(n):le(n)),  &
!                 n1i(1:4,1:4,ls(n):le(n)),  &
!                 n2i(1:4,1:4,ls(n):le(n)),  &
!                 mc(1:4,1:4,ls(n):le(n)),   &
!                 spr(1:3,ls(n):le(n)),   &
!                 aj(ls(n):le(n)),        &
!                 wd(ls(n):le(n)),        &
!                 xnut(ls(n):le(n)),      &
!                 pk(1:4,ls(n):le(n)),    &
!                 rh(1:4,ls(n):le(n)) )
! !!$           end if















