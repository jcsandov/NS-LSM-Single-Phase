subroutine solver_daf ( me, decide_grid_level, decide_recalc_rh , & ! 1
                        decide_calc_pk, il, iu, jl, ju, kl, ku  , & ! 2
                        igp, jgp, kgp                           , & ! 3
                        dc, de, dz                              , & ! 4
                        q                                       , & ! 5
                        qn                                      , & ! 6
                        qnm1                                    , & ! 7
                        csi                                     , & ! 8
                        eta                                     , & ! 9
                        zet                                     , & ! 10
                        aj                                      , & ! 11
                        xnut                                    , & ! 12
                        rh                                      , & ! 14
                        phi                                     , &
                        phi_gradient                            , &
                        h_gradient                              , &
                        x, y ,z                                 , &
                        iteraciontiempo                           &
                      )                       
!     uij )                             ! 31


   ! Solve the system of equations using approximate factorization
   ! implicit solver with local time stepping and residual smoothing

   ! input
   !     decide_grid_level
   !     decide_recalc_rh
   !     q(4,ijk)
   !     qn(4,ijk)
   !     qnm1(4,ijk)
   !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
   !     aj(ijk)
   !     xnut(ijk)
   !     dc, de, dz
   !     dc2, de2, dz2
   !     dcsq, desq, dzsq

   ! output
   !     xnut(ijk)
   !     rh(4,ijk)
   !     
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
   use global_param
   use global_app
   use global_mpi
   use checksum

   use global_debug
   use global_lsm, only : rhoLSM , muLSM , thetainc , FrLSM        , &
                          WeLSM , drho , deltaf , epslsm , dissLSM , &
                          heaviside , BigPhi

   use global_obstacle!, only : csi_obs, eta_obs 

   use pascal_tdma

   !use SolitaryWave

   implicit none

   integer :: me

   integer :: il                 ! i lower bound
   integer :: jl                 ! j lower bound
   integer :: kl                 ! k lower bound
   integer :: iu                 ! i upper bound
   integer :: ju                 ! j upper bound
   integer :: ku                 ! k upper bound

   integer :: i_mysta            ! first interior i node
   integer :: j_mysta            ! first interior j node
   integer :: k_mysta            ! first interior k node
   integer :: i_myend            ! last interior i node
   integer :: j_myend            ! last interior j node
   integer :: k_myend            ! last interior k node

   integer :: igp                ! i-direction ghost points
   integer :: jgp                ! j-direction ghost points
   integer :: kgp                ! k-direction ghost points

   real (kind = rdf) :: dc       ! csi-direction grid spacing
   real (kind = rdf) :: de       ! eta-direction grid spacing
   real (kind = rdf) :: dz       ! zeta-direction grid spacing

   ! arrays for rhs & solver routines
   !
   real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: q
   real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(in)    :: qn, qnm1
   real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(out) :: rh
   !real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: pk
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi,eta,zet
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)    , intent(in) :: aj

   !Level set method
   !real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi_n, phi_static
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi !, phi_static
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: phi_gradient
   !real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)    , intent(in) :: h
   real (kind = rdf), dimension(1:2,il:iu,jl:ju,kl:ku), intent(in) :: h_gradient
   !real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: rsign
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: x,y,z
   !real (kind = rdf), dimension(:,:,:,:), allocatable :: gforce,sforce, pinterfaz        

   !Obstaculo
   !real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: wd
   !real (kind = rdf), dimension(:,:,:), allocatable :: identificador 


   ! turbulence stuff
   ! 
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: xnut

   ! Reynolds stresses
   !
   !real (kind = rdf), dimension(1:6,il:iu,jl:ju,kl:ku) :: uij
   !real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) :: rstr

   ! dummy variables
   real (kind = rdf), dimension(:,:,:,:), allocatable   :: qold_af , &
                                                           visc    , &
                                                           ucn_j   , &
                                                           spr     , &
                                                           fv

   real (kind = rdf), dimension(:,:,:,:,:), allocatable ::  mai , &
                                                            n1i , &
                                                            n2i , &
                                                            mc

   real (kind = rdf), dimension(:,:,:), allocatable     ::  dtau , diss
  

   !for debugging
   character (len = 256) :: debugname

   ! control from previous program
   ! 
   
   integer :: decide_grid_level ! actually just contains current grid level
   integer :: decide_calc_pk    ! calculate fine-grid forcing function
   integer :: decide_recalc_rh

   integer :: i, j, k, idend, m 
   integer :: ilsearch, jlsearch, klsearch, iusearch, jusearch, kusearch
   integer :: istage

   ! start using this again;
   ! decide prefix is just confusing
   ! 
   integer :: n
  
   integer :: inout, ibs, ibe, ibse

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

   integer, intent(in) :: iteraciontiempo
   logical :: LinearFlux, ViscousFluxes, NonLinearFluxes, GravitySource , &
              UnsteadyTerm  , NearFreeSurfaceCorrection, DissipationFlux

   DissipationFlux    = .true.

   LinearFlux         = .true.

   ViscousFluxes      = .true.
   NonLinearFluxes    = .true.
   UnsteadyTerm       = .true.

   GravitySource      = .true.

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


   ! ▶ TO DO: check all the subroutines here to apply single-phase
   !        method

   n = decide_grid_level

   ! allocate variables: dummy, flux variables
   allocate ( qold_af   (1:4 , il:iu , jl:ju , kl:ku ), &
              diss      (      il:iu , jl:ju , kl:ku ), &
              visc      (1:3 , il:iu , jl:ju , kl:ku ), &
              ucn_j     (1:3 , il:iu , jl:ju , kl:ku ), &
              fv        (1:3 , il:iu , jl:ju , kl:ku ), &
              dtau      (      il:iu , jl:ju , kl:ku )  &
            )

   ! solver variables initialisation
  
   qold_af   = zero
   visc      = zero 
   diss      = zero  
   ucn_j     = zero
   dtau      = zero
   fv        = zero
   
   rh      = zero

   ! daf, arrays for model matrices
   allocate ( mai ( 1:4,1:4 , il:iu , jl:ju , kl:ku ), &
              n1i ( 1:4,1:4 , il:iu , jl:ju , kl:ku ), &
              n2i ( 1:4,1:4 , il:iu , jl:ju , kl:ku ), &
              mc  ( 1:4,1:4 , il:iu , jl:ju , kl:ku ), &
              spr ( 1:3     , il:iu , jl:ju , kl:ku )  &
            )

   mai = zero ; n1i = zero ; n2i = zero ; mc = zero ; spr = zero


   ! set first and last ** interior ** grid nodes for this process
   ! 
   ! *** note ***
   !
   ! using the 'my' notation because this varies depending on
   ! whether or not the process owns any nodes on the boundary
   ! of the computational domain
   ! 
  
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! processes on the domain boundaries
   ! 
   if (myback == mpi_proc_null)  i_mysta = il + igp + 1
   if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
   if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

   if (myfront == mpi_proc_null) i_myend = iu - igp - 1
   if (myright == mpi_proc_null) j_myend = ju - jgp - 1
   if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

   ! include boundary grid plane for characteristics-based boundary conditions
   !
   idend = i_myend

   if ( myfront == mpi_proc_null .and. n == 1 .and. btype(2,myzone) == 5 ) then
      idend = i_myend + 1
   end if

   ! save solution for iterations; include ghost points because we
   ! don't know if the boundary is included (and needed or not)
   ! because of the characteristics-based boundary conditions or the
   ! fact that we might be on an interior process
   !    

   qold_af = q

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   ! I'm gonna extrapolate the pressure values to the air phase, no matter the zero
   ! pressure condition at the free-surface
   ! call rhs_ghost_fluid_nodes_extp_3d( q(1,:,:,:)  )
   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !
   !   ####   ####  #    # ##### #####    ##   #    #   ##   #####  #   ##   #    # ##### 
   !  #    # #    # ##   #   #   #    #  #  #  #    #  #  #  #    # #  #  #  ##   #   #   
   !  #      #    # # #  #   #   #    # #    # #    # #    # #    # # #    # # #  #   #   
   !  #      #    # #  # #   #   #####  ###### #    # ###### #####  # ###### #  # #   #   
   !  #    # #    # #   ##   #   #   #  #    #  #  #  #    # #   #  # #    # #   ##   #   
   !   ####   ####  #    #   #   #    # #    #   ##   #    # #    # # #    # #    #   #   
   !                                                                                      
   !                                                                                      
   !  #    # ###### #       ####   ####  # ##### #   #                                    
   !  #    # #      #      #    # #    # #   #    # #                                     
   !  #    # #####  #      #    # #      #   #     #                                      
   !  #    # #      #      #    # #      #   #     #                                      
   !   #  #  #      #      #    # #    # #   #     #                                      
   !    ##   ###### ######  ####   ####  #   #     #                                      
   !                                                                                      
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  It computes the ξ, η and ζ contravariant velocities divided by the jacobianc U^j/J.
   !  
   !     U^j      1      /         ∂ξ^j  \
   !    ----- =  --- *  (   u_i * ------  )   
   !      J       J      \         ∂x_i  /
   !                
   !  it returns the array ucn_j with the three components of U^j/J on every node

   call rhs_contra_j ()
  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      

!#ifdef DEBUG
!   call cksum_4d_par ('daf.ucn', ucn_j)
!#endif

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                             
   !   #       ####   ####    ##   #         ##### # #    # ###### 
   !   #      #    # #    #  #  #  #           #   # ##  ## #      
   !   #      #    # #      #    # #           #   # # ## # #####  
   !   #      #    # #      ###### #           #   # #    # #      
   !   #      #    # #    # #    # #           #   # #    # #      
   !   ######  ####   ####  #    # ######      #   # #    # ###### 
   !                                                               
   !                                                               
   !    ####  ##### ###### #####                                   
   !   #        #   #      #    #                                  
   !    ####    #   #####  #    #                                  
   !        #   #   #      #####                                   
   !   #    #   #   #      #                                       
   !    ####    #   ###### #                                       
   !                                                                                      
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  it computes the local time step for the pseudo-time time-marching solving algorithm
   !  
   !              /           CFL                             V N                   \
   !   Δτ =  min (   --------------------- , -------------------------------------   )   
   !              \    max(ρ^1, ρ^2, ρ^3)      (1/Re + νt) * max(g^11, g^22, g^33)  /
   !                
   !  Where ρ^j is the spectral radii of the matrix A. ρ^j is computed as ρ^j = sqrt(g^jj)
   !
   !  it returns the array dtau with local pseudo - time on every node

   call rhs_daf_dtau()

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      

!#ifdef DEBUG
!   call cksum_3d_par ('daf.dtau', dtau)
!#endif
   
   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
   !    
   ! ########  ##     ##  ######                                                                                 
   ! ##     ## ##     ## ##    ##                                                                                
   ! ##     ## ##     ## ##                                                                                      
   ! ########  #########  ######                                                                                 
   ! ##   ##   ##     ##       ##                                                                                
   ! ##    ##  ##     ## ##    ##                                                                                
   ! ##     ## ##     ##  ######                                                                                                                                                                                 
   !                                                                                                             
   !                                                                                                             
   !  ######   #######  ##    ##  ######  ######## ########  ##     ##  ######  ######## ####  #######  ##    ## 
   ! ##    ## ##     ## ###   ## ##    ##    ##    ##     ## ##     ## ##    ##    ##     ##  ##     ## ###   ## 
   ! ##       ##     ## ####  ## ##          ##    ##     ## ##     ## ##          ##     ##  ##     ## ####  ## 
   ! ##       ##     ## ## ## ##  ######     ##    ########  ##     ## ##          ##     ##  ##     ## ## ## ## 
   ! ##       ##     ## ##  ####       ##    ##    ##   ##   ##     ## ##          ##     ##  ##     ## ##  #### 
   ! ##    ## ##     ## ##   ### ##    ##    ##    ##    ##  ##     ## ##    ##    ##     ##  ##     ## ##   ### 
   !  ######   #######  ##    ##  ######     ##    ##     ##  #######   ######     ##    ####  #######  ##    ## 
   ! 
   ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
   

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                             
   !                                             
   !   #    # #  ####   ####   ####  #    #  ####  
   !   #    # # #      #    # #    # #    # #      
   !   #    # #  ####  #      #    # #    #  ####  
   !   #    # #      # #      #    # #    #      # 
   !    #  #  # #    # #    # #    # #    # #    # 
   !     ##   #  ####   ####   ####   ####   ####  
   !                                               
   !                                               
   !   ###### #      #    # #    # ######  ####    
   !   #      #      #    #  #  #  #      #        
   !   #####  #      #    #   ##   #####   ####    
   !   #      #      #    #   ##   #           #   
   !   #      #      #    #  #  #  #      #    #   
   !   #      ######  ####  #    # ######  ####    
   !                                                                                                                                   
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  it computes the ξ, η and ζ components of the curvilinear viscous fluxes
   !
   !            _                                                                            _                                                               
   !           |                                    0                                         |
   !           |                                                                              |
   !           |    ∂     /   1   /   1           \    /   ∂u_1                 ∂ξ^j    \  \  |
   !           | ------- ( - --- (  ----- + νsgs   )* (   ------ * g^mj + R_m1 * ------  )  ) |
   !           |  ∂ξ^j    \   J   \   Re          /    \   ∂ξ^m                 ∂x_m    /  /  |    
   !           |                                                                              |
   !           |    ∂     /   1   /   1           \    /   ∂u_2                 ∂ξ^j    \  \  |
   !   visc =  | ------- ( - --- (  ----- + νsgs   )* (   ------ * g^mj + R_m2 * ------  )  ) |  
   !           |  ∂ξ^j    \   J   \   Re          /    \   ∂ξ^m                 ∂x_m    /  /  |
   !           |                                                                              |
   !           |    ∂     /   1   /   1           \    /   ∂u_3                 ∂ξ^j    \  \  |
   !           | ------- ( - --- (  ----- + νsgs   )* (   ------ * g^mj + R_m3 * ------  )  ) |
   !           |  ∂ξ^j    \   J   \   Re          /    \   ∂ξ^m                 ∂x_m    /  /  |
   !           |                                                                              |
   !           ⋅-                                                                            -⋅
   !
   !  Where the tensor R is the velocity gradient tensor
   !
   !  it returns the array visc with the three components of the viscous fluxes at every 
   !  node

   ! calculate expensive terms only during first stage
   ! viscous terms
   !
   
   if ( ViscousFluxes ) call rhs_viscous ()

   
   ! write(debugname, fmt ='(a,i6.6)') 'viscX',itc
   ! call outputD3_real( visc(1,:,:,:) , debugname )  

   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      

!#ifdef DEBUG
!   call cksum_4d_par ('daf.visc', visc)
!#endif

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                                  
   !                                                                 
   !      ##   #####  ##### # ###### #  ####  #   ##   #             
   !     #  #  #    #   #   # #      # #    # #  #  #  #             
   !    #    # #    #   #   # #####  # #      # #    # #             
   !    ###### #####    #   # #      # #      # ###### #             
   !    #    # #   #    #   # #      # #    # # #    # #             
   !    #    # #    #   #   # #      #  ####  # #    # ######        
   !                                                                 
   !                                                                 
   !    #####  #  ####   ####  # #####    ##   ##### #  ####  #    # 
   !    #    # # #      #      # #    #  #  #    #   # #    # ##   # 
   !    #    # #  ####   ####  # #    # #    #   #   # #    # # #  # 
   !    #    # #      #      # # #####  ######   #   # #    # #  # # 
   !    #    # # #    # #    # # #      #    #   #   # #    # #   ## 
   !    #####  #  ####   ####  # #      #    #   #   #  ####  #    # 
   !                                                                    !                                                                                                                                                                                                     
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   ! 
   ! Computes the third-order fourth difference artificial dissipation of Sotiropoulos and 
   ! Abdallah (1992), added to eliminate the odd-even decoupling of the pressure field due 
   ! to the nonstaggered mesh layout   
   !             _                                              _                                                               
   !            |                                                |
   !            |                                                |
   !            |                                                |
   !            |          ∂     /  g^11 * ∆τ    ∂^2 p  \        |
   !            |       ------- (  -----------  -------  )       |
   !            |        ∂ξ^2    \      J        ∂ξ^2   /        |    
   !            |                                                |
   !            |          ∂     /  g^22 * ∆τ    ∂^2 p  \        |
   !            |    +  ------- (  -----------  -------  )       |
   !            |        ∂η^2    \      J        ∂η^2   /        |
   !            |                                                |
   !            |                                                |
   !            |          ∂     /  g^33 * ∆τ    ∂^2 p  \        |
   !   diss  =  |    +  ------- (  -----------  -------  )       |
   !            |        ∂ζ^2    \      J        ∂ζ^2   /        |
   !            |                                                |
   !            |                                                |
   !            |                                                |
   !            |                                                |
   !            |                                                |
   !            |                                                |
   !            |                     0                          |  
   !            |                                                |
   !            |                     0                          |
   !            |                                                |
   !            |                     0                          |
   !            |                                                |
   !            ⋅-                                              -⋅

   !if (quick) then
   if ( DissipationFlux ) call rhs_diss_p ()

!#ifdef DEBUG
!   call cksum_4d_par ('daf.diss', diss)
!#endif


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                                  
   !                                                                        
   !    #    # ###### #       ####   ####  # ##### #   #                    
   !    #    # #      #      #    # #    # #   #    # #                     
   !    #    # #####  #      #    # #      #   #     #                      
   !    #    # #      #      #    # #      #   #     #                      
   !     #  #  #      #      #    # #    # #   #     #                      
   !      ##   ###### ######  ####   ####  #   #     #                      
   !                                                                        
   !                                                                        
   !    #####  # #    # ###### #####   ####  ###### #    #  ####  ######    
   !    #    # # #    # #      #    # #    # #      ##   # #    # #         
   !    #    # # #    # #####  #    # #      #####  # #  # #      #####     
   !    #    # # #    # #      #####  #  ### #      #  # # #      #         
   !    #    # #  #  #  #      #   #  #    # #      #   ## #    # #         
   !    #####  #   ##   ###### #    #  ####  ###### #    #  ####  ######    
   !                                                                        
   !                     ##                                                 
   !                    #  #                                                
   !                     ##                                                 
   !                    ###                                                 
   !                   #   # #                                              
   !                   #    #                                               
   !                    ###  #                                              
   !                                                                        
   !                                                                        
   !    #####  #####  ######  ####   ####  #    # #####  ######             
   !    #    # #    # #      #      #      #    # #    # #                  
   !    #    # #    # #####   ####   ####  #    # #    # #####              
   !    #####  #####  #           #      # #    # #####  #                  
   !    #      #   #  #      #    # #    # #    # #   #  #                  
   !    #      #    # ######  ####   ####   ####  #    # ######             
   !                                                                        
   !                                                                        
   !     ####  #####    ##   #####  # ###### #    # #####                   
   !    #    # #    #  #  #  #    # # #      ##   #   #                     
   !    #      #    # #    # #    # # #####  # #  #   #                     
   !    #  ### #####  ###### #    # # #      #  # #   #                     
   !    #    # #   #  #    # #    # # #      #   ##   #                     
   !     ####  #    # #    # #####  # ###### #    #   #                     
   !                                                                           !                                                                                                                                                                                                     
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  it computes the velocity divergence for the first component of the system equations
   !               _                              _                                                               
   !              |                                |
   !              |         ∂     /  U^j  \        |
   !              |      ------- (  -----  )       |
   !              |       ∂ξ^j    \   J   /        |    
   !              |                                |
   !              |       ∂     /  p   ∂ξ^j \      |
   !              |    ------- (  --- -----  )     |
   !              |     ∂ξ^j    \  J   ∂x1  /      |
   !              |                                |
   !   rh = rh +  |                                | 
   !              |      ∂     /  p   ∂ξ^j \       |
   !              |   ------- (  --- -----  )      |
   !              |    ∂ξ^j    \  J   ∂x2  /       |
   !              |                                |      
   !              |                                |
   !              |       ∂     /  p   ∂ξ^j \      | 
   !              |    ------- (  --- -----  )     | 
   !              |     ∂ξ^j    \  J   ∂x3  /      |
   !              ⋅-                              -⋅
   !
   !  Where the tensor R is the velocity gradient tensor
   !
   !  it returns the array rh updated with the velocity divergence in the first component

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      

   if( LinearFlux ) call rhs_flux_sans_convec ()

!#ifdef DEBUG
!     call cksum_4d_par ('daf.flux', rh)
!#endif


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                                       
   !   #    #  ####  #    #             #      # #    # ######   ##   #####  
   !   ##   # #    # ##   #             #      # ##   # #       #  #  #    # 
   !   # #  # #    # # #  #    #####    #      # # #  # #####  #    # #    # 
   !   #  # # #    # #  # #             #      # #  # # #      ###### #####  
   !   #   ## #    # #   ##             #      # #   ## #      #    # #   #  
   !   #    #  ####  #    #             ###### # #    # ###### #    # #    # 
   !                                                                         
   !                                                                         
   !   ###### #      #    # #    # ######  ####                              
   !   #      #      #    #  #  #  #      #                                  
   !   #####  #      #    #   ##   #####   ####                              
   !   #      #      #    #   ##   #           #                             
   !   #      #      #    #  #  #  #      #    #                             
   !   #      ######  ####  #    # ######  ####                              
   !                                                                            !                                                                                                                                   
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  it computes the ξ, η and ζ components of the curvilinear viscous fluxes
   !
   !               _                         _                                                               
   !              |                           |
   !              |            0              |
   !              |                           |
   !              |    ∂     /       U^j  \   |
   !              | ------- ( u_1 * -----  )  |
   !              |  ∂ξ^j    \        J   /   |    
   !              |                           |
   !              |    ∂     /       U^j  \   |
   !   rh =  rh + | ------- ( u_2 * -----  )  |  
   !              |  ∂ξ^j    \        J   /   |
   !              |                           |
   !              |    ∂     /       U^j  \   |
   !              | ------- ( u_3 * -----  )  |
   !              |  ∂ξ^j    \        J   /   |
   !              |                           |
   !              ⋅-                         -⋅
   !
   !
   !  it returns the array rh updated with the non-linear momentum flux

   if ( NonLinearFluxes ) then
   
      if ( n == 1 ) then
         call rhs_convec_quick (2)
      else
         call rhs_convec_quick (1)
      end if

   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      

!#ifdef DEBUG
!   call cksum_4d_par ('daf.quick', rh)
!#endif

   !else
   !   call rhs_diss_matrix ()
   !   call rhs_flux ()

   !#ifdef DEBUG
   !     call cksum_4d_par ('daf.matrix', rh)
   !#endif

   !end if


  ! compute non-linear closure terms
  ! 
!   if ( nlinc ) then
!      if ( craft .and. n == 1 ) then
!         call rhs_ns_nlin_craft ()
        
! #ifdef DEBUG
!         call cksum_4d_par ('daf.uij', uij)
! #endif

!      end if
!      call rhs_ns_rstress ()

! #ifdef DEBUG
!      call cksum_4d_par ('daf.rstr', rstr)
! #endif
  
!   end if
  

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                      
   !                                                                                  
   !   #    # #    #  ####  ##### ######   ##   #####  #   #      ##   #    # #####   
   !   #    # ##   # #        #   #       #  #  #    #  # #      #  #  ##   # #    #  
   !   #    # # #  #  ####    #   #####  #    # #    #   #      #    # # #  # #    #  
   !   #    # #  # #      #   #   #      ###### #    #   #      ###### #  # # #    #  
   !   #    # #   ## #    #   #   #      #    # #    #   #      #    # #   ## #    #  
   !    ####  #    #  ####    #   ###### #    # #####    #      #    # #    # #####   
   !                                                                                  
   !                                                                                  
   !    ####   ####  #    # #####   ####  ######    ##### ###### #####  #    #  ####  
   !   #      #    # #    # #    # #    # #           #   #      #    # ##  ## #      
   !    ####  #    # #    # #    # #      #####       #   #####  #    # # ## #  ####  
   !        # #    # #    # #####  #      #           #   #      #####  #    #      # 
   !   #    # #    # #    # #   #  #    # #           #   #      #   #  #    # #    # 
   !    ####   ####   ####  #    #  ####  ######      #   ###### #    # #    #  ####  
   !                                                                                                           
   !                                                                            !                                                                                                                                   
   ! *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *                                                                                         
   !
   !  it computes the ξ, η and ζ components of the curvilinear viscous fluxes
   !
   !    rh = rh + diss + visc
   !               _                                              _                                                               
   !              |                                                |
   !              |                       0                        |
   !              |                                                |
   !              |   3u^{l,n+1} - 4u^n + u^{n-1}                  |
   !              | ------------------------------                 |
   !              |           2 * Δt * J                           |    
   !              |                                                |
   !              |   3v^{l,n+1} - 4v^n + v^{n-1}                  |
   !   rh =  rh + | ------------------------------                 |  
   !              |           2 * Δt * J                           |
   !              |                                                |
   !              |   3w^{l,n+1} - 4w^n + w^{n-1}          1       |
   !              | ------------------------------ +  -----------  |
   !              |           2 * Δt * J                J * Fr^2   |
   !              |                                                |
   !              ⋅-                                              -⋅
   !
   !  it returns the array rh updated with the unsteady and source terms
   
   if ( UnsteadyTerm  ) call rhs_unst_visc_diss ()
   if ( GravitySource ) call rhs_gravity_source_term ()

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    

   ! Putting everything together
   ! I clean the rhs values from the air phase
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   ! i_myend
      rh(1,i,j,k) = rsign( phi(i,j,k) ) * ( rh(1,i,j,k) + diss(  i,j,k) )
      rh(2,i,j,k) = rsign( phi(i,j,k) ) * ( rh(2,i,j,k) - visc(1,i,j,k) )
      rh(3,i,j,k) = rsign( phi(i,j,k) ) * ( rh(3,i,j,k) - visc(2,i,j,k) )
      rh(4,i,j,k) = rsign( phi(i,j,k) ) * ( rh(4,i,j,k) - visc(3,i,j,k) )
   end do
   end do
   end do
   
   ! RHS extrapolation to the first air layer

   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update again
   call rhs_exchng3_4d( rh )

   !write(debugname, fmt ='(a,i6.6)') 'rsign_solverdaf',iteraciontiempo
   !call outputD3_real( rsign(:,:,:) , debugname ) 


   ! add forcing term on coarse grids
   !
   !if ( n > 1 ) then
   !   if (decide_calc_pk == 1) then
   !      call rhs_pk (1)
   !   else
   !      call rhs_pk (0)
   !   end if
   !end if


   ! rh is replaced in rhs_diag_solver by brh
   if ( myback  == mpi_proc_null ) rh(:,i_mysta-1,:,:) = zero
   if ( myfront == mpi_proc_null ) rh(:,i_myend+1,:,:) = zero
   if ( myleft  == mpi_proc_null ) rh(:,:,j_mysta-1,:) = zero
   if ( myright == mpi_proc_null ) rh(:,:,j_myend+1,:) = zero
   if ( mydown  == mpi_proc_null ) rh(:,:,:,k_mysta-1) = zero
   if ( myup    == mpi_proc_null ) rh(:,:,:,k_myend+1) = zero

   if (nblk /= 0) then
   
      do nb = 1, nblk
         rh( : , li_blk_ia(n,nb) : li_blk_ib(n,nb) , &
                 li_blk_ja(n,nb) : li_blk_jb(n,nb) , &
                 li_blk_ka(n,nb) : li_blk_kb(n,nb)     ) = zero
      end do
   
   end if

  
   ! calculate model matrices
   call rhs_modal_matrices()

   
   ! implicit diagonal solver
   call rhs_diag_solver ()

   
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   ! i_myend
      
         q(1,i,j,k) = qold_af(1,i,j,k) + rsign( phi(i,j,k) ) * rh(1,i,j,k)
         q(2,i,j,k) = qold_af(2,i,j,k) + rsign( phi(i,j,k) ) * rh(2,i,j,k)
         q(3,i,j,k) = qold_af(3,i,j,k) + rsign( phi(i,j,k) ) * rh(3,i,j,k)
         q(4,i,j,k) = qold_af(4,i,j,k) + rsign( phi(i,j,k) ) * rh(4,i,j,k)
   
   end do
   end do
   end do


   if (n == 1) call bcond_fm ( il, iu, jl, ju, kl, ku, igp, jgp, kgp, dc, de, dz, &
                               q , eta , aj, x,y,z , xnut , phi )


   call rhs_exchng3_4d (q)   !tiene que llegar con gp actualizados



   ! deallocate dummy variables
   deallocate ( qold_af, diss, visc, ucn_j, dtau, fv )
   deallocate ( mai, n1i, n2i, mc, spr )

contains

   include 'rhs_contra_j.F90'
   include 'rhs_daf_dtau.F90'
   include 'rhs_viscous.F90'
   include 'rhs_diss_p.F90'
   include 'rhs_flux_sans_convec.F90'
   include 'rhs_convec_quick.F90'
   !include 'rhs_pk.F90'
   include 'rhs_exchng3_3d.F90'
   include 'rhs_exchng3_4d.F90'
   !include 'rhs_ghost_fluid_nodes_extp_3d.F90'
   include 'rhs_ghost_fluid_nodes_extp_4d.F90'
   include 'rhs_modal_matrices.F90'
   include 'rhs_diag_solver.F90'
   include 'rhs_unst_visc_diss.F90'
   include 'rhs_gravity_source_term.F90'
   !include 'near_free_surface_q_update.F90'
   include 'PhiGradientVector.F90'
   include 'VelocityGradientTensor.F90'

!  include 'rhs_flux.F90'
!  include 'rhs_diss_quick.F90'
!  include 'rhs_diss_matrix.F90'
!  include 'rhs_ns_nlin_craft.F90'
!  include 'rhs_ns_rstress.F90'
!  include 'rhs_ns_unst_visc_diss.F90'

end subroutine solver_daf


