subroutine levelsetmethod( il,iu                   , &
                           jl,ju                   , &
                           kl,ku                   , &
                           igp, jgp, kgp           , &
                           dc, de, dz              , &
                           q                       , &
                           csi                     , &
                           eta                     , &
                           zet                     , &
                           aj                      , &
                           phi, phi_n              , &
                           x, y, z                 , &
                           iteraciontiempo   )
   use global_param
   use global_app
   use global_mpi
   use checksum

   use global_debug
   use global_lsm, only : IGitermax, IGtimestep,RNitermax,RNtimestep, RNfreq                       , & 
                          heaviside, deltaf, btype_lsm, numiter_lsm, delti_lsm,isnarrow,narrowcoef , &
                          orderAD, bc_extrapolation_order, call_levelsetmethod                     , &
                          call_reinitialisation, hybrid_reinitialisation, TotalVolumeComputation   , &
                          OrderLSAdvectionBoundaries, phi_outputiter                               , & 
                          OrderReinitialisationBoundaries , ENOBCReinitialisation                  , &
                          GlobalMassCorrection , BigPhi , limit_ghost_velocities, epslsm

   use global_obstacle, only : obstacle_lsm_ad, boundary_obstacle, is_obstacle, act_obstacle_ad, &
                               act_obstacle_rn,li_obs_ia, li_obs_ib, li_obs_ja, li_obs_jb, li_obs_ka,&
                               li_obs_kb

   use AdvectionMethods

   implicit none

   !argumentos recibidos !asumir que los gp vienen actualizados!
   integer, intent(in) :: il,iu,jl,ju,kl,ku
   integer, intent(in) :: igp,jgp,kgp
   real (kind = rdf), intent(in) :: dc,de,dz
   real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: q
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi , eta, zet
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: aj
   !real (kind = rdf), intent(in) :: epslsm
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku),intent(in) :: x,y,z
   !real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: rsign
   !real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) , intent(inout) :: phi_gradient ! ∂phi/∂x_j

   !argumentos de entrada y salida

   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in)    :: phi_n
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: phi

   !local
   real (kind = rdf), dimension(:,:,:), allocatable ::  phicorrected

   integer, dimension(:,:,:), allocatable :: InterfaceNodesID 
   integer, dimension(:,:,:), allocatable :: AdvectionNodes 

   ! ___________________________________________________________________________________
   ! 
   !    WATER VOLUMES FOR GEOMETRIC REINITIALISATION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   !
   ! Volumes Preadvection
   real (kind = rdf) :: Vol_KTet_PreAdv               , &
                        Vol_2ndTet_PreAdv             , & 
                        Vol_BulkCells_PreAdv 
   
   ! Total Local Volulme Preadvection
   real (kind = rdf) :: TotVol_PreAdv_Local
   
   ! Total Global Volulme Preadvection
   real (kind = rdf) :: TotVol_PreAdv

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Volumes after advection 
   real (kind = rdf) :: Vol_KTet_PostAdv_Local       , &
                        Vol_2ndTet_PostAdv_Local     , &
                        Vol_BulkCells_PostAdv_Local 


   ! Global Volumes after advection after mpi_allreduce
   real (kind = rdf) :: Vol_KTet_PostAdv           , &
                        Vol_2ndTet_PostAdv         , &
                        Vol_BulkCells_PostAdv 

   ! Total Local Volulme Preadvection
   real (kind = rdf) :: TotVol_PostAdv_Local
   
   ! Total Global Volulme Preadvection
   real (kind = rdf) :: TotVol_PostAdv


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Volumes after reinitialisation 

   real (kind = rdf) :: Vol_KTet_PostReini_Local        , &
                        Vol_2ndTet_PostReini_Local      , &
                        Vol_BulkCells_PostReini_Local 
   

   real (kind = rdf) :: Vol_KTet_PostReini              , &
                        Vol_2ndTet_PostReini            , &
                        Vol_BulkCells_PostReini 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real (kind = rdf) :: ConvergenceVolume
   
   ! ///////////////////////////////////////////////////////////////////////////////////


   real (kind = rdf), dimension(:,:,:,:), allocatable :: ucn_j

   integer :: RK
   integer :: i,  j,  k, iter_lsm

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend
              
   integer :: TriangulationCaseAux
   integer, save :: TriangulationBitFlag 

   ! - - - -  DEBUG VARIABLES - - - - - - - - - - - - - - - - 
   
   character(len = 256) :: debugname
   integer :: iteraciontiempo
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Initialisation

   allocate ( InterfaceNodesID( il:iu , jl:ju , kl:ku ) , &
                AdvectionNodes( il:iu , jl:ju , kl:ku )  )

   allocate ( phicorrected( il:iu , jl:ju , kl:ku ) )
   phicorrected = zero


   TriangulationCaseAux = 1 

   !========================================================================================= 
   !                                                                            
   !    ####  #       ####  #####    ##   #         #    #   ##    ####   ####  
   !   #    # #      #    # #    #  #  #  #         ##  ##  #  #  #      #      
   !   #      #      #    # #####  #    # #         # ## # #    #  ####   ####  
   !   #  ### #      #    # #    # ###### #         #    # ######      #      # 
   !   #    # #      #    # #    # #    # #         #    # #    # #    # #    # 
   !    ####  ######  ####  #####  #    # ######    #    # #    #  ####   ####  
   !                                                                            
   !                                                                            
   !    ####   ####  #    # #####  #    # #####   ##   ##### #  ####  #    #    
   !   #    # #    # ##  ## #    # #    #   #    #  #    #   # #    # ##   #    
   !   #      #    # # ## # #    # #    #   #   #    #   #   # #    # # #  #    
   !   #      #    # #    # #####  #    #   #   ######   #   # #    # #  # #    
   !   #    # #    # #    # #      #    #   #   #    #   #   # #    # #   ##    
   !    ####   ####  #    # #       ####    #   #    #   #   #  ####  #    #    
   !                                                                               
   !========================================================================================= 

   
   ! This command swaps TriangulationBitFlag between 0 and 1 circulary every iteration
   !TriangulationBitFlag = xor( TriangulationBitFlag, 1 )

   call TotalWaterVolume ( phi                                        , &             
                           Vol_KTet_PreAdv                            , &
                           Vol_2ndTet_PreAdv                          , &
                           Vol_BulkCells_PreAdv                       , &
                           AdvectionNodes                             , &
                           WriteGlobalVolume = .true.                 , &    
                           TriangulationCase = TriangulationCaseAux     &
                         )

   TotVol_PreAdv_Local =   Vol_KTet_PreAdv        + &
                           Vol_2ndTet_PreAdv      + &
                           Vol_BulkCells_PreAdv


   ! Global volume over processors by reduction
   call mpi_allreduce( TotVol_PreAdv_Local, TotVol_PreAdv , &
                       1 , MPI_REAL_TYPE , mpi_sum , mpi_comm_world, ierr )

   if ( myid == root ) call WriteTotalWaterVolume( TotVol_PreAdv )


   !========================================================================================= 
   !                                                                                         
   !    ##   #####  #    # ######  ####  ##### #  ####  #    #     ####  ##### ###### #####  
   !   #  #  #    # #    # #      #    #   #   # #    # ##   #    #        #   #      #    # 
   !  #    # #    # #    # #####  #        #   # #    # # #  #     ####    #   #####  #    # 
   !  ###### #    # #    # #      #        #   # #    # #  # #         #   #   #      #####  
   !  #    # #    #  #  #  #      #    #   #   # #    # #   ##    #    #   #   #      #      
   !  #    # #####    ##   ######  ####    #   #  ####  #    #     ####    #   ###### #      
   !                                                                                         
   !========================================================================================= 
   ! ϕ-advection:                             !
   ! -------------                            !
   !                                          !
   !      ∂ϕ    U^j    ∂ϕ                     !  
   !      -- + ----- ------ = 0               !
   !      ∂t     J    ∂ξ^j                    !
   !                                          !
   ! - - - - - - - - - - - - - - - - - - - - -*
   
   if ( ( .not. call_reinitialisation   ) .or. &
        ( .not. hybrid_reinitialisation ) ) AdvectionNodes = 1
   
   if ( call_levelsetmethod ) then
         
      ! U^j/J for ϕ-advection  
      allocate (ucn_j(1:3,il:iu,jl:ju,kl:ku))
      ucn_j = zero
      
      call rhs_contra_j()
   
      ! I use phicorrected as auxiliar array to store phi_n each iteration
      phicorrected = phi_n

      ! Level-Set advection iterations
      do iter_lsm = 1, numiter_lsm
     
         !call narrowband()
         call ad_RKTVD3( phi , phicorrected , AdvectionNodes )
         phicorrected = phi     !gp de phi vienen ya actualizados de rutina ad_RKTVD3
     
      end do
   
      deallocate (ucn_j)
   
   end if
   
   ! If GlobalMassCorrection is false, ConvergenceVolume is a dummy variable
   ConvergenceVolume = TotVol_PreAdv
   
   if ( GlobalMassCorrection ) then

      call TotalWaterVolume( phi                          , &             
                             Vol_KTet_PostAdv_Local       , &
                             Vol_2ndTet_PostAdv_Local     , &
                             Vol_BulkCells_PostAdv_Local      )
      
      
      ! Get the global volume corresponding to K tetrahedrons and 2nd
      ! tetrahedrons by reducing over all the processors
      call mpi_allreduce( Vol_KTet_PostAdv_Local, Vol_KTet_PostAdv , &
                          1 , MPI_REAL_TYPE , mpi_sum , mpi_comm_world, ierr )
      
      call mpi_allreduce( Vol_2ndTet_PostAdv_Local, Vol_2ndTet_PostAdv , &
                          1 , MPI_REAL_TYPE , mpi_sum , mpi_comm_world, ierr )
         
      ConvergenceVolume =   TotVol_PreAdv - ( Vol_KTet_PostAdv + Vol_2ndTet_PostAdv ) 

   end if
   
   !========================================================================================= 
   !                                                                                      
   !  #####  ###### # #    # # ##### #   ##   #      #  ####    ##   ##### #  ####  #    #    
   !  #    # #      # ##   # #   #   #  #  #  #      # #       #  #    #   # #    # ##   #    
   !  #    # #####  # # #  # #   #   # #    # #      #  ####  #    #   #   # #    # # #  #    
   !  #####  #      # #  # # #   #   # ###### #      #      # ######   #   # #    # #  # #    
   !  #   #  #      # #   ## #   #   # #    # #      # #    # #    #   #   # #    # #   ##    
   !  #    # ###### # #    # #   #   # #    # ###### #  ####  #    #   #   #  ####  #    #    
   !                                                                                          
   !                                                                                          
   !   ####  ##### ###### #####                                                               
   !  #        #   #      #    #                                                              
   !   ####    #   #####  #    #                                                              
   !       #   #   #      #####                                                               
   !  #    #   #   #      #                                                                   
   !   ####    #   ###### #                                                                   
   !                                                                                      
   !========================================================================================= 
   
   !write(debugname, fmt ='(a,i6.6)') 'phi_postadvection',iteraciontiempo
   !call outputD3_real( phi(:,:,:) , debugname ) 

   if ( call_reinitialisation .and. mod( ntime , RNfreq ) == 0 ) then
             
      ! Flag variable for interface nodes initialisation
      InterfaceNodesID = 0
      
      !---------------------------------------------------------------------------------
      ! When the hybrid reinitialisation is activated. The ϕ function is updated
      ! at the nodes next to the free-surface and the rest remain unchanged.
      !---------------------------------------------------------------------------------
   
      if (hybrid_reinitialisation) then
         
         where ( phi < -eps_sims ) phicorrected = -9.9_rdf
         where ( phi >  eps_sims ) phicorrected =  9.9_rdf
            
         call geometric_reinitialisation( phi , phicorrected , &
                                          ConvergenceVolume      , & 
                                          InterfaceNodesID       , & 
                                          TriangulationCase = TriangulationCaseAux    )
   
      else
         ! Sussman reinitialisation
           
         call reinitialisation_test2 ( phi, phicorrected, InterfaceNodesID )
         
      end if

      ! Apply Boundary Contions to the reinitialised ϕ 
      call bcond_lsm( phicorrected )    

      ! ϕ is finally updated after the reinitialisation step         
      phi = phicorrected
      
      ! Update the ghost nodes with the reinitialised phi values
      call rhs_exchng3_3d( phi )
      
   end if
   
   deallocate ( phicorrected )
   deallocate ( InterfaceNodesID , AdvectionNodes )

   contains
     
   include 'bcond_lsm.F90'
   include 'bcond_fm.F90'
   include 'bcond_lsm_obstacle.F90'
   include 'rhs_contra_j.F90'
   include 'ad_RKTVD3.F90'
   !include 'calc_RH_AD_ENO3.F90'
   include 'calc_RH_AD_WENO3.F90'
   include 'calc_RH_AD.F90'
   !include 'calc_RH_AD_obstacle.F90'
   include 'rhs_exchng3_3d.F90'
   include 'rhs_exchng3_4d.F90'
   !include 'narrowband.F90'
   !include 'reinitialization.F90'
   include 'reinitialization_benchmark.F90'
   include 'reinitialisation_test2.F90'
   include 'calc_NormGrad_d.F90'
   include 'calc_RH_RN.F90'
   include 'calc_RH_RN_2.F90'
   include 'calc_RH_RN_obstacle.F90'
   include 'calc_sign_ini.F90'
   include 'calc_sign_ini_2.F90'
   include 'calc_sign_ini_obstacle.F90'
   include 'calc_lambda_RN.F90'
   include 'testfilter.F90'
   !include 'p_correction_post_advection.F90'
   !include 'pressure_extrapolation.F90'
   !include 'get_phi_gradient.F90'
   ! geometric reinitialisation routines
   
   include 'geometric_reinitialisation.F90'
   include 'GetTetrahedraList.F90'
   include 'GetIsosurfaceAndDistances.F90'
   include 'PiecewiseConstantFunction.F90'
   include 'OrthogonalProjection.F90'
   include 'NeighboursPhiCorrection.F90'
   include 'WriteTriangulation.F90'
   include 'WriteTotalWaterVolume.F90'
   include 'TotalWaterVolume.F90'

end subroutine levelsetmethod
