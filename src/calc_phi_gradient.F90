subroutine calc_phi_gradient( il,iu          ,&
                              jl,ju          ,&
                              kl,ku          ,&
                              igp, jgp, kgp  ,&
                              dc, de, dz     ,&
                              csi            ,&
                              eta            ,&
                              zet            ,&
                              phi            ,&
                              phi_gradient    &
                            )

   use global_mpi
   use global_app

   use AdvectionMethods

   implicit none

   ! Input variables

   integer, intent(in) :: il,iu,jl,ju,kl,ku ! external nodes
   integer, intent(in) :: igp, jgp, kgp ! # of ghostpoint at local processor 
   real (kind = rdf), intent(in) :: dc,de,dz ! dcsi, deta, dzet
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi , eta, zet ! metrics
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi

   ! Output variable
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) , intent(inout) :: phi_gradient

   !local

   !indexes
   integer  :: i_mysta , i_myend
   integer  :: j_mysta , j_myend
   integer  :: k_mysta , k_myend

   integer  :: ista , iend
   integer  :: jsta , jend
   integer  :: ksta , kend
   
   integer :: i , j , k 
   
   ! phi gradient at internal points
   real (kind = rdf) :: dphi_dcsi , dphi_deta , dphi_dzet
   real (kind = rdf) :: dphi_dx   , dphi_dy   , dphi_dz
   real (kind = rdf) :: dc2 , de2 , dz2 ! ( 1/(2Δξ) , 1/(2Δη) , 1/(2Δζ) )
   real (kind = rdf) :: dummy
   
   integer :: iBiasDirection, jBiasDirection, kBiasDirection

   !Nodos incluyendo el borde
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp
         
   ! processes on the domain boundaries

   if ( myback  == mpi_proc_null )  i_mysta = il + igp + 1
   if ( myleft  == mpi_proc_null )  j_mysta = jl + jgp + 1
   if ( mydown  == mpi_proc_null )  k_mysta = kl + kgp + 1
   
   if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1
   if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1
   if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1


   ! Physical boundaries of the processor
   
   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp
   
   dc2 = one_half * dc
   de2 = one_half * de
   dz2 = one_half * dz

   dummy = zero

   ! variable initialisation
   phi_gradient = zero

   do k = ksta+1, kend-1 !k_mysta, k_myend
   do j = jsta+1, jend-1 !j_mysta, j_myend
   do i = ista+1, iend-1 !i_mysta, i_myend
     
      ! central difference for phi curvilinear derivatives
      dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) ) ! ∂ϕ/∂ξ 
      dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) ) ! ∂ϕ/∂η
      dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) ) ! ∂ϕ/∂ζ
    
      ! ∂ϕ/∂x
      phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                 + dphi_deta * eta(1,i,j,k) &
                                 + dphi_dzet * zet(1,i,j,k) 
            
      ! ∂ϕ/∂y
      phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                 + dphi_deta * eta(2,i,j,k) &
                                 + dphi_dzet * zet(2,i,j,k) 
            
      ! ∂ϕ/∂z
      phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                 + dphi_deta * eta(3,i,j,k) &
                                 + dphi_dzet * zet(3,i,j,k) 
   
   end do
   end do
   end do
   
   
   ! phi gradient at the boundaries
   ! One-sided derivatives at the whole domain borders
   ! expresions obtained from section 3.7.1 (Implementation of
   ! Boundary conditions using internal grid points) of Ferziger's
   ! book
   
   if ( myback == mpi_proc_null ) then

      !i global = 1
      i = ista

      do k = k_mysta , k_myend
      do j = j_mysta , j_myend
            
         ! one sided derivative
         dphi_dcsi = dc2 * (  - one   * phi(i+2,j,k) &
                              + four  * phi(i+1,j,k) &
                              - three * phi(i  ,j,k) )
         
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      
      end do
      end do
   
   end if
   
   if ( myfront == mpi_proc_null ) then

      !i global = im
      i = iend
   
      do k = k_mysta , k_myend
      do j = j_mysta , j_myend
            
         ! one-sided derivative
         dphi_dcsi = dc2 * (    one   * phi(i-2,j,k) &
                              - four  * phi(i-1,j,k) &
                              + three * phi(i  ,j,k) )
            
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      end do
      end do
   
   end if
   
   
   if ( myleft == mpi_proc_null ) then
   
      ! j global  = 1
      j = jsta
   
      do k = k_mysta , k_myend
      do i = i_mysta , i_myend
            
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) )
         
         ! one sided derivative
         dphi_deta = de2 * (  - one   * phi(i,j+2,k) &
                              + four  * phi(i,j+1,k) &
                              - three * phi(i,j  ,k) )
            
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      end do
      end do

   end if
   
   
   
   if ( myright == mpi_proc_null ) then   
   
      ! j global  = jm
      j = jend
      
      do k = k_mysta , k_myend
      do i = i_mysta , i_myend
            
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) )
         
         ! one sided derivative
         dphi_deta = de2 * (    one   * phi(i,j-2,k) &
                              - four  * phi(i,j-1,k) &
                              + three * phi(i,j  ,k) )
            
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      end do
      end do
      
   end if
   
   
   if ( mydown == mpi_proc_null ) then
   
      ! k global = 1
      k = ksta
      
      do j = j_mysta , j_myend
      do i = i_mysta , i_myend
            
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) )
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
            
         ! one-sided derivative
         dphi_dzet = dz2 * ( - one   * phi(i,j,k+2) &
                           & + four  * phi(i,j,k+1) &
                           & - three * phi(i,j,k  ) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      end do
      end do
   
   
   end if
      
   if ( myup == mpi_proc_null ) then

      ! k global = km
      k = kend
   
      do j = j_mysta , j_myend
      do i = i_mysta , i_myend
            
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) )
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
            
         ! one-sided derivative
         dphi_dzet = dz2 * (    one   * phi(i,j,k-2) &
                              - four  * phi(i,j,k-1) &
                              + three * phi(i,j,k  ) )
      
         ! ∂ϕ/∂x
         phi_gradient(1,i,j,k) =      dphi_dcsi * csi(1,i,j,k) &
                                    + dphi_deta * eta(1,i,j,k) &
                                    + dphi_dzet * zet(1,i,j,k) 
      
         ! ∂ϕ/∂y
         phi_gradient(2,i,j,k) =      dphi_dcsi * csi(2,i,j,k) &
                                    + dphi_deta * eta(2,i,j,k) &
                                    + dphi_dzet * zet(2,i,j,k) 
            
         ! ∂ϕ/∂z
         phi_gradient(3,i,j,k) =      dphi_dcsi * csi(3,i,j,k) &
                                    + dphi_deta * eta(3,i,j,k) &
                                    + dphi_dzet * zet(3,i,j,k) 
      
      end do
      end do
   
   end if



   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   !
   !    
   ! ####### ######   #####  #######  #####  
   ! #       #     # #       #       #      
   ! #####   #     # #  #### #####    #####  
   ! #       #     # #     # #             # 
   ! ####### ######   #####  #######  #####  
   !                                         
   !
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   !
   !
   !            i1,jm,km ---E11--- im,jm,km    
   !                /|(8)           /|(7)                            
   !        E12----/ |             /--------E10                                   
   !            i1,j1,km---E9---im,j1,km                   
   !              |(5)           |(6)|                                   
   !              | E8       E3  |   E7                                
   !             E5  |       |   E6  |                                
   !              | i1,jm,k1-----|-im,jm,k1           
   !          E4--|-/(4)         |  /(3)                 
   !              |/             | / -------E2
   ! ζ,k      i1,j1,k1---E1----im,j1,k1
   ! ^   η,j     (1)            (2)   
   ! |  7             
   ! | /
   ! |/
   ! +-----------> ξ,i
   
   ! Edge 1  : 1 - 2 --> i free ; j = 1  , k = 1
   ! Edge 2  : 2 - 3 --> j free ; i = im , k = 1
   ! Edge 3  : 3 - 4 --> i free ; j = jm , k = 1
   ! Edge 4  : 4 - 1 --> j free ; i = 1  , k = 1
   !
   ! Edge 5  : 1 - 5 --> k free ; i = 1  , j = 1
   ! Edge 6  : 2 - 6 --> k free ; i = im , j = 1
   ! Edge 7  : 3 - 7 --> k free ; i = im , j = jm
   ! Edge 8  : 4 - 8 --> k free ; i = 1  , j = jm
   !
   ! Edge 9  : 5 - 6 --> i free ; j = 1  , k = km
   ! Edge 10 : 6 - 7 --> j free ; i = im , k = km
   ! Edge 11 : 7 - 8 --> i free ; j = jm , k = km
   ! Edge 12 : 8 - 5 --> j free ; i = 1  , k = km
   !
   ! i free : E1 , E3 , E9  , E11
   ! j free : E2 , E4 , E10 , E12
   ! k free : E5 , E6 , E7  , E8
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   
   !**************************
   !
   ! i - free edges
   !
   !**************************
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 1
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myleft == mpi_proc_null .and. mydown == mpi_proc_null)  then
   
      j = jsta !jl + jgp
      k = ksta !kl + kgp
   
      jBiasDirection = 1
      kBiasDirection = 1
   
      !do i = iRN_sta , iRN_end
      do i = i_mysta , i_myend
         
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) ) ! ∂ϕ/∂ξ 
      
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                           )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 3
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myright == mpi_proc_null .and. mydown == mpi_proc_null)  then
   
      j = jend ! ju - jgp
      k = ksta ! kl + kgp
   
      jBiasDirection = -1
      kBiasDirection =  1
   
      !do i = iRN_sta , iRN_end
      do i = i_mysta , i_myend
            
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) ) ! ∂ϕ/∂ξ 
      
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         ! |grad(phi)| 
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 9
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myleft == mpi_proc_null .and. myup == mpi_proc_null)  then
   
      j = jsta ! jl + jgp
      k = kend ! ku - kgp
   
      jBiasDirection =  1
      kBiasDirection = -1
   
   !   do i = iRN_sta , iRN_end
      do i = i_mysta , i_myend
          
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) ) ! ∂ϕ/∂ξ 
      
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 11
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myright == mpi_proc_null .and. myup == mpi_proc_null)  then
   
      j = jend ! ju - jgp
      k = kend ! ku - kgp
   
      jBiasDirection = -1
      kBiasDirection = -1
   
   !   do i = iRN_sta , iRN_end
      do i = i_mysta , i_myend
         
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dcsi = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) ) ! ∂ϕ/∂ξ 
      
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   !**************************
   !
   ! j - free edges
   !
   !**************************
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 2
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myfront == mpi_proc_null .and. mydown == mpi_proc_null)  then
   
      i = iend ! iu - igp
      k = ksta ! kl + kgp
   
      iBiasDirection = -1
      kBiasDirection =  1
   
   !   do j = jRN_sta , jRN_end
      do j = j_mysta , j_myend
      
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
      
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 4
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myback == mpi_proc_null .and. mydown == mpi_proc_null)  then
   
      i = ista ! il + igp
      k = ksta ! kl + kgp
   
      iBiasDirection =  1
      kBiasDirection =  1
   
   !   do j = jRN_sta , jRN_end
      do j = j_mysta , j_myend
      
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
      
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 10
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myfront == mpi_proc_null .and. myup == mpi_proc_null)  then
   
      i = iend ! iu - igp
      k = kend ! ku - kgp
   
      iBiasDirection = -1
      kBiasDirection = -1
   
      do j = j_mysta , j_myend
      
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
      
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 12
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myback == mpi_proc_null .and. myup == mpi_proc_null)  then
   
      i = ista ! il + igp
      k = kend ! ku - kgp
   
      iBiasDirection =  1
      kBiasDirection = -1
   
      do j = j_mysta , j_myend
      
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_deta = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )
      
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                            phi( i , j , k + kBiasDirection * 1 ) , &
                                            phi( i , j , k + kBiasDirection * 2 ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            kBiasDirection                          &
                                          )
   
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   !**************************
   !
   ! k-free edges
   !
   !**************************
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 5
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myback == mpi_proc_null .and. myleft == mpi_proc_null)  then
   
      i = ista ! il + igp
      j = jsta ! jl + jgp
   
      iBiasDirection =  1
      jBiasDirection =  1
   
      do k = k_mysta , k_myend
         
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 6
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null)  then
   
      i = iend ! iu - igp
      j = jsta ! jl + jgp
   
      iBiasDirection = -1
      jBiasDirection =  1
   
   !   do k = kRN_sta , kRN_end
      do k = k_mysta , k_myend
         
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
     
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 7
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myfront == mpi_proc_null .and. myright == mpi_proc_null)  then
   
      i = iend ! iu - igp
      j = jend ! ju - jgp
   
      iBiasDirection = -1
      jBiasDirection = -1
   
   !   do k = kRN_sta , kRN_end
      do k = k_mysta , k_myend
         
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   ! EDGE 8
   ! - - - - - - - - - - - - - - - - - - - - - - - - 
   
   if ( myback == mpi_proc_null .and. myright == mpi_proc_null)  then
   
      i = ista ! il + igp
      j = jend ! ju - jgp
   
      iBiasDirection =  1
      jBiasDirection = -1
   
      do k = k_mysta , k_myend
   
         !--------------------------------------------------------------------
         ! ξ - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                            phi( i + iBiasDirection * 1 , j , k ) , &
                                            phi( i + iBiasDirection * 2 , j , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            iBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! η - DIRECTION
         !--------------------------------------------------------------------
   
         dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                            phi( i , j + jBiasDirection * 1 , k ) , &
                                            phi( i , j + jBiasDirection * 2 , k ) , &
                                            dummy                                 , &
                                            2                                     , &
                                            jBiasDirection                          &
                                          )
   
         !--------------------------------------------------------------------
         ! ζ - DIRECTION
         !--------------------------------------------------------------------
      
         dphi_dzet = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )
      
         ! Calculo del gradiente en coordenadas curvilineas
      
         dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                   eta( 1 , i , j , k ) * dphi_deta + &
                   zet( 1 , i , j , k ) * dphi_dzet
      
         dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                   eta( 2 , i , j , k ) * dphi_deta + &
                   zet( 2 , i , j , k ) * dphi_dzet
      
         dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                   eta( 3 , i , j , k ) * dphi_deta + &
                   zet( 3 , i , j , k ) * dphi_dzet
      
         phi_gradient(1,i,j,k) = dphi_dx
         phi_gradient(2,i,j,k) = dphi_dy
         phi_gradient(3,i,j,k) = dphi_dz
   
      end do
   
   end if
   
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   !
   ! #     # ####### ######  ####### ###  #####  #######  #####  
   ! #     # #       #     #    #     #  #       #       #       
   ! #     # #####   ######     #     #  #       #####    #####  
   !  #   #  #       #   #      #     #  #       #             # 
   !    #    ####### #    #     #    ###  #####  #######  #####  
   !                                                            
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   
   
      !    i,j+1,k+1 ----i+1,j+1,k+1    
      !     /|(8)           /|(7)                            
      !    / |             / |                                  
      ! i,j,k+1-------i+1,j,k+1                   
      !   |(5)           |(6)|                                   
      !   |  |           |   |                                
      !   |  |           |   |                                
      !   |  i,j+1,k-----|-i+1,j+1,k           
      !   | /(4)         |  /(3)                 
      !   |/             | /
      ! i,j,k-----------i+1,j,k
      !  (1)               (2)
   
   
   
   ! VERTEX 1
   
   if ( myback == mpi_proc_null .and. &
        myleft == mpi_proc_null .and. &
        mydown == mpi_proc_null           )  then
   
      i = ista ! il + igp
      j = jsta ! jl + jgp
      k = ksta ! kl + kgp
   
      iBiasDirection =  1
      jBiasDirection =  1
      kBiasDirection =  1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   end if
   
   
   ! VERTEX 2
   
   if ( myfront == mpi_proc_null .and. &
        myleft  == mpi_proc_null .and. &
        mydown  == mpi_proc_null          )  then
   
      
      i = iend ! iu - igp
      j = jsta ! jl + jgp
      k = ksta ! kl + kgp
   
      iBiasDirection = -1
      jBiasDirection =  1
      kBiasDirection =  1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   end if
   
   ! VERTEX 3
   
   if ( myfront == mpi_proc_null .and. &
        myright == mpi_proc_null .and. &
        mydown  == mpi_proc_null          )  then
   
      i = iend ! iu - igp
      j = jend ! ju - jgp
      k = ksta ! kl + kgp
   
      iBiasDirection = -1
      jBiasDirection = -1
      kBiasDirection =  1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )

      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   
   end if
   
   ! VERTEX 4
   
   if ( myback  == mpi_proc_null .and. & 
        myright == mpi_proc_null .and. &
        mydown  == mpi_proc_null          )  then
   
      i = ista ! il + igp
      j = jend ! ju - jgp
      k = ksta ! kl + kgp
   
      iBiasDirection =  1
      jBiasDirection = -1
      kBiasDirection =  1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   end if
   
   
   ! VERTEX 5
   
   if ( myback == mpi_proc_null .and. &
        myleft == mpi_proc_null .and. &
        myup   == mpi_proc_null          )  then
   
      i = ista ! il + igp
      j = jsta ! jl + jgp
      k = kend ! ku - kgp
   
      iBiasDirection =  1
      jBiasDirection =  1
      kBiasDirection = -1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
      
   
   end if
   
   ! VERTEX 6
   
   if ( myfront == mpi_proc_null .and. &
        myleft  == mpi_proc_null .and. &
        myup    == mpi_proc_null          )  then
      
      i = iend ! iu - igp
      j = jsta ! jl + jgp
      k = kend ! ku - kgp
   
      iBiasDirection = -1
      jBiasDirection =  1
      kBiasDirection = -1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   
   end if
   
   ! VERTEX 7
   
   if ( myfront == mpi_proc_null .and. &
        myright == mpi_proc_null .and. &
        myup    == mpi_proc_null          )  then
   
      i = iend ! iu - igp
      j = jend ! ju - jgp
      k = kend ! ku - kgp
      
      iBiasDirection = -1
      jBiasDirection = -1
      kBiasDirection = -1
   
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   
   end if
   
   ! VERTEX 8
   
   if ( myback  == mpi_proc_null .and. &
        myright == mpi_proc_null .and. &
        myup    == mpi_proc_null          )  then
   
      i = ista ! il + igp
      j = jend ! ju - jgp
      k = kend ! ku - kgp
   
      iBiasDirection =  1
      jBiasDirection = -1
      kBiasDirection = -1
      
      !--------------------------------------------------------------------
      ! ξ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dcsi = dc * BiasedDerivative( phi( i + iBiasDirection * 0 , j , k ) , &
                                         phi( i + iBiasDirection * 1 , j , k ) , &
                                         phi( i + iBiasDirection * 2 , j , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         iBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! η - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_deta = de * BiasedDerivative( phi( i , j + jBiasDirection * 0 , k ) , &
                                         phi( i , j + jBiasDirection * 1 , k ) , &
                                         phi( i , j + jBiasDirection * 2 , k ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         jBiasDirection                          &
                                       )
   
      !--------------------------------------------------------------------
      ! ζ - DIRECTION
      !--------------------------------------------------------------------
   
      dphi_dzet = dz * BiasedDerivative( phi( i , j , k + kBiasDirection * 0 ) , &
                                         phi( i , j , k + kBiasDirection * 1 ) , &
                                         phi( i , j , k + kBiasDirection * 2 ) , &
                                         dummy                                 , &
                                         2                                     , &
                                         kBiasDirection                          &
                                       )
   
      ! Calculo del gradiente en coordenadas curvilineas
   
      dphi_dx = csi( 1 , i , j , k ) * dphi_dcsi + &
                eta( 1 , i , j , k ) * dphi_deta + &
                zet( 1 , i , j , k ) * dphi_dzet
   
      dphi_dy = csi( 2 , i , j , k ) * dphi_dcsi + &
                eta( 2 , i , j , k ) * dphi_deta + &
                zet( 2 , i , j , k ) * dphi_dzet
   
      dphi_dz = csi( 3 , i , j , k ) * dphi_dcsi + &
                eta( 3 , i , j , k ) * dphi_deta + &
                zet( 3 , i , j , k ) * dphi_dzet
   
      phi_gradient(1,i,j,k) = dphi_dx
      phi_gradient(2,i,j,k) = dphi_dy
      phi_gradient(3,i,j,k) = dphi_dz
   
   end if
   
   ! Update the ghost nodes
   call rhs_exchng3_4d( phi_gradient )

   contains

   include 'rhs_exchng3_4d.F90'

end subroutine calc_phi_gradient