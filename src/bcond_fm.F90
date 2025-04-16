
subroutine bcond_fm ( il    ,  iu        ,  &       
                      jl    ,  ju        ,  &       
                      kl    ,  ku        ,  &       
                      igp   ,               &
                      jgp   ,               & 
                      kgp   ,               & 
                      dc    ,  de  ,  dz ,  &
                      q     ,               & 
                      csi   ,               &                                   
                      eta   ,               &                                   
                      zet   ,               &
                      aj    ,               &                                   
                      x     ,  y  ,   z  ,  &                                   
                      xnut  ,               &          
                      rsign              ,  &
                      phi                   &
                     )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! general version !Kim's curved duct (full 3d geometry)
   ! 
   ! fine-mesh
   ! boundary conditions
   ! 
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   use global_app
   use global_mpi
   use global_param

   implicit none
  
   ! extents
   ! 
   integer :: il
   integer :: jl
   integer :: kl
   integer :: iu
   integer :: ju
   integer :: ku

   ! ghost points
   ! 
   integer :: igp
   integer :: jgp
   integer :: kgp

   ! solution vector
   ! 
   real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: q
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(inout) ::  xnut
   real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi,eta,zet
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)    , intent(in) :: aj
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)    , intent(in) :: x,y,z
   real (kind = rdf) :: dc , de , dz

   ! From bcond_lsm
   real (kind = rdf), parameter :: itp1stOrd_a = one , itp1stOrd_b = zero
   real (kind = rdf), parameter :: itp2ndOrd_a = four/three, itp2ndOrd_b = -one/three
   real (kind = rdf) :: itpa , itpb ! first and second nodes
   real (kind = rdf) :: exsign

   ! solid wall flag
   logical, parameter :: solid_wall_blanking = .true.

   ! flag variable to identify air or water phase
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) ::  rsign

   ! flag variable to identify air or water phase
   real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) ::  phi

   ! fixed pressure
   ! 
   real (kind = rdf) :: pfix

   real (kind = rdf) :: paux_extp, paux_vert, phase_extp, phase_vert
   real (kind = rdf) :: coeff1, coeff2
   real (kind = rdf) :: paux

   integer :: i, j, k, b , n
  
   integer :: i_mysta
   integer :: j_mysta
   integer :: k_mysta

   integer :: i_myend
   integer :: j_myend
   integer :: k_myend


   ! pressure extrapolation variables:

   real (kind = rdf) :: c1 , c2 , c3      ! coeff of each component of the extrapolation
   real (kind = rdf) :: g0 , g1           ! g22 at the wall and 1 node away from it
   real (kind = rdf) :: ret0 , ret1       ! (1/Re + xnut) at the wall and 1 node away from it
   real (kind = rdf) :: f1aux , f2aux     ! internal derivative of the viscous flux
   real (kind = rdf) :: fv1star , fv2star ! external derivative of the viscous flux
   real (kind = rdf) :: ds1               ! distance to the first node



   ! make definition of i_mysta etc., consistent with
   ! all other uses, e.g., solver_daf & rhs_kw_bcond
   !
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! processes on the domain boundaries
   ! 
   !if (myback == mpi_proc_null)  i_mysta = il + igp + 1
   !if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
   !if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

   !if (myfront == mpi_proc_null) i_myend = iu - igp - 1
   !if (myright == mpi_proc_null) j_myend = ju - jgp - 1
   !if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

   !------------------------------------------------------------------------------------------
   !
   !                                   ξ - DIRECTION
   !
   !------------------------------------------------------------------------------------------
   !
   !              INITIAL boundary in ξ-direction (i = 1, b = 1)
   !
   !------------------------------------------------------------------------------------------

   coeff1 = two
   coeff2 = -one

   if (myback == mpi_proc_null) then
      b = 1
      i = i_mysta

      csi1: select case( btype( b , myzone ) )
      
      case(-1) ! -> dummy option to skip bcs

      case(0) ! -> interface

      case(1:3) ! wall, symmetric plane & freestream
      
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
         
            !if( rsign(i,j,k) > one_half ) then

               q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i+1,j,k) + &
                              sb(1:4,b) * q(1:4,i+2,j,k)
            
               ! Pressure extrapolation to the wall
               
               ! q(1,i,j,k) = two * q(1,i+1,j,k) - q(1,i+2,j,k)

            !end if
         
         end do
         end do
     
      case(4) ! inflow
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
            
            !if( rsign(i,j,k)  > one_half ) then
   
               ! I commented this for debbuging the straight channel case
               ! uncomment if corresponds
               q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i+1,j,k) - &
                                   rat(j,k,b)  * q(1,i+2,j,k) 
           
               ! I commented this for the bump case. This is to
               ! keep the initial hydrostatic pressure distribution at
               ! the entrance and the constant inflow .Uncomment in any other 
               ! cases

               ! q(1:4,i,j,k) = two * q(1:4,i+1,j,k) - q(1:4,i+2,j,k)

               !q(1,i,j,k) = sa(1,b) * q(1,i+1,j,k) + &
               !             sb(1,b) * q(1,i+2,j,k)
      
               !q(2:4,i,j,k) = q(2:4,i_myend,j,k) 
            
            !end if
         
         end do
         end do
     
      case(5) ! exit (MOC)
      case(6) ! periodic condition
        !do k = k_mysta, k_myend
        !do j = j_mysta, j_myend
        !   q(1:4,i,j,k) = q(1:4,i_myend,j,k)
        !end do
        !end do
      case(9) !slip condition
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

            !if( rsign(i,j,k) > one_half ) then

               ! q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i+1,j,k) + &
               !               sb(1:4,b) * q(1:4,i+2,j,k)
               !
               
               !q(1:4,i,j,k) = two * q(1:4,i+1,j,k) - &
               !               one * q(1:4,i+2,j,k)

               !q(2:4,i,j,k) = (one + rat(j,k,b)) * q(2:4,i+1,j,k) - &
               !                      rat(j,k,b)  * q(2:4,i+2,j,k)    

               phase_extp =   rsign(i+1,j,k) * rsign(i+2,j,k) 

               phase_vert =   rsign(i,j,max( k-1 , k_mysta )) & 
                            * rsign(i,j,max( k-2 , k_mysta )) 

               ! Slip-walls --> zero shear --> ∂ / ∂ξ = 0
               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i+1,j,k) + &
                              sb(2:4,b) * q(2:4,i+2,j,k)    

               !paux_extp = (one + rat(j,k,b)) * q(1,i+1,j,k) - &
               !                   rat(j,k,b)  * q(1,i+2,j,k)    

               ! ∂p / ∂ξ = 0 at the boundary
               paux_extp = sa(1,b) * q(1,i+1,j,k) + &
                           sb(1,b) * q(1,i+2,j,k)    


               paux_vert =   two * q(1,i,j,max( k-1 , k_mysta ))  &
                           - one * q(1,i,j,max( k-2 , k_mysta ))  

               ! If the extrapolation and the vertical stencil have both
               ! air nodes (rsign = 0), then this formulation leads to p=0
               paux =                         phase_extp * paux_extp &
                       + (one - phase_extp) * phase_vert * paux_vert

               ! If the node is in the air-phase, I don't change it because
               ! it was extrapolated in the pressure_extrapolation()
               ! routine
               q(1,i,j,k) = ( one - rsign(i,j,k) ) * q(1,i,j,k) + &
                                    rsign(i,j,k)   * paux

               ! non-penetration BC
               q(2,i,j,k) = zero

            !end if

         end do
         end do          
      
      end select csi1

   end if
   
   !------------------------------------------------------------------------------------------
   !
   !              FINAL boundary in ξ-direction (i = im, b = 2)
   !
   !------------------------------------------------------------------------------------------

   if (myfront == mpi_proc_null) then
      
      b = 2
      i = i_myend

      csi2: select case( btype( b , myzone ) )

      case(-1) ! -> dummy option to skip bcs
      
      case(0) ! -> interface
      
      case(1:3) ! wall, symmetric plane & freestream
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
      
            !if( rsign(i,j,k)> one_half ) then

               q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i-1,j,k) + &
                              sb(1:4,b) * q(1:4,i-2,j,k)
            
               ! q(1,i,j,k) = two * q(1,i-1,j,k) - q(1,i-2,j,k)

           ! end if

         end do
         end do
      
      case(4) ! inflow
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
               
            if( rsign(i,j,k) > one_half ) then
               
               !q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i-1,j,k) - &
               !                    rat(j,k,b)  * q(1,i-2,j,k)
               !q(1,i,j,k) = sa(1,b) * q(1,i-1,j,k) + &
               !             sb(1,b) * q(1,i-2,j,k)
            
            end if

         end do
         end do
      
      case(5) ! exit (MOC)
         !j = j_mysta
         !do k = k_mysta, k_myend
         !   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
         !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
         !                  sb(2:4,b) * q(2:4,i-2,j,k)
         !end do

         !j = j_myend
         !do k = k_mysta, k_myend
         !   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
         !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
         !                  sb(2:4,b) * q(2:4,i-2,j,k)
         !end do

      case(6) ! periodic condition

         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

            if( rsign(i,j,k) > one_half ) then

               q(1:4,i,j,k) = q(1:4,i_mysta,j,k)
            
            end if
         
         end do
         end do
      
      case(8) ! free surface / outflow
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
            
            !if( rsign(i,j,k)> one_half ) then
               
               !q(1:4,i,j,k) = q(1:4,i-1,j,k) !sa(1:4,b)*q(1:4,i-1,j,k) + sb(1:4,b)*q(1:4,i-2,j,k)               
               
               q(1,i,j,k) = two * q(1,i-1,j,k) - q(1,i-2,j,k)

               ! I commented this for the bump case. This is to
               ! keep the initial hydrostatic pressure distribution at
               ! the outlet .Uncomment in other cases

               q(2:4,i,j,k) = two * q(2:4,i-1,j,k) - q(2:4,i-2,j,k)
               
            !end if
            
         end do
         end do

      case(9) !slip condition   
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
            
            !if( rsign(i,j,k) > one_half ) then

               !q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i-1,j,k) + &
               !               sb(1:4,b) * q(1:4,i-2,j,k)
               
               !q(2:4,i,j,k) = (one + rat(j,k,b)) * q(2:4,i-1,j,k) - &
               !                      rat(j,k,b)  * q(2:4,i-2,j,k)    
               phase_extp =   rsign(i-1,j,k) * rsign(i-2,j,k) 

               phase_vert =   rsign(i,j,max( k-1 , k_mysta )) & 
                            * rsign(i,j,max( k-2 , k_mysta )) 


               ! Slip-walls --> zero shear --> ∂ / ∂ξ = 0
               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
                              sb(2:4,b) * q(2:4,i-2,j,k)    

               !paux_extp = (one + rat(j,k,b)) * q(1,i-1,j,k) - &
               !                   rat(j,k,b)  * q(1,i-2,j,k)    

               ! ∂p / ∂ξ = 0 at the boundary
               paux_extp = sa(1,b) * q(1,i-1,j,k) + &
                           sb(1,b) * q(1,i-2,j,k)    

               paux_vert =   two * q(1,i,j,max( k-1 , k_mysta ))  &
                           - one * q(1,i,j,max( k-2 , k_mysta ))  

               ! If the extrapolation and the vertical stencil have both
               ! air nodes (rsign = 0), then this formulation leads to p=0
               paux =                         phase_extp * paux_extp &
                       + (one - phase_extp) * phase_vert * paux_vert

               ! If the node is in the air-phase, I don't change it because
               ! it was extrapolated in the pressure_extrapolation()
               ! routine
               q(1,i,j,k) = ( one - rsign(i,j,k) ) * q(1,i,j,k) + &
                                    rsign(i,j,k)   * paux

               ! non-penetration BC
               q(2,i,j,k) = zero
            
            !end if
            
         end do
         end do   

      end select csi2
   
   end if ! if(myfront == mpi_proc_null)
     
   !------------------------------------------------------------------------------------------
   !
   !                                   η - DIRECTION
   !
   !------------------------------------------------------------------------------------------
   !
   !            INITIAL boundary in η-direction (j = 1, b = 3)
   !
   !------------------------------------------------------------------------------------------

   if (myleft == mpi_proc_null) then
      
      b = 3
      j = j_mysta
      c1 = four / three
      c2 = -one / three
      c3 = -two / three

      eta1: select case( btype( b , myzone ) )

      case(-1) ! -> dummy option to skip bcs

      case(0) ! -> interface
      
      case(1:3) ! wall, symmetric plane & freestream
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend
            
            !if( rsign(i,j,k) > one_half ) then      
               
               !call p_no_slip_wall_bcond(b)

               ! g0 = g22 at the wall (0 nodes away from the wall)
               g0 = eta(1,i,j,k)**two + eta(2,i,j,k)**two 

               ! g1 = g22 1 node away from the wall
               g1 = eta(1,i,j+1,k)**two + eta(2,i,j+1,k)**two 

               ! ret0 = ( 1/Re + xnut ) at the wall
               ret0 = one / ren + xnut(i,j,k) 

               ! ret1 = ( 1/Re + xnut ) 1 node away from the wall
               ret1 = one / ren + xnut(i,j+1,k) 

               ! u-derivatives
               f1aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * ( four * q(2,i,j+1,k) - q(2,i,j+2,k) ) 
               f2aux = de * g1 / ( two * aj(i,j+1,k) ) * ret1 * q(2,i,j+2,k) 

               fv1star = eta(1,i,j,k) / sqrt( g0 ) * aj(i,j,k) * ( f2aux - f1aux )

               ! v-derivatives
               f1aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * ( four * q(3,i,j+1,k) - q(3,i,j+2,k) ) 
               f2aux = de * g1 / ( two * aj(i,j+1,k) ) * ret1 * q(3,i,j+2,k) 

               fv2star = eta(2,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )
              
 
               ! distance to the first node
               ds1 = norm2( (/x(i,j+1,k)-x(i,j,k) , y(i,j+1,k)-y(i,j,k) , z(i,j+1,k)-z(i,j,k)/) )


               q(1,i,j,k) = c1 * q(1,i,j+1,k) + c2 * q(1,i,j+2,k) + ds1 * c3 * ( fv1star + fv2star )


               !q(1,i,j,k) = coeff1 * q(1,i,j+1,k) + &
               !             coeff2 * q(1,i,j+2,k)
               
               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i,j+1,k) + &
                              sb(2:4,b) * q(2:4,i,j+2,k)

               ! q(1,i,j,k) = two * q(1,i,j+1,k) - q(1,i,j+2,k)
            
            !end if
            
         end do
         end do

      case(4) ! inflow
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend
            
            !if( rsign(i,j,k) > one_half ) then

               q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j+1,k) - &
                                   rat(i,k,b)  * q(1,i,j+2,k)
               
               !q(1,i,j,k) = sa(1,b) * q(1,i,j+1,k) + &
               !             sb(1,b) * q(1,i,j+2,k)
            
            !end if
         
         end do
         end do

      case(5) ! exit (MOC)
      
      case(6) ! periodic condition
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend

            !if( rsign(i,j,k) > one_half ) then

               q(1:4,i,j,k) = q(1:4,i,j_myend,k)
            
            !end if

         end do
         end do

      case(9) !slip condition
         
         do k = k_mysta, k_myend
         do i = i_mysta, i_myend
            
            !if( rsign(i,j,k) > one_half ) then
            
               !q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j+1,k) + &
               !               sb(1:4,b) * q(1:4,i,j+2,k)
               
               !q(2:4,i,j,k) = (one + rat(j,k,b)) * q(2:4,i,j+1,k) - &
               !                      rat(j,k,b)  * q(2:4,i,j+2,k)    

               phase_extp =   rsign(i,j+1,k) * rsign(i,j+2,k) 

               phase_vert =   rsign(i,j,max( k-1 , k_mysta )) & 
                            * rsign(i,j,max( k-2 , k_mysta )) 

               ! Slip-walls --> zero shear --> ∂ / ∂η = 0
               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i,j+1,k) + &
                              sb(2:4,b) * q(2:4,i,j+2,k)    

               !paux_extp = (one + rat(j,k,b)) * q(1,i,j+1,k) - &
               !                   rat(j,k,b)  * q(1,i,j+2,k)    

               ! ∂p / ∂η = 0 at the boundary
               paux_extp = sa(1,b) * q(1,i,j+1,k) + &
                           sb(1,b) * q(1,i,j+2,k)    

               paux_vert =   two * q(1,i,j,max( k-1 , k_mysta ))  &
                           - one * q(1,i,j,max( k-2 , k_mysta ))  

               ! If the extrapolation and the vertical stencil have both
               ! air nodes (rsign = 0), then this formulation leads to p=0
               paux =                         phase_extp * paux_extp &
                       + (one - phase_extp) * phase_vert * paux_vert

               ! If the node is in the air-phase, I don't change it because
               ! it was extrapolated in the pressure_extrapolation()
               ! routine
               q(1,i,j,k) = ( one - rsign(i,j,k) ) * q(1,i,j,k) + &
                                    rsign(i,j,k)   * paux

               ! non-penetration BC
               q(3,i,j,k) = zero
            
            !end if
         
         end do
         end do
      
      end select eta1   

   end if

   !------------------------------------------------------------------------------------------
   !
   !            FINAL boundary in η-direction (j = jm, b = 4)
   !
   !------------------------------------------------------------------------------------------

   if (myright == mpi_proc_null) then

      b = 4
      j = j_myend
      c1 = four / three
      c2 = -one / three
      c3 =  two / three

      eta2: select case( btype( b , myzone ) )
 
      case(-1) ! -> dummy option to skip bcs
     
      case(0) ! -> interface
      
      case(1:3) ! wall, symmetric plane & freestream
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend
            
            !if( rsign(i,j,k) > one_half ) then
  
               !call p_no_slip_wall_bcond(b)
               ! g0 = g22 at the wall (0 nodes away from the wall)
               g0 = eta(1,i,j,k)**two + eta(2,i,j,k)**two 

               ! g1 = g22 1 node away from the wall
               g1 = eta(1,i,j-1,k)**two + eta(2,i,j-1,k)**two 

               ! ret0 = ( 1/Re + xnut ) at the wall
               ret0 = one / ren + xnut(i,j,k) 

               ! ret1 = ( 1/Re + xnut ) 1 node away from the wall
               ret1 = one / ren + xnut(i,j-1,k) 

               ! u-derivatives
               f1aux = de * g1 / ( two * aj(i,j-1,k) ) * ret1 * (-q(2,i,j-2,k) ) 
               f2aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * (-four * q(2,i,j-1,k) + q(2,i,j-2,k) ) 

               fv1star = eta(1,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

               ! v-derivatives
               f1aux = de * g1 / ( two * aj(i,j-1,k) ) * ret1 * (-q(3,i,j-2,k) ) 
               f2aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * (-four * q(3,i,j-1,k) + q(3,i,j-2,k) ) 

               fv2star = eta(2,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

               ! distance to the first node
               ds1 = norm2( (/x(i,j,k)-x(i,j-1,k) , y(i,j,k)-y(i,j-1,k) , z(i,j,k)-z(i,j-1,k)/) )

               ! pressure extrapolation
               q(1,i,j,k) = c1 * q(1,i,j-1,k) + c2 * q(1,i,j-2,k) + ds1 * c3 * ( fv1star + fv2star )

               !q(1,i,j,k) = coeff1 * q(1,i,j-1,k) + &
               !             coeff2 * q(1,i,j-2,k)
               
               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i,j-1,k) + &
                              sb(2:4,b) * q(2:4,i,j-2,k)

               ! q(1,i,j,k) = two * q(1,i,j-1,k) - q(1,i,j-2,k)
            
            !end if
         
         end do
         end do
      
      case(4) ! inflow
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend

            if( rsign(i,j,k) > one_half ) then            
               
               q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j-1,k) - &
                                   rat(i,k,b)  * q(1,i,j-2,k)
               
               !q(1,i,j,k) = sa(1,b) * q(1,i,j-1,k) + &
               !             sb(1,b) * q(1,i,j-2,k)
            
            end if

         end do
         end do
      
      case(5) ! exit (MOC)
      
      case(6) ! periodic condition
         
         do k = k_mysta , k_myend
         do i = i_mysta , i_myend

            if( rsign(i,j,k) > one_half ) then
               
               q(1:4,i,j,k) = q(1:4,i,j_mysta,k)
            
            end if
            
         end do
         end do
      
      case(9)
         
         do k = k_mysta, k_myend
         do i = i_mysta, i_myend
            
            !if( rsign(i,j,k) > one_half ) then
           
               !q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j-1,k) + &
               !               sb(1:4,b) * q(1:4,i,j-2,k)

               !q(2:4,i,j,k) = (one + rat(j,k,b)) * q(2:4,i,j-1,k) - &
               !                      rat(j,k,b)  * q(2:4,i,j-2,k)    

               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i,j-1,k) + &
                              sb(2:4,b) * q(2:4,i,j-2,k)    

               phase_extp =   rsign(i,j-1,k) * rsign(i,j-2,k) 

               phase_vert =   rsign(i,j,max( k-1 , k_mysta )) & 
                            * rsign(i,j,max( k-2 , k_mysta )) 


               !paux_extp = (one + rat(j,k,b)) * q(1,i,j-1,k) - &
               !                   rat(j,k,b)  * q(1,i,j-2,k)    

               ! ∂p / ∂η = 0 at the boundary
               paux_extp = sa(1,b) * q(1,i,j-1,k) + &
                           sb(1,b) * q(1,i,j-2,k)    

               paux_vert =   two * q(1,i,j,max( k-1 , k_mysta ))  &
                           - one * q(1,i,j,max( k-2 , k_mysta ))  

               ! If the extrapolation and the vertical stencil have both
               ! air nodes (rsign = 0), then this formulation leads to p=0
               paux =                         phase_extp * paux_extp &
                       + (one - phase_extp) * phase_vert * paux_vert

               ! If the node is in the air-phase, I don't change it because
               ! it was extrapolated in the pressure_extrapolation()
               ! routine
               q(1,i,j,k) = ( one - rsign(i,j,k) ) * q(1,i,j,k) + &
                                    rsign(i,j,k)   * paux

               ! non-penetration BC
               q(3,i,j,k) = zero
            
            !end if
         
         end do
         end do
      
      end select eta2

   end if

   !------------------------------------------------------------------------------------------
   !
   !                                   ζ - DIRECTION
   !
   !------------------------------------------------------------------------------------------
   !
   !            INITIAL boundary in ζ-direction (k = 1, b = 5)
   !
   !------------------------------------------------------------------------------------------

   if (mydown == mpi_proc_null) then
      b = 5
      k = k_mysta

      zet1: select case( btype( b , myzone ) )

      case(-1) ! -> dummy option to skip bcs
      
      case(0) ! -> interface
      
      case(1:3) ! wall, symmetric plane & freestream
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend
               
            !if( rsign(i,j,k) > one_half ) then

               q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k+1) + &
                              sb(1:4,b) * q(1:4,i,j,k+2)

               !q(1,i,j,k) = two * q(1,i,j,k+1) - q(1,i,j,k+2)
            
            !end if
         
         end do
         end do
      
      case(4) ! inflow
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend

            if( rsign(i,j,k) > one_half ) then
            
               q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k+1) - &
                                   rat(i,j,b)  * q(1,i,j,k+2)
            
               !q(1,i,j,k) = sa(1,b) * q(1,i,j,k+1) + &
               !             sb(1,b) * q(1,i,j,k+2)
            
            end if

         end do
         end do
      
      case(5) ! exit (MOC)
      
      case(6) ! periodic condition
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend
               
            if( rsign(i,j,k) > one_half ) then
               
               q(1:4,i,j,k) = q(1:4,i,j,k_myend)
            
            end if
         
         end do
         end do
      
      case(9) !slip condition
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend

            !if( rsign(i,j,k) > one_half ) then
               
               !q(2:4,i,j,k) = (one + rat(j,k,b)) * q(2:4,i,j,k+1) - &
               !                      rat(j,k,b)  * q(2:4,i,j,k+2)    

               q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i,j,k+1) + &
                              sb(2:4,b) * q(2:4,i,j,k+2)    

               phase_extp =   rsign(i,j,k+1) * rsign(i,j,k+2) 

               phase_vert =   rsign(i,j,max( k-1 , k_mysta )) & 
                            * rsign(i,j,max( k-2 , k_mysta )) 


               !paux_extp = (one + rat(j,k,b)) * q(1,i,j,k+1) - &
               !                   rat(j,k,b)  * q(1,i,j,k+2)    

               ! ∂p / ∂ζ = 0 at the boundary
               paux_extp = sa(1,b) * q(1,i,j,k+1) + &
                           sb(1,b) * q(1,i,j,k+2)    

               paux_vert =   two * q(1,i,j,max( k+1 , k_mysta ))  &
                           - one * q(1,i,j,max( k+2 , k_mysta ))  

               ! If the extrapolation and the vertical stencil have both
               ! air nodes (rsign = 0), then this formulation leads to p=0
               paux =                         phase_extp * paux_extp &
                       + (one - phase_extp) * phase_vert * paux_vert

               ! If the node is in the air-phase, I don't change it because
               ! it was extrapolated in the pressure_extrapolation()
               ! routine
               q(1,i,j,k) = ( one - rsign(i,j,k) ) * q(1,i,j,k) + &
                                    rsign(i,j,k)   * paux

               ! non-penetration BC
               !q(3,i,j,k) = zero
               
               q(4,i,j,k) = zero
            
            !end if
   
         end do
         end do
      end select zet1

   end if

   !------------------------------------------------------------------------------------------
   !
   !             FINAL boundary in ζ-direction (k = km, b = 6)
   !
   !------------------------------------------------------------------------------------------

   if (myup == mpi_proc_null) then
      
      b = 6
      k = k_myend

      zet2: select case( btype( b , myzone ) )

      case(-1) ! -> dummy option to skip bcs
      
      case(0) ! -> interface
      
      case(1:3) ! wall, symmetric plane & freestream
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend
            
            !if( rsign(i,j,k) > one_half ) then
            
               q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k-1) + &
                              sb(1:4,b) * q(1:4,i,j,k-2)
            
               !q(1,i,j,k) = two * q(1,i,j,k-1) - q(1,i,j,k-2)

            !end if
         
         end do
         end do
      
      case(4) ! inflow
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend
               
            if( rsign(i,j,k) > one_half ) then
               
               q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k-1) - &
                                   rat(i,j,b)  * q(1,i,j,k-2)

               !q(1,i,j,k) = sa(1,b) * q(1,i,j,k-1) + &
               !             sb(1,b) * q(1,i,j,k-2)
            
            end if
         
         end do
         end do
     
     case(5) ! exit (MOC)
     
     case(6) ! periodic condition
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend
               
            if( rsign(i,j,k) > one_half ) then
               
               q(1:4,i,j,k) = q(1:4,i,j,k_mysta)
            
            end if

         end do
         end do

      case(9) !slip condition
         
         do j = j_mysta , j_myend
         do i = i_mysta , i_myend

            !if( rsign(i,j,k) > one_half ) then
               
               q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k-1) + &
                              sb(1:4,b) * q(1:4,i,j,k-2)

               !q(1:4,i,j,k) = two * q(1:4,i,j,k-1) - &
               !               one * q(1:4,i,j,k-2)
               
               !q(1,i,j,k) =    three * q(1,i,j,k-1)  &
               !              - three * q(1,i,j,k-2)  &
               !              + one   * q(1,i,j,k-3) 


               q(4,i,j,k) = zero
            
            !end if
         
         end do
         end do 
      end select zet2
   
   end if


   ! For the moment, I'm not gonna remove the pressure reference, because the
   ! code is meant to be working with the absolute non-dimensional pressure.
   ! Some correction to this approach may need to be included as the hydrostatic
   ! pressure is explicitly incorporated in the initial condition in the case
   ! of the sloshin tank 

   ! probe default pressure value
   !
   !if (myid == pfix_proc) pfix = q( 1 , local_ifix , local_jfix , local_kfix )

   ! distribute fixed pressure to all processes
   ! 
   !call mpi_bcast( pfix, 1, MPI_REAL_TYPE, pfix_proc, mpi_comm_world, ierr )

   !do k = k_mysta , k_myend
   !do j = j_mysta , j_myend
   !do i = i_mysta , i_myend

   !   if( rsign(i,j,k) > one_half ) then
         
   !      q(1,i,j,k) = q(1,i,j,k) - pfix
      
   !   end if

   !end do
   !end do
   !end do

   !contains

   !include 'p_no_slip_wall_bcond.F90'

   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   !
   !     ######  #          #    #     # #    # ### #     #  #####  
   !     #     # #         # #   ##    # #   #   #  ##    # #     # 
   !     #     # #        #   #  # #   # #  #    #  # #   # #       
   !     ######  #       #     # #  #  # ###     #  #  #  # #  #### 
   !     #     # #       ####### #   # # #  #    #  #   # # #     # 
   !     #     # #       #     # #    ## #   #   #  #    ## #     # 
   !     ######  ####### #     # #     # #    # ### #     #  #####  
   !
   ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   n = 1

   if (nblk /= 0) then

      do nb = 1, nblk
         
         ! Let's zero all the variables in the interior region of the blanking zone
         do k = max( k_mysta , li_blk_ka(n,nb)+1 ) , min(k_myend , li_blk_kb(n,nb)-1 )
         do j = max( j_mysta , li_blk_ja(n,nb)+1 ) , min(j_myend , li_blk_jb(n,nb)-1 )
         do i = max( i_mysta , li_blk_ia(n,nb)+1 ) , min(i_myend , li_blk_ib(n,nb)-1 )

            q(1:4,i,j,k) = zero
            xnut(i,j,k)  = zero

         end do
         end do
         end do

         i = li_blk_ia(n,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 ) then
               
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
            do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j,k) )- five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i-1,j,k) + itpb * q(1,i-2,j,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i-1,j,k) - q(2:4,i-2,j,k)    
               xnut(i,j,k)  = two * xnut(i-1,j,k)  - xnut(i-2,j,k)    
               ! No penetration condition
               q(2,i,j,k)   = zero    
            
               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do
            end do

         end if
      
         i = li_blk_ib(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 ) then
               
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
            do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i+1,j,k) + itpb * q(1,i+2,j,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i+1,j,k) - q(2:4,i+2,j,k)    
               xnut(i,j,k)  = two * xnut(i+1,j,k)  - xnut(i+2,j,k)    
               ! No penetration condition
               q(2,i,j,k)   = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero
            
            end do
            end do

         end if

         j = li_blk_ja(n,nb)

         if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then
               
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
            do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j-2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i,j-1,k) + itpb * q(1,i,j-2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i,j-1,k) - q(2:4,i,j-2,k)    
               xnut(i,j,k)  = two * xnut(i,j-1,k)  - xnut(i,j-2,k)    
               ! No penetration condition
               q(3,i,j,k)   = zero    
            
               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do
            end do

         end if
      
         j = li_blk_jb(n,nb)

         if ( blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then
               
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
            do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j+2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i,j+1,k) + itpb * q(1,i,j+2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i,j+1,k) - q(2:4,i,j+2,k)    
               xnut(i,j,k)  = two * xnut(i,j+1,k)  - xnut(i,j+2,k)    
               ! No penetration condition
               q(3,i,j,k)   = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero
            
            end do
            end do

         end if

         !------------------------------------------------------------------------------
         ! Corners
         !------------------------------------------------------------------------------

         i = li_blk_ia(n,nb)
         j = li_blk_ja(n,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j-2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i-1,j-1,k) + itpb * q(1,i-2,j-2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i-1,j-1,k) - q(2:4,i-2,j-2,k)    
               xnut(i,j,k)  = two * xnut(i-1,j-1,k)  - xnut(i-2,j-2,k)    
               
               ! No penetration condition
               q(2,i,j,k) = zero    
               q(3,i,j,k) = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do

         end if


         i = li_blk_ib(n,nb)
         j = li_blk_ja(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j-2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i+1,j-1,k) + itpb * q(1,i+2,j-2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i+1,j-1,k) - q(2:4,i+2,j-2,k)    
               xnut(i,j,k)  = two * xnut(i+1,j-1,k)  - xnut(i+2,j-2,k)    
               
               ! No penetration condition
               q(2,i,j,k) = zero    
               q(3,i,j,k) = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do

         end if

         i = li_blk_ib(n,nb)
         j = li_blk_jb(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j+2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i+1,j+1,k) + itpb * q(1,i+2,j+2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i+1,j+1,k) - q(2:4,i+2,j+2,k)    
               xnut(i,j,k)  = two * xnut(i+1,j+1,k)  - xnut(i+2,j+2,k)    
               
               ! No penetration condition
               q(2,i,j,k) = zero    
               q(3,i,j,k) = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do

         end if


         i = li_blk_ia(n,nb)
         j = li_blk_jb(n,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j+2,k) )-five ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               ! Zero pressure gradient
               q(1,i,j,k)   = itpa * q(1,i-1,j+1,k) + itpb * q(1,i-2,j+2,k)    

               ! Extrapolated velocities and xnut (considering uniform grid)
               q(2:4,i,j,k) = two * q(2:4,i-1,j+1,k) - q(2:4,i-2,j+2,k)    
               xnut(i,j,k)  = two * xnut(i-1,j+1,k)  - xnut(i-2,j+2,k)    
               
               ! No penetration condition
               q(2,i,j,k) = zero    
               q(3,i,j,k) = zero    

               if ( non_slip_wall_blanking ) q(2:4,i,j,k) = zero

            end do

         end if

         !k = li_blk_ka(n,nb)
         !
         !if ( blktype(5,nb,myzone) == 0 .and. k > k_mysta+1 ) then
         !      
         !   do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 
         !   do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 
         !
         !      ! Zero pressure gradient
         !      q(1,i,j,k)   = four/three * q(1,i,j,k-1) - one/three * q(1,i,j,k-2)    
         !
         !
         !      ! Extrapolated velocities and xnut (considering uniform grid)
         !      q(2:4,i,j,k) = two * q(2:4,i,j,k-1) - q(2:4,i,j,k-2)    
         !      xnut(i,j,k)  = two * xnut(i,j,k-1)  - xnut(i,j,k-2)    
         !      ! No penetration condition
         !      q(4,i,j,k)   = zero    
         !   
         !   end do
         !   end do
         !
         !end if
         !
         !k = li_blk_kb(n,nb)
         !
         !if ( blktype(6,nb,myzone) == 0 .and. k < k_myend-1 ) then
         !      
         !   do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 
         !   do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 
         !
         !      ! Zero pressure gradient
         !      q(1,i,j,k)   = four/three * q(1,i,j,k+1) - one/three * q(1,i,j,k+2)    
         !
         !
         !      ! Extrapolated velocities and xnut (considering uniform grid)
         !      q(2:4,i,j,k) = two * q(2:4,i,j,k+1) - q(2:4,i,j,k+2)    
         !      xnut(i,j,k)  = two * xnut(i,j,k+1)  - xnut(i,j,k+2)    
         !      ! No penetration condition
         !      q(4,i,j,k)   = zero    
         !   
         !   end do
         !   end do
         !
         !end if

      end do

   end if




end subroutine bcond_fm

