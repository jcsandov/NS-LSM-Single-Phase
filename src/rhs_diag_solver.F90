
subroutine rhs_diag_solver ()

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! implicit, diagonal approximate factorization solver
   !
   ! input
   !     rh(4,ijk)
   !     dtau(ijk)
   !     csi(ijk), eta(ijk), zet(ijk)
   !     aj(ijk)
   !     dc, de, dz
   !     dc2, de2, dz2
   !
   ! output
   !     rh(4,ijk)
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   implicit none
   
   integer :: info
   integer :: ista , iend , jsta , jend , ksta , kend

   ! dimension lhs and temporary rhs
   real (kind = rdf), dimension(:,:), allocatable :: aw, ap, ae, aq

   ! boundary condition array
   real (kind = rdf), dimension(:,:,:,:), allocatable :: brh

   ! local dummy variables
   real (kind = rdf), dimension(4) :: lambda, r, rp
   real (kind = rdf), dimension(4,4) :: a, s

   real (kind = rdf) :: dc2, de2, dz2
   real (kind = rdf) :: dcsq, desq, dzsq
   real (kind = rdf) :: dtc, dds
   real (kind = rdf) :: tmp
   real (kind = rdf) :: alpha, alpha_mult

   integer :: ilength, jlength, klength
   integer :: va, vb

   integer :: var

   logical :: BlankingFlag

   type( ptdma_plan_single ) :: plan_csi , plan_eta , plan_zet

   ! Processor boundaries
   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp


   ! eigenvalue multipliers
   !
   lambda(1) = zero ; lambda(2) = zero ; lambda(3) = one ; lambda(4) = -one

   ! grid spacings
   !
   dc2  = pt5*dc ; de2  = pt5*de ; dz2  = pt5*dz
   dcsq = dc*dc  ; desq = de*de  ; dzsq = dz*dz

   ! compute rh = (-dt*S^(-1)*R)
   !
   alpha_mult = e_source * one_pt_five / delti ! 3/(2Δt)
   s = zero

   ! before: rh = R̂ 
   ! after : rh = -Δτ * E^-1 * J * R̂ 
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      ! α = 1 + 3Δτ / 2Δt  
      alpha  = one + alpha_mult * dtau(i,j,k)

      ! s  = E^-1 = diag[ β , 1/α , 1/α , 1/α ]
      s(1,1) = beta
      s(2,2) = one/alpha
      s(3,3) = one/alpha
      s(4,4) = one/alpha

      ! r  = -Δτ * J * R 
      ! rp = -Δτ * E^-1 * J * R = s * R 
      r(:)        = -dtau(i,j,k) * aj(i,j,k) * rh(:,i,j,k)
      rp(:)       =  matmul(s,r)
      rh(:,i,j,k) =  rp(:)

   end do
   end do
   end do


   !==========================================================================================
   !
   !                                   ξ - SWEEP
   !
   !==========================================================================================
  
   va = i_mysta
   vb = idend
  
   ilength = vb-va+1
  
   if ( allocated(aw) ) deallocate ( aw , ap , ae , aq )
   allocate ( aw( va:vb ,1:4 ) , ap( va:vb , 1:4 ) , ae( va:vb , 1:4 ) , aq( va:vb , 1:4 ) )
   
   aw = zero ; ap = one ; ae = zero ; aq = zero

   ! save rh for (semi-explicit) bc treatment
   ! 
   
   ! brh = boundary right-hand side

   if (n == 1) then

      if ( allocated( brh ) ) deallocate (brh)
      
      !allocate ( brh( 1:4 , 1:2 , j_mysta:j_myend , k_mysta:k_myend ) ) 
      allocate ( brh( 1:4 , 1:4 , j_mysta:j_myend , k_mysta:k_myend ) ) 
      
      brh = zero

      if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

            brh(:,1,j,k) = sa(1:4,1) * rh(:,i_mysta  ,j,k) + & ! i = 2
                           sb(1:4,1) * rh(:,i_mysta+1,j,k)     ! i = 3
 
         end do
         end do
      
      end if

      if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
         
            brh(:,2,j,k) = sa(1:4,2) * rh(:,i_myend  ,j,k) + & ! i = im-1
                           sb(1:4,2) * rh(:,i_myend-1,j,k)     ! i = im-2
         
         end do
         end do
      
      end if

      ! Blanking rhs extrapolation
      if ( nblk /=0 ) then
         
         do nb = 1 , nblk
   
            i = li_blk_ia(1,nb)
      
            if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 ) then
                     
               do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
               do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 
      
                  rh(1,i,j,k)   = four/three * rh(1  ,i-1,j,k) - & ! i = ib_ini - 1
                                  one /three * rh(1  ,i-2,j,k)     ! i = ib_ini - 2
                  
                  rh(2:4,i,j,k) =        two * rh(2:4,i-1,j,k) - & ! i = ib_ini - 1
                                         one * rh(2:4,i-2,j,k)     ! i = ib_ini - 2
      
                  if ( non_slip_wall_blanking ) rh(2:4,i,j,k) = zero
   
               end do
               end do
      
            end if
      
            i = li_blk_ib(1,nb)
      
            if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 ) then
                     
               do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
               do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 
      
                  rh(1,i,j,k)   = four/three * rh(1  ,i+1,j,k) - & ! i = ib_end + 1
                                  one /three * rh(1  ,i+2,j,k)     ! i = ib_end + 2
                  
                  rh(2:4,i,j,k) =        two * rh(2:4,i+1,j,k) - & ! i = ib_end + 1
                                         one * rh(2:4,i+2,j,k)     ! i = ib_end + 2
      
                  if ( non_slip_wall_blanking ) rh(2:4,i,j,k) = zero
   
               end do
               end do
      
            end if
      
         end do
      
      end if

   end if ! if (n == 1)

   ! compute Ma^(-1) * (R)
   !

   ! before rh = -Δτ * E^-1 * J * Ř
   ! after  rh = M1^-1  * ( -Δτ * E^-1 * J * Ř )


   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend
   
      a(:,:)      = mai(:,:,i,j,k)
      r(:)        = rh(:,i,j,k)

      ! rp = M1^-1 * ( -Δτ * J * inv(E) * R ) 
      rp(:)       = matmul(a, r)
      rh(:,i,j,k) = rp(:)
   
   end do
   end do
   end do


   if (n == 1) then

      if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
        
         do k = k_mysta, k_myend
         do j = j_mysta, j_myend
                    
            a(:,:)       = mai(:,:,i_mysta-1,j,k)
            r(:)         = brh(:,1,j,k)
            rp(:)        = matmul(a, r)
            brh(:,1,j,k) = rp(:)
         
        end do
        end do

      end if

      if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

            a(:,:)        = mai(:,:,i_myend+1,j,k)
            r(:)          = brh(:,2,j,k)
            rp(:)         = matmul(a, r)
            brh(:,2,j,k)  = rp(:)

         end do
         end do

      end if
  
   end if
  
   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)

   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

   ! csi-sweep each jk line
   ! 
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

      ! Create PaScaL plan
      ! Each proc creates its plan in each direction
      call PaScaL_TDMA_plan_single_create( plan_csi               , &
                                           comm1d_csi(1)%myrank   , &
                                           comm1d_csi(1)%nprocs   , &
                                           comm1d_csi(1)%mpi_comm , &
                                           0                        &
                                          )

      do i = i_mysta , i_myend ! don't include boundary (idend)
                               ! need special treatment for boundary
      
         ! left-hand side
   
         ! center pt (i)            
         !
         dtc = one
         
         if ( dtau(i,j,k) == zero ) dtc = zero
         
         tmp     = dtau(i,j,k) *  ep(1,1,myzone) * spr(1,i,j,k) * dcsq
         
         ap(i,1) = one + two * tmp
         ap(i,2) = ap(i,1)
         ap(i,3) = ap(i,1)
         ap(i,4) = ap(i,1)
   
         ! left pt (i-1)
         !
         dds     =   dtc * dtau(i-1,j,k)  * dc2 * spr(1,i-1,j,k)
         aw(i,1) = - tmp
         aw(i,2) = - tmp
         aw(i,3) = - dds * lambda(3) - tmp
         aw(i,4) = - dds * lambda(4) - tmp
   
         ! right pt (i+1)
         !
         dds     =   dtc * dtau(i+1,j,k) * dc2 * spr(1,i+1,j,k)
         ae(i,1) = - tmp
         ae(i,2) = - tmp
         ae(i,3) =   dds * lambda(3) - tmp
         ae(i,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !
         aq(i,1:4) = rh(1:4,i,j,k)
      
      end do ! do i = i_mysta , i_myend

      ! semi-explicit bc treatment
      ! I don't set ae and aw to 0 at the boundaries, because
      ! those indexes are just not included in the matrix system

      if (n == 1) then

         ! back-side bc
         ! 
         if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
            
            i = i_mysta
            aq(i,1:4)= aq(i,1:4) - aw(i,1:4) * brh(1:4,1,j,k)

         end if

         ! front-side bc (adjust equations so we can solve this plane)
         ! 
         if (myfront == mpi_proc_null .and. btype(2,myzone) /= 5) then
            
            i = i_myend
            aq(i,1:4)= aq(i,1:4) - ae(i,1:4) * brh(1:4,2,j,k)

         end if
        
         !no corregire esto pues nunca deberia entrar, ya que no uso condicion 5  
         if ( myfront == mpi_proc_null .and. btype(2,myzone) == 5) then
         
            dds = dtau(idend,j,k)*dc*spr(1,idend,j,k)
               
            ap(idend,1) = one
            ap(idend,2) = one
            ap(idend,3) = one+dds*lambda(3)
            ap(idend,4) = one

            dds = dtau(idend-1,j,k)*dc*spr(1,idend-1,j,k)
               
            aw(idend,1) =  zero
            aw(idend,2) =  zero
            aw(idend,3) = -dds * lambda(3)
            aw(idend,4) =  zero
           
            aq(idend,1:4) = rh(1:4,idend,j,k)
         
         end if

         ! Blanking correction
         if ( nblk/=0 ) then

            do nb = 1 , nblk

               if ( j >= li_blk_ja(1,nb) .and. j<=li_blk_jb(1,nb) ) then
                  
                  if ( li_blk_ia(1,nb) > i_mysta ) then
                  
                     i = li_blk_ia(1,nb) - 1 
                     aq(i,1:4) = aq(i,1:4) - ae(i,1:4) * rh(1:4,i+1,j,k)

                  end if
                  
                  if ( li_blk_ib(1,nb) < i_myend ) then
               
                     i = li_blk_ib(1,nb) + 1 
                     aq(i,1:4) = aq(i,1:4) - aw(i,1:4) * rh(1:4,i-1,j,k)
               
                  end if   
               
               end if
         
            end do
          
         end if ! nblk/=0  

         ! The blanking node is one node away from the exterior nodes limit
         ! NOTE: Poner else if
         if ( nblke/=0 ) then

            do nb = 1 , nblke

               if ( j >= le_blk_ja(1,nb) .and. j<=le_blk_jb(1,nb) ) then
                  
                  if ( le_blk_ia(1,nb) == iend-1 ) then

                     i = le_blk_ia(1,nb) - 1 
                     aq(i,1:4) = aq(i,1:4) - ae(i,1:4) * rh(1:4,i+1,j,k)

                  end if

                  if ( le_blk_ib(1,nb) == ista+1 ) then

                     i = le_blk_ib(1,nb) + 1 
                     aq(i,1:4) = aq(i,1:4) - aw(i,1:4) * rh(1:4,i-1,j,k)
               
                  end if   
               
               end if
         
            end do
         
         end if ! nblke/=0 

      end if ! if (n == 1)


      ! Free-surface correction
      do i = i_mysta , i_myend

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i-1,j,k) ) < one_half ) then

            aq(i,1:4) = aq(i,1:4) - aw(i,1:4) * rh(1:4,i-1,j,k)

         end if

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i+1,j,k) ) < one_half ) then
            
            aq(i,1:4) = aq(i,1:4) - ae(i,1:4) * rh(1:4,i+1,j,k)
         
         end if

      end do

      do i = i_mysta , i_myend

         aw(i,1:4) = rsign( phi(i-1,j,k) ) * rsign( phi(i,j,k) ) * aw(i,1:4)

         aw(i,1:4) = rsign( phi(i,j,k) ) * aw(i,1:4)
         ap(i,1:4) = rsign( phi(i,j,k) ) * ap(i,1:4)  
         ae(i,1:4) = rsign( phi(i,j,k) ) * ae(i,1:4)  
         
         ae(i,1:4) = rsign( phi(i+1,j,k) ) * rsign( phi(i,j,k) ) * ae(i,1:4)  

      end do


      ! blanking area
      !
      if (nblk /= 0) then
         
         do nb = 1, nblk
            
            if ( k >= li_blk_ka(n,nb) .and. k <= li_blk_kb(n,nb) .and. &
                 j >= li_blk_ja(n,nb) .and. j <= li_blk_jb(n,nb) ) then
               
               aw( max(li_blk_ia(n,nb),va) : min(li_blk_ib(n,nb),vb) , 1:4 ) = zero
               ap( max(li_blk_ia(n,nb),va) : min(li_blk_ib(n,nb),vb) , 1:4 ) = zero
               ae( max(li_blk_ia(n,nb),va) : min(li_blk_ib(n,nb),vb) , 1:4 ) = zero
            
               if (li_blk_ia(n,nb) > i_mysta) ae( max(li_blk_ia(n,nb)-1,va) , 1:4 )=zero
               if (li_blk_ib(n,nb) < i_myend) aw( min(li_blk_ib(n,nb)+1,vb) , 1:4 )=zero
            
            end if

         end do
      
      end if

      ! solve linear tridiagonal eqns

      ! p
      call PaScaL_TDMA_single_solve( plan_csi , &
                                     aw(:,1)  , &
                                     ap(:,1)  , &
                                     ae(:,1)  , &
                                     aq(:,1)  , &
                                     ilength    &
                                    )

      ! u
      call PaScaL_TDMA_single_solve( plan_csi , &
                                     aw(:,2)  , &
                                     ap(:,2)  , &
                                     ae(:,2)  , &
                                     aq(:,2)  , &
                                     ilength    &
                                    )

      ! v
      call PaScaL_TDMA_single_solve( plan_csi , &
                                     aw(:,3)  , &
                                     ap(:,3)  , &
                                     ae(:,3)  , &
                                     aq(:,3)  , &
                                     ilength    &
                                    )

      ! w
      call PaScaL_TDMA_single_solve( plan_csi , &
                                     aw(:,4)  , &
                                     ap(:,4)  , &
                                     ae(:,4)  , &
                                     aq(:,4)  , &
                                     ilength    &
                                    )

      ! put right-hand side back
      !
      do i = i_mysta , idend
         rh(1:4,i,j,k) = aq(i,1:4)
      end do

      call PaScaL_TDMA_plan_single_destroy( plan_csi )

   end do ! j = j_mysta , j_myend
   end do ! k = k_mysta , k_myend


   ! I extrapolate the rhs again
   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )


   !==========================================================================================
   !
   !                                   η - SWEEP
   !
   !==========================================================================================
  
   va = j_mysta
   vb = j_myend
   
   jlength = vb-va+1
   
   if ( allocated( aw ) ) deallocate ( aw , ap , ae , aq )
   allocate ( aw(va:vb,1:4) , ap(va:vb,1:4) , ae(va:vb,1:4) , aq(va:vb,1:4) )
   
   aw = zero ; ap = one ; ae = zero ; aq = zero

   !
   ! save rh for semi-explicit bc treatment
   !
   ! apply bc on fine grid only
   !
   
   if (n == 1) then

      if ( allocated( brh ) ) deallocate ( brh )
      
      !allocate ( brh( 1:4,1:2 , i_mysta:idend , k_mysta:k_myend ) ) 
      allocate ( brh( 1:4,1:4 , i_mysta:idend , k_mysta:k_myend ) ) 
     
      brh=zero

      if ( myleft == mpi_proc_null .and. btype(3,myzone) /= 5 ) then
         
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend

            ! the rh vector here, comes overwritten from the csi-sweep
            brh(:,1,i,k) = sa(1:4,3) * rh(:,i,j_mysta,k) + sb(1:4,3) * rh(:,i,j_mysta+1,k)

         end do
         end do
      
      end if
     
      if ( myright == mpi_proc_null .and. btype(4,myzone) /= 5 ) then
      
         do k = k_mysta , k_myend
         do i = i_mysta , idend !i_myend
            
            brh(:,2,i,k) = sa(1:4,4) * rh(:,i,j_myend,k) + sb(1:4,4) * rh(:,i,j_myend-1,k)

         end do
         end do
      
      end if

      ! Blanking rhs extrapolation
      if (nblk /=0) then
   
         do nb = 1 , nblk
         
            j = li_blk_ja(1,nb)
      
            if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then
                     
               do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
               do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 
      
                  rh(1,i,j,k)   = four/three * rh(1  ,i,j-1,k) - & ! j = jb_ini - 1
                                  one /three * rh(1  ,i,j-2,k)     ! j = jb_ini - 2
                  
                  rh(2:4,i,j,k) =        two * rh(2:4,i,j-1,k) - & ! j = jb_ini - 1
                                         one * rh(2:4,i,j-2,k)     ! j = jb_ini - 2
      
                  if ( non_slip_wall_blanking ) rh(2:4,i,j,k) = zero
      
               end do
               end do
      
            end if
      
            j = li_blk_jb(1,nb)
      
            if ( blktype(4,nb,myzone) == 0 .and. j < j_mysta-1 ) then
                     
      
               do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
               do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) ) 
      
                  rh(1,i,j,k)   = four/three * rh(1  ,i,j+1,k) - & ! j = jb_end + 1
                                  one /three * rh(1  ,i,j+2,k)     ! j = jb_end + 2
                  
                  rh(2:4,i,j,k) =        two * rh(2:4,i,j+1,k) - & ! j = jb_end + 1
                                         one * rh(2:4,i,j+2,k)     ! j = jb_end + 2
      
                  if ( non_slip_wall_blanking ) rh(2:4,i,j,k) = zero
      
               end do
               end do
      
            end if

         end do

      end if   
   
   end if ! if (n == 1)

   ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh
   !
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      r(:)        = rh(:,i,j,k)
      a(:,:)      = n1i(:,:,i,j,k)
      rp(:)       = matmul(a, r)
      rh(:,i,j,k) = rp(:)
      
   end do
   end do
   end do

   ! semi-explicit bc treatment
   ! 
   if (n == 1) then   

      if (myleft == mpi_proc_null .and. btype(3,myzone) /= 5) then
      
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend
                     
            r(:)         = brh(:,1,i,k)
            a(:,:)       = n1i(:,:,i,1,k)
            rp(:)        = matmul(a, r)
            brh(:,1,i,k) = rp(:)
            
         end do
         end do
      
      end if
     
      if (myright == mpi_proc_null .and. btype(4,myzone) /= 5) then
  
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend
                        
            r(:)   = brh(:,2,i,k)
            a(:,:) = n1i(:,:,i,j_myend+1,k)
            rp(:)  = matmul(a, r)
            brh(:,2,i,k) = rp(:)
            
         end do
         end do
  
      end if
  
   end if

   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

   ! eta-sweep each ik line
   !
   do k = k_mysta , k_myend
   do i = i_mysta , idend   !i_myend

      ! Create PaScaL plan
      ! Each proc creates its plan in each direction
      call PaScaL_TDMA_plan_single_create( plan_eta               , &
                                           comm1d_eta(1)%myrank   , &
                                           comm1d_eta(1)%nprocs   , &
                                           comm1d_eta(1)%mpi_comm , &
                                           0                        &
                                          )

      do j = j_mysta , j_myend

         ! left-hand side
   
         ! center pt (j)            
         !
         
         dtc = one
            
         if ( dtau(i,j,k) == zero ) dtc = zero
            
         tmp     = dtau(i,j,k)  * ep(1,2,myzone) * spr(2,i,j,k) * desq
         
         ap(j,1) = one + two * tmp
         ap(j,2) = ap(j,1)
         ap(j,3) = ap(j,1)
         ap(j,4) = ap(j,1)
   
         ! left pt (j-1)
         !
         dds     =   dtc * dtau(i,j-1,k) * de2 * spr(2,i,j-1,k)
         
         aw(j,1) = - tmp
         aw(j,2) = - tmp
         aw(j,3) = - dds * lambda(3) - tmp
         aw(j,4) = - dds * lambda(4) - tmp
   
         ! right pt (j+1)
         !
         dds     =   dtc * dtau(i,j+1,k) * de2 * spr(2,i,j+1,k)
         
         ae(j,1) = - tmp
         ae(j,2) = - tmp
         ae(j,3) =   dds * lambda(3) - tmp
         ae(j,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !
         aq(j,1:4) = rh(1:4,i,j,k)

      end do

      ! semi-explicit bc treatment
      !
      if (n == 1) then

         ! left-side bc
         ! 
         if ( myleft == mpi_proc_null .and. btype(3,myzone) /= 5 ) then
            
            j = j_mysta

            aq(j,1:4)= aq(j,1:4) - aw(j,1:4) * brh(1:4,1,i,k)
         
         end if

         ! right-side bc
         ! 
         if ( myright == mpi_proc_null .and. btype(4,myzone) /= 5 ) then
            
            j = j_myend

            aq(j,1:4)= aq(j,1:4) - ae(j,1:4) * brh(1:4,2,i,k)
         
         end if

         ! Blanking correction

         if ( nblk /= 0 ) then

            do nb = 1 , nblk

               if ( i >= li_blk_ia(1,nb) .and. i<=li_blk_ib(1,nb) ) then
                  
                  if ( li_blk_ja(1,nb) > j_mysta ) then
                  
                     j = li_blk_ja(1,nb) - 1 
                     aq(j,1:4) = aq(j,1:4) - ae(j,1:4) * rh(1:4,i,j+1,k)
                  
                  end if

                  if ( li_blk_jb(1,nb) < j_myend ) then
               
                     j = li_blk_jb(1,nb) + 1 
                     aq(j,1:4) = aq(j,1:4) - aw(j,1:4) * rh(1:4,i,j-1,k)
               
                  end if
               
               end if
         
            end do
         
         end if ! nblk /= 0

         ! The blanking node is one node away from the exterior nodes limit         
         if ( nblke /= 0 ) then

            do nb = 1 , nblke

               if ( i >= le_blk_ia(1,nb) .and. i<=le_blk_ib(1,nb) ) then
                  
                  if ( le_blk_ja(1,nb) == jend-1 ) then
                  
                     j = le_blk_ja(1,nb) - 1 
                     aq(j,1:4) = aq(j,1:4) - ae(j,1:4) * rh(1:4,i,j+1,k)
                  
                  end if

                  if ( le_blk_jb(1,nb) == jsta+1 ) then
               
                     j = le_blk_jb(1,nb) + 1 
                     aq(j,1:4) = aq(j,1:4) - aw(j,1:4) * rh(1:4,i,j-1,k)
               
                  end if
               
               end if
         
            end do
         
         end if ! nblke /= 0

      end if ! if (n == 1)


      ! Free-surface correction

      do j = j_mysta , j_myend

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i,j-1,k) ) < one_half ) then

            aq(j,1:4) = aq(j,1:4) - aw(j,1:4) * rh(1:4,i,j-1,k)

         end if

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i,j+1,k) ) < one_half ) then
            
            aq(j,1:4) = aq(j,1:4) - ae(j,1:4) * rh(1:4,i,j+1,k)
         
         end if

      end do

      do j = j_mysta , j_myend

         aw(j,1:4) = rsign( phi(i,j-1,k) ) * rsign( phi(i,j,k) ) * aw(j,1:4)

         aw(j,1:4) = rsign( phi(i,j,k) ) * aw(j,1:4)
         ap(j,1:4) = rsign( phi(i,j,k) ) * ap(j,1:4)  
         ae(j,1:4) = rsign( phi(i,j,k) ) * ae(j,1:4)  
         
         ae(j,1:4) = rsign( phi(i,j+1,k) ) * rsign( phi(i,j,k) ) * ae(j,1:4)  

      end do


      ! blanking area
      !
      if (nblk /= 0) then
      
         do nb = 1, nblk
      
            if ( li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
                 li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
                     
                  aw( max(li_blk_ja(n,nb),va) : min(li_blk_jb(n,nb),vb) , 1:4 ) = zero
                  ap( max(li_blk_ja(n,nb),va) : min(li_blk_jb(n,nb),vb) , 1:4 ) = zero
                  ae( max(li_blk_ja(n,nb),va) : min(li_blk_jb(n,nb),vb) , 1:4 ) = zero
                  
                  if ( li_blk_ja(n,nb) > j_mysta ) ae( max(li_blk_ja(n,nb)-1,va) , 1:4 ) = zero
                  if ( li_blk_jb(n,nb) < j_myend ) aw( min(li_blk_jb(n,nb)+1,vb) , 1:4 ) = zero
            
            end if
        
         end do
      
      end if

      ! solve linear tridiagonal eqns

      ! p 
      call PaScaL_TDMA_single_solve( plan_eta , &
                                     aw(:,1)  , &
                                     ap(:,1)  , &
                                     ae(:,1)  , &
                                     aq(:,1)  , &
                                     jlength    &
                                    )

      ! u 
      call PaScaL_TDMA_single_solve( plan_eta , &
                                     aw(:,2)  , &
                                     ap(:,2)  , &
                                     ae(:,2)  , &
                                     aq(:,2)  , &
                                     jlength    &
                                   )

      ! v 
      call PaScaL_TDMA_single_solve( plan_eta , &
                                     aw(:,3)  , &
                                     ap(:,3)  , &
                                     ae(:,3)  , &
                                     aq(:,3)  , &
                                     jlength    &
                                   )

      ! w 
      call PaScaL_TDMA_single_solve( plan_eta , &
                                     aw(:,4)  , &
                                     ap(:,4)  , &
                                     ae(:,4)  , &
                                     aq(:,4)  , &
                                     jlength    &
                                   )

      ! put right-hand side back
      !
      do j = j_mysta, j_myend
         rh(1:4,i,j,k) = aq(j,1:4)
      end do

      call PaScaL_TDMA_plan_single_destroy( plan_eta )

   end do !k = k_mysta, k_myend 
   end do !i = i_mysta, idend  


   ! I extrapolate the rhs again
   !
   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)
   
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

   !==========================================================================================
   !
   !                                   ζ - SWEEP
   !
   !==========================================================================================
   
   va = k_mysta
   vb = k_myend
   
   klength = vb-va+1
   
   if ( allocated(aw) ) deallocate (aw,ap,ae,aq)
   allocate ( aw(va:vb,1:4) , ap(va:vb,1:4) , ae(va:vb,1:4) , aq(va:vb,1:4) )
   
   aw = zero ; ap = one ; ae = zero ; aq = zero

   ! save rh for semi-explicit bc treatment
   !
   ! apply bc on fine grid only
   !
   
   if (n == 1) then

      if ( allocated( brh ) ) deallocate ( brh )
      
      !allocate ( brh( 1:4 , 1:2 , i_mysta:idend , j_mysta:j_myend ) ) 
      allocate ( brh( 1:4 , 1:4 , i_mysta:idend , j_mysta:j_myend ) ) 
      
      brh=zero

      if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
      
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
         
            brh(:,1,i,j)=sa(1:4,5)*rh(:,i,j,k_mysta)+sb(1:4,5)*rh(:,i,j,k_mysta+1)

         end do
         end do
      
      end if
     
      if ( myup == mpi_proc_null .and. btype(6,myzone) /= 6 ) then
   
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
                  
            brh(:,2,i,j) = sa(1:4,6)*rh(:,i,j,k_myend) + sb(1:4,6)*rh(:,i,j,k_myend-1)

         end do
         end do
   
      end if
   
   end if ! (n == 1) 

   ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh
   !

   do k = k_mysta, k_myend
   do j = j_mysta, j_myend
   do i = i_mysta, idend    !i_myend

      r(:)        = rh(:,i,j,k)
      a(:,:)      = n2i(:,:,i,j,k)
      rp(:)       = matmul(a, r)
      rh(:,i,j,k) = rp(:)

   end do
   end do
   end do

   ! semi-explicit bc treatment
   ! 
   if (n == 1) then   

      if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
      
         do j = j_mysta , j_myend
         do i = i_mysta , idend     !i_myend

            r(:)         = brh(:,1,i,j)
            a(:,:)       = n2i(:,:,i,j,1)
            rp(:)        = matmul(a, r)
            brh(:,1,i,j) = rp(:)
            
         end do
         end do
      
      end if
     
      if (myup == mpi_proc_null .and. btype(6,myzone) /= 5) then
         
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
         
            r(:)   = brh(:,2,i,j)
            a(:,:) = n2i(:,:,i,j,k_myend+1)
            rp(:)  = matmul(a, r)
            brh(:,2,i,j) = rp(:)
                     
         end do
         end do
      
      end if

   end if

   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

   ! zet-sweep each ik line
   !
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      ! Each proc creates its plan in each direction
      call PaScaL_TDMA_plan_single_create( plan_zet               , &
                                           comm1d_zet(1)%myrank   , &
                                           comm1d_zet(1)%nprocs   , &
                                           comm1d_zet(1)%mpi_comm , &
                                           0                        &
                                          )

      do k = k_mysta, k_myend

         ! left-hand side
   
         ! center pt (k)            
         !
         dtc = one

         if ( dtau(i,j,k) == zero ) dtc = zero
            
         tmp     = dtau(i,j,k)  * ep(1,3,myzone) * spr(3,i,j,k) * dzsq
            
         ap(k,1) = one + two * tmp
         ap(k,2) = ap(k,1)
         ap(k,3) = ap(k,1)
         ap(k,4) = ap(k,1)
   
         ! left pt (k-1)
         !
         dds     = dtc * dtau(i,j,k-1) *dz2 * spr(3,i,j,k-1)

         aw(k,1) = - tmp
         aw(k,2) = - tmp
         aw(k,3) = - dds * lambda(3) - tmp
         aw(k,4) = - dds * lambda(4) - tmp
   
         ! right pt (k+1)
         !
            
         dds     = dtc * dtau(i,j,k+1) *dz2 * spr(3,i,j,k+1)
            
         ae(k,1) = - tmp
         ae(k,2) = - tmp
         ae(k,3) =   dds * lambda(3) - tmp
         ae(k,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !
         aq(k,1:4) = rh(1:4,i,j,k)
      
      end do

      ! semi-explicit bc treatment
      !
      if (n == 1) then

         ! left-side bc
         ! 
         if ( mydown == mpi_proc_null .and. btype(5,myzone) /= 5 ) then
            
            k = k_mysta
         
            aq(k,1:4) = aq(k,1:4) - aw(k,1:4) * brh(1:4,1,i,j)
         
         end if

         ! right-side bc
         ! 
         if ( myup == mpi_proc_null .and. btype(6,myzone) /= 5 ) then
            
            k = k_myend
            aq(k,1:4) = aq(k,1:4) - ae(k,1:4) * brh(1:4,2,i,j)
         
         end if

      end if

      ! Free-surface correction

      do k = k_mysta , k_myend

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i,j,k-1) ) < one_half ) then

            aq(k,1:4) = aq(k,1:4) - aw(k,1:4) * rh(1:4,i,j,k-1)

         end if

         if ( rsign( phi(i,j,k) ) > one_half .and. rsign( phi(i,j,k+1) ) < one_half ) then
            
            aq(k,1:4) = aq(k,1:4) - ae(k,1:4) * rh(1:4,i,j,k+1)
         
         end if

      end do

      do k = k_mysta , k_myend

         aw(k,1:4) = rsign( phi(i,j,k-1) ) * rsign( phi(i,j,k) ) * aw(k,1:4)

         aw(k,1:4) = rsign( phi(i,j,k) ) * aw(k,1:4)
         ap(k,1:4) = rsign( phi(i,j,k) ) * ap(k,1:4)  
         ae(k,1:4) = rsign( phi(i,j,k) ) * ae(k,1:4)  
         
         ae(k,1:4) = rsign( phi(i,j,k+1) ) * rsign( phi(i,j,k) ) * ae(k,1:4)  

      end do


      ! blanking area
      !
      if (nblk /= 0) then
         
         do nb = 1, nblk
            
            if ( li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) .and. &
                 li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
            
               aw( max(li_blk_ka(n,nb),va) : min(li_blk_kb(n,nb),vb) , 1:4 ) = zero
               ap( max(li_blk_ka(n,nb),va) : min(li_blk_kb(n,nb),vb) , 1:4 ) = zero
               ae( max(li_blk_ka(n,nb),va) : min(li_blk_kb(n,nb),vb) , 1:4 ) = zero
            
               if (li_blk_ka(n,nb) > k_mysta) ae( max(li_blk_ka(n,nb)-1,va) , 1:4 ) = zero
               if (li_blk_kb(n,nb) < k_myend) aw( min(li_blk_kb(n,nb)+1,vb) , 1:4 ) = zero
            
            end if
         
         end do
      
      end if

      ! solve linear tridiagonal eqns
      !
      klength = k_myend - k_mysta + 1
           
      ! p
      call PaScaL_TDMA_single_solve( plan_zet , &
                                     aw(:,1)  , &
                                     ap(:,1)  , &
                                     ae(:,1)  , &
                                     aq(:,1)  , &
                                     klength    &
                                    )
      ! u
      call PaScaL_TDMA_single_solve( plan_zet , &
                                     aw(:,2)  , &
                                     ap(:,2)  , &
                                     ae(:,2)  , &
                                     aq(:,2)  , &
                                     klength    &
                                    )

      ! v
      call PaScaL_TDMA_single_solve( plan_zet , &
                                     aw(:,3)  , &
                                     ap(:,3)  , &
                                     ae(:,3)  , &
                                     aq(:,3)  , &
                                     klength    &
                                    )
      ! w
      call PaScaL_TDMA_single_solve( plan_zet , &
                                     aw(:,4)  , &
                                     ap(:,4)  , &
                                     ae(:,4)  , &
                                     aq(:,4)  , &
                                     klength    &
                                    )

      ! put right-hand side back
      !
      do k = k_mysta, k_myend
        rh(1:4,i,j,k) = aq(k,1:4)
      end do

      call PaScaL_TDMA_plan_single_destroy( plan_zet )

   end do
   end do

   ! I extrapolate the rhs again
   !
   !rh(1,:,:,:) = rsign(:,:,:) * rh(1,:,:,:)
   !rh(2,:,:,:) = rsign(:,:,:) * rh(2,:,:,:)
   !rh(3,:,:,:) = rsign(:,:,:) * rh(3,:,:,:)
   !rh(4,:,:,:) = rsign(:,:,:) * rh(4,:,:,:)
   !
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )


   if ( allocated( aw  ) ) deallocate ( aw,ap,ae,aq )
   if ( allocated( brh ) ) deallocate ( brh         )

   !==============================
   !
   ! Final Update
   !
   !==============================
   ! dQ = Mc * dQ^***

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend !i_myend
           
      a(:,:)      = mc(:,:,i,j,k)
      r(:)        = rh(:,i,j,k)
      rp(:)       = matmul(a,r)
      rh(:,i,j,k) = rp(:)
      
   end do
   end do
   end do

   ! Ghost nodes update
   call rhs_exchng3_4d( rh )
   ! ghost fluid nodes extrapolation
   call rhs_ghost_fluid_nodes_extp_4d( rh , InteriorNodesOnly = .true. )
   ! Ghost nodes update
   call rhs_exchng3_4d( rh )

end subroutine rhs_diag_solver


