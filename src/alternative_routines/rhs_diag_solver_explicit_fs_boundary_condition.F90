
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

   use InterpolationMethods, only: Invert4x4Matrix

   implicit none
   
   integer :: info

   ! dimension lhs and temporary rhs
   real (kind = rdf), dimension(:,:), allocatable :: aw, ap, ae, aq

   ! boundary condition array
   real (kind = rdf), dimension(:,:,:,:), allocatable :: brh

   ! local dummy variables
   real (kind = rdf), dimension(4) :: lambda, r, rp
   real (kind = rdf), dimension(4,4) :: Mlocal, Einv

   real (kind = rdf) :: dc2, de2, dz2
   real (kind = rdf) :: dcsq, desq, dzsq
   real (kind = rdf) :: dtc, dds
   real (kind = rdf) :: tmp
   real (kind = rdf) :: alpha, alpha_mult

   real (kind = rdf), parameter :: tol = 1.0E-14 ! tolerance for interpolation near the free surface

   integer :: ilength, jlength, klength
   integer :: va, vb

   ! Near free surface extrapolation

   real(kind = rdf) :: InterpCoeff 
   real(kind = rdf) :: Sa_fs , Sb_fs
   real(kind = rdf) , dimension(4)   :: ae_fs , aw_fs
   real (kind = rdf), dimension(4,4) :: Minv_fs , Maux1, Maux2
   real(kind = rdf) , dimension(4)   :: rh_w, rh_p, rh_e, rp_fs

   ! TO DO: figure out how to avoid unecessary computation in this subroutine
   ! using rsing

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
   Einv = zero


   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      !alpha =one+alpha_mult*dtau(i,j,k) * rhoLSM(phi_n(i,j,k))
         
      ! α = 1 + 3Δτ / 2Δt  

      alpha  = one + alpha_mult * dtau(i,j,k)

      ! s  = E^-1 = diag[ β , 1/α , 1/α , 1/α ]

      Einv(1,1) = beta
      Einv(2,2) = one/alpha
      Einv(3,3) = one/alpha
      Einv(4,4) = one/alpha

      ! r  = -Δτ * J * R 
      ! rp = -Δτ * E^-1 * J * R  

      r (1:4)       = -dtau(i,j,k) * aj(i,j,k) * rh(:,i,j,k)
      rp(1:4)       =  matmul( Einv , r )
      rh(1:4,i,j,k) =  rp(1:4) ! rh = -Δτ * E^-1 * J * R 

   end do
   end do
   end do

   !==============================
   !
   ! CSI SWEEP
   !
   !==============================
  
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
      
      !allocate ( brh( 1:4 , 1:4 , j_mysta:j_myend , k_mysta:k_myend ) ) 
      allocate ( brh( 1:4 , 1:2 , j_mysta:j_myend , k_mysta:k_myend ) ) 
      
      brh = zero

      if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

!            if( rsign(i_mysta,j,k) > one_half .and. rsign(i_mysta+1,j,k) > one_half) then

               ! extrapolation of the rhs to the i = 1 boundary
               ! rh_{1,j,k} = sa * rh_{2,j,k} + sb * rh_{3,j,k}
               brh(:,1,j,k) = sa(:,1) * rh(:,i_mysta,j,k) + sb(:,1) * rh(:,i_mysta+1,j,k)
         
!            end if

         end do
         end do
      
      end if

      if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend
         
!            if( rsign(i_myend,j,k) > one_half .and. rsign(i_myend-1,j,k) > one_half ) then

               ! extrapolation of the rhs to the i = im boundary
               ! rh_{im,j,k} = sa * rh_{im-1,j,k} + sb * rh_{im-2,j,k}
               ! remember: i_myend = im-1
               brh(:,2,j,k) = sa(:,2) * rh(:,i_myend,j,k) + sb(:,2) * rh(:,i_myend-1,j,k)
         
!            end if

         end do
         end do
      
      end if
   
   end if ! if (n == 1)


   !aw, ap, and ae coefficients

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

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
         
      end do
   
   end do 
   end do


   ! compute Ma^(-1) * (R) over interior nodes

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend
      
      Mlocal(1:4,1:4) = mai(1:4,1:4,i,j,k)
      r(1:4)          = rh (1:4    ,i,j,k)

      ! rp = M1^-1 * ( -Δτ * J * inv(E) * R ) 
      rp(1:4)       = matmul(Mlocal, r)
      rh(1:4,i,j,k) = rp(1:4)

   end do
   end do
   end do


   ! Near free-surface correction

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

      do i = i_mysta , idend  
   
         ! NEAR FREE SURFACE LAYER BOUNDARY CONDITIONS
   
         ! i node within the WATER phase, i+1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i+1,j,k) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the ξ direction between i and i+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i+1,j,k) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i+1,j,k) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            ae_fs  (1:4)     = Sa_fs * ae(i,1:4)          + Sb_fs * ae(i-1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * mai(1:4,1:4,i,j,k) + Sb_fs * mai(1:4,1:4,i-1,j,k)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly         
   
            Maux1(1:4,1:4) = Invert4x4Matrix( mai(1:4,1:4,i  ,j,k) )
            Maux2(1:4,1:4) = Invert4x4Matrix( mai(1:4,1:4,i-1,j,k) )
   
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i  ,j,k) )
            rh_w(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i-1,j,k) )
   
            ! Extrapolated rhs
   
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p(1:4) + Sb_fs * rh_w(1:4) )
   
            ! Updated RHS
            aq(i,1:4) = rh(1:4,i,j,k) - ae_fs(1:4) * rp_fs(1:4)
   
         end if
   
         ! i node within the WATER phase, i-1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i-1,j,k) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the ξ direction between i and i+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i-1,j,k) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i-1,j,k) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            aw_fs(1:4)       = Sa_fs * aw(i,1:4)          + Sb_fs * aw(i+1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * mai(1:4,1:4,i,j,k) + Sb_fs * mai(1:4,1:4,i+1,j,k)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly
   
            Maux1(1:4,1:4) = Invert4x4Matrix( mai(1:4,1:4,i  ,j,k) )
            Maux2(1:4,1:4) = Invert4x4Matrix( mai(1:4,1:4,i+1,j,k) )
   
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i  ,j,k) )
            rh_e(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i+1,j,k) )
   
            ! Extrapolated rhs
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p + Sb_fs * rh_e )
   
            ! Updated RHS
            aq(i,1:4) = rh(1:4,i,j,k) - aw_fs(1:4) * rp_fs(1:4)
   
         end if
   
      end do

   end do
   end do

   ! compute Ma^(-1) * (R) at the boundaries
   if (n == 1) then

      if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
        
         do k = k_mysta, k_myend
         do j = j_mysta, j_myend
           
!            if( rsign(i_mysta-1,j,k) > one_half ) then
         
               Mlocal(:,:)  = mai(:,:,i_mysta-1,j,k) ! Ma^{-1} at i = 1
               r(:)         = brh(:,1,j,k)
               rp(:)        = matmul( Mlocal , r)
               brh(:,1,j,k) = rp(:)
         
!            end if  

        end do
        end do

      end if

      if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
         
         do k = k_mysta , k_myend
         do j = j_mysta , j_myend

!            if( rsign(i_myend+1,j,k) > one_half ) then

               Mlocal(:,:)   = mai(:,:,i_myend+1,j,k) ! Ma^{-1} at i = im
               r(:)          = brh(:,2,j,k)
               rp(:)         = matmul( Mlocal , r )
               brh(:,2,j,k)  = rp(:)

!            end if

         end do
         end do

      end if
  
   end if
  
   ! csi-sweep each jk line
   ! 
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend

      do i = i_mysta , i_myend ! don't include boundary (idend)
                               ! need special treatment for boundary
   
         ! left-hand side
   
         ! center pt (i)            
         !
         !dtc = one
         !
         !if ( dtau(i,j,k) == zero ) dtc = zero
         !
         !tmp     = dtau(i,j,k) *  ep(1,1,myzone) * spr(1,i,j,k) * dcsq
         !
         !ap(i,1) = one + two * tmp
         !ap(i,2) = ap(i,1)
         !ap(i,3) = ap(i,1)
         !ap(i,4) = ap(i,1)
   !
         !! left pt (i-1)
         !!
         !dds     =   dtc * dtau(i-1,j,k)  * dc2 * spr(1,i-1,j,k)
         !aw(i,1) = - tmp
         !aw(i,2) = - tmp
         !aw(i,3) = - dds * lambda(3) - tmp
         !aw(i,4) = - dds * lambda(4) - tmp
   !
         !! right pt (i+1)
         !!
         !dds     =   dtc * dtau(i+1,j,k) * dc2 * spr(1,i+1,j,k)
         !ae(i,1) = - tmp
         !ae(i,2) = - tmp
         !ae(i,3) =   dds * lambda(3) - tmp
         !ae(i,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !

         if ( ( rsign(i,j,k) > one_half .and. rsign(i+1,j,k) < one_half ) .or. & 
              ( rsign(i,j,k) > one_half .and. rsign(i-1,j,k) < one_half )        ) then

            cycle

         else
         
            aq(i,1:4) = rh(1:4,i,j,k)
         
         end if

      end do

      ! semi-explicit bc treatment
      !
      if (n == 1) then

         ! back-side bc
         ! 
         if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
            
            i = i_mysta

            !if( rsign(i,j,k) > one_half ) 
            aq(i,1:4)= aq(i,1:4) - aw(i,1:4) * brh(1:4,1,j,k)

         end if

         ! front-side bc (adjust equations so we can solve this plane)
         ! 
         if (myfront == mpi_proc_null .and. btype(2,myzone) /= 5) then
            
            i = i_myend
            
            !if( rsign(i,j,k) > one_half ) 
            aq(i,1:4)= aq(i,1:4) - ae(i,1:4) * brh(1:4,2,j,k)

         end if
        
         !no corregire esto pues nunca deberia entrar, ya que no uso condicion 5  
         if ( myfront == mpi_proc_null .and. btype(2,myzone) == 5) then
         
!            if( rsign(idend,j,k) > one_half ) then

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
         
!            end if
        
         end if

      end if ! if (n == 1)

      if ( rsign(i,j,k) > one_half .and. rsign(i+1,j,k) < one_half ) ae(i,1:4) = zero
      if ( rsign(i,j,k) > one_half .and. rsign(i-1,j,k) < one_half ) aw(i,1:4) = zero

      ! blanking area
      !
      if (nblk /= 0) then
         
         do nb = 1, nblk
            
            if ( li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
                 li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) ) then
               
               aw( li_blk_ia(n,nb) : li_blk_ib(n,nb) , 1:4 ) = zero
               ap( li_blk_ia(n,nb) : li_blk_ib(n,nb) , 1:4 ) = zero
               ae( li_blk_ia(n,nb) : li_blk_ib(n,nb) , 1:4 ) = zero
            
               if (li_blk_ia(n,nb) > i_mysta) ae( li_blk_ia(n,nb)-1 , 1:4 )=zero
               if (li_blk_ib(n,nb) < i_myend) aw( li_blk_ib(n,nb)+1 , 1:4 )=zero
            
            end if

         end do
      
      end if

      ! solve linear tridiagonal eqns
      !
      !call sgtsv(ilength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
      !           ilength, info)
      !call sgtsv(ilength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
      !           ilength, info)
      !call sgtsv(ilength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
      !           ilength, info)
      !call sgtsv(ilength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
      !           ilength, info)

      ! TO DO: maybe we can make these array operations more computationally efficient

      !aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign_aux( va+1:vb , j , k ) 
      !aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign_aux( va+1:vb , j , k ) 
      !aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign_aux( va+1:vb , j , k ) 
      !aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign_aux( va+1:vb , j , k ) 

      !ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign_aux( va:vb   , j , k ) 
      !ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign_aux( va:vb   , j , k ) 
      !ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign_aux( va:vb   , j , k ) 
      !ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign_aux( va:vb   , j , k ) 

      !ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign_aux( va:vb-1 , j , k ) 
      !ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign_aux( va:vb-1 , j , k ) 
      !ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign_aux( va:vb-1 , j , k ) 
      !ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign_aux( va:vb-1 , j , k ) 

      !aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign_aux( va:vb   , j , k ) 
      !aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign_aux( va:vb   , j , k ) 
      !aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign_aux( va:vb   , j , k ) 
      !aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign_aux( va:vb   , j , k ) 


      aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign( va+1:vb , j , k ) 
      aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign( va+1:vb , j , k ) 
      aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign( va+1:vb , j , k ) 
      aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign( va+1:vb , j , k ) 

      ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign( va:vb   , j , k ) 
      ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign( va:vb   , j , k ) 
      ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign( va:vb   , j , k ) 
      ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign( va:vb   , j , k ) 

      ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign( va:vb-1 , j , k ) 
      ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign( va:vb-1 , j , k ) 
      ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign( va:vb-1 , j , k ) 
      ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign( va:vb-1 , j , k ) 

      aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign( va:vb   , j , k ) 
      aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign( va:vb   , j , k ) 
      aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign( va:vb   , j , k ) 
      aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign( va:vb   , j , k ) 


      ! p
      call dgtsv( ilength , 1                                                               , &
                  aw( va+1:vb , 1 ) , ap( va:vb , 1 ) , ae( va:vb-1 , 1 ) , aq( va:vb , 1 ) , &
                  ilength , info                                                              &
                )
      ! u

      call dgtsv( ilength , 1                                                               , &
                  aw( va+1:vb , 2 ) , ap( va:vb , 2 ) , ae( va:vb-1 , 2 ) , aq( va:vb , 2 ) , &
                  ilength , info                                                              &
                )
      ! v
      call dgtsv( ilength , 1                                                               , &
                  aw( va+1:vb , 3 ) , ap( va:vb , 3 ) , ae( va:vb-1 , 3 ) , aq( va:vb , 3 ) , &
                  ilength , info                                                              &
                )
      ! w
      call dgtsv( ilength , 1                                                               , &
                  aw( va+1:vb , 4 ) , ap( va:vb , 4 ) , ae( va:vb-1 , 4 ) , aq( va:vb , 4 ) , &
                  ilength , info                                                              &
                )

      ! put right-hand side back
      !
      do i = i_mysta , idend

         !if( rsign(i,j,k) > one_half ) 
         rh(1:4,i,j,k) = aq(i,1:4)
            
      end do

   end do ! j = j_mysta , j_myend
   end do ! k = k_mysta , k_myend

   !==============================
   !
   ! ETA SWEEP
   !
   !==============================
  
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
      
!      allocate ( brh( 1:4,1:4 , i_mysta:idend , k_mysta:k_myend ) ) 
      allocate ( brh( 1:4,1:2 , i_mysta:idend , k_mysta:k_myend ) ) 
     
      brh=zero

      if ( myleft == mpi_proc_null .and. btype(3,myzone) /= 5 ) then
         
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend

!            if( rsign(i,j_mysta,k) > one_half .and. rsign(i,j_mysta+1,k) > one_half ) then

               ! the rh vector here, comes overwritten from the csi-sweep

               brh(:,1,i,k) = sa(:,3) * rh(:,i,j_mysta,k) + sb(:,3) * rh(:,i,j_mysta+1,k)

!            end if

         end do
         end do
      
      end if
     
      if ( myright == mpi_proc_null .and. btype(4,myzone) /= 5 ) then
      
         do k = k_mysta , k_myend
         do i = i_mysta , idend !i_myend
            
!            if( rsign(i,j_myend,k) > one_half .and. rsign(i,j_myend-1,k) > one_half ) then

               brh(:,2,i,k) = sa(:,4) * rh(:,i,j_myend,k) + sb(:,4) * rh(:,i,j_myend-1,k)
            
!            end if

         end do
         end do
      
      end if
   
   end if ! if (n == 1)

   !aw, ap, and ae coefficients

   do k = k_mysta, k_myend
   do i = i_mysta, idend   !i_myend

      do j = j_mysta, j_myend

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

      end do

   end do
   end do


   ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh on interior nodes
   !
   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      Mlocal(1:4,1:4) = n1i(1:4,1:4,i,j,k)
      r(1:4)          = rh (1:4,    i,j,k)

      rp(1:4)         = matmul( Mlocal , r )
      rh(1:4,i,j,k)   = rp(1:4)
      
   end do
   end do
   end do

   ! Near free-surface correction

   do k = k_mysta , k_myend
   do i = i_mysta , idend  

      do j = j_mysta , j_myend
   
         ! NEAR FREE SURFACE LAYER BOUNDARY CONDITIONS
   
         ! j node within the WATER phase, j+1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i,j+1,k) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the η direction between j and j+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i,j+1,k) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i,j+1,k) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            ae_fs(1:4)       = Sa_fs * ae(j,1:4)          + Sb_fs * ae(j-1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * n1i(1:4,1:4,i,j,k) + Sb_fs * n1i(1:4,1:4,i,j-1,k)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly
            Maux1(1:4,1:4) = Invert4x4Matrix( n1i(1:4,1:4,i,j  ,k) )
            Maux2(1:4,1:4) = Invert4x4Matrix( n1i(1:4,1:4,i,j-1,k) )
            
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i,j  ,k) ) 
            rh_w(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i,j-1,k) ) 
   
            ! Extrapolated rhs
   
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p(1:4) + Sb_fs * rh_w(1:4) )
   
            ! Updated RHS
            aq(j,1:4) = rh(1:4,i,j,k) - ae_fs(1:4) * rp_fs(1:4)
   
         end if
   
         ! j node within the WATER phase, j-1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i,j-1,k) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the η direction between j and j+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i,j-1,k) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i,j-1,k) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            aw_fs(1:4)       = Sa_fs * aw(j,1:4)          + Sb_fs * aw(j+1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * n1i(1:4,1:4,i,j,k) + Sb_fs * n1i(1:4,1:4,i,j+1,k)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly
   
            Maux1(1:4,1:4) = Invert4x4Matrix( n1i(1:4,1:4,i,j  ,k) )
            Maux2(1:4,1:4) = Invert4x4Matrix( n1i(1:4,1:4,i,j+1,k) )
            
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i,j  ,k) ) 
            rh_e(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i,j+1,k) ) 
   
            ! Extrapolated rhs
   
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p(1:4) + Sb_fs * rh_e(1:4) )
   
            ! Updated RHS
            aq(j,1:4) = rh(1:4,i,j,k) - aw_fs(1:4) * rp_fs(1:4)
   
         end if
   
      end do
   
   end do
   end do


   ! semi-explicit bc treatment
   ! 
   if (n == 1) then   

      if (myleft == mpi_proc_null .and. btype(3,myzone) /= 5) then
      
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend
            
            Mlocal(:,:)  = n1i(:,:,i,1,k)
            r(:)         = brh(:,1,i,k)

            rp(:)        = matmul( Mlocal , r)
            brh(:,1,i,k) = rp(:)
                     
         end do
         end do
      
      end if
     
      if (myright == mpi_proc_null .and. btype(4,myzone) /= 5) then
  
         do k = k_mysta, k_myend
         do i = i_mysta, idend !i_myend
                        
            r(:)         = brh(:,2,i,k)
            Mlocal(:,:)  = n1i(:,:,i,j_myend+1,k)
            rp(:)        = matmul( Mlocal , r )
            brh(:,2,i,k) = rp(:)
            
         end do
         end do
  
      end if
  
   end if

   ! eta-sweep each ik line
   !
   do k = k_mysta, k_myend
   do i = i_mysta, idend   !i_myend

      do j = j_mysta, j_myend

         ! left-hand side
   
         ! center pt (j)            
         !
         
         !dtc = one
         
         !if ( dtau(i,j,k) == zero ) dtc = zero
         
         !tmp     = dtau(i,j,k)  * ep(1,2,myzone) * spr(2,i,j,k) * desq
         
         !ap(j,1) = one + two * tmp
         !ap(j,2) = ap(j,1)
         !ap(j,3) = ap(j,1)
         !ap(j,4) = ap(j,1)
   
         ! left pt (j-1)
         !
         !dds     =   dtc * dtau(i,j-1,k) * de2 * spr(2,i,j-1,k)
         
         !aw(j,1) = - tmp
         !aw(j,2) = - tmp
         !aw(j,3) = - dds * lambda(3) - tmp
         !aw(j,4) = - dds * lambda(4) - tmp
   
         ! right pt (j+1)
         !
         !dds     =   dtc * dtau(i,j+1,k) * de2 * spr(2,i,j+1,k)
         
         !ae(j,1) = - tmp
         !ae(j,2) = - tmp
         !ae(j,3) =   dds * lambda(3) - tmp
         !ae(j,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !

         if ( ( rsign(i,j,k) > one_half .and. rsign(i,j+1,k) < one_half ) .or. &
              ( rsign(i,j,k) > one_half .and. rsign(i,j-1,k) < one_half )         ) then

            cycle

         else

            aq(j,1:4) = rh(1:4,i,j,k)

         end if

      end do

      ! semi-explicit bc treatment
      !
      if (n == 1) then

         ! left-side bc
         ! 
         if ( myleft == mpi_proc_null .and. btype(3,myzone) /= 5 ) then
            
            j = j_mysta

            !if( rsign(i,j,k) > one_half ) 
            aq(j,1:4)= aq(j,1:4) - aw(j,1:4) * brh(1:4,1,i,k)
         
         end if

         ! right-side bc
         ! 
         if ( myright == mpi_proc_null .and. btype(4,myzone) /= 5 ) then
            
            j = j_myend

            !if( rsign(i,j,k) > one_half ) 
            aq(j,1:4)= aq(j,1:4) - ae(j,1:4) * brh(1:4,2,i,k)
         
         end if

      end if

      if ( rsign(i,j,k) > one_half .and. rsign(i,j+1,k) < one_half ) ae(j,1:4) = zero
      if ( rsign(i,j,k) > one_half .and. rsign(i,j-1,k) < one_half ) aw(j,1:4) = zero      

      ! blanking area
      !
      if (nblk /= 0) then
      
         do nb = 1, nblk
      
            if ( li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
                 li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
                     
                  aw( li_blk_ja(n,nb) : li_blk_jb(n,nb) , 1:4 ) = zero
                  ap( li_blk_ja(n,nb) : li_blk_jb(n,nb) , 1:4 ) = zero
                  ae( li_blk_ja(n,nb) : li_blk_jb(n,nb) , 1:4 ) = zero
                  
                  if (li_blk_ja(n,nb) > j_mysta) ae( li_blk_ja(n,nb)-1 , 1:4 ) = zero
                  if (li_blk_jb(n,nb) < j_myend) aw( li_blk_jb(n,nb)+1 , 1:4 ) = zero
            
            end if
        
         end do
      
      end if

      ! solve linear tridiagonal eqns
      !
      !call sgtsv(jlength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
      !           jlength, info)
      !call sgtsv(jlength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
      !           jlength, info)
      !call sgtsv(jlength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
      !           jlength, info)
      !call sgtsv(jlength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
      !           jlength, info)


      !aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign_aux( i , va+1:vb , k ) 
      !aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign_aux( i , va+1:vb , k ) 
      !aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign_aux( i , va+1:vb , k ) 
      !aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign_aux( i , va+1:vb , k ) 

      !ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign_aux( i , va:vb   , k ) 
      !ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign_aux( i , va:vb   , k ) 
      !ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign_aux( i , va:vb   , k ) 
      !ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign_aux( i , va:vb   , k ) 

      !ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign_aux( i , va:vb-1 , k ) 
      !ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign_aux( i , va:vb-1 , k ) 
      !ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign_aux( i , va:vb-1 , k ) 
      !ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign_aux( i , va:vb-1 , k ) 

      !aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign_aux( i , va:vb   , k ) 
      !aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign_aux( i , va:vb   , k ) 
      !aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign_aux( i , va:vb   , k ) 
      !aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign_aux( i , va:vb   , k ) 



      aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign( i , va+1:vb , k ) 
      aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign( i , va+1:vb , k ) 
      aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign( i , va+1:vb , k ) 
      aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign( i , va+1:vb , k ) 

      ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign( i , va:vb   , k ) 
      ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign( i , va:vb   , k ) 
      ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign( i , va:vb   , k ) 
      ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign( i , va:vb   , k ) 

      ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign( i , va:vb-1 , k ) 
      ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign( i , va:vb-1 , k ) 
      ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign( i , va:vb-1 , k ) 
      ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign( i , va:vb-1 , k ) 

      aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign( i , va:vb   , k ) 
      aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign( i , va:vb   , k ) 
      aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign( i , va:vb   , k ) 
      aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign( i , va:vb   , k ) 


      ! p 

      call dgtsv( jlength , 1                                                               , &
                  aw( va+1:vb , 1 ) , ap( va:vb , 1 ) , ae( va:vb-1 , 1 ) , aq( va:vb , 1 ) , &
                  jlength , info                                                              &
                )

      ! u 
      call dgtsv( jlength , 1                                                               , &
                  aw( va+1:vb , 2 ) , ap( va:vb , 2 ) , ae( va:vb-1 , 2 ) , aq( va:vb , 2 ) , &
                  jlength , info                                                              &
                )
      ! v 
      call dgtsv( jlength , 1                                                               , &
                  aw( va+1:vb , 3 ) , ap( va:vb , 3 ) , ae( va:vb-1 , 3 ) , aq( va:vb , 3 ) , &
                  jlength , info                                                              &
                )
      ! w 
      call dgtsv( jlength , 1                                                               , &
                  aw( va+1:vb , 4 ) , ap( va:vb , 4 ) , ae( va:vb-1 , 4 ) , aq( va:vb , 4 ) , &
                  jlength , info                                                              &
                )

      ! put right-hand side back
      !
      do j = j_mysta, j_myend

         !if( rsign(i,j,k) > one_half ) 
         rh(1:4,i,j,k) = aq(j,1:4)
      
      end do

   end do !k = k_mysta, k_myend 
   end do !i = i_mysta, idend  

   !==============================
   !
   ! ZET SWEEP
   !
   !==============================
   
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
      
!      allocate ( brh( 1:4 , 1:4 , i_mysta:idend , j_mysta:j_myend ) ) 
      allocate ( brh( 1:4 , 1:2 , i_mysta:idend , j_mysta:j_myend ) ) 
      
      brh=zero

      if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
      
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
         
!            if( rsign(i,j,k_mysta) > one_half .and. rsign(i,j,k_mysta+1) > one_half) then

               brh(:,1,i,j)=sa(:,5)*rh(:,i,j,k_mysta)+sb(:,5)*rh(:,i,j,k_mysta+1)
      
!            end if

         end do
         end do
      
      end if
     
      if ( myup == mpi_proc_null .and. btype(6,myzone) /= 6 ) then
   
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
   
!            if( rsign(i,j,k_myend) > one_half .and. rsign(i,j,k_myend-1) > one_half) then
               
               brh(:,2,i,j)=sa(:,6)*rh(:,i,j,k_myend)+sb(:,6)*rh(:,i,j,k_myend-1)

!            end if

         end do
         end do
   
      end if
   
   end if ! (n == 1) 


   !aw, ap, and ae coefficients


   ! zet-sweep each ik line
   !
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

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
   
      end do

   end do
   end do


   ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh
   !

   do k = k_mysta, k_myend
   do j = j_mysta, j_myend
   do i = i_mysta, idend    !i_myend

!      if( rsign(i,j,k) > one_half ) then

         Mlocal(:,:) = n2i(:,:,i,j,k)
         r(:)        = rh(:,i,j,k)

         rp(:)       = matmul( Mlocal , r )
         rh(:,i,j,k) = rp(:)
      
!      end if

   end do
   end do
   end do


   ! Near free-surface correction

   do j = j_mysta , j_myend
   do i = i_mysta , idend  

      do k = k_mysta , k_myend
   
         ! NEAR FREE SURFACE LAYER BOUNDARY CONDITIONS
   
         ! i node within the WATER phase, i+1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i,j,k+1) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the ξ direction between i and i+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i,j,k+1) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i,j,k+1) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            ae_fs(1:4)       = Sa_fs * ae(k,1:4)          + Sb_fs * ae(k-1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * n2i(1:4,1:4,i,j,k) + Sb_fs * n2i(1:4,1:4,i,j,k-1)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly
   
            Maux1(1:4,1:4) = Invert4x4Matrix( n2i(1:4,1:4,i,j  ,k) )
            Maux2(1:4,1:4) = Invert4x4Matrix( n2i(1:4,1:4,i,j-1,k) )
   
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i,j,k  ) ) 
            rh_w(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i,j,k-1) ) 
   
            ! Extrapolated rhs
   
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p(1:4) + Sb_fs * rh_w(1:4) )
   
            ! Updated RHS
            aq(k,1:4) = rh(1:4,i,j,k) - ae_fs(1:4) * rp_fs(1:4)
   
         end if
   
         ! i node within the WATER phase, i-1 node within the AIR phase
         if ( rsign(i,j,k) > one_half .and. rsign(i,j,k-1) < one_half ) then
   
            ! Interpolation coefficient based on the location of the free surface
            ! along the ξ direction between i and i+1 nodes
   
            InterpCoeff = zero
   
            if ( abs ( phi(i,j,k) - phi(i,j,k-1) ) > tol ) then   
               InterpCoeff = abs( phi(i,j,k) / ( phi(i,j,k) - phi(i,j,k-1) ) )
            end if
   
            ! Extrapolation coefficients: from {i-1, i} to the free surface (fs)
            Sa_fs = ( one + InterpCoeff )
            Sb_fs = -InterpCoeff
   
            ! Extrapolated variables: matrix coefficients and inverse modal matrix
            aw_fs(1:4)       = Sa_fs * aw(k,1:4)          + Sb_fs * aw(k+1,1:4)
            Minv_fs(1:4,1:4) = Sa_fs * n2i(1:4,1:4,i,j,k) + Sb_fs * n2i(1:4,1:4,i,j,k+1)
   
            ! I have to revert the matrix multiplication to extrapolate the rhs properly
   
            Maux1(1:4,1:4) = Invert4x4Matrix( n2i(1:4,1:4,i,j,k  ) )
            Maux2(1:4,1:4) = Invert4x4Matrix( n2i(1:4,1:4,i,j,k+1) )         
            
            rh_p(1:4) = matmul( Maux1(1:4,1:4) , rh(1:4,i,j,k  ) )
            rh_e(1:4) = matmul( Maux2(1:4,1:4) , rh(1:4,i,j,k+1) )
   
            ! Extrapolated rhs
   
            rp_fs(1:4) = matmul ( Minv_fs , Sa_fs * rh_p(1:4) + Sb_fs * rh_e(1:4) )
   
            ! Updated RHS
            aq(k,1:4) = rh(1:4,i,j,k) - aw_fs(1:4) * rp_fs(1:4)
   
         end if
   
      end do

   end do
   end do


   ! semi-explicit bc treatment
   ! 
   if (n == 1) then   

      if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
      
         do j = j_mysta , j_myend
         do i = i_mysta , idend     !i_myend
      
!            if( rsign(i,j,1) > one_half ) then

               r(:)         = brh(:,1,i,j)
               Mlocal(:,:)  = n2i(:,:,i,j,1)
               rp(:)        = matmul( Mlocal , r )
               brh(:,1,i,j) = rp(:)
            
!            end if

         end do
         end do
      
      end if
     
      if (myup == mpi_proc_null .and. btype(6,myzone) /= 5) then
         
         do j = j_mysta, j_myend
         do i = i_mysta, idend !i_myend
         
!            if( rsign(i,j,1) > one_half ) then

               r(:)         = brh(:,2,i,j)
               Mlocal(:,:)  = n2i(:,:,i,j,k_myend+1)
               rp(:)        = matmul( Mlocal , r )
               brh(:,2,i,j) = rp(:)
            
!            end if
         
         end do
         end do
      
      end if

   end if

   ! zet-sweep each ik line
   !
   do j = j_mysta , j_myend
   do i = i_mysta , idend   !i_myend

      do k = k_mysta, k_myend

         ! left-hand side
   
         ! center pt (k)            
         !
         !dtc = one

         !if ( dtau(i,j,k) == zero ) dtc = zero
         
         !tmp     = dtau(i,j,k)  * ep(1,3,myzone) * spr(3,i,j,k) * dzsq
         
         !ap(k,1) = one + two * tmp
         !ap(k,2) = ap(k,1)
         !ap(k,3) = ap(k,1)
         !ap(k,4) = ap(k,1)
   
         ! left pt (k-1)
         !
         !dds     = dtc * dtau(i,j,k-1) *dz2 * spr(3,i,j,k-1)

         !aw(k,1) = - tmp
         !aw(k,2) = - tmp
         !aw(k,3) = - dds * lambda(3) - tmp
         !aw(k,4) = - dds * lambda(4) - tmp
   
         ! right pt (k+1)
         !
         
         !dds     = dtc * dtau(i,j,k+1) *dz2 * spr(3,i,j,k+1)
         
         !ae(k,1) = - tmp
         !ae(k,2) = - tmp
         !ae(k,3) =   dds * lambda(3) - tmp
         !ae(k,4) =   dds * lambda(4) - tmp
   
         ! right-hand side
         !

         if ( ( rsign(i,j,k) > one_half .and. rsign(i,j,k+1) < one_half ) .or. &
              ( rsign(i,j,k) > one_half .and. rsign(i,j,k-1) < one_half )        ) then

            cycle ! aq is updated above in the code

         else

            aq(k,1:4) = rh(1:4,i,j,k)

         end if

      end do

      ! semi-explicit bc treatment
      !
      if (n == 1) then

         ! left-side bc
         ! 
         if ( mydown == mpi_proc_null .and. btype(5,myzone) /= 5 ) then
            
            k = k_mysta
         
            !if( rsign(i,j,k) > one_half ) 
            aq(k,1:4) = aq(k,1:4) - aw(k,1:4) * brh(1:4,1,i,j)
         
         end if

         ! right-side bc
         ! 
         if ( myup == mpi_proc_null .and. btype(6,myzone) /= 5 ) then
            
            k = k_myend
            
            !if( rsign(i,j,k) > one_half ) 
            aq(k,1:4) = aq(k,1:4) - ae(k,1:4) * brh(1:4,2,i,j)
         
         end if

      end if

      if ( rsign(i,j,k) > one_half .and. rsign(i,j,k+1) < one_half ) ae(k,1:4) = zero
      if ( rsign(i,j,k) > one_half .and. rsign(i,j,k-1) < one_half ) aw(k,1:4) = zero      

      ! blanking area
      !
      if (nblk /= 0) then
         
         do nb = 1, nblk
            
            if ( li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) .and. &
                 li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
            
               aw(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
               ap(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
               ae(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
            
               if (li_blk_ka(n,nb) > k_mysta) ae(li_blk_ka(n,nb)-1,1:4) = zero
               if (li_blk_kb(n,nb) < k_myend) aw(li_blk_kb(n,nb)+1,1:4) = zero
            
            end if
         
         end do
      
      end if

      ! solve linear tridiagonal eqns
      !
      klength = k_myend - k_mysta + 1
     
      !call sgtsv(klength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
      !           klength, info)
      !call sgtsv(klength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
      !           klength, info)
      !call sgtsv(klength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
      !           klength, info)
      !call sgtsv(klength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
      !           klength, info)


      !aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign_aux( i, j , va+1:vb ) 
      !aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign_aux( i, j , va+1:vb ) 
      !aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign_aux( i, j , va+1:vb ) 
      !aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign_aux( i, j , va+1:vb ) 

      !ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign_aux( i, j , va:vb   ) 
      !ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign_aux( i, j , va:vb   ) 
      !ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign_aux( i, j , va:vb   ) 
      !ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign_aux( i, j , va:vb   ) 

      !ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign_aux( i, j , va:vb-1 ) 
      !ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign_aux( i, j , va:vb-1 ) 
      !ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign_aux( i, j , va:vb-1 ) 
      !ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign_aux( i, j , va:vb-1 ) 

      !aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign_aux( i, j , va:vb   ) 
      !aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign_aux( i, j , va:vb   ) 
      !aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign_aux( i, j , va:vb   ) 
      !aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign_aux( i, j , va:vb   ) 


      aw( va+1:vb , 1 ) = aw( va+1:vb , 1 ) * rsign( i, j , va+1:vb ) 
      aw( va+1:vb , 2 ) = aw( va+1:vb , 2 ) * rsign( i, j , va+1:vb ) 
      aw( va+1:vb , 3 ) = aw( va+1:vb , 3 ) * rsign( i, j , va+1:vb ) 
      aw( va+1:vb , 4 ) = aw( va+1:vb , 4 ) * rsign( i, j , va+1:vb ) 

      ap( va:vb   , 1 ) = ap( va:vb , 1   ) * rsign( i, j , va:vb   ) 
      ap( va:vb   , 2 ) = ap( va:vb , 2   ) * rsign( i, j , va:vb   ) 
      ap( va:vb   , 3 ) = ap( va:vb , 3   ) * rsign( i, j , va:vb   ) 
      ap( va:vb   , 4 ) = ap( va:vb , 4   ) * rsign( i, j , va:vb   ) 

      ae( va:vb-1 , 1 ) = ae( va:vb-1 , 1 ) * rsign( i, j , va:vb-1 ) 
      ae( va:vb-1 , 2 ) = ae( va:vb-1 , 2 ) * rsign( i, j , va:vb-1 ) 
      ae( va:vb-1 , 3 ) = ae( va:vb-1 , 3 ) * rsign( i, j , va:vb-1 ) 
      ae( va:vb-1 , 4 ) = ae( va:vb-1 , 4 ) * rsign( i, j , va:vb-1 ) 

      aq( va:vb   , 1 ) = aq( va:vb , 1   ) * rsign( i, j , va:vb   ) 
      aq( va:vb   , 2 ) = aq( va:vb , 2   ) * rsign( i, j , va:vb   ) 
      aq( va:vb   , 3 ) = aq( va:vb , 3   ) * rsign( i, j , va:vb   ) 
      aq( va:vb   , 4 ) = aq( va:vb , 4   ) * rsign( i, j , va:vb   ) 

      ! p
      call dgtsv( klength , 1                                                               , & 
                  aw( va+1:vb , 1 ) , ap( va:vb , 1 ) , ae( va:vb-1 , 1 ) , aq( va:vb , 1 ) , &
                  klength , info                                                              &
                 )
      ! u
      call dgtsv( klength , 1                                                               , & 
                  aw( va+1:vb , 2 ) , ap( va:vb , 2 ) , ae( va:vb-1 , 2 ) , aq( va:vb , 2 ) , &
                  klength , info                                                              &
                 )
      ! v
      call dgtsv( klength , 1                                                               , & 
                  aw( va+1:vb , 3 ) , ap( va:vb , 3 ) , ae( va:vb-1 , 3 ) , aq( va:vb , 3 ) , &
                  klength , info                                                              &
                 )
      ! w
      call dgtsv( klength , 1                                                               , & 
                  aw( va+1:vb , 4 ) , ap( va:vb , 4 ) , ae( va:vb-1 , 4 ) , aq( va:vb , 4 ) , &
                  klength , info                                                              &
                 )

      ! put right-hand side back
      !
      do k = k_mysta, k_myend
        !if( rsign(i,j,k) > one_half ) 
        rh(1:4,i,j,k) = aq(k,1:4)
      end do

   end do
   end do

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
   
!      if( rsign(i,j,k) > one_half ) then
        
         Mlocal(:,:) = mc(:,:,i,j,k)
         r(:)        = rh(:,i,j,k)

         rp(:)       = matmul( Mlocal , r )
         rh(:,i,j,k) = rp(:)
      
!      end if

   end do
   end do
   end do


end subroutine rhs_diag_solver


