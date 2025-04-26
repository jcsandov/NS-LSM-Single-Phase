subroutine rhs_ghost_fluid_nodes_extp_4d( var , InteriorNodesOnly )

   use AdvectionMethods
   use InterpolationMethods

   implicit none

   ! variable to extrapolate
   ! 
   real (kind = rdf), dimension (:,il:,jl:,kl:), intent(inout) :: var
   logical, intent(in), optional :: InteriorNodesOnly

   ! local variables
   real (kind = rdf), dimension(:,:) , allocatable :: bwvec
   real (kind = rdf), dimension(1:3,1:3)  :: Gmatrix , Ginv
   logical :: OK_FLAG
   real (kind = rdf) :: search_radius
   real (kind = rdf) :: wlocal
   real (kind = rdf) :: phase

   real(kind=rdf) :: var_grad(3) , dr(3) , drnearest(3)
   real(kind=rdf) :: drmin
   real(kind=rdf) :: dx , dy , dz
   integer        :: inearest , jnearest , knearest

   real(kind=rdf) :: alpha1 , alpha2
   real(kind=rdf) :: num , denom 
   real(kind=rdf), parameter :: tol = 1E-6 

   ! neighbourhood phase counter
   integer :: nWaterNodesExtp

   ! stencil indexes  
   integer :: i1,j1,k1,i2,j2,k2 

   ! array bounds
   integer :: ml, mu

   ! Interior nodes
   integer :: i_mysta, j_mysta, k_mysta
   integer :: i_myend, j_myend, k_myend

   ! physical boundaries indexes
   integer :: ista , iend , jsta , jend , ksta , kend ! 
  
   ! Loop boundaries
   integer :: ilbound, iubound, jlbound, jubound, klbound, kubound

   ! Dummy indexes
   integer :: n,i,j,k,ii,jj,kk,v,g

   ! Number of variables
   integer :: nvars

   ! Blanking flag
   logical :: BlankingFlag

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

   ! Counter for averaging
   real(kind = rdf) :: cont
   real(kind = rdf) :: rhavg1, rhavg2 , rhavg3 , rhavg4
   real(kind = rdf) :: drh_dx(4) , drh_dy(4) , drh_dz(4)

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

   ! Interior nodes
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

   ! Physical boundaries
   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp

   ! amount of variables in the first index
   ml = lbound(var,1)
   mu = ubound(var,1)

   nvars = mu-ml+1

   allocate( bwvec( 1:3 , 1:nvars ) )

   ! Loop boundaries
   ilbound = i_mysta - 1 ; iubound = i_myend + 1 
   jlbound = j_mysta - 1 ; jubound = j_myend + 1 
   klbound = k_mysta - 1 ; kubound = k_myend + 1 

   ! Search restriction to interior nodes only
   if ( present ( InteriorNodesOnly ) ) then

      if ( InteriorNodesOnly ) then

         ! Loop boundaries
         ilbound = i_mysta ; iubound = i_myend 
         jlbound = j_mysta ; jubound = j_myend 
         klbound = k_mysta ; kubound = k_myend 
  
         if ( myback  == mpi_proc_null )  ista = i_mysta 
         if ( myleft  == mpi_proc_null )  jsta = j_mysta 
         if ( mydown  == mpi_proc_null )  ksta = k_mysta 
    
         if ( myfront == mpi_proc_null )  iend = i_myend
         if ( myright == mpi_proc_null )  jend = j_myend
         if ( myup    == mpi_proc_null )  kend = k_myend

      end if

   end if

   do g = 1 , gfmnodes_sdaf

!   do k = klbound , kubound
!   do j = jlbound , jubound
!   do i = ilbound , iubound
!
!      ! search radius for extrapolation
!      !search_radius = three * aj(i,j,k)**( -one_third )
!      !search_radius = ten * aj(i,j,k)**( -one_third )
!
!      ! If I'm in the water phase or too far from the fs, I just cycle
!      if ( rsign( phi(i,j,k) ) > one_half .or. abs( phi(i,j,k) ) > BigPhi ) cycle 
!
!      BlankingFlag = .false.
!
!      if ( nblk /= 0 ) then
!
!         do nb = 1,nblk
!
!            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. &  
!                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
!                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
!  
!               BlankingFlag = .true.
!  
!            end if
!    
!         end do  
!    
!      end if
!
!      if ( BlankingFlag ) cycle ! cycle i,j,k
!
!      ! Check if there are water or near water nodes in the neighbourhood 
!      nWaterNodesExtp = 0
!
!      inearest = i
!      jnearest = j
!      knearest = k
!
!      drmin = ten * ten
!
!
!      do kk = max( k-1 , ksta ) , min( k+1 , kend )
!      do jj = max( j-1 , jsta ) , min( j+1 , jend )
!      do ii = max( i-1 , ista ) , min( i+1 , iend )
!      
!         if ( ii == i .and. jj == j .and. kk == k ) cycle 
!
!         BlankingFlag = .false.
!
!         if ( nblk /= 0 ) then
!            
!            do nb = 1,nblk
!
!               ! I consider the nodes where the rh values extrapolated 
!               ! were to the obstacles 
!               if ( ii > li_blk_ia(1,nb) .and. ii < li_blk_ib(1,nb) .and. &  
!                    jj > li_blk_ja(1,nb) .and. jj < li_blk_jb(1,nb) .and. &
!                    kk > li_blk_ka(1,nb) .and. kk < li_blk_kb(1,nb) ) then
!                  
!                  BlankingFlag = .true.
!               
!               end if
!            
!            end do
!         
!         end if
!
!         if ( BlankingFlag ) cycle ! cycle ii,jj,kk
!
!         ! ii,jj,kk is in the water phase
!         if ( rsign( phi(ii,jj,kk) ) > one_half ) then
!            
!            nWaterNodesExtp = 1
!                    
!            dr = (/ x(i,j,k) - x(ii,jj,kk) , &
!                    y(i,j,k) - y(ii,jj,kk) , & 
!                    z(i,j,k) - z(ii,jj,kk) /)
!
!            if ( norm2(dr) < drmin ) then
!               drmin     = norm2(dr)
!               drnearest = dr
!               inearest  = ii
!               jnearest  = jj
!               knearest  = kk
!            end if
!
!         end if
!
!      end do
!      end do
!      end do
!
!      if ( nWaterNodesExtp == 0 ) cycle

      i = gfmnodes_list_sdaf( g , 1 )
      j = gfmnodes_list_sdaf( g , 2 )
      k = gfmnodes_list_sdaf( g , 3 )

      inearest = gfmnodes_list_sdaf( g , 4 )
      jnearest = gfmnodes_list_sdaf( g , 5 )
      knearest = gfmnodes_list_sdaf( g , 6 )

      drnearest = (/ x(i,j,k) - x(inearest,jnearest,knearest) , &
                     y(i,j,k) - y(inearest,jnearest,knearest) , & 
                     z(i,j,k) - z(inearest,jnearest,knearest) /)

      Gmatrix = zero
      bwvec   = zero
      Ginv    = zero
    
      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      cont   = zero
      rhavg1 = zero
      rhavg2 = zero
      rhavg3 = zero
      rhavg4 = zero
      !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

      ! Now I calculate the variable gradients at (inearest,jnearest,knearest)
      ! using a WLSM
      do kk = max( knearest-3 , ksta ) , min( knearest+3 , kend )
      do jj = max( jnearest-3 , jsta ) , min( jnearest+3 , jend )
      do ii = max( inearest-3 , ista ) , min( inearest+3 , iend )
      
         !if ( rsign ( phi(ii,jj,kk) ) < one_half          .or. & ! ii,jj,kk in the air phase
         !     norm2( phi_gradient(1:3,ii,jj,kk) ) > two  .or. & ! ∇ϕ not defined
         !     (ii == inearest .and. jj == jnearest .and. kk == knearest ) ) cycle

         if ( rsign ( phi(ii,jj,kk) ) < one_half          .or. & ! ii,jj,kk in the air phase
              (ii == inearest .and. jj == jnearest .and. kk == knearest ) ) cycle

         BlankingFlag = .false.

         if ( nblk /= 0 ) then

            do nb = 1,nblk

               ! I consider the nodes where the rh values extrapolated 
               ! were to the obstacles 
               if ( ii > li_blk_ia(1,nb) .and. ii < li_blk_ib(1,nb) .and. &  
                    jj > li_blk_ja(1,nb) .and. jj < li_blk_jb(1,nb) .and. &
                    kk > li_blk_ka(1,nb) .and. kk < li_blk_kb(1,nb) ) then

                  BlankingFlag = .true.
            
               end if
            
            end do

         else if ( nblke /= 0 ) then

            do nb = 1,nblke

               ! I consider the nodes where the rh values extrapolated 
               ! were to the obstacles 
               if ( ii > le_blk_ia(1,nb) .and. ii < le_blk_ib(1,nb) .and. &  
                    jj > le_blk_ja(1,nb) .and. jj < le_blk_jb(1,nb) .and. &
                    kk > le_blk_ka(1,nb) .and. kk < le_blk_kb(1,nb) ) then

                  BlankingFlag = .true.
            
               end if
            
            end do

         end if

         if ( BlankingFlag ) cycle ! cycle ii,jj,kk

         dx = x( inearest , jnearest , knearest ) - x(ii,jj,kk)
         dy = y( inearest , jnearest , knearest ) - y(ii,jj,kk)
         dz = z( inearest , jnearest , knearest ) - z(ii,jj,kk)

         dr = (/dx,dy,dz/)

         num   = abs( dot_product( dr/norm2(dr) , phi_gradient(:,inearest,jnearest,knearest) ) )!**(three/two) 
         denom = one !norm2(dr)**(three/two) + eps_denom

         wlocal = num / denom

         !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
         cont   = cont   + one 
         rhavg1 = rhavg1 + var(1,ii,jj,kk)
         rhavg2 = rhavg2 + var(2,ii,jj,kk)
         rhavg3 = rhavg3 + var(3,ii,jj,kk)
         rhavg4 = rhavg4 + var(4,ii,jj,kk)
         !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

         Gmatrix(1,1) = Gmatrix(1,1) + ( dx*dx ) * ( wlocal**two ) 
         Gmatrix(1,2) = Gmatrix(1,2) + ( dx*dy ) * ( wlocal**two ) 
         Gmatrix(1,3) = Gmatrix(1,3) + ( dx*dz ) * ( wlocal**two ) 

         Gmatrix(2,1) = Gmatrix(2,1) + ( dy*dx ) * ( wlocal**two ) 
         Gmatrix(2,2) = Gmatrix(2,2) + ( dy*dy ) * ( wlocal**two ) 
         Gmatrix(2,3) = Gmatrix(2,3) + ( dy*dz ) * ( wlocal**two ) 

         Gmatrix(3,1) = Gmatrix(3,1) + ( dz*dx ) * ( wlocal**two ) 
         Gmatrix(3,2) = Gmatrix(3,2) + ( dz*dy ) * ( wlocal**two ) 
         Gmatrix(3,3) = Gmatrix(3,3) + ( dz*dz ) * ( wlocal**two ) 

         bwvec(1,1:nvars) = bwvec(1,1:nvars) + ( wlocal**two ) * dx *   &
                            ( var(1:nvars,inearest,jnearest,knearest) - &
                              var(1:nvars,ii      ,jj      ,kk      )      )  

         bwvec(2,1:nvars) = bwvec(2,1:nvars) + ( wlocal**two ) * dy *   &
                            ( var(1:nvars,inearest,jnearest,knearest) - &
                              var(1:nvars,ii      ,jj      ,kk      )      )  

         bwvec(3,1:nvars) = bwvec(3,1:nvars) + ( wlocal**two ) * dz *   &
                            ( var(1:nvars,inearest,jnearest,knearest) - &
                              var(1:nvars,ii      ,jj      ,kk      )      )  
         

      end do
      end do
      end do 

      rhavg1 = ( rhavg1 + var(1,inearest,jnearest,knearest) ) / ( cont + 1 )
      rhavg2 = ( rhavg2 + var(2,inearest,jnearest,knearest) ) / ( cont + 1 )
      rhavg3 = ( rhavg3 + var(3,inearest,jnearest,knearest) ) / ( cont + 1 )
      rhavg4 = ( rhavg4 + var(4,inearest,jnearest,knearest) ) / ( cont + 1 )

      ! Gmatrix inversion
      call M33INV( Gmatrix , Ginv , OK_FLAG )

      if ( OK_FLAG ) then
         
         do v = 1,nvars

            var_grad = matmul( Ginv , bwvec(1:3,v) )
            var(v,i,j,k) = var(v,inearest,jnearest,knearest) + dot_product(var_grad , drnearest)

            drh_dx(v) = var_grad(1)
            drh_dy(v) = var_grad(2)
            drh_dz(v) = var_grad(3)

         end do

      else
         
         print *,'Problem inverting the matrix in rhs_ghost_fluid_nodes_extp_4d at i,j,k: ',i,j,k               
         var(:,i,j,k) = var(:,inearest,jnearest,knearest)

      end if

!   end do
!   end do
!   end do
 
   end do 

   deallocate( bwvec )

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

end subroutine rhs_ghost_fluid_nodes_extp_4d



