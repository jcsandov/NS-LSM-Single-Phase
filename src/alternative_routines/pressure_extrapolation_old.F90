subroutine pressure_extrapolation_old ( )

    real ( kind = rdf ) :: pExtp

   ! local variables

   integer :: i,j,k
   integer :: ii, jj, kk
   integer :: isearch, jsearch, ksearch
   integer :: iic, jic, kic
   integer :: ioc, joc, koc
   integer :: iaux, jaux, kaux

   real(kind=rdf) :: PhiWater, pWater
   real(kind=rdf) :: xIntpNode, yIntpNode, zIntpNode
   real(kind=rdf) :: dWater, dAir, dTotal
   real(kind=rdf) :: pAux, LocalWeight, TotalWeight
   real(kind=rdf) :: dx   
   real(kind=rdf) :: extp_threshold

   integer :: nWaterNodesExtp

   extp_threshold  = 1.0E-08


   !do k = k_mysta, k_myend
   !do j = j_mysta, j_myend
   !do i = i_mysta, i_myend
            
   do k = k_mysta-1, k_myend+1
   do j = j_mysta-1, j_myend+1
   do i = i_mysta-1, i_myend+1

      ! ----------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
      ! ----------------------------------------------------------------------------
      ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
      ! phase (rsign = 1), then this if is true. This condition identifies air nodes
      ! next to the free-surface
      ! ----------------------------------------------------------------------------

      ! If the node is closer to the free-surface than the machine-precision   
      if ( abs( phi(i,j,k) ) < eps_sims ) then

         q(1,i,j,k) = zero
         cycle

      end if

      !write(*,'(A,I0,A,3(I0,A))') &
      !'myid = ', myid , ', ijk = ',i,',',j,',',k,' '

      iaux = i
      jaux = j
      kaux = k

      if ( iaux <= i_mysta-1 )  iaux = iaux + 1
      if ( jaux <= j_mysta-1 )  jaux = jaux + 1
      if ( kaux <= k_mysta-1 )  kaux = kaux + 1
      
      if ( iaux >= i_myend+1 )  iaux = iaux - 1
      if ( jaux >= j_myend+1 )  jaux = jaux - 1
      if ( kaux >= k_myend+1 )  kaux = kaux - 1


      dx = aj(iaux,jaux,kaux)**( -one_third )
      nWaterNodesExtp = 0

      if ( rsign(i,j,k) < one_half .and. abs( phi(i,j,k) ) < radius_lsqm * dx) then ! Air-phase
      
         ! Replace this loop for an any or a searchloop with an exit
         ! as nWaterNeighbours is not important itself   
         
         searchloop2 : do kk = max( k-1 , k_mysta-1 ) , min( k+1 , k_myend+1 )
                       do jj = max( j-1 , j_mysta-1 ) , min( j+1 , j_myend+1 )
                       do ii = max( i-1 , i_mysta-1 ) , min( i+1 , i_myend+1 )
      
            if ( phi (ii,jj,kk) > -extp_threshold ) then

            !if ( rsign (ii,jj,kk) > one_half ) then
            
               nWaterNodesExtp = 1
            
               exit searchloop2
            
            end if

         end do
         end do
         end do searchloop2
   
      end if


      if( nWaterNodesExtp > 0 ) then

         ! variables initialisation
      
         pAux        = zero
         TotalWeight = zero

      
         do isearch = -1,1
         do jsearch = -1,1
         do ksearch = -1,1
      

            ! Central node is skipped
            if ( isearch /= 0 .or. jsearch /= 0 .or. ksearch /= 0 ) then
      
               iic = isearch
               jic = jsearch
               kic = ksearch
      
               ! Trucated search area near the boundaries
               !if (myback  == mpi_proc_null .and. i + isearch <= i_mysta-1 )  iic = 0
               !if (myleft  == mpi_proc_null .and. j + jsearch <= j_mysta-1 )  jic = 0
               !if (mydown  == mpi_proc_null .and. k + ksearch <= k_mysta-1 )  kic = 0
               !
               !if (myfront == mpi_proc_null .and. i + isearch >= i_myend+1 )  iic = 0
               !if (myright == mpi_proc_null .and. j + jsearch >= j_myend+1 )  jic = 0
               !if (myup    == mpi_proc_null .and. k + ksearch >= k_myend+1 )  kic = 0
      

               if ( i + isearch <= i_mysta-1 )  iic = 0
               if ( j + jsearch <= j_mysta-1 )  jic = 0
               if ( k + ksearch <= k_mysta-1 )  kic = 0
               
               if ( i + isearch >= i_myend+1 )  iic = 0
               if ( j + jsearch >= j_myend+1 )  jic = 0
               if ( k + ksearch >= k_myend+1 )  kic = 0
      
               ! I update the inner cell indexes with the actual indexes in the domain
               ! ic = inner cell
      
               iic = i + iic
               jic = j + jic
               kic = k + kic

               ! If the neighbour is in the water phase, or closer to the free-surface
               ! than a defined theshold, then, the extrapolation is performed
      
               if ( phi (iic,jic,kic) > extp_threshold  ) then
      
               !if ( phi (iic,jic,kic) > -eps_sims ) then

                  ! Indexes outer cell
                  ioc = iic + isearch
                  joc = jic + jsearch
                  koc = kic + ksearch

                  ! Trucated search area near the boundaries
                  !if (myback  == mpi_proc_null .and. ioc <= i_mysta-1 )  ioc = iic
                  !if (myleft  == mpi_proc_null .and. joc <= j_mysta-1 )  joc = jic
                  !if (mydown  == mpi_proc_null .and. koc <= k_mysta-1 )  koc = kic
                  
                  !if (myfront == mpi_proc_null .and. ioc >= i_myend+1 )  ioc = iic
                  !if (myright == mpi_proc_null .and. joc >= j_myend+1 )  joc = jic
                  !if (myup    == mpi_proc_null .and. koc >= k_myend+1 )  koc = kic

                  if ( ioc <= i_mysta-1 )  ioc = iic
                  if ( joc <= j_mysta-1 )  joc = jic
                  if ( koc <= k_mysta-1 )  koc = kic
                  
                  if ( ioc >= i_myend+1 )  ioc = iic
                  if ( joc >= j_myend+1 )  joc = jic
                  if ( koc >= k_myend+1 )  koc = kic

                  !write(*,'(A,I0,A,3(I0,A),A,3(I0,A))') &
                  !'myid = ', myid , ' - ijk inner cell = ',iic,',',jic,',',kic,' ', ' - ijk outter cell = ',ioc,',',joc,',',koc,' ' 
         
                  PhiWater = one_half * ( phi( iic , jic , kic )  + phi( ioc , joc , koc ) )
      
                  !if ( abs( phi( iic , jic , kic ) ) < eps_sims ) q(1 , iic , jic , kic ) = zero
                  !if ( abs( phi( ioc , joc , koc ) ) < eps_sims ) q(1 , ioc , joc , koc ) = zero

                  pWater   = one_half * ( q(1 , iic , jic , kic ) + q(1 , ioc , joc , koc ) )
      
                  !if ( rsign( iic , jic , kic ) < one_half .or. rsign( ioc , joc , koc ) < one_half .or. PhiWater < eps_sims ) cycle
      
                  xIntpNode = one_half * (   x( iic , jic , kic ) + x( ioc , joc , koc ) )
                  yIntpNode = one_half * (   y( iic , jic , kic ) + y( ioc , joc , koc ) )
                  zIntpNode = one_half * (   z( iic , jic , kic ) + z( ioc , joc , koc ) )
      
                  ! Distance between node i,j,k and the current node where the extrapolation 
                  ! is performed
      
                  dTotal = norm2( (/ x(i,j,k) - xIntpNode , &
                                     y(i,j,k) - yIntpNode , &
                                     z(i,j,k) - zIntpNode    /) )

                  ! Avoid division by zero
                  if ( abs( PhiWater ) < eps_sims) cycle

                  dWater = dTotal / ( one + abs( phi(i,j,k) ) / abs( PhiWater ) )
                  dAir   = dTotal - dWater
      
                  LocalWeight = one / dTotal**two

                  pAux        = pAux - LocalWeight * ( dAir / dWater * pWater )
                  TotalWeight = TotalWeight + LocalWeight

      
               else if ( abs( phi (iic,jic,kic) ) <= extp_threshold  ) then

                  ! I skip the node from the inner cell

                  ! Indexes outer cell
                  ioc = iic + isearch
                  joc = jic + jsearch
                  koc = kic + ksearch

                  ! Trucated search area near the boundaries
                  !if (myback  == mpi_proc_null .and. ioc <= i_mysta-1 )  ioc = iic
                  !if (myleft  == mpi_proc_null .and. joc <= j_mysta-1 )  joc = jic
                  !if (mydown  == mpi_proc_null .and. koc <= k_mysta-1 )  koc = kic
                  
                  !if (myfront == mpi_proc_null .and. ioc >= i_myend+1 )  ioc = iic
                  !if (myright == mpi_proc_null .and. joc >= j_myend+1 )  joc = jic
                  !if (myup    == mpi_proc_null .and. koc >= k_myend+1 )  koc = kic

                  if ( ioc <= i_mysta-1 )  ioc = iic
                  if ( joc <= j_mysta-1 )  joc = jic
                  if ( koc <= k_mysta-1 )  koc = kic
                  
                  if ( ioc >= i_myend+1 )  ioc = iic
                  if ( joc >= j_myend+1 )  joc = jic
                  if ( koc >= k_myend+1 )  koc = kic

                  if ( rsign(ioc,joc,koc) > one_half ) then

                     PhiWater = phi( ioc , joc , koc )
                     pWater   = q(1 , ioc , joc , koc )

                     dTotal = norm2( (/ x(i,j,k) - x(ioc,joc,koc) , &
                                        y(i,j,k) - y(ioc,joc,koc) , &
                                        z(i,j,k) - z(ioc,joc,koc)    /) )


                     ! Avoid division by zero
                     if ( abs( PhiWater ) < eps_sims ) cycle

                     dWater = dTotal / ( one + abs( phi(i,j,k) ) / abs( PhiWater ) )
                     dAir   = dTotal - dWater
      
                     LocalWeight = one / dTotal**two

                     pAux        = pAux - LocalWeight * ( dAir / dWater * pWater )
                     TotalWeight = TotalWeight + LocalWeight


                  end if

               end if
            
            end if
      
      
         end do
         end do
         end do
      

         if ( TotalWeight > zero ) q(1,i,j,k) = pAux / TotalWeight
   
      end if
   
   end do
   end do
   end do


end subroutine pressure_extrapolation_old