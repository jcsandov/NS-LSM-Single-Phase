subroutine get_gfm_sdaf_list()

   ! Interior nodes
   integer :: i_mysta , j_mysta , k_mysta
   integer :: i_myend , j_myend , k_myend

   ! physical boundaries indexes
   integer :: ista , iend , jsta , jend , ksta , kend ! 
  
   ! Loop boundaries
   integer :: ilbound, iubound, jlbound, jubound, klbound, kubound

   ! Nearest water node to the current air node
   integer :: inearest , jnearest , knearest

   ! Dummy indexes
   integer :: i,j,k,ii,jj,kk,cont

   ! neighbourhood phase counter
   integer :: nWaterNodesExtp

   logical :: BlankingFlag
   real(kind=rdf) :: drmin
   real(kind=rdf) :: dr(3)

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

   ilbound = i_mysta ; iubound = i_myend 
   jlbound = j_mysta ; jubound = j_myend 
   klbound = k_mysta ; kubound = k_myend 
  
   if ( myback  == mpi_proc_null )  ista = i_mysta 
   if ( myleft  == mpi_proc_null )  jsta = j_mysta 
   if ( mydown  == mpi_proc_null )  ksta = k_mysta 
    
   if ( myfront == mpi_proc_null )  iend = i_myend
   if ( myright == mpi_proc_null )  jend = j_myend
   if ( myup    == mpi_proc_null )  kend = k_myend

   gfmnodes_sdaf = 0
   if ( allocated( gfmnodes_list_sdaf ) ) deallocate ( gfmnodes_list_sdaf )

   do k = klbound , kubound
   do j = jlbound , jubound
   do i = ilbound , iubound

      ! search radius for extrapolation
      !search_radius = three * aj(i,j,k)**( -one_third )
      !search_radius = ten * aj(i,j,k)**( -one_third )

      ! If I'm in the water phase or too far from the fs, I just cycle
      if ( rsign( phi(i,j,k) ) > one_half .or. abs( phi(i,j,k) ) > BigPhi ) cycle 

      BlankingFlag = .false.

      if ( nblk /= 0 ) then

         do nb = 1,nblk

            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. &  
                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
  
               BlankingFlag = .true.
  
            end if
    
         end do  
    
      end if

      if ( BlankingFlag ) cycle ! cycle i,j,k

      ! Check if there are water or near water nodes in the neighbourhood 
      nWaterNodesExtp = 0

      inearest = i
      jnearest = j
      knearest = k

      drmin = ten * ten

      do kk = max( k-1 , ksta ) , min( k+1 , kend )
      do jj = max( j-1 , jsta ) , min( j+1 , jend )
      do ii = max( i-1 , ista ) , min( i+1 , iend )
      
         if ( ii == i .and. jj == j .and. kk == k ) cycle 

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
         
         end if

         if ( BlankingFlag ) cycle ! cycle ii,jj,kk

         ! ii,jj,kk is in the water phase
         if ( rsign( phi(ii,jj,kk) ) > one_half ) then
            
            nWaterNodesExtp = 1
                    
            dr = (/ x(i,j,k) - x(ii,jj,kk) , &
                    y(i,j,k) - y(ii,jj,kk) , & 
                    z(i,j,k) - z(ii,jj,kk) /)

            if ( norm2(dr) < drmin ) then
               drmin     = norm2(dr)
               inearest  = ii
               jnearest  = jj
               knearest  = kk
            end if

         end if

      end do
      end do
      end do

      if ( nWaterNodesExtp > 0 ) gfmnodes_sdaf = gfmnodes_sdaf + 1

   end do
   end do
   end do

   allocate ( gfmnodes_list_sdaf( gfmnodes_sdaf , 6 ) )

   cont = 1 

   do k = klbound , kubound
   do j = jlbound , jubound
   do i = ilbound , iubound

      ! search radius for extrapolation
      !search_radius = three * aj(i,j,k)**( -one_third )
      !search_radius = ten * aj(i,j,k)**( -one_third )

      ! If I'm in the water phase or too far from the fs, I just cycle
      if ( rsign( phi(i,j,k) ) > one_half .or. abs( phi(i,j,k) ) > BigPhi ) cycle 

      BlankingFlag = .false.

      if ( nblk /= 0 ) then

         do nb = 1,nblk

            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. &  
                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
  
               BlankingFlag = .true.
  
            end if
    
         end do  
    
      end if

      if ( BlankingFlag ) cycle ! cycle i,j,k

      ! Check if there are water or near water nodes in the neighbourhood 
      nWaterNodesExtp = 0

      drmin = ten * ten

      do kk = max( k-1 , ksta ) , min( k+1 , kend )
      do jj = max( j-1 , jsta ) , min( j+1 , jend )
      do ii = max( i-1 , ista ) , min( i+1 , iend )
      
         if ( ii == i .and. jj == j .and. kk == k ) cycle 

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
         
         end if

         if ( BlankingFlag ) cycle ! cycle ii,jj,kk

         ! ii,jj,kk is in the water phase
         if ( rsign( phi(ii,jj,kk) ) > one_half ) then
            
            nWaterNodesExtp = 1
                    
            dr = (/ x(i,j,k) - x(ii,jj,kk) , &
                    y(i,j,k) - y(ii,jj,kk) , & 
                    z(i,j,k) - z(ii,jj,kk) /)

            if ( norm2(dr) < drmin ) then
               drmin     = norm2(dr)
               inearest  = ii
               jnearest  = jj
               knearest  = kk
            end if

         end if

      end do
      end do
      end do

      if ( nWaterNodesExtp > 0 ) then

         ! Point to be extrapolated in the air phase
         gfmnodes_list_sdaf(cont,1) = i
         gfmnodes_list_sdaf(cont,2) = j
         gfmnodes_list_sdaf(cont,3) = k

         ! Nearest point in the water phase
         gfmnodes_list_sdaf(cont,4) = inearest
         gfmnodes_list_sdaf(cont,5) = jnearest
         gfmnodes_list_sdaf(cont,6) = knearest

         cont = cont + 1

      end if

   end do
   end do
   end do


end subroutine get_gfm_sdaf_list