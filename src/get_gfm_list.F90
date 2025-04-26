subroutine get_gfm_list()

   gfmnodes = 0
   if ( allocated( gfmnodes_list ) ) deallocate ( gfmnodes_list )

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , i_myend
  
      ! water-phase or nodes too far from the free surface
      if ( rsign( phi(i,j,k) ) > one_half .or. abs( phi(i,j,k) ) > BigPhi ) cycle
  
      BlankingFlag = .false.
   
      if ( nblk/=0 ) then
         
         do nb = 1,nblk
        
            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. & 
                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
  
               BlankingFlag = .true.
        
            end if
  
         end do  
    
      end if
        
      ! ----------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
      ! ----------------------------------------------------------------------------
      ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
      ! phase (rsign = 1), then this if is true. This condition identifies air nodes
      ! next to the free-surface
      ! ----------------------------------------------------------------------------
  
      ! Nodes inside the blanking region
      if ( BlankingFlag ) cycle 
  
      nWaterNeighbours = zero
  
      search_water_neighbour1 : &   
      do kk = k-1 , k+1
      do jj = j-1 , j+1
      do ii = i-1 , i+1
        
         !nWaterNeighbours = nWaterNeighbours + rsign(ii,jj,kk)
           
         if ( rsign( phi(ii,jj,kk) ) > one_half ) then
  
            nWaterNeighbours = one 
            exit search_water_neighbour1
  
         end if
  
      end do
      end do
      end do search_water_neighbour1
  
      if ( nWaterNeighbours < one_half ) cycle  
    
      ! Update the gfmnodes list
      gfmnodes = gfmnodes + 1
    
   end do
   end do
   end do 

   allocate ( gfmnodes_list( gfmnodes , 3 ) )

   cont = 1

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , i_myend
  
      ! water-phase or nodes too far from the free surface
      if ( rsign( phi(i,j,k) ) > one_half .or. abs( phi(i,j,k) ) > BigPhi ) cycle
  
      BlankingFlag = .false.
  
      if ( nblk/=0 ) then
         do nb = 1,nblk
        
            if ( i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) .and. & 
                 j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) .and. &
                 k >= li_blk_ka(1,nb) .and. k <= li_blk_kb(1,nb) ) then
  
               BlankingFlag = .true.
        
            end if
  
         end do  
    
      end if
        
      ! ----------------------------------------------------------------------------
      ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
      ! ----------------------------------------------------------------------------
      ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
      ! phase (rsign = 1), then this if is true. This condition identifies air nodes
      ! next to the free-surface
      ! ----------------------------------------------------------------------------
  
      ! Nodes inside the blanking region
      if ( BlankingFlag ) cycle 
  
      nWaterNeighbours = zero
  
      search_water_neighbour2 : &   
      do kk = k-1 , k+1
      do jj = j-1 , j+1
      do ii = i-1 , i+1
        
      !nWaterNeighbours = nWaterNeighbours + rsign(ii,jj,kk)
           
         if ( rsign( phi(ii,jj,kk) ) > one_half ) then
  
            nWaterNeighbours = one 
            exit search_water_neighbour2
  
         end if
  
      end do
      end do
      end do search_water_neighbour2
  
      if ( nWaterNeighbours < one_half ) cycle  
    
      ! Update the gfmnodes list
      gfmnodes_list(cont,1) = i
      gfmnodes_list(cont,2) = j
      gfmnodes_list(cont,3) = k
    
      cont = cont + 1

   end do
   end do
   end do   
  
end subroutine get_gfm_list