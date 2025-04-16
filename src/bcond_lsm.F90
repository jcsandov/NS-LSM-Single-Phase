
  subroutine bcond_lsm (phi)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !extrapolo bordes a la funcion phi, del level set
  ! - - - - - - - - - - - -

  implicit none
   
  real (kind = rdf) :: inter1, inter2
  real (kind = rdf) :: inter1Blk, inter2Blk

  real (kind = rdf), parameter :: itp1stOrd_a = one , itp1stOrd_b = zero
  real (kind = rdf), parameter :: itp2ndOrd_a = four/three, itp2ndOrd_b = -one/three
  real (kind = rdf) :: exsign

  real (kind = rdf) :: itpa , itpb ! first and second nodes

  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(inout) :: phi

  integer :: i_mysta,i_myend
  integer :: j_mysta,j_myend
  integer :: k_mysta,k_myend

  integer :: i,j,k, b , n

  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! UPWIND extrapolation is the default option to
  ! extrapolate phi

  !inter1 =  4.0_rdf/3.0_rdf
  !inter2 = -1.0_rdf/3.0_rdf
  
  ! If user enters 0, the code uses CD extrapolation

  if (bc_extrapolation_order == 0) then

    inter1 = 2.0_rdf
    inter2 = -1.0_rdf

  end if

  ! case = 0 : absorption
  ! case = 1 : extrapolation
  ! case = 2 : advection (this routines do nothing in that case)

   if(myback == mpi_proc_null) then
    
      b = 1
      i = i_mysta

      csi1: select case( btype_lsm ( b , myzone ) )
    
         case(0)
            do k = k_mysta , k_myend
            do j = j_mysta , j_myend                                
               phi(i,j,k) = phi(i+1,j,k)
            end do
            end do 

         case(1)
            do k = k_mysta , k_myend
            do j = j_mysta , j_myend

               ! If |ϕ|< Bigϕ ----> exsign = 1. exsign = 0, otherwise
               exsign = ( one + sign( one , abs( phi(i+2,j,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i+1,j,k) + &
                            itpb * phi(i+2,j,k)
          end do
          end do
        
        case(2)
  
      end select csi1

   end if        


   if(myfront == mpi_proc_null) then
  
      b = 2
      i = i_myend

      csi2: select case(btype_lsm(b,myzone))
       
         case(0)
            
            do k = k_mysta,k_myend
            do j = j_mysta, j_myend
            
               phi(i,j,k) = phi(i-1,j,k)
                                            
            end do
            end do
   
         case(1)
            
            do k = k_mysta , k_myend
            do j = j_mysta , j_myend
               
               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i-1,j,k) + &
                            itpb * phi(i-2,j,k)
                                            
            end do
            end do
   
         case(2)
              
            do k = k_mysta,k_myend
            do j = j_mysta, j_myend
               
               phi(i,j,k) = two * phi(i-1,j,k) - &
                            one * phi(i-2,j,k)
                                            
            end do
            end do


      end select csi2
       
   end if       
  
   if(myleft == mpi_proc_null) then
  
      b = 3
      j = j_mysta

      eta1: select case(btype_lsm(b,myzone))
    
         case(0)
            do k = k_mysta , k_myend                  
            do i = i_mysta , i_myend
               phi(i,j,k) = phi(i,j+1,k) 
            end do
            end do
 
         case(1)
            do k = k_mysta , k_myend                  
            do i = i_mysta , i_myend

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j+2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j+1,k) + &
                            itpb * phi(i,j+2,k)

            end do
            end do
         
         case(2)
    
      end select eta1
   
   end if    

   if(myright == mpi_proc_null) then 
  
      b = 4
      j = j_myend
    
      eta2: select case(btype_lsm(b,myzone))
    
         case(0)
            do k = k_mysta , k_myend                  
            do i = i_mysta , i_myend
               phi(i,j,k) = phi(i,j-1,k) 
            end do
            end do
 
         case(1)
            do k = k_mysta , k_myend                  
            do i = i_mysta , i_myend

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j-2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j-1,k) + &
                            itpb * phi(i,j-2,k)

            end do
            end do
         
         case(2)
    
      end select eta2

   end if

   if(mydown == mpi_proc_null) then
    
      b = 5
      k = k_mysta

      zet1: select case(btype_lsm(b,myzone))
    
         case(0)
      
            do j = j_mysta , j_myend                  
            do i = i_mysta , i_myend
               phi(i,j,k) = phi(i,j,k + 1)
                                      
            end do
            end do

         case(1)

            do j = j_mysta , j_myend                  
            do i = i_mysta , i_myend
          
               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j,k+2) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j,k+1) + &
                            itpb * phi(i,j,k+2)
                              
            end do
            end do
    
         case(2)
    
      end select zet1
   
   end if   

   if(myup == mpi_proc_null) then
  
      b = 6
      k = k_myend

      zet2: select case(btype_lsm(b,myzone))
    
         case(0)
      
            do j = j_mysta , j_myend                  
            do i = i_mysta , i_myend
               phi(i,j,k) = phi(i,j,k - 1)
                                      
            end do
            end do

         case(1)

            do j = j_mysta , j_myend                  
            do i = i_mysta , i_myend
          
               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j,k-2) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j,k-1) + &
                            itpb * phi(i,j,k-2)
                                
            end do
            end do
    
         case(2)

      end select zet2
    
   end if   
 
   !=====================================================================================

   ! blanking area
   !

   n = 1

   if (nblk /= 0) then
    
      do nb = 1, nblk

         do k = max( k_mysta , li_blk_ka(n,nb)+2 ) , min(k_myend , li_blk_kb(n,nb)-2 )
         do j = max( j_mysta , li_blk_ja(n,nb)+2 ) , min(j_myend , li_blk_jb(n,nb)-2 )
         do i = max( i_mysta , li_blk_ia(n,nb)+2 ) , min(i_myend , li_blk_ib(n,nb)-2 )
            phi(i,j,k) = zero
         end do
         end do
         end do
  
         !=============================================================================== 
         ! ξ - direction
         !=============================================================================== 
         !
         ! Front side of an obstacle
         !
         !             i = li_blk_ib(1,nb)
         !             |##############
         !             |##          ##
         !   o----o----o## OBSTACLE ##
         !  i-2  i-1   |##          ##     
         !             |##############
         !   

         i = li_blk_ia(n,nb)
        
         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 ) then
                
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 
            do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i-1,j,k) + &
                            itpb * phi(i-2,j,k)

               phi(i+1,j,k) = phi(i,j,k) 

            end do
            end do

         end if    
      
         ! 
         ! Rear side of an obstacle
         !  
         !                i = li_blk_ib(1,nb)
         !  ##############|        
         !  ##          ##|        
         !  ## OBSTACLE ##o----o----o
         !  ##          ##|   i+1  i+2       
         !  ##############|
         !   
      
         i = li_blk_ib(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 ) then
                
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) )
            do j = max( j_mysta , li_blk_ja(n,nb) ) , min( j_myend , li_blk_jb(n,nb) )

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i+1,j,k) + &
                            itpb * phi(i+2,j,k)

               phi(i-1,j,k) = phi(i,j,k) 
        
            end do
            end do
       
         end if
       
         !=============================================================================== 
         ! η - direction
         !=============================================================================== 
              
         ! left side of an obstacle
         j = li_blk_ja(n,nb)
       
         if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then
                    
            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) )
            do i = max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) )

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j-2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j-1,k) + &
                            itpb * phi(i,j-2,k)

               phi(i,j+1,k) = phi(i,j,k) 
        
            end do
            end do     
        
         end if    
       
         j = li_blk_jb(n,nb)
        
         ! right side of an obstacle
         if ( blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then
                
            do k =max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) )
            do i =max( i_mysta , li_blk_ia(n,nb) ) , min( i_myend , li_blk_ib(n,nb) )

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i,j+2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i,j+1,k) + &
                            itpb * phi(i,j+2,k)

               phi(i,j-1,k) = phi(i,j,k) 

            end do
            end do     
        
         end if

         !=============================================================================== 
         ! Corners
         !=============================================================================== 

         i = li_blk_ia(n,nb)
         j = li_blk_ja(n,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j-2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i-1,j-1,k) + &
                            itpb * phi(i-2,j-2,k)

               phi(i+1,j+1,k) = phi(i,j,k) 

            end do

         end if

         i = li_blk_ib(n,nb)
         j = li_blk_ja(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j-2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i+1,j-1,k) + &
                            itpb * phi(i+2,j-2,k)

               phi(i-1,j+1,k) = phi(i,j,k) 

            end do

         end if

         i = li_blk_ib(n,nb)
         j = li_blk_jb(n,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i+2,j+2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i+1,j+1,k) + &
                            itpb * phi(i+2,j+2,k)

               phi(i-1,j-1,k) = phi(i,j,k) 

            end do

         end if

         i = li_blk_ia(n,nb)
         j = li_blk_jb(n,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(n,nb) ) , min( k_myend , li_blk_kb(n,nb) ) 

               ! either 1 or 0
               exsign = ( one + sign( one , abs( phi(i-2,j+2,k) )- BigPhi ) )/two

               itpa = itp2ndOrd_a + exsign * ( itp1stOrd_a - itp2ndOrd_a )
               itpb = itp2ndOrd_b + exsign * ( itp1stOrd_b - itp2ndOrd_b )

               phi(i,j,k) = itpa * phi(i-1,j+1,k) + &
                            itpb * phi(i-2,j+2,k)

               phi(i+1,j-1,k) = phi(i,j,k) 

            end do

         end if


        !=============================================================================== 
        ! ζ - direction
        !=============================================================================== 

        !k = li_blk_kb(1,nb)
        
        !if (k <= k_myend) then
                
          !do j = li_blk_ja(1,nb) , li_blk_jb(1,nb)
          !do i = li_blk_ia(1,nb) , li_blk_ib(1,nb)

          !  phi(i,j,k) = inter1 * phi(i,j,min(k+1,k_myend)) + &
          !               inter2 * phi(i,j,min(k+2,k_myend))
          !end do
          !end do     
        
        !end if  
       
        !k = li_blk_ka(1,nb)
        
        !if (k >= k_mysta) then
                
          !do j = li_blk_ja(1,nb) , li_blk_jb(1,nb)
          !do i = li_blk_ia(1,nb) , li_blk_ib(1,nb)

            !phi(i,j,k) = inter1 * phi(i,j,max(k-1,k_mysta)) + &
            !             inter2 * phi(i,j,max(k-2,k_mysta))
          !end do
          !end do

        !end if       
                   
      end do
   
   end if





end subroutine bcond_lsm


