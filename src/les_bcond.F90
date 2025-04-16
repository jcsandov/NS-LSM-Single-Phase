
subroutine les_bcond ()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! boundary conditions for the eddy viscosity in des
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  integer :: i, j, k, b
  
  integer :: i_mysta
  integer :: j_mysta
  integer :: k_mysta

  integer :: i_myend
  integer :: j_myend
  integer :: k_myend

  real (kind = rdf) :: reinol

  real (kind = rdf) :: x1,y1,z1,x2,y2,z2,xb,yb,zb
  real (kind = rdf) :: alpha1 , alpha2


  ! make definition of i_mysta etc., consistent with
  ! all other uses, e.g., solver_daf & rhs_kw_bcond
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp


  reinol = one / ren

  ! boundary in csi-direction (i = 1)
  !
  if ( myback == mpi_proc_null ) then
      
      b = 1
      i = i_mysta
     
     csi1: select case( btype ( b , myzone ) )
     
     case(0)
     
     case(1) ! solid
        do k = k_mysta , k_myend
        do j = j_mysta , j_myend
           xnut(i,j,k) = zero
        end do
        end do
 
     case(2:3)
        do k = k_mysta , k_myend
        do j = j_mysta , j_myend
           xnut(i,j,k) = sa(1,b) * xnut(i+1,j,k) + &
                         sb(1,b) * xnut(i+2,j,k)
        end do
        end do
     case(4) ! inlet
     case(5)
        do k = k_mysta , k_myend
        do j = j_mysta , j_myend
           xnut(i,j,k) = (one + rat(j,k,b)) * xnut(i+1,j,k) - &
                                rat(j,k,b)  * xnut(i+2,j,k)
        end do
        end do
     case(6)
        do k = k_mysta , k_myend
        do j = j_mysta , j_myend
           xnut(i,j,k) = xnut(i+1,j,k)
        end do
        end do

     case(9)
        do k = k_mysta , k_myend
        do j = j_mysta , j_myend

           xnut(i,j,k) = two * xnut(i+1,j,k) - &
                         one * xnut(i+2,j,k)
        
        end do
        end do

     end select csi1
  end if

  ! boundary in csi-direction (i = imax)
  ! 
  if ( myfront == mpi_proc_null ) then
     
     b = 2
     i = i_myend
     
     csi2: select case( btype ( b , myzone ) )
     
     case(0)
     case(1) ! solid
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           xnut(i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           xnut(i,j,k) = sa(1,b) * xnut(i-1,j,k) + &
                         sb(1,b) * xnut(i-2,j,k)
        end do
        end do
     case(4)
     case(5)
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           xnut(i,j,k) = (one + rat(j,k,b))* xnut(i-1,j,k) - &
                                rat(j,k,b) * xnut(i-2,j,k)
        end do
        end do
     case(6)
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           xnut(i,j,k) = xnut(i-1,j,k)
        end do
        end do
     case(8)
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           xnut(i,j,k) = xnut(i-1,j,k)
        end do
        end do

     case(9)
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend

           xnut(i,j,k) = two * xnut(i-1,j,k) - & 
                         one * xnut(i-2,j,k)
        end do
        end do

     end select csi2
  end if

  ! boundary in eta-direction (j = 1)
  !
  if ( myleft == mpi_proc_null ) then

      b = 3
      j = j_mysta

     eta1: select case(btype(b,myzone))

     case(0)
     case(1) ! solid
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = sa(1,b) * xnut(i,j+1,k) + &
                         sb(1,b) * xnut(i,j+2,k)
        end do
        end do
     case(4)
     case(5)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = (one + rat(i,k,b))* xnut(i,j+1,k) - &
                                rat(i,k,b) * xnut(i,j+2,k)
        end do
        end do
     case(6)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = xnut(i,j+1,k)
        end do
        end do

     case(9)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend

           xnut(i,j,k) = two * xnut(i,j+1,k) - &
                         one * xnut(i,j+2,k)        
        
        end do
        end do

     end select eta1
  end if

  ! boundary in eta-direction (j = jmax)
  !
  if ( myright == mpi_proc_null ) then

      b = 4
      j = j_myend
   
     eta2: select case(btype(b,myzone))
   
     case(0)
     case(1) ! solid
        do k = k_mysta , k_myend
        do i = i_mysta , i_myend
           xnut(i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = sa(1,b) * xnut(i,j-1,k) + &
                         sb(1,b) * xnut(i,j-2,k)
        end do
        end do
     case(4)
     case(5)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = (one + rat(i,k,b))* xnut(i,j-1,k) - &
                                rat(i,k,b) * xnut(i,j-2,k)
        end do
        end do
     case(6)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = xnut(i,j-1,k)
        end do
        end do

     case(9)
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend

           xnut(i,j,k) = two * xnut(i,j-1,k) - &
                         one * xnut(i,j-2,k)
        end do
        end do


     end select eta2
  end if

  ! boundary in zet-direction (k = 1)
  !
  if ( mydown == mpi_proc_null ) then

      b = 5
      k = k_mysta
   
     zet1: select case( btype ( b , myzone ) )
   
     case(0)
     case(1) ! solid
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = sa(1,b) * xnut(i,j,k+1) + &
                         sb(1,b) * xnut(i,j,k+2)
        end do
        end do
     case(4)
     case(5)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = (one + rat(i,j,b))* xnut(i,j,k+1) - &
                                rat(i,j,b) * xnut(i,j,k+2)
        end do
        end do
     case(6)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = xnut(i,j,k_myend)
        end do
        end do

     case(9)

        do j = j_mysta, j_myend
        do i = i_mysta, i_myend

           xnut(i,j,k) = two * xnut(i,j,k+1) - &
                         one * xnut(i,j,k+2)

        end do
        end do

     end select zet1
  end if

  ! boundary in zet-direction (k = kmax)
  !
  if ( myup == mpi_proc_null ) then

      b = 6
      k = k_myend
   
     zet2: select case( btype( b , myzone ) )
   
     case(0)
     case(1) ! solid
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = sa(1,b) * xnut(i,j,k-1) + &
                         sb(1,b) * xnut(i,j,k-2)
        end do
        end do
     case(4)
     case(5)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = (one + rat(i,j,b))* xnut(i,j,k-1) - &
                                rat(i,j,b) * xnut(i,j,k-2)
        end do
        end do
     case(6)
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           xnut(i,j,k) = xnut(i,j,k_mysta)
        end do
        end do

     case(9)

        do j = j_mysta, j_myend
        do i = i_mysta, i_myend

           xnut(i,j,k) = two * xnut(i,j,k-1) - &
                         one * xnut(i,j,k-2)

        end do
        end do


     end select zet2
  end if


   if (nblk /= 0) then

      do nb = 1, nblk
         
         ! Let's zero all the variables in the interior region of the blanking zone
         do k = max( k_mysta , li_blk_ka(1,nb)+1 ) , min(k_myend , li_blk_kb(1,nb)-1 )
         do j = max( j_mysta , li_blk_ja(1,nb)+1 ) , min(j_myend , li_blk_jb(1,nb)-1 )
         do i = max( i_mysta , li_blk_ia(1,nb)+1 ) , min(i_myend , li_blk_ib(1,nb)-1 )

            xnut(i,j,k)  = zero

         end do
         end do
         end do

         i = li_blk_ia(1,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 ) then
               
            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
            do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i-1,j,k) - xnut(i-2,j,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i-1,j,k) ; y1 = y(i-1,j,k) ; z1 = z(i-1,j,k) ;
               x2 = x(i-2,j,k) ; y2 = y(i-2,j,k) ; z2 = z(i-2,j,k) ;
               xb = x(i  ,j,k) ; yb = y(i  ,j,k) ; zb = z(i  ,j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i-1,j,k) + alpha2 * xnut(i-2,j,k)

            end do
            end do

         end if
      
         i = li_blk_ib(1,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 ) then
               
            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
            do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i+1,j,k)  - xnut(i+2,j,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i+1,j,k) ; y1 = y(i+1,j,k) ; z1 = z(i+1,j,k) ;
               x2 = x(i+2,j,k) ; y2 = y(i+2,j,k) ; z2 = z(i+2,j,k) ;
               xb = x(i  ,j,k) ; yb = y(i  ,j,k) ; zb = z(i  ,j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i+1,j,k) + alpha2 * xnut(i+2,j,k)
            
            end do
            end do

         end if

         j = li_blk_ja(1,nb)

         if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then
               
            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
            do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i,j-1,k)  - xnut(i,j-2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i,j-1,k) ; y1 = y(i,j-1,k) ; z1 = z(i,j-1,k) ;
               x2 = x(i,j-2,k) ; y2 = y(i,j-2,k) ; z2 = z(i,j-2,k) ;
               xb = x(i  ,j,k) ; yb = y(i  ,j,k) ; zb = z(i  ,j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i,j-1,k) + alpha2 * xnut(i,j-2,k)
            
            end do
            end do

         end if
      
         j = li_blk_jb(1,nb)

         if ( blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then
               
            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
            do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i,j+1,k)  - xnut(i,j+2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i,j+1,k) ; y1 = y(i,j+1,k) ; z1 = z(i,j+1,k) ;
               x2 = x(i,j+2,k) ; y2 = y(i,j+2,k) ; z2 = z(i,j+2,k) ;
               xb = x(i  ,j,k) ; yb = y(i  ,j,k) ; zb = z(i  ,j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i,j+1,k) + alpha2 * xnut(i,j+2,k)
            
            end do
            end do

         end if

         ! Corners

         i = li_blk_ia(1,nb)
         j = li_blk_ja(1,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i-1,j-1,k)  - xnut(i-2,j-2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i-1,j-1,k) ; y1 = y(i-1,j-1,k) ; z1 = z(i-1,j-1,k) ;
               x2 = x(i-2,j-2,k) ; y1 = y(i-2,j-2,k) ; z1 = z(i-2,j-2,k) ;
               xb = x(i  ,  j,k) ; yb = y(i  ,  j,k) ; zb = z(i  ,  j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i-1,j-1,k) + alpha2 * xnut(i-2,j-2,k)

            end do

         end if


         i = li_blk_ib(1,nb)
         j = li_blk_ja(1,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i+1,j-1,k)  - xnut(i+2,j-2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i+1,j-1,k) ; y1 = y(i+1,j-1,k) ; z1 = z(i+1,j-1,k) ;
               x2 = x(i+2,j-2,k) ; y1 = y(i+2,j-2,k) ; z1 = z(i+2,j-2,k) ;
               xb = x(i  ,  j,k) ; yb = y(i  ,  j,k) ; zb = z(i  ,  j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i+1,j-1,k) + alpha2 * xnut(i+2,j-2,k)
               
            end do

         end if

         i = li_blk_ib(1,nb)
         j = li_blk_jb(1,nb)

         if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i+1,j+1,k)  - xnut(i+2,j+2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i+1,j+1,k) ; y1 = y(i+1,j+1,k) ; z1 = z(i+1,j+1,k) ;
               x2 = x(i+2,j+2,k) ; y1 = y(i+2,j+2,k) ; z1 = z(i+2,j+2,k) ;
               xb = x(i  ,  j,k) ; yb = y(i  ,  j,k) ; zb = z(i  ,  j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i+1,j+1,k) + alpha2 * xnut(i+2,j+2,k)
               
            end do

         end if


         i = li_blk_ia(1,nb)
         j = li_blk_jb(1,nb)

         if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
              blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

            do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               ! xnut(i,j,k)  = two * xnut(i-1,j+1,k)  - xnut(i-2,j+2,k)    

               ! Extrapolated xnut for arbitrary grid
               x1 = x(i-1,j+1,k) ; y1 = y(i-1,j+1,k) ; z1 = z(i-1,j+1,k) ;
               x2 = x(i-2,j+2,k) ; y1 = y(i-2,j+2,k) ; z1 = z(i-2,j+2,k) ;
               xb = x(i  ,  j,k) ; yb = y(i  ,  j,k) ; zb = z(i  ,  j,k) ;

               call GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                        x2 , y2 , z2 ,     &
                                                        xb , yb , zb ,     &
                                                        alpha1 , alpha2    &
                                                      )
            
               xnut(i,j,k) = alpha1 * xnut(i-1,j+1,k) + alpha2 * xnut(i-2,j+2,k)
               
            end do

         end if

         k = li_blk_ka(1,nb)

         if ( blktype(5,nb,myzone) == 0 .and. k > k_mysta+1 ) then
               
            do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) ) 
            do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               xnut(i,j,k)  = two * xnut(i,j,k-1)  - xnut(i,j,k-2)    
            
            end do
            end do

         end if
      
         k = li_blk_kb(1,nb)

         if ( blktype(6,nb,myzone) == 0 .and. k < k_myend-1 ) then
               
            do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) ) 
            do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 

               ! Extrapolated velocities and xnut (considering uniform grid)
               xnut(i,j,k)  = two * xnut(i,j,k+1)  - xnut(i,j,k+2)    
            
            end do
            end do

         end if

      end do

   end if


  !xnut
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    if ( xnut(i,j,k) + reinol < zero ) then
      xnut(i,j,k) = - reinol
    end if


  end do
  end do
  end do



end subroutine les_bcond


