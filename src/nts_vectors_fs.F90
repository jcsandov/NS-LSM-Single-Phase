

subroutine nts_vectors_fs( i,j,k, xs,ys,zs,nvec, tvec, svec, xnut_fs, PointWithinCell)
   
   use InterpolationMethods

   implicit none

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! input - output variables
   !
   ! values at (0,0,0) mean values at (i,j,k), so values at (-1,0,1) imply values at
   ! (i-1,j,k+1). This structure is to avoid big arrays communication between program
   ! units when it is no needed.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   integer, intent(in) :: i,j,k
   real (kind = rdf) , intent(in) :: xs,ys,zs
   real (kind = rdf), dimension(3), intent(out) :: nvec, tvec, svec
   real (kind = rdf) , intent(out) :: xnut_fs
   logical :: PointWithinCell 

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! local variables

   real (kind = rdf) :: phi_gradient_norm

   ! local variables interpolation

   integer                                      :: iv1 , jv1 , kv1
   integer                                      :: iv  , jv  , kv 
   real ( kind = rdf ) , dimension(8,3)         :: CellVerticesCoordinates
   real ( kind = rdf ) , dimension(6)           :: DistanceToFaces

   integer :: CellVertexLoop

   ! Phi gradietnt
   real ( kind = rdf ) , dimension(8)           :: dphidx_CellArray , dphidy_CellArray , dphidz_CellArray
   ! Tangential vectors
   real ( kind = rdf ) , dimension(8)           :: t1_CellArray , t2_CellArray , t3_CellArray
   real ( kind = rdf ) , dimension(8)           :: s1_CellArray , s2_CellArray , s3_CellArray

   ! xnut
   real ( kind = rdf ) , dimension(8)           :: xnut_CellArray

   integer , parameter                          :: nvars = 10 ! 9
   real ( kind = rdf ) , dimension( nvars , 8 ) :: VarsToInterpolate_CellArray
   real ( kind = rdf ) , dimension( nvars )     :: nts_interpolated

   integer, parameter :: iOffset(1:8) = (/0,1,1,0,0,1,1,0/) ! i - index offset
   integer, parameter :: jOffset(1:8) = (/0,0,1,1,0,0,1,1/) ! j - index offset
   integer, parameter :: kOffset(1:8) = (/0,0,0,0,1,1,1,1/) ! k - index offset
   logical :: BlankingFlag 
   logical :: Blk_ifront , Blk_iback , Blk_jleft , Blk_jright  


   PointWithinCell = .false.
   BlankingFlag    = .false.

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! normal vector interpolation
   
   ! I look over the eight octants aroun i,j,k
   searchloop_extp : &
   do iv1 = i-1 , i 
   do jv1 = j-1 , j 
   do kv1 = k-1 , k 

      ! 
      !     i,j+1,k+1 ----i+1,j+1,k+1    
      !      /|(8)           /|(7)                            
      !     / |             / |                                  
      !  i,j,k+1-------i+1,j,k+1                   
      !    |(5)           |(6)|                                   
      !    |  |           |   |                                
      !    |  |           |   |                                
      !    |  i,j+1,k-----|-i+1,j+1,k           
      !    | /(4)         |  /(3)                 
      !    |/             | /
      !  i,j,k-----------i+1,j,k
      !   (1)               (2)
      ! 
      
      ! Coordinates of the eight vertices of the cell
      CellVerticesCoordinates(1,:) = (/ x( iv1   , jv1   , kv1   ) , &
                                        y( iv1   , jv1   , kv1   ) , &
                                        z( iv1   , jv1   , kv1   ) /)

      CellVerticesCoordinates(2,:) = (/ x( iv1+1 , jv1   , kv1   ) , &
                                        y( iv1+1 , jv1   , kv1   ) , &
                                        z( iv1+1 , jv1   , kv1   ) /)

      CellVerticesCoordinates(3,:) = (/ x( iv1+1 , jv1+1 , kv1   ) , &
                                        y( iv1+1 , jv1+1 , kv1   ) , &
                                        z( iv1+1 , jv1+1 , kv1   ) /)

      CellVerticesCoordinates(4,:) = (/ x( iv1   , jv1+1 , kv1   ) , &
                                        y( iv1   , jv1+1 , kv1   ) , &
                                        z( iv1   , jv1+1 , kv1   ) /)
   
      CellVerticesCoordinates(5,:) = (/ x( iv1   , jv1   , kv1+1 ) , &
                                        y( iv1   , jv1   , kv1+1 ) , &
                                        z( iv1   , jv1   , kv1+1 ) /)

      CellVerticesCoordinates(6,:) = (/ x( iv1+1 , jv1   , kv1+1 ) , &
                                        y( iv1+1 , jv1   , kv1+1 ) , &
                                        z( iv1+1 , jv1   , kv1+1 ) /)

      CellVerticesCoordinates(7,:) = (/ x( iv1+1 , jv1+1 , kv1+1 ) , &
                                        y( iv1+1 , jv1+1 , kv1+1 ) , &
                                        z( iv1+1 , jv1+1 , kv1+1 ) /)

      CellVerticesCoordinates(8,:) = (/ x( iv1   , jv1+1 , kv1+1 ) , &
                                        y( iv1   , jv1+1 , kv1+1 ) , &
                                        z( iv1   , jv1+1 , kv1+1 ) /)

      ! verification if the fs lies on the current cell

      call PointWithinCellCheck ( CellVerticesCoordinates , &
                                  (/xs,ys,zs/)            , &
                                  PointWithinCell         , &
                                  DistanceToFaces           &
                                )

      if ( PointWithinCell ) exit searchloop_extp

      
   end do
   end do
   end do searchloop_extp

   if ( .not. PointWithinCell ) return

   ! I came out this searchloop with iv1, jv1 and kv1 which tells me the v1 node
   ! of the cell where the free surface is located. I construct my 8 vertices cell
   ! starting at that iv1,jv1,kv1

   ! ---------------------------------------------------------------------------
   ! ϕ GRADIENT FOR NORMAL VECTOR CALCULATION
   ! ---------------------------------------------------------------------------


   do CellVertexLoop = 1,8
             
      ! Current vertex
      iv = iv1 + iOffset( CellVertexLoop )
      jv = jv1 + jOffset( CellVertexLoop )
      kv = kv1 + kOffset( CellVertexLoop )

      dphidx_CellArray( CellVertexLoop ) = phi_gradient( 1 , iv , jv , kv )
      dphidy_CellArray( CellVertexLoop ) = phi_gradient( 2 , iv , jv , kv )
      dphidz_CellArray( CellVertexLoop ) = phi_gradient( 3 , iv , jv , kv )

      ! xnut
      xnut_CellArray( CellVertexLoop ) = xnut( iv , jv , kv )

      ! I calculate the tangential vectors away from the walls 

      BlankingFlag = .false.

      if ( nblke /= 0 ) then
        do nb = 1,nblke
    
          if ( iv >= le_blk_ia(1,nb) .and. iv <= le_blk_ib(1,nb) .and. & 
               jv >= le_blk_ja(1,nb) .and. jv <= le_blk_jb(1,nb) .and. &
               kv >= le_blk_ka(1,nb) .and. kv <= le_blk_kb(1,nb) ) then
    
            BlankingFlag = .true.
          
          end if

        end do
      
      end if

      if( iv /= ista .and. iv /= iend  .and. &
          jv /= jsta .and. jv /= jend  .and. &
          kv /= ksta .and. kv /= kend  .and. &
          .not. BlankingFlag ) then

        ! This call needs a further ±1 stencil to calculate the Hessian matrix
        call tangential_vectors_fs( iv , jv , kv , tvec , svec )

      end if

      t1_CellArray( CellVertexLoop ) = tvec(1)
      t2_CellArray( CellVertexLoop ) = tvec(2)
      t3_CellArray( CellVertexLoop ) = tvec(3)

      s1_CellArray( CellVertexLoop ) = svec(1)
      s2_CellArray( CellVertexLoop ) = svec(2)
      s3_CellArray( CellVertexLoop ) = svec(3)

   end do

   Blk_ifront = .false.    
   Blk_iback  = .false. 
   Blk_jleft  = .false. 
   Blk_jright = .false. 


   if (nblk /= 0) then
      
      do nb = 1,nblk

         ! Blanking stencil biasing
         if( iv1+1 == li_blk_ia(1,nb) .and. &
             j     >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) Blk_ifront  = .true.

         if( iv1   == li_blk_ib(1,nb) .and. &
             j     >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) Blk_iback   = .true.

         if( jv1+1 == li_blk_ja(1,nb) .and. &
             i     >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) Blk_jright  = .true.

         if( jv1   == li_blk_jb(1,nb) .and. &
             i     >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) Blk_jleft   = .true.

      end do

   end if

   ! In case the blanking region is somewhere within the ghost nodes
   if (nblke /= 0) then
      
      do nb = 1,nblke

         ! Blanking stencil biasing
         if( iv1+1 == le_blk_ia(1,nb) .and. &
             j     >= le_blk_ja(1,nb) .and. j <= le_blk_jb(1,nb) ) Blk_ifront  = .true.

         if( iv1   == le_blk_ib(1,nb) .and. &
             j     >= le_blk_ja(1,nb) .and. j <= le_blk_jb(1,nb) ) Blk_iback   = .true.

         if( jv1+1 == le_blk_ja(1,nb) .and. &
             i     >= le_blk_ia(1,nb) .and. i <= le_blk_ib(1,nb) ) Blk_jright  = .true.

         if( jv1   == le_blk_jb(1,nb) .and. &
             i     >= le_blk_ia(1,nb) .and. i <= le_blk_ib(1,nb) ) Blk_jleft   = .true.

      end do

   end if

   if ( iv1 == ista .or. Blk_iback ) then
    
      t1_CellArray(1) = t1_CellArray(2) ; t2_CellArray(1) = t2_CellArray(2) ; t3_CellArray(1) = t3_CellArray(2)
      t1_CellArray(4) = t1_CellArray(3) ; t2_CellArray(4) = t2_CellArray(3) ; t3_CellArray(4) = t3_CellArray(3)
      t1_CellArray(5) = t1_CellArray(6) ; t2_CellArray(5) = t2_CellArray(6) ; t3_CellArray(5) = t3_CellArray(6)
      t1_CellArray(8) = t1_CellArray(7) ; t2_CellArray(8) = t2_CellArray(7) ; t3_CellArray(8) = t3_CellArray(7)

      s1_CellArray(1) = s1_CellArray(2) ; s2_CellArray(1) = s2_CellArray(2) ; s3_CellArray(1) = s3_CellArray(2)
      s1_CellArray(4) = s1_CellArray(3) ; s2_CellArray(4) = s2_CellArray(3) ; s3_CellArray(4) = s3_CellArray(3)
      s1_CellArray(5) = s1_CellArray(6) ; s2_CellArray(5) = s2_CellArray(6) ; s3_CellArray(5) = s3_CellArray(6)
      s1_CellArray(8) = s1_CellArray(7) ; s2_CellArray(8) = s2_CellArray(7) ; s3_CellArray(8) = s3_CellArray(7)

      xnut_CellArray(1) = xnut_CellArray(2)
      xnut_CellArray(4) = xnut_CellArray(3)
      xnut_CellArray(5) = xnut_CellArray(6)
      xnut_CellArray(8) = xnut_CellArray(7)

    end if

   if ( iv1+1 == iend .or. Blk_ifront ) then
    
      t1_CellArray(2) = t1_CellArray(1) ; t2_CellArray(2) = t2_CellArray(1) ; t3_CellArray(2) = t3_CellArray(1)
      t1_CellArray(3) = t1_CellArray(4) ; t2_CellArray(3) = t2_CellArray(4) ; t3_CellArray(3) = t3_CellArray(4)
      t1_CellArray(6) = t1_CellArray(5) ; t2_CellArray(6) = t2_CellArray(5) ; t3_CellArray(6) = t3_CellArray(5)
      t1_CellArray(7) = t1_CellArray(8) ; t2_CellArray(7) = t2_CellArray(8) ; t3_CellArray(7) = t3_CellArray(8)

      s1_CellArray(2) = s1_CellArray(1) ; s2_CellArray(2) = s2_CellArray(1) ; s3_CellArray(2) = s3_CellArray(1)
      s1_CellArray(3) = s1_CellArray(4) ; s2_CellArray(3) = s2_CellArray(4) ; s3_CellArray(3) = s3_CellArray(4)
      s1_CellArray(6) = s1_CellArray(5) ; s2_CellArray(6) = s2_CellArray(5) ; s3_CellArray(6) = s3_CellArray(5)
      s1_CellArray(7) = s1_CellArray(8) ; s2_CellArray(7) = s2_CellArray(8) ; s3_CellArray(7) = s3_CellArray(8)

      xnut_CellArray(2) = xnut_CellArray(1) 
      xnut_CellArray(3) = xnut_CellArray(4) 
      xnut_CellArray(6) = xnut_CellArray(5) 
      xnut_CellArray(7) = xnut_CellArray(8) 

   end if

   if ( jv1 == jsta .or. Blk_jleft ) then
    
      t1_CellArray(1) = t1_CellArray(4) ; t2_CellArray(1) = t2_CellArray(4) ; t3_CellArray(1) = t3_CellArray(4)
      t1_CellArray(2) = t1_CellArray(3) ; t2_CellArray(2) = t2_CellArray(3) ; t3_CellArray(2) = t3_CellArray(3)
      t1_CellArray(5) = t1_CellArray(8) ; t2_CellArray(5) = t2_CellArray(8) ; t3_CellArray(5) = t3_CellArray(8)
      t1_CellArray(6) = t1_CellArray(7) ; t2_CellArray(6) = t2_CellArray(7) ; t3_CellArray(6) = t3_CellArray(7)

      s1_CellArray(1) = s1_CellArray(4) ; s2_CellArray(1) = s2_CellArray(4) ; s3_CellArray(1) = s3_CellArray(4)
      s1_CellArray(2) = s1_CellArray(3) ; s2_CellArray(2) = s2_CellArray(3) ; s3_CellArray(2) = s3_CellArray(3)
      s1_CellArray(5) = s1_CellArray(8) ; s2_CellArray(5) = s2_CellArray(8) ; s3_CellArray(5) = s3_CellArray(8)
      s1_CellArray(6) = s1_CellArray(7) ; s2_CellArray(6) = s2_CellArray(7) ; s3_CellArray(6) = s3_CellArray(7)

      xnut_CellArray(1) = xnut_CellArray(4) 
      xnut_CellArray(2) = xnut_CellArray(3) 
      xnut_CellArray(5) = xnut_CellArray(8) 
      xnut_CellArray(6) = xnut_CellArray(7) 

   end if

   if ( jv1+1 == jend .or. Blk_jright ) then

      t1_CellArray(4) = t1_CellArray(1) ; t2_CellArray(4) = t2_CellArray(1) ; t3_CellArray(4) = t3_CellArray(1)
      t1_CellArray(3) = t1_CellArray(2) ; t2_CellArray(3) = t2_CellArray(2) ; t3_CellArray(3) = t3_CellArray(2)
      t1_CellArray(8) = t1_CellArray(5) ; t2_CellArray(8) = t2_CellArray(5) ; t3_CellArray(8) = t3_CellArray(5)
      t1_CellArray(7) = t1_CellArray(6) ; t2_CellArray(7) = t2_CellArray(6) ; t3_CellArray(7) = t3_CellArray(6)

      s1_CellArray(4) = s1_CellArray(1) ; s2_CellArray(4) = s2_CellArray(1) ; s3_CellArray(4) = s3_CellArray(1)
      s1_CellArray(3) = s1_CellArray(2) ; s2_CellArray(3) = s2_CellArray(2) ; s3_CellArray(3) = s3_CellArray(2)
      s1_CellArray(8) = s1_CellArray(5) ; s2_CellArray(8) = s2_CellArray(5) ; s3_CellArray(8) = s3_CellArray(5)
      s1_CellArray(7) = s1_CellArray(6) ; s2_CellArray(7) = s2_CellArray(6) ; s3_CellArray(7) = s3_CellArray(6)    

      xnut_CellArray(4) = xnut_CellArray(1)
      xnut_CellArray(3) = xnut_CellArray(2)
      xnut_CellArray(8) = xnut_CellArray(5)
      xnut_CellArray(7) = xnut_CellArray(6)    

   end if

   if ( kv1 == ksta ) then
    
      t1_CellArray(1) = t1_CellArray(5) ; t2_CellArray(1) = t2_CellArray(5) ; t3_CellArray(1) = t3_CellArray(5)
      t1_CellArray(2) = t1_CellArray(6) ; t2_CellArray(2) = t2_CellArray(6) ; t3_CellArray(2) = t3_CellArray(6)
      t1_CellArray(3) = t1_CellArray(7) ; t2_CellArray(3) = t2_CellArray(7) ; t3_CellArray(3) = t3_CellArray(7)
      t1_CellArray(4) = t1_CellArray(8) ; t2_CellArray(4) = t2_CellArray(8) ; t3_CellArray(4) = t3_CellArray(8)

      s1_CellArray(1) = s1_CellArray(5) ; s2_CellArray(1) = s2_CellArray(5) ; s3_CellArray(1) = s3_CellArray(5)
      s1_CellArray(2) = s1_CellArray(6) ; s2_CellArray(2) = s2_CellArray(6) ; s3_CellArray(2) = s3_CellArray(6)
      s1_CellArray(3) = s1_CellArray(7) ; s2_CellArray(3) = s2_CellArray(7) ; s3_CellArray(3) = s3_CellArray(7)
      s1_CellArray(4) = s1_CellArray(8) ; s2_CellArray(4) = s2_CellArray(8) ; s3_CellArray(4) = s3_CellArray(8)

      xnut_CellArray(1) = xnut_CellArray(5)
      xnut_CellArray(2) = xnut_CellArray(6)
      xnut_CellArray(3) = xnut_CellArray(7)
      xnut_CellArray(4) = xnut_CellArray(8)
   
   end if

   if ( kv1+1 == kend ) then
    
      t1_CellArray(5) = t1_CellArray(1) ; t2_CellArray(5) = t2_CellArray(1) ; t3_CellArray(5) = t3_CellArray(1)
      t1_CellArray(6) = t1_CellArray(2) ; t2_CellArray(6) = t2_CellArray(2) ; t3_CellArray(6) = t3_CellArray(2)
      t1_CellArray(7) = t1_CellArray(3) ; t2_CellArray(7) = t2_CellArray(3) ; t3_CellArray(7) = t3_CellArray(3)
      t1_CellArray(8) = t1_CellArray(4) ; t2_CellArray(8) = t2_CellArray(4) ; t3_CellArray(8) = t3_CellArray(4)

      s1_CellArray(5) = s1_CellArray(1) ; s2_CellArray(5) = s2_CellArray(1) ; s3_CellArray(5) = s3_CellArray(1)
      s1_CellArray(6) = s1_CellArray(2) ; s2_CellArray(6) = s2_CellArray(2) ; s3_CellArray(6) = s3_CellArray(2)
      s1_CellArray(7) = s1_CellArray(3) ; s2_CellArray(7) = s2_CellArray(3) ; s3_CellArray(7) = s3_CellArray(3)
      s1_CellArray(8) = s1_CellArray(4) ; s2_CellArray(8) = s2_CellArray(4) ; s3_CellArray(8) = s3_CellArray(4)

      xnut_CellArray(5) = xnut_CellArray(1)
      xnut_CellArray(6) = xnut_CellArray(2)
      xnut_CellArray(7) = xnut_CellArray(3)
      xnut_CellArray(8) = xnut_CellArray(4)

   end if

   ! minus sign to align it with the normal vector
   VarsToInterpolate_CellArray( 1  , 1:8 ) =  - dphidx_CellArray 
   VarsToInterpolate_CellArray( 2  , 1:8 ) =  - dphidy_CellArray
   VarsToInterpolate_CellArray( 3  , 1:8 ) =  - dphidz_CellArray

   VarsToInterpolate_CellArray( 4  , 1:8 ) =  t1_CellArray 
   VarsToInterpolate_CellArray( 5  , 1:8 ) =  t2_CellArray
   VarsToInterpolate_CellArray( 6  , 1:8 ) =  t3_CellArray

   VarsToInterpolate_CellArray( 7  , 1:8 ) =  s1_CellArray 
   VarsToInterpolate_CellArray( 8  , 1:8 ) =  s2_CellArray
   VarsToInterpolate_CellArray( 9  , 1:8 ) =  s3_CellArray
   VarsToInterpolate_CellArray( 10 , 1:8 ) =  xnut_CellArray

   call TrilinearInterpolation( CellVerticesCoordinates     , &
                                (/xs,ys,zs/)                , &
                                DistanceToFaces             , &
                                nvars                       , &
                                VarsToInterpolate_CellArray , &
                                nts_interpolated              &
                              )   

   nvec = nts_interpolated(1:3) / norm2 ( nts_interpolated (1:3) )
   tvec = nts_interpolated(4:6) / norm2 ( nts_interpolated (4:6) )
   svec = nts_interpolated(7:9) / norm2 ( nts_interpolated (7:9) )

   xnut_fs = nts_interpolated(10)

end subroutine nts_vectors_fs