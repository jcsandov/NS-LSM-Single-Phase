module InterpolationMethods
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Set of routines to perform interpolation computations 
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
   ! Jorge Sandoval, UoE/PUC. Edinburgh, August the 4th, 2023.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

   use precision
   
   implicit none

   contains

   subroutine TrilinearInterpolation ( CellVerticesCoordinates     , &
                                       point                       , &
                                       DistanceToFaces             , &
                                       nvars                       , &
                                       VarsToInterpolate_CellArray , &
                                       VarsInterpolated              &
                                     )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! input-variables
      ! ------------------
      ! CellVerticesCoordinates       : (x,y,z) coordinates of the eight nodes that form the 
      !                                 cell
      ! point                         : (x,y,z) coordinates of the point to be interpolated
      ! DistanceToFaces               : distance between local i1, im, j1, km, k1, km face 
      !                                 and the point
      ! nvars                         : Number of vars to interpolate
      ! VarsToInterpolate_CellArray   : Array with the value of the variables to perform the
      !                                 interpolation at each of the eight nodes
      !
      ! ------------------
      ! output-variables
      ! ------------------
      !
      ! VarsInterpolated              : Array with the nvars interpolated variables at the 
      !                                 point 
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
      implicit none

      real (kind = rdf), intent(in)  :: CellVerticesCoordinates(8,3), point(3), DistanceToFaces(6) 
      integer , intent(in) :: nvars

      real( kind = rdf ), dimension(1:nvars,8) :: VarsToInterpolate_CellArray

      real( kind = rdf ), dimension(nvars) , intent(out) :: VarsInterpolated

      ! Local variables

      real (kind = rdf) :: coeff(3)
      real (kind = rdf) :: LagrangePolynomial(8)   
      integer :: VarLoop, CellVertexLoop

      ! local variables for Newton-Rhapson subroutine
   
      real ( kind = rdf ) :: p(3),v000(3),v100(3),v010(3),v001(3),v101(3), & 
                             v011(3), v110(3),v111(3)
   
      integer :: nmax,MAXITS,nzones
      parameter(nmax=2000000,MAXITS=200)
      real ( kind = rdf ) :: rst(3)      
      logical :: check      
      integer :: m, its,itss
      
      !-----------------------------------------------------------------------------------------
      ! The COMMON block, a piece of shared memory in the computer, is another method for 
      ! passing information between program units. Data stored in a COMMON block is not passed 
      ! between program units via argument lists, but through the COMMON statement near the 
      ! beginning of each program unit. 
      ! For further information: https://www.obliquity.com/computer/fortran/common.html
   
      common/intp3D/v000,v100,v010,v001,v101,v011,v110,v111,p
   
      !-----------------------------------------------------------------------------------------
      
      ! Variable assignation

      do m=1,3

         v000(m) = CellVerticesCoordinates(1,m)
         v100(m) = CellVerticesCoordinates(2,m)
         v010(m) = CellVerticesCoordinates(4,m)
         v001(m) = CellVerticesCoordinates(5,m)
         v101(m) = CellVerticesCoordinates(6,m)
         v011(m) = CellVerticesCoordinates(8,m)
         v110(m) = CellVerticesCoordinates(3,m)
         v111(m) = CellVerticesCoordinates(7,m)
         p(m)    = point(m)
      
      end do

      ! initial guest NR method
      do m=1,3
         rst(m) = DistanceToFaces(2*m-1) / ( DistanceToFaces(2*m-1) + DistanceToFaces(2*m) + eps_sims)
      end do

      call newt(rst,3,check,MAXITS,its,itss)
      
      if( its < MAXITS ) then
         do m=1,3
            coeff(m) = rst(m)
         end do
      else
      
      ! When Newton-Rhapson doesn't converge, the coefficiente is asigned as the portion of 
      ! the total distance (denominator is the total distance between opposite cell faces)
   
         do m = 1,3
            coeff(m) = DistanceToFaces(2*m-1) / ( DistanceToFaces(2*m-1) + DistanceToFaces(2*m) )
         end do
   
      end if ! its

      ! Polynomials defined in Delandmeter et al., (2019)
   
      LagrangePolynomial(1) = (one-coeff(1))  *  (one-coeff(2))   *  (one-coeff(3))
      LagrangePolynomial(2) = coeff(1)        *  (one-coeff(2))   *  (one-coeff(3))
      LagrangePolynomial(3) = coeff(1)        *  coeff(2)         *  (one-coeff(3))
      LagrangePolynomial(4) = (one-coeff(1))  *  coeff(2)         *  (one-coeff(3))
      LagrangePolynomial(5) = (one-coeff(1))  *  (one-coeff(2))   *  coeff(3)
      LagrangePolynomial(6) = coeff(1)        *  (one-coeff(2))   *  coeff(3)
      LagrangePolynomial(7) = coeff(1)        *  coeff(2)         *  coeff(3)
      LagrangePolynomial(8) = (one-coeff(1))  *  coeff(2)         *  coeff(3)

      VarsInterpolated = zero

      do VarLoop = 1 , nvars
         do CellVertexLoop = 1,8

            VarsInterpolated(VarLoop) =    VarsInterpolated( VarLoop )                             & 
                                         + LagrangePolynomial( CellVertexLoop )                    &
                                         * VarsToInterpolate_CellArray( VarLoop , CellVertexLoop )
         
         end do
      end do


   end subroutine TrilinearInterpolation


   subroutine BilinearInterpolation( VerticesCoordinates                   , &
                                     point                                 , &
                                     nvars                                 , &
                                     VarsToInterpolate_CuadrilateralArray  , &
                                     VarsInterpolated                      , &
                                     PrintResults                            &
                                    ) 
      
      implicit none

      ! Input-output variables
      real( kind = rdf ), dimension(4,3) , intent(in) :: VerticesCoordinates
      real( kind = rdf ), dimension(3)   , intent(in) :: point
      integer                            , intent(in) :: nvars
      real( kind = rdf ), dimension(1:nvars,4) , intent(in) :: VarsToInterpolate_CuadrilateralArray
      real( kind = rdf ), dimension(1:nvars)   , intent(out) :: VarsInterpolated
      logical , optional, intent(in) :: PrintResults


      ! local variables
      real( kind = rdf ), dimension(3)   :: eu, ev
      real( kind = rdf ), dimension(2)   :: p2D
      real( kind = rdf ), dimension(4,2) :: VerticesCoordinates2D

      real (kind = rdf) :: coeff(2)
      real (kind = rdf) :: LagrangePolynomial(4)   
      integer :: VarLoop, VertexLoop


      real ( kind = rdf ) :: p(2),v00(2),v10(2),v11(2),v01(2)
   
      integer :: nmax,MAXITS,nzones
      parameter(nmax=2000000,MAXITS=200)
      real ( kind = rdf ) :: rst(2)      
      logical :: check      
      integer :: m, its,itss
      
      !-----------------------------------------------------------------------------------------
      ! The COMMON block, a piece of shared memory in the computer, is another method for 
      ! passing information between program units. Data stored in a COMMON block is not passed 
      ! between program units via argument lists, but through the COMMON statement near the 
      ! beginning of each program unit. 
      ! For further information: https://www.obliquity.com/computer/fortran/common.html
   
      common/intp2D/v00,v10,v11,v01,p


      ! Transformation of the problem to 2D auxiliary coordinates

      eu = VerticesCoordinates(2,:) - VerticesCoordinates(1,:) 
      ev = VerticesCoordinates(4,:) - VerticesCoordinates(1,:) 

      ! Normalisation
      eu = eu / norm2( eu )
      ev = ev / norm2( ev )

      ! Projection to the 2D plane
      p2D(1) = dot_product( point - VerticesCoordinates(1,:), eu )
      p2D(2) = dot_product( point - VerticesCoordinates(1,:), ev )

      VerticesCoordinates2D = zero

      ! Second vertex in 2D
      VerticesCoordinates2D(2,1) = dot_product( VerticesCoordinates(2,:) - VerticesCoordinates(1,:) , eu ) 
      VerticesCoordinates2D(2,2) = zero

      ! Third vertex in 2D
      VerticesCoordinates2D(3,1) = dot_product( VerticesCoordinates(3,:) - VerticesCoordinates(1,:) , eu )
      VerticesCoordinates2D(3,2) = dot_product( VerticesCoordinates(3,:) - VerticesCoordinates(1,:) , ev )

      ! Fourth vertex in 2S
      VerticesCoordinates2D(4,1) = zero                                        
      VerticesCoordinates2D(4,2) = dot_product( VerticesCoordinates(4,:) - VerticesCoordinates(1,:) , ev )

      ! Variable assignation

      do m=1,2
         v00(m) = VerticesCoordinates2D(1,m)
         v10(m) = VerticesCoordinates2D(2,m)
         v11(m) = VerticesCoordinates2D(3,m)
         v01(m) = VerticesCoordinates2D(4,m)
         p(m)   = p2D(m)
      end do

      ! initial guest NR method
      rst(1) = ( p(1) - v00(1) ) / ( v10(1) - v00(1) )
      rst(2) = ( p(2) - v00(2) ) / ( v01(2) - v00(2) )

      ! Newton Rhapson Algorithm to solve the non-linear system
      ! to obtain the interpolation coefficients
      
      call newt(rst,2,check,MAXITS,its,itss)
      
      if( its < MAXITS ) then
         do m=1,2
            coeff(m) = rst(m)
         end do

      else
         coeff(1) = ( p(1) - v00(1) ) / ( v10(1) - v00(1) )
         coeff(2) = ( p(2) - v00(2) ) / ( v11(2) - v00(1) )
      end if

      ! Polynomials defined in Delandmeter et al., (2019)
   
      LagrangePolynomial(1) = ( one - coeff(1) )  *  ( one - coeff(2) ) 
      LagrangePolynomial(2) =         coeff(1)    *          coeff(2)  
      LagrangePolynomial(3) =         coeff(1)    *  ( one - coeff(2) )       
      LagrangePolynomial(4) = ( one - coeff(1) )  *          coeff(2)        

      VarsInterpolated = zero
      

      do VarLoop = 1 , nvars
         do VertexLoop = 1,4

            VarsInterpolated(VarLoop) =    VarsInterpolated( VarLoop )                                 & 
                                         + LagrangePolynomial( VertexLoop )                            &
                                         * VarsToInterpolate_CuadrilateralArray( VarLoop , VertexLoop )
         
         end do
      end do


      if ( present( PrintResults ) ) then

         if( PrintResults ) then

            print *, 'PRESSURE AT THE CORNERS '
            print *, 'P at V1 = ', VarsToInterpolate_CuadrilateralArray( 1 , 1 )
            print *, 'P at V2 = ', VarsToInterpolate_CuadrilateralArray( 1 , 2 )
            print *, 'P at V3 = ', VarsToInterpolate_CuadrilateralArray( 1 , 3 )
            print *, 'P at V4 = ', VarsToInterpolate_CuadrilateralArray( 1 , 4 )
            print *, 'P interpolated = ', VarsInterpolated(1)
            print *, ' '
   
            print *, 'u-VEL AT THE CORNERS '
            print *, 'u at V1 = ', VarsToInterpolate_CuadrilateralArray( 2 , 1 )
            print *, 'u at V2 = ', VarsToInterpolate_CuadrilateralArray( 2 , 2 )
            print *, 'u at V3 = ', VarsToInterpolate_CuadrilateralArray( 2 , 3 )
            print *, 'u at V4 = ', VarsToInterpolate_CuadrilateralArray( 2 , 4 )
            print *, 'u interpolated = ', VarsInterpolated(2)
            print *, ' '
   
            print *, 'v-VEL AT THE CORNERS '
            print *, 'v at V1 = ', VarsToInterpolate_CuadrilateralArray( 3 , 1 )
            print *, 'v at V2 = ', VarsToInterpolate_CuadrilateralArray( 3 , 2 )
            print *, 'v at V3 = ', VarsToInterpolate_CuadrilateralArray( 3 , 3 )
            print *, 'v at V4 = ', VarsToInterpolate_CuadrilateralArray( 3 , 4 )
            print *, 'v interpolated = ', VarsInterpolated(3)
            print *, ' '
   
            print *, 'w-VEL AT THE CORNERS '
            print *, 'w at V1 = ', VarsToInterpolate_CuadrilateralArray( 4 , 1 )
            print *, 'w at V2 = ', VarsToInterpolate_CuadrilateralArray( 4 , 2 )
            print *, 'w at V3 = ', VarsToInterpolate_CuadrilateralArray( 4 , 3 )
            print *, 'w at V4 = ', VarsToInterpolate_CuadrilateralArray( 4 , 4 )
            print *, 'w interpolated = ', VarsInterpolated(4)
            print *, ' '

            print *, 'INTERPOLATION COEFFICIENTS '
            print *, 'coeff 1 = ', coeff(1)
            print *, 'coeff 2 = ', coeff(2)
            print *, ' '

         end if

      end if

   end subroutine BilinearInterpolation


   subroutine PointWithinCellCheck( CellVerticesCoordinates     , &
                                    point                       , &
                                    PointWithinCell             , &
                                    DistanceToFaces               &             
                                   )
   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! this subroutine determines if a point is inside a cell using their nodal 
      ! coordinates
      !
      ! it calls vector subroutine which compute the distance between cell faces to the
      ! point using a cross-product criterion
      !
      ! For details, reader is referred to
      !
      ! Lackey, T. C. (2004). Numerical Investigation of Chaotic Advection in 
      ! Three-Dimensional Experimentally Realizable Rotating Flows Numerical Investigation 
      ! of Chaotic Advection in Three-Dimensional Experimentally Realizable Rotating Flows. 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
      !   NODE DISTRIBUTION OF THE COMPUTATIONAL CELL
      !
      !    i,j+1,k+1 ----i+1,j+1,k+1    
      !     /|(8)           /|(7)                            
      !    / |             / |                                  
      ! i,j,k+1-------i+1,j,k+1                   
      !   |(5)           |(6)|                                   
      !   |  |           |   |                                
      !   |  |           |   |                                
      !   |  i,j+1,k-----|-i+1,j+1,k           
      !   | /(4)         |  /(3)                 
      !   |/             | /
      ! i,j,k-----------i+1,j,k
      !  (1)               (2)


   
      implicit none

      real ( kind = rdf ) :: CellVerticesCoordinates(8,3) , point(3), dis(6)
      real(kind = rdf), optional, intent(out) :: DistanceToFaces(6)
      logical :: PointWithinCell
      integer :: i,p1(6),p2(6),p3(6),p4(6)
     
      PointWithinCell = .false.
   
      ! local k1 face
      p1(5)=1
      p2(5)=4
      p3(5)=3
      p4(5)=2
   
      ! local j1 face
      p1(3)=1
      p2(3)=2
      p3(3)=6
      p4(3)=5
   
      ! local i1 face
      p1(1)=1
      p2(1)=5
      p3(1)=8
      p4(1)=4
   
      ! local im face
      p1(2)=2
      p2(2)=3
      p3(2)=7
      p4(2)=6
   
      ! local km face
      p1(6)=5
      p2(6)=6
      p3(6)=7
      p4(6)=8
   
      ! local jm face
      p1(4)=4
      p2(4)=8
      p3(4)=7
      p4(4)=3
   
      dis = zero

      do i=1,6

         call vector( CellVerticesCoordinates(p1(i),:) , &
                      CellVerticesCoordinates(p2(i),:) , &
                      CellVerticesCoordinates(p3(i),:) , &
                      CellVerticesCoordinates(p4(i),:) , &
                      point , dis(i) )

         if( abs(dis(i) ) < 1.0E-08 ) dis(i) = zero
   
         if( dis(i)  < zero ) then
            PointWithinCell = .false.
            exit
         else
            PointWithinCell = .true.
         end if
   
      end do
   
      if ( present( DistanceToFaces ) ) DistanceToFaces = dis

   end subroutine PointWithinCellCheck
   

   subroutine LineQuadrilateralIntersectionTest( VerticesCoordinates      , &
                                                 PositionVector           , &
                                                 DirectionVector          , &
                                                 IsIntersection           , &
                                                 IntersectionPoint             )
     
      implicit none
     
      ! Inputs:
     
      real( kind = rdf ), dimension(4,3) , intent(in) :: VerticesCoordinates ! Coordinates of the quadrilateral's vertices
      real( kind = rdf ), dimension(3)   , intent(in) :: PositionVector      ! Coordinates of the starting point of the line
      real( kind = rdf ), dimension(3)   , intent(in) :: DirectionVector     ! Coordinates of the ending point of the line
     
      ! Outputs:
      logical, intent(out) :: IsIntersection  ! True if the line intersects the quadrilateral, false otherwise
      real( kind = rdf ), dimension(3), intent(out) :: IntersectionPoint  ! Intersection point (output only if is_intersection is true)

      ! local variables

      real( kind = rdf ) :: t, u, v

      ! Triangle 1: V1, V2, V3
      ! Triangle 2: V3, V4, V1


      IsIntersection = .false.

      call IntersectionMollerTrumbore( PositionVector                           , &
                                       DirectionVector / norm2(DirectionVector) , &
                                       VerticesCoordinates(1,:)                 , &
                                       VerticesCoordinates(2,:)                 , &
                                       VerticesCoordinates(3,:)                 , &
                                       t, u, v                                    &
                                      )
      
      !print * , ' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
      !write(*, '(A, 3(F15.12, A))') "Position  Vector : ", PositionVector(1), " ", PositionVector(2), " ", PositionVector(3)
      !write(*, '(A, 3(F15.12, A))') "Direction Vector : ", DirectionVector(1) , " ", DirectionVector(2), " ", DirectionVector(3)

      if ( t >= zero ) then

         IsIntersection    = .true. 
         IntersectionPoint = PositionVector + t * DirectionVector / norm2(DirectionVector)


         !print *, ' '
         !print *, 'Inside Quadrilateral check'
         !write(*, '(A, (F15.12, A))') "t = ", t
         !write(*, '(A, 3(F15.12, A))') "IntersectionPoint = ", IntersectionPoint(1), " ", IntersectionPoint(2), " ", IntersectionPoint(3)
         !print *, ' '


         return

      else

         !print *, ' '
         !write(*, '(A)') "El checking de Moller Trumbore fallo "
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(1,:) : ", VerticesCoordinates(1,1) , " ", VerticesCoordinates(1,2), " ", VerticesCoordinates(1,3)
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(2,:) : ", VerticesCoordinates(2,1) , " ", VerticesCoordinates(2,2), " ", VerticesCoordinates(2,3)
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(3,:) : ", VerticesCoordinates(3,1) , " ", VerticesCoordinates(3,2), " ", VerticesCoordinates(3,3)

         !print *, 't = ', t


      end if

      call IntersectionMollerTrumbore( PositionVector                           , &
                                       DirectionVector / norm2(DirectionVector) , &
                                       VerticesCoordinates(3,:)                 , &
                                       VerticesCoordinates(4,:)                 , &
                                       VerticesCoordinates(1,:)                 , &
                                       t, u, v                                    &
                                      )

      if ( t < zero ) then
         
         !print *, ' '
         !write(*, '(A)') "El checking de Moller Trumbore fallo "
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(3,:) : ", VerticesCoordinates(3,1) , " ", VerticesCoordinates(3,2), " ", VerticesCoordinates(3,3)
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(4,:) : ", VerticesCoordinates(4,1) , " ", VerticesCoordinates(4,2), " ", VerticesCoordinates(4,3)
         !write(*, '(A, 3(F15.12, A))') "VerticesCoordinates(1,:) : ", VerticesCoordinates(1,1) , " ", VerticesCoordinates(1,2), " ", VerticesCoordinates(1,3)
         
         return

      else
      
         IsIntersection    = .true. 
         IntersectionPoint = PositionVector + t * DirectionVector / norm2(DirectionVector)

         !print *, ' '
         !print *, 'Inside Quadrilateral check'
         !write(*, '(A, (F15.12, A))') "t = ", t
         !write(*, '(A, 3(F15.12, A))') "IntersectionPoint = ", IntersectionPoint(1), " ", IntersectionPoint(2), " ", IntersectionPoint(3)
         !print *, ' '
                    
         return
      
      end if


   end subroutine LineQuadrilateralIntersectionTest


   subroutine LineQuadrilateralIntersectionTest_old( VerticesCoordinates      , &
                                                     PositionVector           , &
                                                     DirectionVector          , &
                                                     IsIntersection           , &
                                                     IntersectionPoint             )
     
      implicit none
     
      ! Inputs:
     
      real( kind = rdf ), dimension(4,3) , intent(in) :: VerticesCoordinates ! Coordinates of the quadrilateral's vertices
      real( kind = rdf ), dimension(3)   , intent(in) :: PositionVector   ! Coordinates of the starting point of the line
      real( kind = rdf ), dimension(3)   , intent(in) :: DirectionVector     ! Coordinates of the ending point of the line
     
      ! Outputs:
      logical :: IsIntersection  ! True if the line intersects the quadrilateral, false otherwise
      real( kind = rdf ), dimension(3), intent(out) :: IntersectionPoint  ! Intersection point (output only if is_intersection is true)

      real( kind = rdf ), dimension(3) :: p
      real( kind = rdf ) :: det, u, v
      integer :: i, j
   
      ! Initialize outputs
      IsIntersection = .false.
   
      ! Iterate over the three pairs of adjacent vertices of the quadrilateral
      do i = 1, 3

         ! Define two edges of the quadrilateral (vertex i and vertex i+1)
         j = i + 1
         if (j == 5) j = 1

         p = VerticesCoordinates(j,:) - VerticesCoordinates(i,:)
   
         ! Calculate the determinant of the two matrices
         det = dot_product( cross_product( DirectionVector, p )                 , &
                            VerticesCoordinates(1,:) - VerticesCoordinates(i,:)       )
   
         ! If the determinant is close to zero, the line is parallel to the plane of the quadrilateral
         if (abs(det) < 1e-10) cycle
   
         ! Calculate the parameter values u and v
         u = dot_product( cross_product( VerticesCoordinates(1,:) - PositionVector, p )              , &
                          VerticesCoordinates(1,:) - VerticesCoordinates(i,:)                      ) / det

         v = dot_product( cross_product(DirectionVector, VerticesCoordinates(1,:) - PositionVector ) , &
                          VerticesCoordinates(1,:) - VerticesCoordinates(i,:)                      ) / det
   
         ! Check if the intersection point lies within the quadrilateral
         if (u >= eps_sims .and. u <= one .and. v >= eps_sims .and. v <= one) then
      
            ! Calculate the intersection point
            IntersectionPoint = PositionVector + u * DirectionVector
            IsIntersection = .true.
            return
      
         end if
      
      end do
   
   end subroutine LineQuadrilateralIntersectionTest_old


   subroutine vector(p1,p2,p3,p4,point,dis)
     
      implicit none
     
      real ( kind = rdf ) :: dis,p1(3),p2(3),p3(3),p4(3),point(3)
      real ( kind = rdf ) :: diag13(3),diag24(3),FaceNormal(3),FaceCentroid2Point(3)
   

      !             (3)
      !            / |
      !           /  |  
      !         (2)  |
      !          | C |      point  
      !          | *-------->*           
      ! OUTSIDE  |   |                 INSIDE   
      !   THE    |   |                  THE    
      !  CELL    |  /(4)                CELL
      !          | /        
      !          |/         
      !          (1)    

      diag13 = p3-p1
      diag24 = p4-p2

      ! This FaceNormal vector points to the interior of the cell
      FaceNormal = cross_product( diag24 , diag13 )

      ! vector between face centroid and point to be
      ! interpolated
   
      FaceCentroid2Point = point - (p1+p2+p3+p4)/four

      ! projection of the Facenormal onto FaceCentroid2Point
      ! if the projection is negative. The point is outside the cell
      dis = dot_product( FaceNormal , FaceCentroid2Point )
      
      ! possible max(dis) = distance between the centroid of a cell and
      !                     the opposite vertex of the cell

      if( dis >= zero ) dis = norm2( FaceCentroid2Point )

   end subroutine vector


!   subroutine vector_old(p1,p2,p3,p4,point,dis)
!     
!      implicit none
!     
!      real ( kind = rdf ) ::dis,p1(3),p2(3),p3(3),p4(3),point(3)
!      real ( kind = rdf ) ::a1v,a2v,a3v,b1v,b2v,b3v,c1v,c2v,c3v,ab1,ab2,ab3
!   
!      a1v = p4(1) - p2(1)
!      a2v = p4(2) - p2(2)
!      a3v = p4(3) - p2(3)
!      
!      b1v = p3(1) - p1(1)
!      b2v = p3(2) - p1(2)
!      b3v = p3(3) - p1(3)
!   
!      ! distance between face centroid and point to be
!      ! interpolated
!   
!      c1v = point(1) - ( p1(1) + p2(1) + p3(1) + p4(1) ) / four
!      c2v = point(2) - ( p1(2) + p2(2) + p3(2) + p4(2) ) / four
!      c3v = point(3) - ( p1(3) + p2(3) + p3(3) + p4(3) ) / four
!
!      ! av x bv
!
!      ab1 = a2v*b3v - a3v*b2v
!      ab2 = a3v*b1v - a1v*b3v
!      ab3 = a1v*b2v - a2v*b1v
!
!      dis = c1v*ab1 + c2v*ab2 + c3v*ab3
!   
!   end subroutine vector_old


   function OctantIdentification( rijk, r_ip1 , r_jp1 , r_kp1 , rvec ) result(OctantID)
   
      implicit none

      real( kind = rdf ), dimension(3), intent(in) :: rijk, r_ip1 , r_jp1 , r_kp1 , rvec
      real( kind = rdf ), dimension(3) :: ei, ej, ek
      integer :: CaseFlag
      integer :: OctantID

      ei = ( r_ip1 - rijk ) / norm2( r_ip1 - rijk )
      ej = ( r_jp1 - rijk ) / norm2( r_jp1 - rijk )
      ek = ( r_kp1 - rijk ) / norm2( r_kp1 - rijk )
      
      !print *, ' '
      !print *, ' ei = ', ei
      !print *, ' ej = ', ej
      !print *, ' ek = ', ek
      !print *, 'rvec = ', rvec/norm2(rvec)
      
      ! unset bit --> positive projection
      !   set bit --> negative projection
      
      !---------------------------------------*
      !  Bits array   |  Integer  |   Octant  | 
      !---------------|-----------|-----------|
      !      000      |     0     |     1     | 
      !      001      |     1     |     2     | 
      !      011      |     3     |     3     |  
      !      010      |     2     |     4     |   
      !---------------------------------------*
      !      100      |     4     |     5     |   
      !      101      |     5     |     6     | 
      !      111      |     7     |     7     |  
      !      110      |     6     |     8     |
      !---------------------------------------*
      
      CaseFlag = 0
      
      if ( dot_product (rvec,ei) < zero ) CaseFlag = ibset(CaseFlag,0) 
      if ( dot_product (rvec,ej) < zero ) CaseFlag = ibset(CaseFlag,1) 
      if ( dot_product (rvec,ek) < zero ) CaseFlag = ibset(CaseFlag,2) 
      
      select case( CaseFlag )
      
         case(0) 
            OctantID = 1
         case(1) 
            OctantID = 2
         case(3) 
            OctantID = 3
         case(2) 
            OctantID = 4
      
         case(4) 
            OctantID = 5
         case(5) 
            OctantID = 6
         case(7) 
            OctantID = 7
         case(6) 
            OctantID = 8
      
      end select
   
   end function OctantIdentification
   
   
   function GetOctantFacesIndexes (OctantID) result(Faces)
   
      implicit none

   
      integer, intent(in) :: OctantID
      !-----------------------------------------------------
      !        FaceID   VertexID   i,j,k index shift (±1) 
      !           \        |        /  
      !            \       |       /
      !             \      |      /
      !              \     |     /
      !               \    |    /
      !                \   |   /
      !                 \  |  /
      integer, dimension(3,4,3) :: Faces
      !-----------------------------------------------------
      
      ! Local variables
      integer, dimension(8,3) :: OctantDiagonalShift
      integer, dimension(3)   :: DiagonalShiftBase
      integer, dimension(3)   :: base1, base2        ! Face vector bases
      
      
      ! Look up table of the shift indexes for the diagonals of every octant
      OctantDiagonalShift(1,:) = (/  1 ,  1 ,  1 /)
      OctantDiagonalShift(2,:) = (/ -1 ,  1 ,  1 /)
      OctantDiagonalShift(3,:) = (/ -1 , -1 ,  1 /)
      OctantDiagonalShift(4,:) = (/  1 , -1 ,  1 /)
      OctantDiagonalShift(5,:) = (/  1 ,  1 , -1 /)
      OctantDiagonalShift(6,:) = (/ -1 ,  1 , -1 /)
      OctantDiagonalShift(7,:) = (/ -1 , -1 , -1 /)
      OctantDiagonalShift(8,:) = (/  1 , -1 , -1 /)
      
      ! It gives you the vector base from the opposite diagonal pointing
      ! i,j,k point
   
      DiagonalShiftBase = - OctantDiagonalShift( OctantID , :)
      
      !-----------------------------------------------------------------------!
      !                FACE 1: Defined by i,j base vectors                    !
      !-----------------------------------------------------------------------!
      
      base1 = (/1,0,0/)
      base2 = (/0,1,0/)
      
      ! Face 1 - V1
      Faces(1,1,:) = OctantDiagonalShift(OctantID,:)
      
      ! Face 1 - V2
      Faces(1,2,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(1) * base1
      
      ! Face 1 - V3
      Faces(1,3,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(1) * base1 + DiagonalShiftBase(2) * base2 
      
      ! Face 1 - V4
      Faces(1,4,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(2)*base2 
      
      !-----------------------------------------------------------------------!
      !                FACE 2: Defined by j,k base vectors                    !
      !-----------------------------------------------------------------------!
      
      base1 = (/0,1,0/)
      base2 = (/0,0,1/)
      
      ! Face 2 - V1
      Faces(2,1,:) = OctantDiagonalShift(OctantID,:)
      
      ! Face 2 - V2
      Faces(2,2,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(2) * base1
      
      ! Face 2 - V3
      Faces(2,3,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(2) * base1 + DiagonalShiftBase(3) * base2 
      
      ! Face 2 - V4
      Faces(2,4,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(3) * base2 
      
      
      !-----------------------------------------------------------------------!
      !                FACE 3: Defined by k,i base vectors                    !
      !-----------------------------------------------------------------------!
      
      base1 = (/0,0,1/)
      base2 = (/1,0,0/)
      
      ! Face 3 - V1
      Faces(3,1,:) = OctantDiagonalShift(OctantID,:)
      
      ! Face 3 - V2
      Faces(3,2,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(3) * base1
      
      ! Face 3 - V3
      Faces(3,3,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(3) * base1 + DiagonalShiftBase(1) * base2 
      
      ! Face 3 - V4
      Faces(3,4,:) = OctantDiagonalShift(OctantID,:) + &
                     DiagonalShiftBase(1) * base2 
      
   
   end function GetOctantFacesIndexes



   subroutine IntersectionMollerTrumbore(orig, dir, vert0, vert1, vert2, t, u, v)

      ! This algorithm is fully described in:
      ! Möller, T., & Trumbore, B. (2005). Fast, minimum storage ray/triangle intersection. 
      ! In ACM SIGGRAPH 2005 Courses (pp. 7-es).
      ! https://dl.acm.org/doi/pdf/10.1145/1198555.1198746

      real( kind = rdf ), intent(in) :: orig(3), dir(3), vert0(3), vert1(3), vert2(3)
      real( kind = rdf ), intent(out) :: t, u, v
      
      real( kind = rdf ) :: edge1(3), edge2(3), tvec(3), pvec(3), qvec(3)
      real( kind = rdf ) :: det, inv_det
      
      real( kind = rdf ), parameter :: TOL = 0.000001_rdf

      ! orig is the origin of the ray in 3D space.
      ! dir is the normalized direction of the ray.
      ! vert0, vert1, and vert2 are the vertices of the triangle.
      ! t is the output parameter for the distance from the ray origin to the intersection point.
      ! u and v are the output parameters for the barycentric coordinates of the intersection point.
    
      ! Calculate edges of the triangle.
      edge1 = vert1 - vert0
      edge2 = vert2 - vert0
      
      ! Calculate cross product of ray direction and edge2.
      pvec = cross_product(dir, edge2)
      
      ! Calculate dot product of edge1 and h.
      det = dot_product(edge1, pvec)
      
      ! Non-culling branch
      if (det > -TOL .and. det < TOL) then
         t = -one
         return
      end if

      inv_det = one / det

      ! Calculate distance from vert0 to ray origin
      tvec = orig - vert0

      ! Calculate U parameter and test bounds
      u = dot_product(tvec, pvec) * inv_det

      if (u < zero .or. u > one) then
        t = -one
        return
      end if

      ! Prepare to test V parameter
      qvec = cross_product( tvec, edge1 )

      ! Calculate V parameter and test bounds
      v = dot_product(dir, qvec) * inv_det
      
      if ( v < zero .or. u + v > one ) then
         t = -one
         return
      end if

      ! Calculate t, ray intersects triangle
      t = dot_product(edge2, qvec) * inv_det

   end subroutine IntersectionMollerTrumbore


   subroutine IntersectionMollerTrumbore_old(orig, dir, vert0, vert1, vert2, t, u, v)

      ! This algorithm is fully described in:
      ! Möller, T., & Trumbore, B. (2005). Fast, minimum storage ray/triangle intersection. 
      ! In ACM SIGGRAPH 2005 Courses (pp. 7-es).
      ! https://dl.acm.org/doi/pdf/10.1145/1198555.1198746

      real( kind = rdf ), intent(in) :: orig(3), dir(3), vert0(3), vert1(3), vert2(3)
      real( kind = rdf ), intent(out) :: t, u, v
      
      real( kind = rdf ) :: edge1(3), edge2(3), h(3), s(3), q(3), a, f
      
      ! orig is the origin of the ray in 3D space.
      ! dir is the normalized direction of the ray.
      ! vert0, vert1, and vert2 are the vertices of the triangle.
      ! t is the output parameter for the distance from the ray origin to the intersection point.
      ! u and v are the output parameters for the barycentric coordinates of the intersection point.
    
      ! Calculate edges of the triangle.
      edge1 = vert1 - vert0
      edge2 = vert2 - vert0
      
      ! Calculate cross product of ray direction and edge2.
      h = cross_product(dir, edge2)
      
      ! Calculate dot product of edge1 and h.
      a = dot_product(edge1, h)
      
      ! Check for nearly parallel ray and triangle.
      if ( a > -eps_sims .and. a < eps_sims ) then
         t = -one ! No intersection
         return
      end if
      
      ! Calculate f, which is the reciprocal of a.
      f = one / a
      
      ! Calculate vector from vertex 0 to ray origin.
      s = orig - vert0
      
      ! Calculate u, the barycentric coordinate of the intersection point.
      u = f * dot_product(s, h)
      
      ! Check if u is outside the valid range [0, 1].
      if ( u < zero .or. u > one ) then
         t = -one ! No intersection
         return
      end if
      
      ! Calculate cross product of s and edge1.
      q = cross_product(s, edge1)
      
      ! Calculate v, the second barycentric coordinate.
      v = f * dot_product(dir, q)
      
      ! Check if v is outside the valid range [0, 1] or if u+v is greater than 1.
      if (v < zero .or. u + v > one) then
         t = -one ! No intersection
         return
      end if
      
      ! Calculate t, the distance along the ray to the intersection point.
      t = f * dot_product(edge2, q)
      
      ! At this point, if t >= 0.0, an intersection has been found.
      ! The intersection point's barycentric coordinates are given by u and v.
      ! If t < 0.0, there is no intersection.

   end subroutine IntersectionMollerTrumbore_old



   !-------------------------------------------------------------------------------
   !
   !       A globally convergent Newton-Raphson method (from Numerical
   !       receipes in Fortran, Vol I, by  W. H. Press, 2nd ed., p379 ) 
   !       Parameter: NP--dimension of the sysytem, MAXITS--maximum 
   !       number of iterations within which convergence shall be 
   !       obtained if it is supposed. 
   !
   !       * newt, fdjac, fmin, lnsrch, ludcmp, lubksb were obtained 
   !         directly from the reference presented above.
   !
   !       * funcv3D and fdjac3D were defined specifically for the equations 
   !         to be solved, but the guidelines and in-out parameters are 
   !         described in the same reference.
   !
   !       * fmin was modified and was converted from a function to a subroutine.
   !         it was because the use of external variables was generating 
   !         compilation issues because of the references (Jorge Sandoval,
   !         Edinburgh, November 8th, 2021)
   !
   !
   !-------------------------------------------------------------------------------
   
   SUBROUTINE newt(x,n,check,MAXITS,its,itss)
     
     !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
     integer :: its, itss
     INTEGER :: n,nn,NP
     LOGICAL :: check
     real(kind = rdf) :: x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
   
     PARAMETER (NP=40,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,   &
          &  STPMX=100.)
     COMMON /newtv/ fvec(NP),nn
     SAVE /newtv/
   !C     USES fdjac,fmin,lnsrch,lubksb,ludcmp
   
     INTEGER :: i,j,indx(NP),MAXITS
     real(kind = rdf) :: d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP), &
                         xold(NP),fmin_result
   
   !  EXTERNAL fmin
   
     nn=n
   
     ! f and fmin_result initialisation
     
     f = zero
     fmin_result = zero
   
     call fmin(x, fmin_result)
   
     f=fmin_result
   
     test=0.
     do i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
     end do
     if(test.lt..01*TOLF) then
        check=.false.
        return
     endif
   
     sum=0.
     do i=1,n    
        sum=sum+x(i)**2
     end do
   
     stpmax=STPMX*max(sqrt(sum),float(n))
   
     do its=1,MAXITS
        !write(*,*) its
        
        ! Jacobian
        if ( n == 2 ) call fdjac2D(n,x,fjac(1:n,1:n))
        if ( n == 3 ) call fdjac3D(n,x,fjac(1:n,1:n))

        ! Function
        if ( n == 2 ) call funcv2D(n,x,fvec)
        if ( n == 3 ) call funcv3D(n,x,fvec)
   
   !     print *,fjac(1:n,1:n)
        do i=1,n
           sum=0.
           do j=1,n
              sum=sum+fjac(j,i)*fvec(j)
           end do
           g(i)=sum
        end do
        do i=1,n
           xold(i)=x(i)
        end do
        fold=f
        do i=1,n
           p(i)=-fvec(i)
        end do
   !     write(*,*) '********'
        call ludcmp(fjac,n,NP,indx,d,itss)
   !     write(*,*) '*******1'
        call lubksb(fjac,n,NP,indx,p)
   !     write(*,*) '*******2'
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
   !     write(*,*) '*******3'
        test=0.
        do i=1,n
           if(abs(fvec(i)).gt.test)test=abs(fvec(i))
        end do
        if(test.lt.TOLF)then
           check=.false.
           return
        endif
        if(check)then
           test=0.
           den=max(f,.5*n)
           do i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.)/den
              if(temp.gt.test)test=temp
           end do
           if(test.lt.TOLMIN)then
              check=.true.
           else
              check=.false.
           endif
           return
        endif
        test=0.
        do i=1,n
           temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
           if(temp.gt.test)test=temp
        end do
        if(test.lt.TOLX)return
     end do
   !c      pause 'MAXITS exceeded in newt'
   END SUBROUTINE newt
   
   
   SUBROUTINE fdjac(n,x,fvec,np,df)
     INTEGER :: n,np,NMAX
     real(kind = rdf) :: df(np,np),fvec(n),x(n),EPS
     PARAMETER (NMAX=40,EPS=1.e-4)
   !    USES funcv3D
     INTEGER i,j
     real(kind = rdf) :: h,temp,f(NMAX)
     do j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv3D(n,x,f)
        x(j)=temp
        do i=1,n
           df(i,j)=(f(i)-fvec(i))/h
        end do
     end do
     return
   END SUBROUTINE fdjac
   
   
   SUBROUTINE fmin(x, fmin_result)
     INTEGER :: n,NP
     real(kind = rdf) , intent(out) :: fmin_result
     real(kind = rdf) :: x(*),fvec
     
     PARAMETER (NP=40)
     COMMON /newtv/ fvec(NP),n
     SAVE /newtv/
   !CU    USES funcv3D
     INTEGER :: i
     real(kind = rdf) :: sum
   
     fmin_result = zero
   
     if ( n == 2 ) call funcv2D(n,x,fvec)
     if ( n == 3 ) call funcv3D(n,x,fvec)

     sum=0.
     do i=1,n
        sum=sum+fvec(i)**2
     end do
     fmin_result=0.5*sum
     return
   END SUBROUTINE fmin
   
   SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
     INTEGER :: n
     LOGICAL :: check
     real(kind = rdf) :: f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF,TOLX, fmin_result
     PARAMETER (ALF=1.e-4,TOLX=1.e-7)
     !EXTERNAL func
   !CU    USES func
     INTEGER :: i
     real(kind = rdf) :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp, &
                         test,tmplam
     check=.false.
     sum=0.
     do i=1,n
        sum=sum+p(i)*p(i)
     end do
     sum=sqrt(sum)
     if(sum.gt.stpmax)then
        do i=1,n
           p(i)=p(i)*stpmax/sum
        end do
     endif
     slope=0.
     do i=1,n
        slope=slope+g(i)*p(i)
     end do
     !if (slope.ge.0.) pause 'roundoff problem in lnsrch'
     test=0.
     do i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
     end do
     alamin=TOLX/test
     alam=1.
   
     alam2=2.
     f2=2.
     fold2=fold
   
   1     continue
     do i=1,n
        x(i)=xold(i)+alam*p(i)
     end do
   
     fmin_result = zero
   
     f = zero
     fmin_result = zero
   
     call fmin(x,fmin_result)
     f = fmin_result
   
     if(alam.lt.alamin)then
        do i=1,n
           x(i)=xold(i)
        end do
        check=.true.
        return
     else if(f.le.fold+ALF*alam*slope)then
        return
     else
        if(alam.eq.1.)then
           tmplam=-slope/(2.*(f-fold-slope))
        else
           rhs1=f-fold-alam*slope
           rhs2=f2-fold-alam2*slope
           a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
           b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
           if(a.eq.0.)then
              tmplam=-slope/(2.*b)
           else
              disc=b*b-3.*a*slope
              if(disc.lt.0.) then
                 tmplam=0.5*alam
              else if (b.le.0.) then
                 tmplam=(-b+sqrt(disc))/(3.*a)
              else
                 tmplam=-slope/(b+sqrt(disc))
              endif
   !           C              tmplam=(-b+sqrt(disc))/(3.*a)
           endif
           if(tmplam.gt..5*alam)tmplam=.5*alam
        endif
     endif
     alam2=alam
     f2=f
   !  fold2=fold
     alam=max(tmplam,.1*alam)
     goto 1
   END SUBROUTINE lnsrch
   
   SUBROUTINE ludcmp(a,n,np,indx,d,itss)
     integer :: itss
     INTEGER :: n,np,indx(n),NMAX
     real(kind = rdf) :: d,a(np,np),TINY
     PARAMETER (NMAX=500,TINY=1.0e-20)
     INTEGER :: i,imax,j,k
     real(kind = rdf) :: aamax,dum,sum,vv(NMAX)
     itss=0
     d=1.
     do i=1,n
        aamax=0.
        do j=1,n
           if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        end do
        !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        if (aamax.eq.0.) itss=1
        if (aamax.eq.0.) return
        vv(i)=1./aamax
     end do
     do j=1,n
        do i=1,j-1
           sum=a(i,j)
           do k=1,i-1
              sum=sum-a(i,k)*a(k,j)
           end do
           a(i,j)=sum
        end do
        aamax=0.
        do i=j,n
           sum=a(i,j)
           do k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           end do
           a(i,j)=sum
           dum=vv(i)*abs(sum)
           if (dum.ge.aamax) then
              imax=i
              aamax=dum
           endif
        end do
        if (j.ne.imax)then
           do k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
           end do
           d=-d
           vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
           dum=1./a(j,j)
           do i=j+1,n
              a(i,j)=a(i,j)*dum
           end do
        endif
     end do
     return
   END SUBROUTINE ludcmp
   
   SUBROUTINE lubksb(a,n,np,indx,b)
     INTEGER :: n,np,indx(n)
     real(kind = rdf) :: a(np,np),b(n)
     INTEGER :: i,ii,j,ll
     real(kind = rdf) :: sum
   
     ii=0
     do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
           do j=ii,i-1
              sum=sum-a(i,j)*b(j)
           end do
        else if (sum.ne.0.) then
           ii=i
        endif
        b(i)=sum
     end do
     do i=n,1,-1
        sum=b(i)
        do j=i+1,n
           sum=sum-a(i,j)*b(j)
        end do
        b(i)=sum/a(i,i)
     end do
     return
   END SUBROUTINE lubksb
   
   SUBROUTINE funcv3D(n,rst,fvec)
     
     !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
   
     integer :: i,n
     real (kind = rdf) ::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
     common/intp3D/v000,v100,v010,v001,v101,v011,v110,v111,p
     real(kind = rdf) :: rst(n),fvec(n)
     
     do i=1,n
        fvec(i) = v000(i) * (1-rst(1)) * (1-rst(2)) * (1-rst(3)) + &
                  v100(i) *     rst(1) * (1-rst(2)) * (1-rst(3)) + &
                  v010(i) * (1-rst(1)) *    rst(2)  * (1-rst(3)) + &
                  v001(i) * (1-rst(1)) * (1-rst(2)) *    rst(3)  + &
                  v101(i) *     rst(1) * (1-rst(2)) *    rst(3)  + &
                  v011(i) * (1-rst(1)) *    rst(2)  *    rst(3)  + &
                  v110(i) *     rst(1) *    rst(2)  * (1-rst(3)) + &
                  v111(i) *     rst(1) *    rst(2)  *    rst(3)  - p(i)
     end do
     
     return

   END SUBROUTINE funcv3D

   SUBROUTINE funcv2D(n,rst,fvec)
     
     !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
   
     integer :: i,n
     real (kind = rdf) ::p(2),v00(2),v10(2),v11(2),v01(2)
     common/intp2D/v00,v10,v11,v01,p
     real(kind = rdf) :: rst(n),fvec(n)
     
     do i=1,n
        fvec(i) =   v00(i) * (1-rst(1)) * (1-rst(2)) &
                  + v10(i) *     rst(1) * (1-rst(2)) &
                  + v11(i) *     rst(1) *    rst(2)  &
                  + v01(i) * (1-rst(1)) *    rst(2)  &
                  - p(i)
     end do

     return

   END SUBROUTINE funcv2D

   
   SUBROUTINE fdjac3D(n,rst,fjac)
   
     !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
     integer :: n
     real (kind = rdf) :: p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
     common/intp3D/v000,v100,v010,v001,v101,v011,v110,v111,p
     real(kind = rdf) :: rst(n),fjac(n,n)
   
     fjac(1,1) = -v000(1) * (1-rst(2)) * (1-rst(3)) + &
                  v100(1) * (1-rst(2)) * (1-rst(3)) - &
                  v010(1) *    rst(2)  * (1-rst(3)) - &
                  v001(1) * (1-rst(2)) *    rst(3)  + &
                  v101(1) * (1-rst(2)) *    rst(3)  - &
                  v011(1) *    rst(2)  *    rst(3)  + &
                  v110(1) *    rst(2)  * (1-rst(3)) + &
                  v111(1) *    rst(2)  *    rst(3)

     fjac(2,1) = -v000(2) * (1-rst(2)) * (1-rst(3)) + &
                  v100(2) * (1-rst(2)) * (1-rst(3)) - &
                  v010(2) *    rst(2)  * (1-rst(3)) - &
                  v001(2) * (1-rst(2)) *    rst(3)  + &
                  v101(2) * (1-rst(2)) *    rst(3)  - &
                  v011(2) *    rst(2)  *    rst(3)  + &
                  v110(2) *    rst(2)  * (1-rst(3)) + &
                  v111(2) *    rst(2)  *    rst(3)

     fjac(3,1) = -v000(3) * (1-rst(2)) * (1-rst(3)) + &
                  v100(3) * (1-rst(2)) * (1-rst(3)) - &
                  v010(3) *    rst(2)  * (1-rst(3)) - &
                  v001(3) * (1-rst(2)) *    rst(3)  + &
                  v101(3) * (1-rst(2)) *    rst(3)  - &
                  v011(3) *    rst(2)  *    rst(3)  + &
                  v110(3) *    rst(2)  * (1-rst(3)) + &
                  v111(3) *    rst(2)  *    rst(3)

     fjac(1,2) = -v000(1) * (1-rst(1)) * (1-rst(3)) - &
                  v100(1) *    rst(1)  * (1-rst(3)) + &
                  v010(1) * (1-rst(1)) * (1-rst(3)) - &
                  v001(1) * (1-rst(1)) *    rst(3)  - &
                  v101(1) *    rst(1)  *    rst(3)  + &
                  v011(1) * (1-rst(1)) *    rst(3)  + &
                  v110(1) *    rst(1)  * (1-rst(3)) + &
                  v111(1) *    rst(1)  *    rst(3)

     fjac(2,2) = -v000(2) * (1-rst(1)) * (1-rst(3)) - &
                  v100(2) *    rst(1)  * (1-rst(3)) + &
                  v010(2) * (1-rst(1)) * (1-rst(3)) - &
                  v001(2) * (1-rst(1)) *    rst(3)  - &
                  v101(2) *    rst(1)  *    rst(3)  + &
                  v011(2) * (1-rst(1)) *    rst(3)  + &
                  v110(2) *    rst(1)  * (1-rst(3)) + &
                  v111(2) *    rst(1)  *    rst(3)

     fjac(3,2) = -v000(3) * (1-rst(1)) * (1-rst(3)) - &
                  v100(3) *    rst(1)  * (1-rst(3)) + &
                  v010(3) * (1-rst(1)) * (1-rst(3)) - &
                  v001(3) * (1-rst(1)) *    rst(3)  - &
                  v101(3) *    rst(1)  *    rst(3)  + &
                  v011(3) * (1-rst(1)) *    rst(3)  + &
                  v110(3) *    rst(1)  * (1-rst(3)) + &
                  v111(3) *    rst(1)  *    rst(3)

     fjac(1,3) = -v000(1) * (1-rst(1)) * (1-rst(2)) - &
                  v100(1) *    rst(1)  * (1-rst(2)) - &
                  v010(1) * (1-rst(1)) *    rst(2)  + &
                  v001(1) * (1-rst(1)) * (1-rst(2)) + &
                  v101(1) *    rst(1)  * (1-rst(2)) + &
                  v011(1) * (1-rst(1)) *    rst(2)  - &
                  v110(1) *    rst(1)  *    rst(2)  + &
                  v111(1) *    rst(1)  *    rst(2)

     fjac(2,3) = -v000(2) * (1-rst(1)) * (1-rst(2)) - &
                  v100(2) *    rst(1)  * (1-rst(2)) - &
                  v010(2) * (1-rst(1)) *    rst(2)  + &
                  v001(2) * (1-rst(1)) * (1-rst(2)) + &
                  v101(2) *    rst(1)  * (1-rst(2)) + &
                  v011(2) * (1-rst(1)) *    rst(2)  - &
                  v110(2) *    rst(1)  *    rst(2)  + &
                  v111(2) *    rst(1)  *    rst(2)

     fjac(3,3) = -v000(3) * (1-rst(1)) * (1-rst(2)) - &
                  v100(3) *    rst(1)  * (1-rst(2)) - &
                  v010(3) * (1-rst(1)) *    rst(2)  + &
                  v001(3) * (1-rst(1)) * (1-rst(2)) + &
                  v101(3) *    rst(1)  * (1-rst(2)) + &
                  v011(3) * (1-rst(1)) *    rst(2)  - &
                  v110(3) *    rst(1)  *    rst(2)  + &
                  v111(3) *    rst(1)  *    rst(2)

   END SUBROUTINE fdjac3D


   SUBROUTINE fdjac2D(n,rst,fjac)
   
     !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
     integer :: n
     real (kind = rdf) ::p(2),v00(2),v10(2),v11(2),v01(2)
     common/intp2D/v00,v10,v11,v01,p
     real(kind = rdf) :: rst(n),fjac(n,n)
   

     fjac(1,1) = - v00(1) * (1-rst(2)) &
                 + v10(1) * (1-rst(2)) &
                 + v11(1) *    rst(2)  &
                 - v01(1) *    rst(2) 

     fjac(2,1) = - v00(2) * (1-rst(2)) &
                 + v10(2) * (1-rst(2)) &
                 + v11(2) *    rst(2)  &
                 - v01(2) *    rst(2)  

     fjac(1,2) = - v00(1) * (1-rst(1)) &
                 - v10(1) *    rst(1)  &
                 + v11(1) *    rst(1)  &
                 + v01(1) * (1-rst(1)) 

     fjac(2,2) = - v00(2) * (1-rst(1)) &
                 - v10(2) *    rst(1)  &
                 + v11(2) *    rst(1)  &
                 + v01(2) * (1-rst(1)) 


   END SUBROUTINE fdjac2D

   function cross_product(a, b) result(c)

      real( kind = rdf ), intent(in) :: a(3), b(3)
      real( kind = rdf ) :: c(3)
    
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)

   end function cross_product


   SUBROUTINE M33INV (A, AINV, OK_FLAG)

      !***********************************************************************************************************************************
      !  M33INV  -  Compute the inverse of a 3x3 matrix.
      !
      !  A       = input 3x3 matrix to be inverted
      !  AINV    = output 3x3 inverse of matrix A
      !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
      !***********************************************************************************************************************************

      IMPLICIT NONE

      real(kind = rdf), DIMENSION(3,3), INTENT(IN)  :: A
      real(kind = rdf), DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      real(kind = rdf), PARAMETER :: EPS = 1.0E-14
      real(kind = rdf) :: DET
      real(kind = rdf), DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = zero
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

   END SUBROUTINE M33INV


   function Invert4x4Matrix( A ) result( Ainv )
   
      implicit none
   
      real( kind = rdf ), intent(in)  :: A(4,4)
      real( kind = rdf ) :: Ainv(4,4)
   
      ! local variables
      
      integer :: i,j
      real( kind = rdf ) :: det
   
      Ainv(1,1) = A(2,2) * A(3,3) * A(4,4) - &
                  A(2,2) * A(3,4) * A(4,3) - &
                  A(3,2) * A(2,3) * A(4,4) + &
                  A(3,2) * A(2,4) * A(4,3) + &
                  A(4,2) * A(2,3) * A(3,4) - &
                  A(4,2) * A(2,4) * A(3,3)
   
      Ainv(2,1) =-A(2,1) * A(3,3) * A(4,4) + &
                  A(2,1) * A(3,4) * A(4,3) + &
                  A(3,1) * A(2,3) * A(4,4) - &
                  A(3,1) * A(2,4) * A(4,3) - &
                  A(4,1) * A(2,3) * A(3,4) + &
                  A(4,1) * A(2,4) * A(3,3)
   
      Ainv(3,1) = A(2,1) * A(3,2) * A(4,4) - &
                  A(2,1) * A(3,4) * A(4,2) - &
                  A(3,1) * A(2,2) * A(4,4) + &
                  A(3,1) * A(2,4) * A(4,2) + &
                  A(4,1) * A(2,2) * A(3,4) - &
                  A(4,1) * A(2,4) * A(3,2)
   
      Ainv(4,1) =-A(2,1) * A(3,2) * A(4,3) + &
                  A(2,1) * A(3,3) * A(4,2) + &
                  A(3,1) * A(2,2) * A(4,3) - &
                  A(3,1) * A(2,3) * A(4,2) - &
                  A(4,1) * A(2,2) * A(3,3) + &
                  A(4,1) * A(2,3) * A(3,2)
   
      Ainv(1,2) =-A(1,2) * A(3,3) * A(4,4) + &
                  A(1,2) * A(3,4) * A(4,3) + &
                  A(3,2) * A(1,3) * A(4,4) - &
                  A(3,2) * A(1,4) * A(4,3) - &
                  A(4,2) * A(1,3) * A(3,4) + &
                  A(4,2) * A(1,4) * A(3,3)
   
      Ainv(2,2) = A(1,1) * A(3,3) * A(4,4) - &
                  A(1,1) * A(3,4) * A(4,3) - &
                  A(3,1) * A(1,3) * A(4,4) + &
                  A(3,1) * A(1,4) * A(4,3) + &
                  A(4,1) * A(1,3) * A(3,4) - &
                  A(4,1) * A(1,4) * A(3,3)
   
      Ainv(3,2) =-A(1,1) * A(3,2) * A(4,4) + &
                  A(1,1) * A(3,4) * A(4,2) + &
                  A(3,1) * A(1,2) * A(4,4) - &
                  A(3,1) * A(1,4) * A(4,2) - &
                  A(4,1) * A(1,2) * A(3,4) + &
                  A(4,1) * A(1,4) * A(3,2)
   
      Ainv(4,2) = A(1,1) * A(3,2) * A(4,3) - &
                  A(1,1) * A(3,3) * A(4,2) - &
                  A(3,1) * A(1,2) * A(4,3) + &
                  A(3,1) * A(1,3) * A(4,2) + &
                  A(4,1) * A(1,2) * A(3,3) - &
                  A(4,1) * A(1,3) * A(3,2)
   
      Ainv(1,3) = A(1,2) * A(2,3) * A(4,4) - &
                  A(1,2) * A(2,4) * A(4,3) - &
                  A(2,2) * A(1,3) * A(4,4) + &
                  A(2,2) * A(1,4) * A(4,3) + &
                  A(4,2) * A(1,3) * A(2,4) - &
                  A(4,2) * A(1,4) * A(2,3)
   
      Ainv(2,3) =-A(1,1) * A(2,3) * A(4,4) + &
                  A(1,1) * A(2,4) * A(4,3) + &
                  A(2,1) * A(1,3) * A(4,4) - &
                  A(2,1) * A(1,4) * A(4,3) - &
                  A(4,1) * A(1,3) * A(2,4) + &
                  A(4,1) * A(1,4) * A(2,3)
   
      Ainv(3,3) = A(1,1) * A(2,2) * A(4,4) - &
                  A(1,1) * A(2,4) * A(4,2) - &
                  A(2,1) * A(1,2) * A(4,4) + &
                  A(2,1) * A(1,4) * A(4,2) + &
                  A(4,1) * A(1,2) * A(2,4) - &
                  A(4,1) * A(1,4) * A(2,2)
   
      Ainv(4,3) =-A(1,1) * A(2,2) * A(4,3) + &
                  A(1,1) * A(2,3) * A(4,2) + &
                  A(2,1) * A(1,2) * A(4,3) - &
                  A(2,1) * A(1,3) * A(4,2) - &
                  A(4,1) * A(1,2) * A(2,3) + &
                  A(4,1) * A(1,3) * A(2,2)
   
      Ainv(1,4) =-A(1,2) * A(2,3) * A(3,4) + &
                  A(1,2) * A(2,4) * A(3,3) + &
                  A(2,2) * A(1,3) * A(3,4) - &
                  A(2,2) * A(1,4) * A(3,3) - &
                  A(3,2) * A(1,3) * A(2,4) + &
                  A(3,2) * A(1,4) * A(2,3)
   
      Ainv(2,4) = A(1,1) * A(2,3) * A(3,4) - &
                  A(1,1) * A(2,4) * A(3,3) - &
                  A(2,1) * A(1,3) * A(3,4) + &
                  A(2,1) * A(1,4) * A(3,3) + &
                  A(3,1) * A(1,3) * A(2,4) - &
                  A(3,1) * A(1,4) * A(2,3)
   
      Ainv(3,4) =-A(1,1) * A(2,2) * A(3,4) + &
                  A(1,1) * A(2,4) * A(3,2) + &
                  A(2,1) * A(1,2) * A(3,4) - &
                  A(2,1) * A(1,4) * A(3,2) - &
                  A(3,1) * A(1,2) * A(2,4) + &
                  A(3,1) * A(1,4) * A(2,2)
   
      Ainv(4,4) = A(1,1) * A(2,2) * A(3,3) - &
                  A(1,1) * A(2,3) * A(3,2) - &
                  A(2,1) * A(1,2) * A(3,3) + &
                  A(2,1) * A(1,3) * A(3,2) + &
                  A(3,1) * A(1,2) * A(2,3) - &
                  A(3,1) * A(1,3) * A(2,2)
   
      det = A(1,1) * Ainv(1,1) + A(1,2) * Ainv(2,1) + &
            A(1,3) * Ainv(3,1) + A(1,4) * Ainv(4,1)
   
      if ( abs( det ) < eps_sims ) then
   
         print *, 'Determinant of the matrix equal to ZERO'
         stop
   
      end if
   
      Ainv = Ainv / det
   
   end function Invert4x4Matrix


   function SmoothKernelFunction ( x ) result( y )

      real (kind = rdf) :: x,y


      if ( x > eps_sims ) then
         y = exp( -one / x )
      else
         y = zero
      end if

   end function SmoothKernelFunction

   function SmoothStepFunction_exp ( x ) result( y )

      real (kind = rdf) :: x,y


      if ( x <= eps_sims ) then
         y = zero
      else if ( x > eps_sims .and. x < one ) then
         y =                                       SmoothKernelFunction( x )   /  & 
               ( SmoothKernelFunction( one - x ) + SmoothKernelFunction( x ) )
      else
         y = one
      end if

   end function SmoothStepFunction_exp


   function SmoothStepFunction ( x ) result( y )

      real (kind = rdf) :: x,y

      y = max ( zero , min ( one , x ) )

   end function SmoothStepFunction

end module InterpolationMethods

























