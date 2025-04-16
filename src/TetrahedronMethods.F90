module TetrahedronMethods
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Set of subroutines and functions to perform marching tetrahedra (algorithm to get the
   ! isosurface of a function inside a tetrahedron given its nodal values), for 
   ! computing the distance from the nodes to the isosurface and to obtain the water
   ! volume within a tetrahedron.
   !
   ! All these routines are used for the geometric reinisialisation algorithm proposed in:
   ! 
   ! Ausas, R. F., Dari, E. A., & Buscaglia, G. C. (2011). A geometric mass‐preserving 
   ! redistancing scheme for the level set function. International journal for numerical 
   ! methods in fluids, 65(8), 989-1010.
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
   ! Jorge Sandoval, UoE/PUC. Edinburgh, April 4th, 2022.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

   use precision

   implicit none

   contains

   subroutine MarchingTetrahedron(  vertices_coordinates       , &
                                    phi_tetrahedron            , &
                                    vertices_isosurfaces       , &
                                    ntriangles                 , &
                                    IsosurfaceArea             , &
                                    WaterVolume                , &
                                    Distances                    &
                                  )
   
      implicit none

      ! This subroutine Performs the marching tetrahedrons algorithm
      ! on one tetrahedron, whose verticies are specified in the array
      ! vertices_coordinates.  The function values associated with each vertex 
      ! is given in the phi_tetrahedronn array and the contouring level in the 
      ! real variable level.  

      ! input variables
      real(kind = rdf), dimension(4,3) , intent(in) ::  vertices_coordinates
      real(kind = rdf), dimension (4)  , intent(in) ::  phi_tetrahedron
   
      ! output variables
   
      ! vertices_isosurfaces(triangle, vertex, coordinate component)
      ! for example vertices_isosurfaces(2,1,3) gives me the z-coordinate of the first
      ! vertex of the second triangle
      ! 
      ! ntriangles is the number of isosurfaces extracted using MarchingTetrahedron
   
      real(kind = rdf), dimension (2,3,3), intent(out) :: vertices_isosurfaces
   
      ! ntriangles is a flag variable, where each component tells if there is a triangle in
      ! the slot 1 or 2
      ! if there is a triangle, that position is set to 1. Otherwise, it is set to -1
      
      integer          ,               intent(out) :: ntriangles
      real(kind = rdf) ,               intent(out) :: IsosurfaceArea
      real(kind = rdf) ,               intent(out) :: WaterVolume
      real(kind = rdf) , dimension(4), intent(out) :: Distances
   
      ! Edges vertices connections
      !
      !          3
      !         /|\
      !        / | \
      !       /  2  \
      !      / .´ `. \
      !     /.´     `.\
      !    0´---------`1  
      
      ! order = (/2,1/) writes the matrices row by row (as the arrays are written below)
      integer, dimension(6,2), parameter :: TetrahedronEdgeConnection &
           = reshape( (/ 0,1, &
                         1,2, &
                         2,0, &
                         0,3, &
                         1,3, &
                         2,3 /), shape = (/6,2/), &
                                 order = (/2,1/))
   
      ! Edges numeration
      !
      !          *
      !         /|\
      !        / 5 \
      !       3  |  4
      !      /  .*.  \
      !     /.2´   `1.\
      !    *´----0----`*  
     
      integer, dimension(16,6), parameter :: TetrahedronTriangles &
           = reshape( (/ -1, -1, -1, -1, -1, -1, &
                          0,  3,  2, -1, -1, -1, &
                          0,  1,  4, -1, -1, -1, &
                          1,  4,  2,  2,  4,  3, &
                          1,  2,  5, -1, -1, -1, &
                          0,  3,  5,  0,  5,  1, &
                          0,  2,  5,  0,  5,  4, &
                          5,  4,  3, -1, -1, -1, &
                          3,  4,  5, -1, -1, -1, &
                          4,  5,  0,  5,  2,  0, &
                          1,  5,  0,  5,  3,  0, &
                          5,  2,  1, -1, -1, -1, &
                          3,  4,  2,  2,  4,  1, &
                          4,  1,  0, -1, -1, -1, &
                          2,  3,  0, -1, -1, -1, &
                         -1, -1, -1, -1, -1, -1/) ,shape = (/16,6/), &
                                                   order = (/2,1 /) )
   
      ! In bit representation, TetrahedronEdgeFlags tells (see Edges 
      ! numeration) which edges contain intersection for every case
      ! (bits are set and read from right to left)
      !
      !                    Edges:   6 5 4 3 2 1
      !                            |-|-|-|-|-|-|
      ! TetrahedronEdgeFlags(1)  = |0|0|0|0|0|0|
      ! TetrahedronEdgeFlags(2)  = |0|0|1|1|0|1|
      ! TetrahedronEdgeFlags(3)  = |0|1|0|0|1|1|
      ! TetrahedronEdgeFlags(4)  = |0|1|1|1|1|0| *
      ! TetrahedronEdgeFlags(5)  = |1|0|0|1|1|0|
      ! TetrahedronEdgeFlags(6)  = |1|0|1|0|1|1| *
      ! TetrahedronEdgeFlags(7)  = |1|1|0|1|0|1| *
      ! TetrahedronEdgeFlags(8)  = |1|1|1|0|0|0|
      ! TetrahedronEdgeFlags(9)  = |1|1|1|0|0|0|
      ! TetrahedronEdgeFlags(10) = |1|1|0|1|0|1| *
      ! TetrahedronEdgeFlags(11) = |1|0|1|0|1|1| *
      ! TetrahedronEdgeFlags(12) = |1|0|0|1|1|0|
      ! TetrahedronEdgeFlags(13) = |0|1|1|1|1|0| *
      ! TetrahedronEdgeFlags(14) = |0|1|0|0|1|1|
      ! TetrahedronEdgeFlags(15) = |0|0|1|1|0|1| 
      ! TetrahedronEdgeFlags(16) = |0|0|0|0|0|0|
   
      ! * cases represents the isosurface by means of two triangles
   
      ! 'int(Z'**') is a way to write an integer by setting its bits

      integer, dimension(16), parameter :: TetrahedronEdgeFlags = &
                          (/int(Z'00'),int(Z'0d'),int(Z'13'),int(Z'1e'), &
                            int(Z'26'),int(Z'2b'),int(Z'35'),int(Z'38'), &
                            int(Z'38'),int(Z'35'),int(Z'2b'),int(Z'26'), &
                            int(Z'1e'),int(Z'13'),int(Z'0d'),int(Z'00')/)
   
      ! local variables
      integer :: FlagIndex, EdgeFlag
      integer :: VertexLoop, EdgeLoop, TriangleLoop ! loop counters
      real(kind = rdf), dimension(6,3) :: EdgeIntersection
      integer :: StartVertex, EndVertex, CurrentVertex, cont
      real(kind = rdf) :: s
      real(kind = rdf) :: DeltaPhiDenomInterpolation
         
      integer , dimension(4) :: PhiSign
      integer                :: SumAux, nZeroVertices

      ! output variables initilisation

      vertices_isosurfaces = zero
      ntriangles           = 0
      IsosurfaceArea       = zero
      WaterVolume          = zero
      Distances            =  9.9_rdf

      ! local variables initilisation
      EdgeIntersection      = -9.9_rdf
      s                     = zero

      ! -  - - - - - - - - - - execution - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
      ! First scan the verticies of the tetrahedron testing to see if
      ! the the iso-surface intersects it.
      

      FlagIndex = 0 ! in bit representation FlagIndex = 0000 (the other 28 bits are 0)
                    ! where each bit represents a tetrahedron corner to be potentially
                    ! set to 1 by executing ibset
   
      ! With this flag variable we can check special cases, as when the isosurface crosses an entire
      ! edge or face

      PhiSign = 0

      where ( phi_tetrahedron   >  eps_sims ) PhiSign =  1
      where ( phi_tetrahedron   < -eps_sims ) PhiSign = -1

      if ( all( PhiSign == 1 ) ) then
      
         WaterVolume = GetTetrahedronVolume(vertices_coordinates)      
         return

      else if ( all( PhiSign == -1 ) ) then
      
         WaterVolume = zero
         return

      else if ( any( PhiSign == 0 ) .and. ( all( PhiSign >= 0 ) .or. all( PhiSign <= 0 ) ) ) then

         nZeroVertices = count( PhiSign == 0 ) ! Count how many zero vertices there are in the tetrahedron

         select case ( nZeroVertices )
   
            case(1) ! nZeroVertices == 1 : ONE vertex is zero      --> One of the VERTICES is the isosurface

               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               !                                                     ! 
               !                           v3                        ! 
               !                           .                         !    
               !                          / \                        !  
               !                         /   \                       ! 
               !                        /     \                      ! 
               !                       /       \                     !
               !                      /         \                    !    
               !                  \  /           \                   !    
               !                   \/_____________\                  ! 
               !                v1  \              v2                !  
               !                     \                               !    
               !                  Free surface                       !          
               !                                                     !          
               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - 

               ntriangles = -1
               IsosurfaceArea = zero

               ! Isosurface on VERTEX 1
               if ( PhiSign(1) == 0 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(1,:)

                  WaterVolume = zero
                  if ( PhiSign(2) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! Isosurface on VERTEX 2
               if ( PhiSign(2) == 0 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(2,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! Isosurface on VERTEX 3
               if ( PhiSign(3) == 0 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(3,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! Isosurface on VERTEX 4
               if ( PhiSign(4) == 0 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(4,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(4,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if


            case(2) ! nZeroVertices == 2 : TWO vertices are zero   --> One of the EDGES is the isosurface
   
               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               !                                                     ! 
               !                           v3                        ! 
               !                           .                         !    
               !                          / \                        !  
               !                         /   \                       ! 
               !                        /     \                      ! 
               !                       /       \                     !
               !                      /         \                    !    
               !  Free surface       /           \                   !    
               !          _________ /_____________\ ___________      ! 
               !                v1                 v2                !  
               !                                                     !    
               !                                                     !          
               !                                                     !          
               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - 

               ntriangles = -2
               IsosurfaceArea = zero

               ! The free surface pass by the edge formed by VERTICES 1 and 2 (Edge 1)
               if ( PhiSign(1) == 0 .and. PhiSign(2) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(2,:)

                  WaterVolume = zero
                  if ( PhiSign(3) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! The free surface pass by the edge formed by VERTICES 2 and 3 (Edge 2)
               if ( PhiSign(2) == 0 .and. PhiSign(3) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(3,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! The free surface pass by the edge formed by VERTICES 3 and 1 (Edge 3)
               if ( PhiSign(3) == 0 .and. PhiSign(1) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(1,:)

                  WaterVolume = zero
                  if ( PhiSign(2) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! The free surface pass by the edge formed by VERTICES 3 and 1 (Edge 4)
               if ( PhiSign(1) == 0 .and. PhiSign(4) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(4,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  WaterVolume = zero
                  if ( PhiSign(2) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! The free surface pass by the edge formed by VERTICES 2 and 4 (Edge 5)
               if ( PhiSign(2) == 0 .and. PhiSign(4) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(4,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

               ! The free surface pass by the edge formed by VERTICES 3 and 4 (Edge 6)
               if ( PhiSign(3) == 0 .and. PhiSign(4) == 0) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(4,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

            case(3) ! nZeroVertices == 3 : THREE vertices are zero --> One of the FACES is the isosurface

               ! - - - - - - - - - - - - - - - - - - - - - !
               !                                           !
               !                  *    Free surface        !
               !                 /|\  /                    !  
               !                / * \/                     !
               !               / .../\                     !  
               !              / .../. \                    !
               !             /.........\                   !
               !            *´---------`*                  !  
               !                                           !
               ! - - - - - - - - - - - - - - - - - - - - - !
    
               ntriangles = 1

               ! The face formed by vertices 2,3,4 are the isosurface
               if ( PhiSign(1) == 1 .or. PhiSign(1) == -1 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  IsosurfaceArea = GetTriangleArea(vertices_isosurfaces(1,:,:))

                  WaterVolume = zero
                  if ( PhiSign(1) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if


               ! The face formed by vertices 1,3,4 are the isosurface
               if ( PhiSign(2) == 1 .or. PhiSign(2) == -1 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(3,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  IsosurfaceArea = GetTriangleArea(vertices_isosurfaces(1,:,:))

                  WaterVolume = zero
                  if ( PhiSign(2) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if


               ! The face formed by vertices 1,2,4 are the isosurface
               if ( PhiSign(3) == 1 .or. PhiSign(3) == -1 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(4,:)

                  IsosurfaceArea = GetTriangleArea(vertices_isosurfaces(1,:,:))

                  WaterVolume = zero
                  if ( PhiSign(3) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if


               ! The face formed by vertices 1,2,3 are the isosurface
               if ( PhiSign(4) == 1 .or. PhiSign(4) == -1 ) then

                  vertices_isosurfaces(1,1,:) = vertices_coordinates(1,:)
                  vertices_isosurfaces(1,2,:) = vertices_coordinates(2,:)
                  vertices_isosurfaces(1,3,:) = vertices_coordinates(3,:)

                  IsosurfaceArea = GetTriangleArea(vertices_isosurfaces(1,:,:))

                  WaterVolume = zero
                  if ( PhiSign(4) == 1 ) WaterVolume = GetTetrahedronVolume( vertices_coordinates ) 

                  return

               end if

         end select


      else ! Standar Marching Tetrahedron algorithm

         do VertexLoop = 1,4

            if ( PhiSign( VertexLoop ) == -1 ) then
   
!            if (   phi_tetrahedron(VertexLoop) < -eps_sims .and. &
!                 - phi_tetrahedron(VertexLoop) >  eps_sims           ) then
   
               ! FI means FlagIndex
               !-----------------------------------------------------------
               !  CASE   |   v1  v2  v3  v4  |  FI (value) | FI(binary)   |
               !         |   v0  v1  v2  v3  |             | (read <--)   |
               !-----------------------------------------------------------
               !    1       | + | + | + | + |      0            0000        
               !    2       | - | + | + | + |      1            0001        
               !    3       | + | - | + | + |      2            0010        
               !    4       | - | - | + | + |*     3            0011         
               !    5       | + | + | - | + |      4            0100        
               !    6       | - | + | - | + |*     5            0101         
               !    7       | + | - | - | + |*     6            0110         
               !    8       | - | - | - | + |      7            0111        
               !    9       | + | + | + | - |      8            1000        
               !    10      | - | + | + | - |*     9            1001         
               !    11      | + | - | + | - |*     10           1010         
               !    12      | - | - | + | - |      11           1011        
               !    13      | + | + | - | - |*     12           1100         
               !    14      | - | + | - | - |      13           1101        
               !    15      | + | - | - | - |      14           1110        
               !    16      | - | - | - | - |      15           1111        
               !-----------------------------  ----------------------------
               ! * cases ry means of two tria  epresents the isosurface bngles
      
      
               !                                                      (binary)               (int)
               ! if phi_tetrahedron values are all positive FlagIndex = 0000   -> FlagIndex =  0
               ! if phi_tetrahedron values are all negative FlagIndex = 1111   -> FlagIndex =  15
               
               ! if a phi_tetrahedron value is exactly zero, that corner bit is not set to one 
               ! either (phi_tetrahedron(VertexLoop) < zero .and. -phi_tetrahedron(VertexLoop) > zero
               ! condition)
      
               ! Flag index, as an integer, is also an array of 32 bits.
               ! with this procedure I set to 1 the bit number VertexLoop-1
               ! (0 to 3)
      
               FlagIndex = ibset(FlagIndex,VertexLoop-1)
            
            end if
   
         end do
      
      end if
      
      ! Now we can find the Edges coresponding to FLagIndex, note
      ! FlagIndex goes from 0 to 15, whilst the array is numbered from
      ! 1 to 16.

      ! if FlagIndex == 0 , the tetrahedron is entirely in the water phase
      ! if FlagIndex == 15, the tetrahedron is entirely in the air phase
   
      if (FlagIndex == 0 ) WaterVolume = GetTetrahedronVolume(vertices_coordinates)
      if (FlagIndex == 15) WaterVolume = zero
      
      EdgeFlag = TetrahedronEdgeFlags(FlagIndex+1)
   
      ! If the tetrahedron is entirely inside or outside of the surface,
      ! then there will be no intersections (TetrahedronEdgeFlags(1) or
      ! TetrahedronEdgeFlags(16) )
      
      if (EdgeFlag == 0) return
   
      ! If we have got this far then their is an intersection
      ! between the tetrahedron and the iso-surface.
   
      do EdgeLoop=1,6
   
         if ( btest ( EdgeFlag,EdgeLoop-1 ) ) then
   
            ! If the bit is set to 1 then this edge is intersected by 
            ! the surface. The first job is to locate the coordinate of the
            ! cut, for interpolation of the vertex positions to locate
            ! the coordinates of the cut. 0 < s < 1 is the normalised 
            ! position of the cut along the edge.
            
   
            StartVertex = TetrahedronEdgeConnection(EdgeLoop,1) + 1 ! vertex 1 to 4 
            EndVertex   = TetrahedronEdgeConnection(EdgeLoop,2) + 1 ! vertex 1 to 4 
   
            DeltaPhiDenomInterpolation = phi_tetrahedron(EndVertex)-phi_tetrahedron(StartVertex)
            
            ! linear interpolation coefficient:  0 < s < 1

            if ( abs( DeltaPhiDenomInterpolation ) < eps_sims ) then
               s = zero
               print *, 'Isosurface on more than one corner'
            else
               s  = phi_tetrahedron(EndVertex) / DeltaPhiDenomInterpolation 
            end if

            ! Now store the interpolated coordinates
   
            ! x - coordinates
            EdgeIntersection(EdgeLoop,1) =  vertices_coordinates(EndVertex,1) &
                 + s * (vertices_coordinates(StartVertex,1)-vertices_coordinates(EndVertex,1))
            
            ! y - coordinates
            EdgeIntersection(EdgeLoop,2) =  vertices_coordinates(EndVertex,2) &
                 + s * (vertices_coordinates(StartVertex,2)-vertices_coordinates(EndVertex,2))
   
            ! z - coordinates
            EdgeIntersection(EdgeLoop,3) =  vertices_coordinates(EndVertex,3) &
                 + s * (vertices_coordinates(StartVertex,3)-vertices_coordinates(EndVertex,3))
            
         end if
      end do
   
      ! Finally add any triangles that have been found to the
      ! surface description.  There can be 1, or 2 triangles per
      ! tetrahedron (the 0 - triangles case produce a return when we 
      ! execute: "if (EdgeFlag == 0) return" some lines above)
   
      do TriangleLoop = 1,2
         
         ! If the index of the node of the triangle is -1, the triangle 1 and 2 don't 
         ! exist 
         
         if (TetrahedronTriangles( FlagIndex + 1 , 3*(TriangleLoop-1) + 1) < 0) exit
   
         ntriangles = ntriangles + 1

         do VertexLoop=1,3
   
            ! CurrentVertex is a number that goes from 0 to 5 and gives me the number of the edge
            ! where the currently evaluated vertex lies on. To get the vertex coordinate, I need
            ! to evulate EdgeIntersection(CurrentVertex+1,1:3) (CurrentVertex+1 because 
            ! EdgeIntersection goes from 1 to 6)
   
            !                                       (1:16)    ,    (1:3) or (4:6)
            CurrentVertex = TetrahedronTriangles( FlagIndex+1 , 3*(TriangleLoop-1) + VertexLoop )
   
            ! x - coordinate of the intersection of the vertex VertexLoop of the triangle TriangleLoop
            vertices_isosurfaces(TriangleLoop,VertexLoop,1) = EdgeIntersection(CurrentVertex+1,1) 

            ! y - coordinate of the intersection of the vertex VertexLoop of the triangle TriangleLoop
            vertices_isosurfaces(TriangleLoop,VertexLoop,2) = EdgeIntersection(CurrentVertex+1,2) 

            ! z - coordinate of the intersection of the vertex VertexLoop of the triangle TriangleLoop
            vertices_isosurfaces(TriangleLoop,VertexLoop,3) = EdgeIntersection(CurrentVertex+1,3) 

         end do ! VertexLoop

         ! Add the triangle area to the IsosurfaceArea
         IsosurfaceArea = IsosurfaceArea + &
                 GetTriangleArea(vertices_isosurfaces(TriangleLoop,:,:))

      end do ! TriangleLoop
      
      if (FlagIndex == 0 ) then 
         WaterVolume = GetTetrahedronVolume(vertices_coordinates)
      else if (FlagIndex == 15) then 
         WaterVolume = zero
      else
         WaterVolume = GetWaterVolume( vertices_coordinates , EdgeIntersection , FlagIndex)
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! I never use this distances in practice
      !
      !Distances   = GetNodesIsosurfaceDistances( vertices_coordinates, &
      !                                           vertices_isosurfaces, &
      !                                           ntriangles )   
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   end subroutine MarchingTetrahedron
   
   function GetTetrahedronVolume(vertices_coordinates)
      
      ! This function computes de volume of a tetrahedron, given the coordinates of
      ! its nodes using the determinant formula V = 1/6 * |det((A-D),(B-D),(C-D))|, 
      ! where A, B, C, D are the tetrahedron vertices

      implicit none
   
      real(kind = rdf), dimension(4,3) :: vertices_coordinates
      real(kind = rdf), dimension(4,4) :: AuxMatrix
      real(kind = rdf) :: GetTetrahedronVolume
   
      ! Auxiliary Matrix for computing the tetrahedron volume using
      ! the determinant formula 
   
      AuxMatrix(1,1:3) = vertices_coordinates(1,1:3); AuxMatrix(1,4) = one;
      AuxMatrix(2,1:3) = vertices_coordinates(2,1:3); AuxMatrix(2,4) = one;
      AuxMatrix(3,1:3) = vertices_coordinates(3,1:3); AuxMatrix(3,4) = one;
      AuxMatrix(4,1:3) = vertices_coordinates(4,1:3); AuxMatrix(4,4) = one;
   
      GetTetrahedronVolume = one/six * abs(M44Det(AuxMatrix)) 
   
   end function GetTetrahedronVolume
   

   function GetTriangleArea(vertices_triangle)
      
      ! This function computes de triangle area using the nodes of its coordinates

      implicit none

      real(kind = rdf), dimension(3,3), intent(in) :: vertices_triangle
      real(kind = rdf), dimension(3) :: A, B, C
      real(kind = rdf) :: GetTriangleArea

      A = vertices_triangle(1,:)
      B = vertices_triangle(2,:)
      C = vertices_triangle(3,:)
   
      GetTriangleArea = one_half * norm2(CrossProduct((B-A),(C-A))) 
   
   end function GetTriangleArea


   function GetWaterVolume(vertices_coordinates, EdgeIntersection, FlagIndex)
      
      ! This function returns the water volume of a tetrahedron given its vertex
      ! coordinates, edge isosurface intersections and case number (FlagIndex).
      ! The volume of water is the region where phi>0. This is not included in 
      ! the method as an input because FlagIndex tells what kind of intersection
      ! we have in the tetrahedron.

      ! The cases where the isosurface is composed of two triangles, the method
      ! for computing the volume is subdividing the water volume in new tetrahedra
      ! and adding the subvolumes. The election of that subdivision is arbitrary.

      implicit none
   
      real(kind = rdf), dimension(4,3), intent(in) :: vertices_coordinates
      real(kind = rdf), dimension(6,3), intent(in) :: EdgeIntersection
      integer, intent(in) :: FlagIndex
   
      ! Local Variables
   
      real(kind = rdf), dimension(4,3) :: vertices_coordinates_aux1
      real(kind = rdf), dimension(4,3) :: vertices_coordinates_aux2
      real(kind = rdf), dimension(4,3) :: vertices_coordinates_aux3
   
      integer :: IsosurfaceCase
      real(kind = rdf) :: TetrahedronVolume, AuxVolume
      real(kind = rdf) :: GetWaterVolume
   
      ! Vertices and Edges numeration
      !
      !         v4   
      !        /|\
      !       / 6 \
      !      4  |  5
      !     /   v3  \
      !    /.3´   `2.\
      !  v1´----1----`v2  
   
      IsosurfaceCase = FlagIndex + 1
   
      ! Tetrahedron Volume using its 4 coordinates
      TetrahedronVolume = GetTetrahedronVolume(vertices_coordinates) 

      AuxVolume = zero
   
      select case (IsosurfaceCase)
   
         case(1,16) ! empty or full tetrahedron
   
            select case(IsosurfaceCase)
               case(1)  ! (v1,v2,v3,v4) = ( + , + , + , + )
                  GetWaterVolume = TetrahedronVolume
               case(16) ! (v1,v2,v3,v4) = ( - , - , - , - )
                  GetWaterVolume = zero
            end select 
   
         case(2,15) ! (1 triangle isosurface)
   
            ! Auxiliary Tetrahedron vertices: (v1, intersections edges 1,3,4)
   
            vertices_coordinates_aux1(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(4,:)
   
            AuxVolume = GetTetrahedronVolume(vertices_coordinates_aux1)
   
            select case(IsosurfaceCase)
               case(2)  ! (v1,v2,v3,v4) = ( - , + , + , + )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(15) ! (v1,v2,v3,v4) = ( + , - , - , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(3,14) ! (1 triangle isosurface)
   
            ! Auxiliary Tetrahedron vertices: (v2, intersections edges 1,2,5)
   
            vertices_coordinates_aux1(1,:) = vertices_coordinates(2,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(5,:)
   
            AuxVolume = GetTetrahedronVolume(vertices_coordinates_aux1)
   
            select case(IsosurfaceCase)
               case(3)  ! (v1,v2,v3,v4) = ( + , - , + , + ) 
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(14) ! (v1,v2,v3,v4) = ( - , + , - , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(4,13) ! (2 triangles isosurface)
   
            ! Auxiliary Tetrahedron 1 vertices: (v1, intersections edges 3,4,5)
            vertices_coordinates_aux1(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(4,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(5,:)
   
            ! Auxiliary Tetrahedron 2 vertices: (v1, intersections edges 2,3,5)
   
            vertices_coordinates_aux2(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux2(2,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux2(3,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux2(4,:) = EdgeIntersection(5,:)
   
            ! Auxiliary Tetrahedron 3 vertices: (v1, v2, intersection edges 2,5)
   
            vertices_coordinates_aux3(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux3(2,:) = vertices_coordinates(2,:)
            vertices_coordinates_aux3(3,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux3(4,:) = EdgeIntersection(5,:)
   
            AuxVolume =    GetTetrahedronVolume( vertices_coordinates_aux1 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux2 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux3 )
   
            select case(IsosurfaceCase)
               case(4)  ! (v1,v2,v3,v4) = ( - , - , + , + )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(13) ! (v1,v2,v3,v4) = ( + , + , - , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(5,12) ! (1 triangle isosurface)
   
            ! Auxiliary Tetrahedron vertices: (v3, intersections edges 2,3,6)
   
            vertices_coordinates_aux1(1,:) = vertices_coordinates(3,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(6,:)
   
            AuxVolume = GetTetrahedronVolume(vertices_coordinates_aux1)
   
            select case(IsosurfaceCase)
               case(5)  ! (v1,v2,v3,v4) = ( + , + , - , + )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(12) ! (v1,v2,v3,v4) = ( - , - , + , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(6,11) ! (2 triangles isosurface) 
   
            ! Auxiliary Tetrahedron 1 vertices: (v1, intersections edges 1,2,6)
            vertices_coordinates_aux1(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(6,:)
   
            ! Auxiliary Tetrahedron 2 vertices: (v1, intersections edges 1,4,6)
   
            vertices_coordinates_aux2(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux2(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux2(3,:) = EdgeIntersection(4,:)
            vertices_coordinates_aux2(4,:) = EdgeIntersection(6,:)
   
            ! Auxiliary Tetrahedron 3 vertices: (v1, v3, intersection edges 2,6)
   
            vertices_coordinates_aux3(1,:) = vertices_coordinates(1,:)
            vertices_coordinates_aux3(2,:) = vertices_coordinates(3,:)
            vertices_coordinates_aux3(3,:) = EdgeIntersection(2,:)
            vertices_coordinates_aux3(4,:) = EdgeIntersection(6,:)
   
            AuxVolume =    GetTetrahedronVolume( vertices_coordinates_aux1 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux2 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux3 )
   
            select case(IsosurfaceCase)
               case(6)  ! (v1,v2,v3,v4) = ( - , + , - , + )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(11) ! (v1,v2,v3,v4) = ( + , - , + , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(7,10) ! (2 triangles isosurface)
   
            ! Auxiliary Tetrahedron 1 vertices: (v2, intersections edges 1,3,6)
            vertices_coordinates_aux1(1,:) = vertices_coordinates(2,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(6,:)
   
            ! Auxiliary Tetrahedron 2 vertices: (v2, intersections edges 1,5,6)
   
            vertices_coordinates_aux2(1,:) = vertices_coordinates(2,:)
            vertices_coordinates_aux2(2,:) = EdgeIntersection(1,:)
            vertices_coordinates_aux2(3,:) = EdgeIntersection(5,:)
            vertices_coordinates_aux2(4,:) = EdgeIntersection(6,:)
   
            ! Auxiliary Tetrahedron 3 vertices: (v2, v3, intersection edges 3,6)
   
            vertices_coordinates_aux3(1,:) = vertices_coordinates(2,:)
            vertices_coordinates_aux3(2,:) = vertices_coordinates(3,:)
            vertices_coordinates_aux3(3,:) = EdgeIntersection(3,:)
            vertices_coordinates_aux3(4,:) = EdgeIntersection(6,:)
   
            AuxVolume =    GetTetrahedronVolume( vertices_coordinates_aux1 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux2 ) &
                        +  GetTetrahedronVolume( vertices_coordinates_aux3 )
   
            select case(IsosurfaceCase)
               case(7)  ! (v1,v2,v3,v4) = ( + , - , - , + )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
               case(10) ! (v1,v2,v3,v4) = ( - , + , + , - )
                  GetWaterVolume = AuxVolume
            end select 
   
         case(8,9) ! (1 triangle isosurface)  
   
            ! Auxiliary Tetrahedron vertices: (v4, intersections edges 4,5,6)
   
            vertices_coordinates_aux1(1,:) = vertices_coordinates(4,:)
            vertices_coordinates_aux1(2,:) = EdgeIntersection(4,:)
            vertices_coordinates_aux1(3,:) = EdgeIntersection(5,:)
            vertices_coordinates_aux1(4,:) = EdgeIntersection(6,:)
   
            AuxVolume = GetTetrahedronVolume(vertices_coordinates_aux1)
   
            ! these cases are the only ones where the order 
            ! TetrahedronVolume - AuxVolume / AuxVolume is inverted
            ! compared with the previous ones (it's because the bits 
            ! configuration)
   
            select case(IsosurfaceCase)
               case(8) ! (v1,v2,v3,v4) = ( - , - , - , + )
                  GetWaterVolume = AuxVolume
               case(9) ! (v1,v2,v3,v4) = ( + , + , + , - )
                  GetWaterVolume = TetrahedronVolume - AuxVolume
            end select 
   
      end select

      if ( GetWaterVolume < -eps_sims ) then

         print *, ' - - - - - - - - - - - - - - - - - - - - - -'
         print *, ' case: ', IsosurfaceCase
         print *, ' Negative WaterVolume = ', GetWaterVolume
         print *, ' - - - - - - - - - - - - - - - - - - - - - -'
         stop

      end if

   end function GetWaterVolume
   
   function GetNodesIsosurfaceDistances(  vertices_coordinates, &
                                          vertices_isosurfaces, &
                                          ntriangles )
      
      ! This function computes the distance from every node of a tetrahedron
      ! to the free-surface. 
      ! For doing that it uses the function DistancePointTriange, which 
      ! computing the minimum distance to a triangle in the space.

      implicit none
   
      real(kind = rdf), dimension(4,3), intent(in)  ::  vertices_coordinates
      real(kind = rdf), dimension (2,3,3), intent(in) :: vertices_isosurfaces
      integer, intent(in) :: ntriangles
   
      real(kind = rdf), dimension(4) :: GetNodesIsosurfaceDistances 
   
      ! local variables
      real(kind = rdf), dimension(3) :: CurrentVertexPosition
      real(kind = rdf) :: Distance1, Distance2
      integer :: VertexLoop
   
      do VertexLoop = 1,4
   
         ! Vertex distance initialisation
         Distance1 = 999.9_rdf ! distance from the vertex to the isosurface 1
         Distance2 = 999.9_rdf ! distance from the vertex to the isosurface 2
   
         CurrentVertexPosition(:) = vertices_coordinates(VertexLoop,:)
   
         if( ntriangles == 1 ) then ! triangle 1 is a valid isosurface
            Distance1 = DistancePointTriangle( vertices_isosurfaces(1,:,:) , &
                                               CurrentVertexPosition)

         else if( ntriangles == 2 ) then ! triangle 2 is a valid isosurface

            Distance1 = DistancePointTriangle( vertices_isosurfaces(1,:,:) , &
                                               CurrentVertexPosition)
            Distance2 = DistancePointTriangle( vertices_isosurfaces(2,:,:) , &
                                               CurrentVertexPosition)
         !else
         !   Distance1 = zero
         !   Distance2 = zero
         end if
   
         GetNodesIsosurfaceDistances(VertexLoop) = minval((/Distance1, &
                                                            Distance2/))
      end do ! VertexLoop
   
   end function GetNodesIsosurfaceDistances
   
   
   function DistancePointTriangle(vertices_triangle, Point)

      ! This function computes the distance between an arbitrary point in the 
      ! 3D space and a triangle defined by the coordinates of its vertices.
      ! The minimum distance may be either the normal projection, the distance
      ! to one of the edges or the distance to one of the vertices

      implicit none
   
      real(kind = rdf), dimension(3,3), intent(in) :: vertices_triangle
      real(kind = rdf), dimension(3), intent(in) :: Point
      real(kind = rdf) :: DistancePointTriangle
      logical          :: InsideProjection
      real(kind = rdf) :: ProjectionDistance
      real(kind = rdf), dimension(3) :: ProjectedPoint

      ! local variables
   
      ! A,B,C are the triangle vertices. P is the point which we want to
      ! compute the distance to the triangle
   
      real(kind = rdf), dimension(3) :: A,B,C,P
   
      ! distances between point and vertices
      real(kind = rdf) :: DistPA, DistPB, DistPC
      ! distances between point and edges
      real(kind = rdf) :: DistPEdge1, DistPEdge2, DistPEdge3
   
      A = vertices_triangle(1,:)
      B = vertices_triangle(2,:)
      C = vertices_triangle(3,:)
      P = Point(:)
   
      ! Distances from the point P to the triangles vertices
      DistPA = norm2((P-A))
      DistPB = norm2((P-B))
      DistPC = norm2((P-C))
   
      ! Distances from the point P to the triangles Edges AB, BC, CA
      DistPEdge1 = DistancePointLine(A,B,P)
      DistPEdge2 = DistancePointLine(B,C,P)
      DistPEdge3 = DistancePointLine(C,A,P)
   
      ! Check if the points has a valid orthogonal projection on the triangle
      
      InsideProjection = .false.

      call CheckProjectionPointOnATriangle3( vertices_triangle , Point, &
                                             InsideProjection  , ProjectionDistance, &
                                             ProjectedPoint )
   
      if(InsideProjection) then

         DistancePointTriangle = minval( (/DistPA, DistPB, DistPC, &
                                           DistPEdge1, DistPEdge2, DistPEdge3, &
                                           ProjectionDistance/) )
      else

         DistancePointTriangle = minval( (/DistPA, DistPB, DistPC, &
                                           DistPEdge1, DistPEdge2, DistPEdge3/) )  
      end if
   
   end function DistancePointTriangle
   

   function DistancePointTriangle2(vertices_triangle, P) result(distance)
      
      ! function [dist,PP0] = pointTriangleDistance(TRI,P)
      ! calculate distance between a point and a triangle in 3D
      ! SYNTAX
      !   dist = pointTriangleDistance(TRI,P)
      !   [dist,PP0] = pointTriangleDistance(TRI,P)
      !
      ! DESCRIPTION
      !   Calculate the distance of a given point P from a triangle TRI.
      !   Point P is a row vector of the form 1x3. The triangle is a matrix
      !   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
      !   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
      !   to the triangle TRI.
      !   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
      !   closest point PP0 to P on the triangle TRI.
      !
      ! Author: Gwolyn Fischer
      ! Release: 1.0
      ! Release date: 09/02/02
      ! Release: 1.1 Fixed Bug because of normalization
      ! Release: 1.2 Fixed Bug because of typo in region 5 20101013
      ! Release: 1.3 Fixed Bug because of typo in region 2 20101014
   
      ! Possible extention could be a version tailored not to return the distance
      ! and additionally the closest point, but instead return only the closest
      ! point. Could lead to a small speed gain.
   
      ! Example:
      ! %% The Problem
      ! P0 = [0.5 -0.3 0.5]
      !
      ! P1 = [0 -1 0]
      ! P2 = [1  0 0]
      ! P3 = [0  0 0]
      !
      ! vertices = [P1; P2; P3]
      ! faces = [1 2 3]
      !
      ! %% The Engine
      ! [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0)
      !
      ! %% Visualization
      ! [x,y,z] = sphere(20)
      ! x = dist*x+P0(1)
      ! y = dist*y+P0(2)
      ! z = dist*z+P0(3)
      !
      ! figure
      ! hold all
      ! patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
      ! plot3(P0(1),P0(2),P0(3),'b*')
      ! plot3(PP0(1),PP0(2),PP0(3),'*g')
      ! surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
      ! view(3)
   
      ! The algorithm is based on
      ! "David Eberly, 'Distance Between Point and Triangle in 3D',
      ! Geometric Tools, LLC, (1999)"
      ! http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
      !
      !        ^t
      !  \     |
      !   \reg2|
      !    \   |
      !     \  |
      !      \ |
      !       \|
      !        *P2
      !        |\
      !        | \
      !  reg3  |  \ reg1
      !        |   \
      !        |reg0\
      !        |     \
      !        |      \ P1
      ! -------*-------*------->s
      !        |P0      \
      !  reg4  | reg5    \ reg6
   
   
      implicit none
      real(kind = rdf) :: vertices_triangle(3, 3), P(3), sqrdistance
      real(kind = rdf) :: Bvec(3), E0(3), E1(3), Dvec(3)
      real(kind = rdf) :: a, b, c, d, e, f, det, s, t, invDet
      real(kind = rdf) :: tmp0, tmp1, numer, denom
      real(kind = rdf) :: distance, PP0(3)
   
      ! rewrite triangle in normal form
      
      Bvec  = vertices_triangle(1,:)      !TRI[0, :]
      E0    = vertices_triangle(2,:) - Bvec  !TRI[1, :] - B
      ! E0 = E0/sqrt(sum(E0.^2)); %normalize vector
      E1    = vertices_triangle(3,:) - Bvec  !TRI[2, :] - B
      ! E1 = E1/sqrt(sum(E1.^2)); %normalize vector
      
      Dvec = Bvec - P
      
      a = dot_product( E0    , E0    )
      b = dot_product( E0    , E1    )
      c = dot_product( E1    , E1    )
      d = dot_product( E0    , Dvec  )
      e = dot_product( E1    , Dvec  )
      f = dot_product( Dvec  , Dvec  )
   
      !print "{0} {1} {2} ".format(B,E1,E0)
      det = a * c - b * b
      s   = b * e - c * d
      t   = b * d - a * e
   
      ! Terible tree of conditionals to determine in which region of the diagram
      ! shown above the projection of the point into the triangle-plane lies.
   
   
   
   
      if ( (s + t) <= det ) then
         if ( s < zero ) then         
            if ( t < zero ) then
               ! region4
               if ( d < zero ) then
                  t = zero         
                  if ( -d >= a ) then
                     s = one
                     sqrdistance = a + two * d + f
                  else
                     s = -d / a
                     sqrdistance = d * s + f
                  end if
               else
                  s = zero
                  if ( e >= zero ) then
                     t = zero
                     sqrdistance = f
                  else
                     if ( -e >= c ) then
                        t = one
                        sqrdistance = c + two * e + f
                     else
                        t = -e / c
                        sqrdistance = e * t + f
                        ! of region 4
                     end if
                  end if
               end if
            else
               ! region 3
               s = zero
               if ( e >= zero ) then
                  t = zero
                  sqrdistance = f
               else
                  if ( -e >= c ) then
                     t = one
                     sqrdistance = c + two * e + f
                  else
                     t = -e / c
                     sqrdistance = e * t + f
                     ! of region 3
                  end if
               end if
            end if
         else
            if ( t < zero ) then
               ! region 5
               t = zero
               if ( d >= zero ) then
                  s = zero
                  sqrdistance = f
               else
                  if ( -d >= a ) then
                     s = one
                     sqrdistance = a + two * d + f  ! GF 20101013 fixed typo d*s ->2*d
                  else
                     s = -d / a
                     sqrdistance = d * s + f
                  end if
               end if
            else
               ! region 0
               invDet = one / det
               s = s * invDet
               t = t * invDet
               sqrdistance = s * (a * s + b * t + two * d) + t * (b * s + c * t + two * e) + f
            end if
         end if
      else
         if ( s < zero ) then
            ! region 2
            tmp0 = b + d
            tmp1 = c + e
            if ( tmp1 > tmp0 ) then  ! minimum on edge s+t=1
               numer = tmp1 - tmp0
               denom = a - two * b + c
               if ( numer >= denom ) then
                  s = one
                  t = zero
                  sqrdistance = a + two * d + f  ! GF 20101014 fixed typo 2*b -> 2*d
               else
                  s = numer / denom
                  t = one - s
                  sqrdistance = s * (a * s + b * t + two * d) + t * (b * s + c * t + two * e) + f
               end if
            else  ! minimum on edge s=0
               s = zero
               if ( tmp1 <= zero ) then
                  t = one
                  sqrdistance = c + two * e + f
               else
                  if ( e >= zero ) then
                     t = zero
                     sqrdistance = f
                  else
                     t = -e / c
                     sqrdistance = e * t + f
                     ! of region 2
                  end if
               end if
            end if
         else 
            if ( t < zero ) then
               ! region6
               tmp0 = b + e
               tmp1 = a + d
               if ( tmp1 > tmp0 ) then
                  numer = tmp1 - tmp0
                  denom = a - two * b + c
                  if ( numer >= denom ) then
                     t = one
                     s = zero
                     sqrdistance = c + two * e + f
                  else
                     t = numer / denom
                     s = one - t
                     sqrdistance = s * (a * s + b * t + two * d) + t * (b * s + c * t + two * e) + f
                  end if
               else
                  t = zero
                  if ( tmp1 <= zero ) then
                     s = one
                     sqrdistance = a + two * d + f
                  else
                     if ( d >= zero ) then
                        s = zero
                        sqrdistance = f
                     else
                        s = -d / a
                        sqrdistance = d * s + f
                     end if
                  end if
               end if
            else
               ! region 1
               numer = c + e - b - d
               if ( numer <= zero ) then
                 s = zero
                 t = one
                 sqrdistance = c + two * e + f
               else
                  denom = a - two * b + c
                  if ( numer >= denom ) then
                     s = one
                     t = zero
                     sqrdistance = a + two * d + f
                  else
                     s = numer / denom
                     t = one - s
                     sqrdistance = s * (a * s + b * t + two * d) + t * (b * s + c * t + two * e) + f
                  end if
               end if
            end if
         end if
      end if
   
      ! account for numerical round-off error
      if ( sqrdistance < zero ) sqrdistance = zero
   
      distance = sqrt(sqrdistance)
      PP0 = Bvec + s * E0 + t * E1
   
   end function DistancePointTriangle2


   function DistancePointLine(A,B, Point)
   
      ! Distance between a line, defined by two points A, B and a point in 3D space.
      ! As the formula is defined for an infinite line in the space, I have to check
      ! the angles betwen AB and Point. The angle is checked using the dotproduct 
      ! sign from
   
      ! dotproduct(u,v) = |u|*|v|*cos(theta)
      !
      ! if cos(theta)<0 --> theta>90º
   
      implicit none
   
      real(kind = rdf), dimension(3), intent(in) :: A, B, Point
      real(kind = rdf) :: DistancePointLine
      logical :: PointLineValidProjection
   
      ! local variables
      real(kind = rdf), dimension(3) :: P, AB, AP, BA, BP
   
      P = Point
   
      AB = B-A
      AP = P-A
      BA = A-B
      BP = P-B
   
      PointLineValidProjection = .true.

      ! if points are coincident, then PointLineValidProjection is false
      if ( norm2(AB) < eps_sims ) PointLineValidProjection = .false.

      ! dotproduct(u,v) = |u|*|v|*cos(theta), so if cos(theta)<0 --> theta>90º
      if (Dot_Product(AB, AP) < zero) PointLineValidProjection = .false. 
      if (Dot_Product(BA, BP) < zero) PointLineValidProjection = .false. 
   
      if(PointLineValidProjection) then
         DistancePointLine = norm2(CrossProduct((P-A),(P-B)))/norm2(B-A)
      else
         ! distance to the vertices
         DistancePointLine = minval((/norm2(AP), norm2(BP)/))
      end if
   
   end function DistancePointLine
   
   subroutine CheckProjectionPointOnATriangle( vertices_triangle , &
                                               Point                , &
                                               InsideProjection     , &
                                               ProjectionDistance)
      
      
      ! To check when a point can be projected onto a triangle, we check if
      ! the angle between medians (lines from a vertex to to the midpoint of 
      ! the opposite edge) are <= 90º 
   
      ! The angle is checked using the dotproduct sign from
   
      ! dotproduct(u,v) = |u|*|v|*cos(theta)
      ! if cos(theta)<0 --> theta>90º
   
   
      !         C   
      !        / \         * P
      !       /   \
      !      x     x
      !     /       \
      !    /         \
      !   A ----x---- B  
   
   
      implicit none
   
      real(kind = rdf), dimension(3,3), intent(in) :: vertices_triangle
      real(kind = rdf), dimension(3), intent(in)   :: Point
      logical, intent(out) :: InsideProjection
      real(kind = rdf), intent(out) :: ProjectionDistance
   
      ! local variables
   
      ! A,B,C are the triangle vertices. P is the point which we want to
      ! compute the distance to the triangle
   
      real(kind = rdf), dimension(3) :: A,B,C,P, NormalVector
   
      ! Difference vectors 
      real(kind = rdf), dimension(3) :: AB, AC 
   
      ! Midpoint coordinates 
      real(kind = rdf), dimension(3) :: MidpointCB, MidpointAC, MidpointBA 
      ! Median vectors 
      real(kind = rdf), dimension(3) :: MedianA_CB, MedianB_AC, MedianC_BA 
      ! Midpoints to P vectors 
      real(kind = rdf), dimension(3) :: MidpointCB_P, MidpointAC_P, MidpointBA_P 
      ! dcoeff in the plane quation ax+by+cz+d = 0
      real(kind = rdf) :: dcoeff
   
      ! Triangle points
      A = vertices_triangle(1,:)
      B = vertices_triangle(2,:)
      C = vertices_triangle(3,:)
      P = Point(:)
   
      ! Edges midpoints
      MidpointCB = one_half * (C+B) ! opposite to A
      MidpointAC = one_half * (A+C) ! opposite to B
      MidpointBA = one_half * (B+A) ! opposite to C
   
      MedianA_CB = A - MidpointCB ! Median from A vertex
      MedianB_AC = B - MidpointAC ! Median from C vertex
      MedianC_BA = C - MidpointBA ! Median from B vertex
   
      MidpointCB_P = P-MidpointCB ! Vector from Median from A vertex to P
      MidpointAC_P = P-MidpointAC ! Vector from Median from B vertex to P
      MidpointBA_P = P-MidpointBA ! Vector from Median from A vertex to P
   
      InsideProjection = .False.
      ProjectionDistance = -one 
   
      ! if any of thes dot products are negative, the point projection
      ! lies outside the triangle
   
      if (Dot_Product(MedianA_CB, MidpointCB_P) <= zero) return
      if (Dot_Product(MedianB_AC, MidpointAC_P) <= zero) return
      if (Dot_Product(MedianC_BA, MidpointBA_P) <= zero) return
   
      InsideProjection = .True.
   
      AB = B-A
      AC = C-A
   
      NormalVector = CrossProduct(AB,AC)    ! arbitrary election for the normal to
                                            ! the plane define by the triangle
   
      dcoeff = -Dot_Product(NormalVector,A)  ! dcoeff in the plane equation
                                             ! ax+by+cz+d = 0
   
      ! Point - plane distance formula
   
      ProjectionDistance =    abs( Dot_Product(NormalVector,Point) + dcoeff) &
                                 / norm2(NormalVector)
   
   end subroutine CheckProjectionPointOnATriangle
   
   subroutine CheckProjectionPointOnATriangle2( vertices_triangle    , &
                                                Point                , &
                                                InsideProjection     , &
                                                ProjectionDistance   , &
                                                ProjectedPoint           )
      
      
      ! we check if a point 

      ! The angle is checked using the dotproduct sign from
   
      ! dotproduct(u,v) = |u|*|v|*cos(theta)
      ! if cos(theta)<0 --> theta>90º
   
   
      !         C   
      !        / \         * P
      !       /   \
      !      x     x
      !     /       \
      !    /         \
      !   A ----x---- B  
   
   
      implicit none
   
      real(kind = rdf), dimension(3,3), intent(in) :: vertices_triangle
      real(kind = rdf), dimension(3), intent(in)   :: Point
      logical, intent(out) :: InsideProjection
      real(kind = rdf), intent(out) :: ProjectionDistance
      real(kind = rdf), dimension(3), intent(out) :: ProjectedPoint
   
      ! local variables
   
      ! A,B,C are the triangle vertices. P is the point which we want to
      ! compute the distance to the triangle
   
      real(kind = rdf), dimension(3) :: A,B,C,P, NormalVector
   
      real(kind = rdf), dimension(3) :: CentroidABC, PCentroindVector!, ProjectedPoint
      real(kind = rdf) :: DistancePointPlane, areaABC, alpha, beta, gamma

      ! Triangle points
      A = vertices_triangle(1,:)
      B = vertices_triangle(2,:)
      C = vertices_triangle(3,:)
      
      ! Centroid triangle

      CentroidABC = (A+B+C)/three

      ! Point to project
      P = Point(:)
   
      ! First, we project the point onto the infinite plane that 
      ! contains the triangle. We choose as the origin point, the centroid
      ! of the triangle and the normal as the cross product between edges

      ! unitary normal vector
      NormalVector = CrossProduct(B-A,C-A) / norm2 (CrossProduct(B-A,C-A)) 

      ! vector from P to the centroid of the triangle
      PCentroindVector = P - CentroidABC

      DistancePointPlane = Dot_Product(PCentroindVector,NormalVector)

      ProjectedPoint = P - DistancePointPlane * NormalVector


      ! Now we need to check if the projected point lies into the 
      ! triangle or not

      areaABC = norm2( CrossProduct( (B-A),(C-A) ) ) / two
   
      alpha = norm2( CrossProduct( (ProjectedPoint-B),(ProjectedPoint-C) ) ) / (two*areaABC)
      beta  = norm2( CrossProduct( (ProjectedPoint-C),(ProjectedPoint-A) ) ) / (two*areaABC)
      gamma = one - alpha - beta

      if ( alpha > zero .and. alpha < one .and. & 
           beta  > zero .and. beta  < one .and. & 
           gamma > zero .and. gamma < one          ) then

         InsideProjection = .true.
         ProjectionDistance = DistancePointPlane

      else

         InsideProjection = .false.

      end if

   end subroutine CheckProjectionPointOnATriangle2
   

   subroutine CheckProjectionPointOnATriangle3( vertices_triangle    , &
                                                Point                , &
                                                InsideProjection     , &
                                                ProjectionDistance   , &
                                                ProjectedPoint           )
      
      
      ! we check if a point 

      ! The angle is checked using the dotproduct sign from
   
      ! dotproduct(u,v) = |u|*|v|*cos(theta)
      ! if cos(theta)<0 --> theta>90º
   
   
      !         C   
      !        / \         * P
      !       /   \
      !      x     x
      !     /       \
      !    /         \
      !   A ----x---- B  
   
   
      implicit none
   
      real(kind = rdf), dimension(3,3), intent(in) :: vertices_triangle
      real(kind = rdf), dimension(3), intent(in)   :: Point
      logical, intent(out) :: InsideProjection
      real(kind = rdf), intent(out) :: ProjectionDistance
      real(kind = rdf), dimension(3), intent(out) :: ProjectedPoint
   
      ! local variables
   
      ! A,B,C are the triangle vertices. P is the point which we want to
      ! compute the distance to the triangle
   
      real(kind = rdf), dimension(3) :: P1, P2, P3, P, u, v, w, NormalVector
   
      real(kind = rdf) :: DistancePointPlane, alpha, beta, gamma

      ! Triangle points
      P1 = vertices_triangle(1,:)
      P2 = vertices_triangle(2,:)
      P3 = vertices_triangle(3,:)
      
      ! if vertices_triangle doesn't form a proper triangle, the subroutine just returns

      if (       norm2(P1-P2) <= eps_sims .or. norm2(P1-P3) <= eps_sims &
            .or. norm2(P2-P3) <= eps_sims                                   ) then

         InsideProjection = .false.
         return

      end if

      ! Point to project
      P = Point(:)

      u = P2 - P1
      v = P3 - P1
      w = P - P1

      NormalVector = CrossProduct (u, v)

      alpha =    Dot_Product(CrossProduct(u, w) , NormalVector) &
               / Dot_Product(NormalVector,NormalVector)

      beta  =    Dot_Product(CrossProduct(w, v) , NormalVector) &
               / Dot_Product(NormalVector,NormalVector)

      gamma = one - beta - alpha
   
      ProjectedPoint = alpha * P1 + beta * P2 + gamma * P3

      DistancePointPlane = norm2(P - ProjectedPoint)

      if ( alpha > zero .and. alpha < one .and. & 
           beta  > zero .and. beta  < one .and. & 
           gamma > zero .and. gamma < one          ) then

         InsideProjection   = .true.
         ProjectionDistance = DistancePointPlane

      else

         InsideProjection   = .false.
         ProjectionDistance = 999.9_rdf

      end if

   end subroutine CheckProjectionPointOnATriangle3

   function CrossProduct(v1, v2)
   
     implicit none
     
     real(kind = rdf), dimension(3):: CrossProduct, v1, v2
     
     CrossProduct(1) = v1(2)*v2(3) - v1(3)*v2(2)
     CrossProduct(2) = v1(3)*v2(1) - v1(1)*v2(3)
     CrossProduct(3) = v1(1)*v2(2) - v1(2)*v2(1)
     
   end function CrossProduct
   
   function M44Det (A) result (Det)
   
      ! Determines the determinant of a 4x4 matrix
      ! Obtained from: 
      ! http://web.hku.hk/~gdli/UsefulFiles/matrix/m44det_f90.txt
   
      implicit none
   
      real(kind = rdf), dimension(4,4), intent(in)  :: A
   
      real(kind = rdf) :: Det
   
      ! (ugly, but it works)
      Det =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
             A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
             A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
             A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
             A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   
   
      return
   
   end function M44Det



   function HexahedronVolume( VerticesCoordinates ) result( volume )

      real (kind = rdf ), dimension(8,3) :: VerticesCoordinates
      real(kind = rdf), dimension(3) :: vertA, vertB, vertC, vertD, vertE, vertF, vertG, vertH
      real (kind = rdf ), dimension(3) :: AC, AF, AH, BE, CH, DB, ED, FC, GA, GB, GC, GD, GE, GF, GH, HF
      real (kind = rdf ) :: volume

         !   NODE DISTRIBUTION OF THE COMPUTATIONAL CELL
         !
         !    i,j+1,k+1 ----i+1,j+1,k+1    
         !     /|(H)           /|(G)                            
         !    / |             / |                                  
         ! i,j,k+1-------i+1,j,k+1                   
         !   |(E)           |(F)|                                   
         !   |  |           |   |                                
         !   |  |           |   |                                
         !   |  i,j+1,k-----|-i+1,j+1,k           
         !   | /(D)         |  /(C)                 
         !   |/             | /
         ! i,j,k-----------i+1,j,k
         !  (A)               (B)

         ! 8 vertices coordinates
         vertA(:) = VerticesCoordinates(1,:)
         vertB(:) = VerticesCoordinates(2,:)
         vertC(:) = VerticesCoordinates(3,:)
         vertD(:) = VerticesCoordinates(4,:)
         vertE(:) = VerticesCoordinates(5,:)
         vertF(:) = VerticesCoordinates(6,:)
         vertG(:) = VerticesCoordinates(7,:)
         vertH(:) = VerticesCoordinates(8,:)
   
         ! edges vectors
         AC = vertC - vertA
         AF = vertF - vertA
         AH = vertH - vertA
   
         BE = vertE - vertB
   
         CH = vertH - vertC
   
         DB = vertB - vertD
   
         ED = vertD - vertE
   
         FC = vertC - vertF
   
         GA = vertA - vertG
         GB = vertB - vertG
         GC = vertC - vertG
         GD = vertD - vertG
         GE = vertE - vertG
         GF = vertF - vertG
         GH = vertH - vertG
   
         HF = vertF - vertH
   
         ! volume formula for a general hexahedronn obtained from
         ! Davies, D. E., & Salmond, D. J. (1985). Calculation of the volume of a general  
         ! hexahedron for flow predictions. AIAA journal, 23(6), 954-956.
   
         volume =   -one / twelve * (    Dot_Product( (GA+GB) , ( CrossProduct(DB,AC) ) ) &
                                      +  Dot_Product( (GA+GE) , ( CrossProduct(BE,AF) ) ) &
                                      +  Dot_Product( (GA+GD) , ( CrossProduct(ED,AH) ) ) &
                                      +  Dot_Product( GF      , ( CrossProduct(HF,GE) ) ) &
                                      +  Dot_Product( GH      , ( CrossProduct(CH,GD) ) ) &
                                      +  Dot_Product( GC      , ( CrossProduct(FC,GB) ) ) &
                                    )
 
   end function HexahedronVolume





   function HexahedronVolume2(x,y,z) result( volume )

      real (kind = rdf ) :: x(8), y(8), z(8)
      real (kind = rdf ) :: volume, a, b, c, d, e, f
   
      a = (x(2)-x(1))*(y(3)-y(1))*(z(4)-z(1))
      b = (x(3)-x(1))*(y(4)-y(1))*(z(5)-z(1))
      c = (x(4)-x(1))*(y(5)-y(1))*(z(6)-z(1))
      d = (x(5)-x(1))*(y(6)-y(1))*(z(7)-z(1))
      e = (x(6)-x(1))*(y(7)-y(1))*(z(8)-z(1))
      f = (x(7)-x(1))*(y(8)-y(1))*(z(2)-z(1))
   
      volume = (a+b+c+d+e+f)/6.0d0 - ( (x(2)-x(1))*(y(4)-y(1))*(z(5)-z(1)) + &
                                       (x(3)-x(1))*(y(5)-y(1))*(z(6)-z(1)) + &
                                       (x(4)-x(1))*(y(6)-y(1))*(z(7)-z(1)) + &
                                       (x(5)-x(1))*(y(7)-y(1))*(z(8)-z(1)) + &
                                       (x(6)-x(1))*(y(8)-y(1))*(z(2)-z(1)) + &
                                       (x(7)-x(1))*(y(2)-y(1))*(z(3)-z(1)) )/6.0d0
   
   end function HexahedronVolume2




   subroutine CellAnalysis ( PhiVertices                    , &
                             VerticesCoordinates            , &
                             TriangulationCase              , &
                             PhaseChangeTetrahedronsVolume  , &
                             SecondTetrahedronsVolume         &
                           )

      implicit none


      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !                                      VARIABLES DECLARATION
      ! 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !
      ! This subroutine perfoms the phase change analysis on a 8 vertices cell.
      ! It receives the coordinates of the cell vertices (optional), and the phi
      ! values at those vertices.
      !
      !
      !    i,j+1,k+1 ----i+1,j+1,k+1    
      !     /|(7)           /|(6)                            
      !    / |             / |                                  
      ! i,j,k+1-------i+1,j,k+1                   
      !   |(4)           |(5)|                                   
      !   |  |           |   |                                
      !   |  |           |   |                                
      !   |  i,j+1,k-----|-i+1,j+1,k           
      !   | /(3)         |  /(2)                 
      !   |/             | /
      ! i,j,k-----------i+1,j,k
      !  (0)               (1)
      !
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! INPUT VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
      real(kind = rdf), dimension(0:7), intent(in) :: PhiVertices
      
      ! OPTIONAL INPUTS: 

      integer, optional ,intent(in) :: TriangulationCase
      real(kind = rdf), dimension(0:7,1:3), optional, intent(in) :: VerticesCoordinates
   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OUTPUT VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      character(len = 15) :: CellType ! "Air Cell", "Water Cell" or "Phase-Change Cell"

      integer, dimension (0:7) :: PnodeVertexFlag
      integer, dimension (1:6) :: PhaseChangeTetrahedronList

      ! OPTIONAL OUTPUTS

      real(kind = rdf), optional, intent(out) :: PhaseChangeTetrahedronsVolume
      real(kind = rdf), optional, intent(out) :: SecondTetrahedronsVolume

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! LOCAL VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      integer, dimension (1:6,1:4) :: TetrahedronCellVertices
      logical :: PhaseChangeCell, SecondTetrahedron, PhaseChangeTetrahedron
      real(kind = rdf) :: PhiAux, SignPhiAux
      real(kind = rdf), dimension(4,3) :: VerticesCoordinatesAux
      real(kind = rdf), dimension(4)   :: PhiVerticesAux
      real(kind = rdf) :: WaterVolumeAux

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! MARCHING TETRAHEDRON VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      real(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
      integer                             :: ntriangles
      real(kind = rdf)                    :: IsosurfaceArea
      real(kind = rdf), dimension(4)      :: VerticesDistances
      real(kind = rdf), dimension(4,3)    :: vertices_coordinates_aux

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! COUNTERS
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      integer :: CellVertexLoop , TetrahedronLoop , TetrahedronVertexLoop

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !                                           SUBROUTINE BODY
      ! 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Single phase cells:

      if ( all( PhiVertices  > eps_sims ) ) then

         CellType = 'Water Cell'

         return

      else if ( all( PhiVertices  < -eps_sims ) ) then

         CellType = 'Air Cell'

         return

      else ! Candidate for a phase-changing cell

         if ( present( TriangulationCase ) )  then
      
            select case ( TriangulationCase )
      
               case(1) ! Main diagonal from v0 --> v6
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6
      
               case(2) ! Main diagonal from v1 --> v7
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/1, 6, 2, 7/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/1, 2, 3, 7/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/1, 3, 0, 7/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/1, 0, 4, 7/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/1, 4, 5, 7/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/1, 5, 6, 7/) ! cell vertices of the tetrahedron # 6
            
               case(3) ! Main diagonal from v3 --> v5
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/3, 4, 0, 5/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/3, 0, 1, 5/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/3, 1, 2, 5/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/3, 2, 6, 5/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/3, 6, 7, 5/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/3, 7, 4, 5/) ! cell vertices of the tetrahedron # 6

               case(4)! Main diagonal from v2 --> v4
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/2, 7, 3, 4/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/2, 3, 0, 4/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/2, 0, 1, 4/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/2, 1, 5, 4/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/2, 5, 6, 4/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/2, 6, 7, 4/) ! cell vertices of the tetrahedron # 6
      
            end select
      
         else
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6
      
         end if
   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !     PHASE CHANGE CELL CHECKING
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         PhaseChangeCell = .false.
   
         ! If any vertex has a zero phi value, then, the cell is a phase-changing cell
         if ( any( abs( PhiVertices ) < eps_sims ) ) then

            PhaseChangeCell = .true.
   
         ! Otherwise, I need to explore sign change over the phi values at the cell vertices
         else

            SignPhiAux = sign(one, PhiVertices(0))
            
            ! Loop over all cell vertices      
            SearchPhaseChange: do CellVertexLoop = 1,7
      
               if ( SignPhiAux * PhiVertices( CellVertexLoop ) < zero ) then
                  
                  PhaseChangeCell = .true.
                  exit SearchPhaseChange
               
               end if
      
            end do SearchPhaseChange
         
         end if

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !     PHASE CHANGE TETRAHEDRONS CHECKING
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
         if ( present ( PhaseChangeTetrahedronsVolume ) ) PhaseChangeTetrahedronsVolume = zero
         if ( present ( SecondTetrahedronsVolume      ) ) SecondTetrahedronsVolume = zero

         if ( PhaseChangeCell ) then
   
            PhaseChangeTetrahedronList = (/0,0,0,0,0,0/) ! The six tetrahedrons are set as not-phase change ones
   
            do TetrahedronLoop = 1,6
   
               TetrahedronVertexLoop  = 1   
               PhaseChangeTetrahedron = .false.
               SecondTetrahedron      = .true.     

               PhiAux = PhiVertices( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )
   
               if ( PhiAux <= eps_sims ) SecondTetrahedron = .false.

               if ( abs( PhiAux ) < eps_sims ) then
   
                  PhaseChangeTetrahedron                        = .true.
                  PhaseChangeTetrahedronList( TetrahedronLoop ) = 1
   
               end if
   
               ! SignPhiAux = -1 or 1
               SignPhiAux = sign( one, PhiAux )
   
               ! Loop over the rest of the nodes of the analysed tetrahedron
               do TetrahedronVertexLoop = 2,4
                  
                  PhiAux = PhiVertices( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )
      
                  if ( SignPhiAux * PhiAux < zero .or. abs( PhiAux ) < eps_sims ) then
                     
                     if ( PhiAux <= eps_sims ) SecondTetrahedron = .false.

                     PhaseChangeTetrahedron = .true.
                     PhaseChangeTetrahedronList( TetrahedronLoop ) = 1
   
                  end if
   
               end do ! TetrahedronVertexLoop = 2,4
   
               if ( PhaseChangeTetrahedron ) then
                                    
                  ! Marching Tetrahedra

                  if ( present ( VerticesCoordinates ) .and. present ( PhaseChangeTetrahedronsVolume ) ) then

                     do TetrahedronVertexLoop = 1,4
                        
                        VerticesCoordinatesAux( TetrahedronVertexLoop , : ) = &
                        VerticesCoordinates( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ),:)
                        
                        PhiVerticesAux( TetrahedronVertexLoop ) = &
                        PhiVertices( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )

                     end do

                     call MarchingTetrahedron(  VerticesCoordinatesAux     , &
                                                PhiVerticesAux             , &
                                                VerticesIsosurfaces        , &
                                                ntriangles                 , &
                                                IsosurfaceArea             , &
                                                WaterVolumeAux             , &
                                                VerticesDistances            &
                                             )

                     PhaseChangeTetrahedronsVolume = PhaseChangeTetrahedronsVolume + WaterVolumeAux

                  end if

               end if

               ! Second Neighbour Tetrahedron PhaseChangeTetrahedronList = 2
   
               if ( SecondTetrahedron .and. .not. ( PhaseChangeTetrahedron ) ) then
   
                  PhaseChangeTetrahedronList( TetrahedronLoop ) = 2
                  
                  if ( present ( VerticesCoordinates ) .and. present ( SecondTetrahedronsVolume ) ) then

                     do TetrahedronVertexLoop = 1,4
                        VerticesCoordinatesAux(TetrahedronVertexLoop,:) = &
                        VerticesCoordinates( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ),:)
                     end do

                     SecondTetrahedronsVolume = SecondTetrahedronsVolume + GetTetrahedronVolume( VerticesCoordinatesAux )

                  end if

               end if
   
            end do ! TetrahedronLoop = 1,6
   
      
         end if ! PhaseChangeCell == .true.

      end if ! Candidate for a phase-change cell


   end subroutine CellAnalysis


   subroutine CellPnodeTetrahedronIdentification ( PhiVertices                    , &
                                                   TriangulationCase              , &
                                                   nCellVerticesPnodes            , &
                                                   nPhaseChangeTetrahedrons         &
                                                 )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !                                      VARIABLES DECLARATION
      ! 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !
      ! This subroutine perfoms the phase change analysis on a 8 vertices cell.
      ! It receives the coordinates of the cell vertices (optional), and the phi
      ! values at those vertices.
      !
      !
      !    i,j+1,k+1 ----i+1,j+1,k+1    
      !     /|(7)           /|(6)                            
      !    / |             / |                                  
      ! i,j,k+1-------i+1,j,k+1                   
      !   |(4)           |(5)|                                   
      !   |  |           |   |                                
      !   |  |           |   |                                
      !   |  i,j+1,k-----|-i+1,j+1,k           
      !   | /(3)         |  /(2)                 
      !   |/             | /
      ! i,j,k-----------i+1,j,k
      !  (0)               (1)
      !
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! INPUT VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
      real(kind = rdf), dimension(0:7), intent(in) :: PhiVertices
      
      ! OPTIONAL INPUTS: 

      integer, optional ,intent(in) :: TriangulationCase
   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OUTPUT VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      integer , intent(out) :: nCellVerticesPnodes
      integer , intent(out) :: nPhaseChangeTetrahedrons

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! LOCAL VARIABLES
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      integer, dimension (1:6,1:4) :: TetrahedronCellVertices
      logical :: PhaseChangeCell, SecondTetrahedron, PhaseChangeTetrahedron
      real(kind = rdf) :: PhiAux, SignPhiAux
      real(kind = rdf), dimension(4,3) :: VerticesCoordinatesAux
      real(kind = rdf), dimension(4)   :: PhiVerticesAux
      real(kind = rdf) :: WaterVolumeAux

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! COUNTERS
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      integer :: CellVertexLoop , TetrahedronLoop , TetrahedronVertexLoop

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !
      !                                           SUBROUTINE BODY
      ! 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Single phase cells:

      nCellVerticesPnodes      = 0
      nPhaseChangeTetrahedrons = 0

      ! Water - cell
      if ( all( PhiVertices  > eps_sims ) ) then

         return

      ! Air - cell
      else if ( all( PhiVertices  < -eps_sims ) ) then

         return

      ! Candidate for a phase-changing cell
      else 

         if ( present( TriangulationCase ) )  then
      
            select case ( TriangulationCase )
      
               case(1) ! Main diagonal from v0 --> v6
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6
      
               case(2) ! Main diagonal from v1 --> v7
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/1, 6, 2, 7/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/1, 2, 3, 7/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/1, 3, 0, 7/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/1, 0, 4, 7/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/1, 4, 5, 7/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/1, 5, 6, 7/) ! cell vertices of the tetrahedron # 6
      
               case(3)! Main diagonal from v2 --> v4
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/2, 7, 3, 4/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/2, 3, 0, 4/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/2, 0, 1, 4/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/2, 1, 5, 4/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/2, 5, 6, 4/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/2, 6, 7, 4/) ! cell vertices of the tetrahedron # 6
      
               case(4) ! Main diagonal from v3 --> v5
      
                  TetrahedronCellVertices( 1 , 1:4 ) = (/3, 4, 0, 5/) ! cell vertices of the tetrahedron # 1
                  TetrahedronCellVertices( 2 , 1:4 ) = (/3, 0, 1, 5/) ! cell vertices of the tetrahedron # 2
                  TetrahedronCellVertices( 3 , 1:4 ) = (/3, 1, 2, 5/) ! cell vertices of the tetrahedron # 3
                  TetrahedronCellVertices( 4 , 1:4 ) = (/3, 2, 6, 5/) ! cell vertices of the tetrahedron # 4
                  TetrahedronCellVertices( 5 , 1:4 ) = (/3, 6, 7, 5/) ! cell vertices of the tetrahedron # 5
                  TetrahedronCellVertices( 6 , 1:4 ) = (/3, 7, 4, 5/) ! cell vertices of the tetrahedron # 6
      
            end select
      
         else
      
            TetrahedronCellVertices( 1 , 1:4 ) = (/0, 5, 1, 6/) ! cell vertices of the tetrahedron # 1
            TetrahedronCellVertices( 2 , 1:4 ) = (/0, 1, 2, 6/) ! cell vertices of the tetrahedron # 2
            TetrahedronCellVertices( 3 , 1:4 ) = (/0, 2, 3, 6/) ! cell vertices of the tetrahedron # 3
            TetrahedronCellVertices( 4 , 1:4 ) = (/0, 3, 7, 6/) ! cell vertices of the tetrahedron # 4
            TetrahedronCellVertices( 5 , 1:4 ) = (/0, 7, 4, 6/) ! cell vertices of the tetrahedron # 5
            TetrahedronCellVertices( 6 , 1:4 ) = (/0, 4, 5, 6/) ! cell vertices of the tetrahedron # 6
      
         end if
   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !     PHASE CHANGE CELL CHECKING
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         PhaseChangeCell = .false.
   
         ! If any vertex has a zero phi value, then, the cell is a phase-changing cell
         if ( any( abs( PhiVertices ) < eps_sims ) ) then

            PhaseChangeCell = .true.
   
         ! Otherwise, I need to explore sign change over the phi values at the cell vertices
         else

            SignPhiAux = sign(one, PhiVertices(0))
            
            ! Loop over all cell vertices      
            SearchPhaseChange: do CellVertexLoop = 1,7
      
               if ( SignPhiAux * PhiVertices( CellVertexLoop ) < zero ) then
                  
                  PhaseChangeCell = .true.
                  exit SearchPhaseChange
               
               end if
      
            end do SearchPhaseChange
         
         end if

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !     PHASE CHANGE TETRAHEDRONS CHECKING
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
         if ( PhaseChangeCell ) then
      
            do TetrahedronLoop = 1,6
   
               TetrahedronVertexLoop  = 1   
               PhaseChangeTetrahedron = .false.
               SecondTetrahedron      = .true.     

               PhiAux = PhiVertices( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )
   
               if (       PhiAux <= eps_sims ) SecondTetrahedron     = .false.
               if ( abs( PhiAux ) < eps_sims ) PhaseChangeTetrahedron = .true.
   
               ! SignPhiAux = -1 or 1
               SignPhiAux = sign( one, PhiAux )
   
               ! Loop over the rest of the nodes of the analysed tetrahedron
               do TetrahedronVertexLoop = 2,4
                  
                  PhiAux = PhiVertices( TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )
      
                  if ( SignPhiAux * PhiAux < zero .or. abs( PhiAux ) < eps_sims ) then
                     
                     if ( PhiAux <= eps_sims ) SecondTetrahedron = .false.

                     PhaseChangeTetrahedron = .true.
   
                  end if
   
               end do ! TetrahedronVertexLoop = 2,4
   
               if ( PhaseChangeTetrahedron .and. .not. (SecondTetrahedron) ) then
                  
                  ! Now I tag the cell vertices as Pnodes and update the number of  
                  ! tetrahedrons within the cell using its binary representation to check
                  ! afterwards which tetrahedrons are KTetrahedrons                                    

                  ! it is gonna be setting from the second bit on to keep the numeration of the tetrahedron from 1 to 6
                  ! it means:
                  !
                  !  ...000000 --> ...000010 --> ...000110 etc

                  nPhaseChangeTetrahedrons = ibset( nPhaseChangeTetrahedrons , TetrahedronLoop ) 

                  do TetrahedronVertexLoop = 1,4

                     nCellVerticesPnodes =   ibset( nCellVerticesPnodes , &
                                                    TetrahedronCellVertices ( TetrahedronLoop , TetrahedronVertexLoop ) )

                  end do

               end if
   
            end do ! TetrahedronLoop = 1,6
         
         end if ! PhaseChangeCell == .true.

      end if ! Candidate for a phase-change cell


   end subroutine CellPnodeTetrahedronIdentification   

   function PhaseChangeLimiterCorrection ( value, correction, option ) result (correction_value)

      real ( kind = rdf ) :: value, correction, correction_value
      integer :: option

      select case (option)

         case (0) ! no restriction --> alway returns correction
         
            correction_value = correction

         case (1) ! if phase doesn't change  --> returns correction
                  ! if phase         changes --> returns 0

            correction_value = correction * abs( sign( one , value ) + sign( one , value + correction) ) / two
         
         case (2) ! if phase doesn't change  --> returns correction
                  ! if phase         changes --> returns -value (to make phi = zero)
         
            correction_value = correction * abs( sign( one , value ) + sign( one , value + correction)       ) / two + &
                               value      *    ( sign( one , value ) + sign( one , value + correction) - two ) / two

      end select   

   end function PhaseChangeLimiterCorrection


   function LocalVolumeVariation( etaK , VerticesCoordinates , PhiHast , VolumeIni ) result( DeltaVolume )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! INPUT VARIABLES    |
      ! - - - - - - - - - -
      !

      real( kind = rdf ),                  intent(in) :: etaK, VolumeIni ! Correction and initial volume    
      real( kind = rdf ), dimension(4,3) , intent(in) :: VerticesCoordinates ! Tetrahedron vertices coordinates
      real( kind = rdf ), dimension(4)   , intent(in) :: PhiHast ! ϕh*, Signed distance initial reconstruction

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OUTPUT VARIABLES   |
      ! - - - - - - - - - -
      !

      real( kind = rdf ) :: DeltaVolume

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! LOCAL VARIABLES (MARCHING TETRAHEDRON)   |
      ! - - - - - - - - - - - - - - - - - - - - -
      !

      real( kind = rdf ), dimension(4)      :: PhiAux ! ϕh* + ηK
      real( kind = rdf ), dimension (2,3,3) :: VerticesIsosurfaces
      integer                               :: ntriangles
      real( kind = rdf )                    :: IsosurfaceArea
      real( kind = rdf )                    :: WaterVolumeAux
      real( kind = rdf ), dimension(4)      :: VerticesDistances
      real( kind = rdf ), dimension(4,3)    :: vertices_coordinates_aux

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! ϕ values within the tetrahedron with the correction applied
      PhiAux = PhiHast + etaK

      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 PhiAux                 , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 IsosurfaceArea         , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )


      DeltaVolume = VolumeIni - WaterVolumeAux

   end function LocalVolumeVariation

   subroutine SecantMethodLocalCorrection ( VerticesCoordinates         , &
                                            VolumeIniPhiH               , &
                                            AreaIniPhiH                 , &
                                            PhiHast                     , &
                                            Tolerance                   , &
                                            MaxNumberOfIterations       , &
                                            x2                          , &
                                            IsConverged                 , &
                                            IterationsUntilConvergence    &
                                          )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! INPUT VARIABLES    |
      ! - - - - - - - - - -
      !
      ! ϕh            :    Linear reconstruction of the advected phi distribution      
      ! ϕh* (PhiHast) : Signed distance reconstruction computed as the distance to the free surface reconstructed
      !                 from ϕh                   
      ! VolumeIniPhiH : Water Volume within the tetrahedron from the linear free surface reconstruction obtained
      !                 from ϕh
      !

      real( kind = rdf ), dimension(4,3) , intent(in) :: VerticesCoordinates ! Tetrahedron vertices coordinates
      real( kind = rdf ),                  intent(in) :: VolumeIniPhiH , AreaIniPhiH ! Initial vol and area with (ϕh)
      real( kind = rdf ), dimension(4)   , intent(in) :: PhiHast ! ϕh*, Signed distance function initial reconstruction
      
      real( kind  = rdf ), intent(in) :: Tolerance
      integer,             intent(in) :: MaxNumberOfIterations 
   
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OUTPUT VARIABLES   |
      ! - - - - - - - - - -  
   
      real( kind  = rdf ), intent(out) :: x2
   
      logical, optional, intent(out) :: IsConverged
      integer, optional, intent(out) :: IterationsUntilConvergence

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! LOCAL VARIABLES  |
      ! - - - - - - - - -   

      real( kind = rdf ) :: x0 , x1 , fx0 , fx1 , DeltaVolumeIni
      integer :: iConvergence

      real( kind = rdf ), dimension (2,3,3) :: VerticesIsosurfaces
      integer                               :: ntriangles
      real( kind = rdf )                    :: AreaIsosurfacePhiHast
      real( kind = rdf )                    :: WaterVolumeAux
      real( kind = rdf ), dimension(4)      :: VerticesDistances
      real( kind = rdf ), dimension(4,3)    :: vertices_coordinates_aux

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if ( present ( IsConverged ) ) IsConverged = .false.

      ! Initial guests determination

      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 PhiHast                , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 AreaIsosurfacePhiHast  , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )

      ! ΔV(ϕh, ϕ*) = ΔV(ϕh, ϕSD)   
      DeltaVolumeIni = VolumeIniPhiH - WaterVolumeAux

      ! If the initial ΔV(ϕh, ϕSD) is less than the Tolerance, then I don't need to correct
      ! ϕ in the tetrahedron and I just need to set etaK to 0

      if ( abs( DeltaVolumeIni ) < Tolerance ) then

         x2 = zero ! etaK = 0 ( there's nothing to correct )
         if ( present ( IsConverged                ) ) IsConverged = .true.
         if ( present ( IterationsUntilConvergence ) ) IterationsUntilConvergence = 0 
         
         return

      end if

      ! η0 = -ΔV(ϕh, ϕSD) / S(ϕSD)
      if ( AreaIsosurfacePhiHast > eps_sims ) then
         x0 = DeltaVolumeIni / AreaIsosurfacePhiHast
      else
         x0 = zero
      end if

      ! η1 = -ΔV(ϕh, ϕSD) / S(ϕh)
      if ( AreaIniPhiH > eps_sims ) then
         x1 = DeltaVolumeIni / AreaIniPhiH
      else
         x1 = zero
      end if

      fx0 = LocalVolumeVariation( x0, VerticesCoordinates , PhiHast , VolumeIniPhiH )
      fx1 = LocalVolumeVariation( x1, VerticesCoordinates , PhiHast , VolumeIniPhiH )


      ! If the convergence condition is hold before the loop, I just keep the etaK value 
      ! that minimises the volume difference
      if( abs(fx0) < Tolerance .or. abs(fx1) < Tolerance) then
         
         if ( present ( IsConverged                ) ) IsConverged = .true.
         if ( present ( IterationsUntilConvergence ) ) IterationsUntilConvergence = 0 

         x2 = x0
         if ( abs(fx1) < abs(fx0) ) x2 = x1

         return
         
      end if

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ****************************       CONVERGENCE LOOP        ****************************
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      RootFinderLoop : do iConvergence = 1, MaxNumberOfIterations
      
         if ( abs( fx1 - fx0 ) < eps_sims ) exit RootFinderLoop

         x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)

         x0  = x1
         x1  = x2

         fx0 = fx1 
         fx1 = LocalVolumeVariation( x1, VerticesCoordinates , PhiHast , VolumeIniPhiH )
         
         ! Apply a stopping criterion here (see below)
      
         if ( abs(fx1) < Tolerance ) then

            if ( present ( IsConverged ) ) IsConverged = .true.

            exit RootFinderLoop

         end if

      end do RootFinderLoop

      ! If the resultant volume difference is
      if ( abs(fx1) > DeltaVolumeIni ) x2 = zero

      ! Just to see how many iterations it took to converge 
      if ( present ( IterationsUntilConvergence ) ) IterationsUntilConvergence = iConvergence


   end subroutine SecantMethodLocalCorrection






end module TetrahedronMethods



