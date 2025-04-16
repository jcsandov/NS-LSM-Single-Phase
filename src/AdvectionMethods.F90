module AdvectionMethods
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Set of routines to compute advection terms of conservation laws
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
   ! Jorge Sandoval, UoE/PUC. Edinburgh, February 20th, 2023.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

   use precision
   
   implicit none

   contains


   function GetUPWINDReconstruction(f0, f1, f2, f3, WaveSpeed)
   
      implicit none

      real (kind = rdf) :: f0, f1, f2, f3, WaveSpeed

      ! local variables

      real (kind = rdf) :: d0, d1
      real (kind = rdf) :: epsWENO3
      real (kind = rdf) :: fL, fC, fR
      real (kind = rdf) :: alpha0, alpha1
      real (kind = rdf) :: beta0 , beta1
      real (kind = rdf) :: u0, u1, w0, w1
      real (kind = rdf) :: GetUPWINDReconstruction

      d0       = two / three
      d1       = one / three
      epsWENO3 = 1.0E-06_rdf

      if ( WaveSpeed > zero ) then
         fL = f0
         fC = f1
         fR = f2
      else ! mirror
         fL = f3
         fC = f2
         fR = f1
      end if   

      beta0  = ( fC - fR )**two
      beta1  = ( fL - fC )**two

      alpha0 = d0 / ( beta0 + epsWENO3 )**two
      alpha1 = d1 / ( beta1 + epsWENO3 )**two

      w0 = alpha0 / ( alpha0 + alpha1 )
      w1 = alpha1 / ( alpha0 + alpha1 )

      u0 =  one/two * fC + one   / two * fR
      u1 = -one/two * fL + three / two * fR

      GetUPWINDReconstruction = u0 * w0 + u1 * w1    
   
   end function GetUPWINDReconstruction
   
   
   function GetENO2Reconstruction(f0, fLL, fL, fC, fR, fRR, LeftOutbounded, RightOutbounded)
   
   
      implicit none

      real (kind = rdf) :: f0, fLL, fL, fC, fR, fRR
      logical :: LeftOutbounded , RightOutbounded

      ! local variables

      real (kind = rdf) :: df_dcsi_plus , df_dcsi_minus
      real (kind = rdf) :: a , b , signf0 , Mf
      real (kind = rdf) :: GetENO2Reconstruction


      if ( LeftOutbounded ) then


         !--------------
         ! (∂f/∂ξ)+ 
         !--------------

         a = fR  - ( two * fC ) + fL
         b = fRR - ( two * fR ) + fC

         Mf = a

         if( abs(b) < abs(a) ) Mf = b

         df_dcsi_plus  = fR - fC - one_half * Mf 
         

         !----------------------------------
         ! (∂f/∂ξ)- only one option for Mf
         !----------------------------------

         a = fR - ( two * fC ) + fL

         Mf = a

         df_dcsi_minus = fC - fL + one_half * Mf 

      else if ( RightOutbounded ) then

         !----------------------------------
         ! (∂f/∂ξ)+ only one option for Mf
         !----------------------------------

         a = fR  - ( two * fC ) + fL

         Mf = a

         df_dcsi_plus  = fR - fC - one_half * Mf 
         

         !--------------
         ! (∂f/∂ξ)- 
         !--------------

         a = fR - ( two * fC ) + fL
         b = fC - ( two * fL ) + fLL

         Mf = a

         if( abs(b) < abs(a) ) Mf = b

         df_dcsi_minus = fC - fL + one_half * Mf 

      else ! not outbounded stencil

         !--------------
         ! (∂f/∂ξ)+ 
         !--------------

         a = fR  - ( two * fC ) + fL
         b = fRR - ( two * fR ) + fC

         Mf = a

         if( abs(b) < abs(a) ) Mf = b

         df_dcsi_plus  = fR - fC - one_half * Mf 
         

         !--------------
         ! (∂f/∂ξ)- 
         !--------------

         a = fR - ( two * fC ) + fL
         b = fC - ( two * fL ) + fLL

         Mf = a

         if( abs(b) < abs(a) ) Mf = b

         df_dcsi_minus = fC - fL + one_half * Mf 

      end if


      !--------------
      ! sig(f0)
      !--------------

      if ( abs(f0) > eps_sims) then       
         signf0 = sign( one, f0 )
      else
         signf0 =  zero            
      end if

   
      if      (       signf0 * ( fR - fC ) <  zero               &
                .and. signf0 * ( fC - fL ) < -signf0 * ( fR - fC)   ) then
   
         GetENO2Reconstruction = df_dcsi_plus
   
      else if ( signf0 * ( fC - fL ) >  zero .and. &
                signf0 * ( fR - fC ) > -signf0 * ( fC - fL ) ) then
   
         GetENO2Reconstruction = df_dcsi_minus
   
      else
   
         GetENO2Reconstruction = one_half * ( df_dcsi_plus + df_dcsi_minus )
      
      end if

   end function GetENO2Reconstruction
   
   
   
   function GetWENO3Reconstruction(f0, f1, f2, f3, WaveSpeed)

      implicit none

      real (kind = rdf), intent(in) :: f0, f1, f2, f3, WaveSpeed

      ! local variables

      real (kind = rdf) :: d0, d1
      real (kind = rdf) :: epsWENO3
      real (kind = rdf) :: fL, fC, fR
      real (kind = rdf) :: alpha0, alpha1
      real (kind = rdf) :: beta0 , beta1
      real (kind = rdf) :: u0, u1, w0, w1
      real (kind = rdf) :: GetWENO3Reconstruction

      d0       = two / three
      d1       = one / three
      epsWENO3 = 0.000001_rdf

      if ( WaveSpeed > zero ) then
         fL = f0
         fC = f1
         fR = f2
      else ! mirror
         fL = f3
         fC = f2
         fR = f1
      end if   

      beta0  = ( fC - fR )**two
      beta1  = ( fL - fC )**two

      ! print *, ' fL = ', fL
      ! print *, ' fC = ', fC
      ! print *, ' fR = ', fR
      ! print *, ' '
      ! print *, 'beta1 = ', beta1

      alpha0 = d0 / ( beta0 + epsWENO3 )**two
      alpha1 = d1 / ( beta1 + epsWENO3 )**two

      w0 = alpha0 / ( alpha0 + alpha1 )
      w1 = alpha1 / ( alpha0 + alpha1 )

      u0 = (  one/two ) * fC + ( one   / two ) * fR
      u1 = ( -one/two ) * fL + ( three / two ) * fC

      GetWENO3Reconstruction = u0 * w0 + u1 * w1      

   end function GetWENO3Reconstruction


   function BiasedDerivative(f0, f1, f2, f3, order, BiasDirection, ENOBC)

      ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      !
      ! Biased derivatives formulas for computation of derivatives 
      ! at the boundaries. They were extracted from Ferziger & Peric CFD book,
      ! section 3.7.
      !
      ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ! Input
      ! f0 = f(i), f1 = f(i ± 1), f2 = f(i ± 2), f3 = f(i ± 3) 
      real( kind = rdf ) :: f0, f1, f2, f3
      integer :: order
      integer :: BiasDirection ! 1 --> ForwardBiased, -1 <-- BackwardBiased
      logical, optional :: ENOBC

      ! Output
      real( kind = rdf ) :: BiasedDerivative

      ! local
      real( kind = rdf ) :: as

      ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      as = real( BiasDirection , kind = rdf) ! auxiliary sign

      ! 1st order biased derivative
      if (order == 1) BiasedDerivative = as * ( - one      * f0 & 
                                                + one      * f1   )

      ! 2nd order biased derivative
      if (order == 2) then 

         BiasedDerivative = as * ( - three * f0 & 
                                   + four  * f1 &
                                   - one   * f2   ) / two

         if ( present(ENOBC) ) then
            if ( ENOBC ) then 

               BiasedDerivative = ENO2Boundary( f0, f1, f2, BiasDirection )

            end if
         end if

      end if

      ! 3rd order biased derivative
      if (order == 3) BiasedDerivative = as * ( - eleven   * f0 &
                                                + eighteen * f1 &
                                                - nine     * f2 &
                                                + two      * f3   ) / six
   
   end function BiasedDerivative

   subroutine BCSignedDistanceFunction ( csivecB   , etavecB   , zetvecB    ,  &
                                         dphi_dcsi , dphi_deta , dphi_dzet  ,  & 
                                         boundary                                )

      real (kind = rdf), dimension(3) :: csivecB , etavecB , zetvecB ! Metrics at the boundary
      real (kind = rdf)               :: dphi_dcsi , dphi_deta , dphi_dzet ! gradient components
      integer                         :: boundary ! ( 1:i1, 2:im ) , (3:j1, 4:jm) , (5:k1, 6:km)

      ! local variables
   
      real (kind = rdf)            :: a,b,c,d,e,f,g,h,r
      real (kind = rdf)            :: x,y,z
      real (kind = rdf)            :: aux1, aux2
      real (kind = rdf)            :: SquareRootArgument, NormAux
      real (kind = rdf)            :: solution1 , solution2
      real (kind = rdf), parameter :: EpsSquareRoot = 1.0E-06_rdf
   
      ! variables reassingment for better readibility
      a = csivecB(1) ; d = csivecB(2) ; g = csivecB(3) ; 
      b = etavecB(1) ; e = etavecB(2) ; h = etavecB(3) ; 
      c = zetvecB(1) ; f = zetvecB(2) ; r = zetvecB(3) ; 
   

      !---------------------------------------------------
      ! i - Boundary
      !---------------------------------------------------

      if      ( boundary == 1 .or. boundary == 2 ) then
   
         y = dphi_deta
         z = dphi_dzet
   
         aux1 = a**2 + d**2 + g**2 
         aux2 = y * ( a*b + d*e +g*h ) + z * ( a*c + d*f + g*r )
         NormAux =  ( b*y + c*z )**2 + ( e*y + f*z )**2 + ( h*y + r*z )**2

         if ( NormAux > one ) NormAux = one

         SquareRootArgument = four * ( aux2**2 - aux1 * ( NormAux - one ) )

         if ( abs( SquareRootArgument ) < EpsSquareRoot ) SquareRootArgument = zero
         
         ! solutions of the quadratic equation to obtain ∂Φ / ∂ξ from ǁ ∇Φ ǁ = 1

         solution1 = one / aux1 * (   one/two * sqrt( SquareRootArgument ) - aux2 )
         solution2 = one / aux1 * (  -one/two * sqrt( SquareRootArgument ) - aux2 )

         ! ∂Φ / ∂ξ update impossing ǁ ∇Φ ǁ = 1
         
         if ( abs( dphi_dcsi - solution1 ) < abs( dphi_dcsi - solution2 ) ) then
            dphi_dcsi = solution1
         else
            dphi_dcsi = solution2
         end if
          
      !---------------------------------------------------
      ! j - Boundary
      !---------------------------------------------------
   
      else if ( boundary == 3 .or. boundary == 4 ) then
   
         x = dphi_dcsi
         z = dphi_dzet
   
         aux1 = b**2 + e**2 + h**2
         aux2 = x * ( a*b + d*e + g*h ) + z * ( b*c + e*f + h*r )
         NormAux = ( a*x + c*z )**2 + ( d*x + f*z )**2 + ( g*x + r*z )**2

         if ( NormAux > one ) NormAux = one

         SquareRootArgument = four * ( aux2**2 - aux1 * ( NormAux - one ) )

         if ( abs( SquareRootArgument ) < EpsSquareRoot ) SquareRootArgument = zero

         ! solutions of the quadratic equation to obtain ∂Φ / ∂η from ǁ ∇Φ ǁ = 1

         solution1 = one / aux1 * (   one/two * sqrt( SquareRootArgument ) - aux2 )
         solution2 = one / aux1 * (  -one/two * sqrt( SquareRootArgument ) - aux2 )

         ! ∂Φ / ∂η update impossing ǁ ∇Φ ǁ = 1

         if ( abs( dphi_deta - solution1 ) < abs( dphi_deta - solution2 ) ) then
            dphi_deta = solution1
         else
            dphi_deta = solution2
         end if

   
      !---------------------------------------------------
      ! k - Boundary
      !---------------------------------------------------

      else if ( boundary == 5 .or. boundary == 6 ) then
   
         x = dphi_dcsi
         y = dphi_deta
   
         aux1 = c**2 + f**2 + r**2
         aux2 = x * (a*c + d*f + g*r ) + y * ( b*c + e*f + h*r )
         NormAux = ( a*x + b*y )**2 + ( d*x + e*y )**2 + ( g*x + h*y )**2

         if ( NormAux > one ) NormAux = one

         SquareRootArgument = four * ( aux2**2 - aux1 * ( NormAux - one ) )

         if ( abs( SquareRootArgument ) < EpsSquareRoot ) SquareRootArgument = zero

         ! solutions of the quadratic equation to obtain ∂Φ / ∂ζ from ǁ ∇Φ ǁ = 1

         solution1 = one / aux1 * (   one/two * sqrt( SquareRootArgument ) - aux2 )
         solution2 = one / aux1 * (  -one/two * sqrt( SquareRootArgument ) - aux2 )

         ! ∂Φ / ∂ζ update impossing ǁ ∇Φ ǁ = 1

         if ( abs( dphi_dzet - solution1 ) < abs( dphi_dzet - solution2 ) ) then
            dphi_dzet = solution1
         else
            dphi_dzet = solution2
         end if
      
      end if

   end subroutine BCSignedDistanceFunction



   function ThreePointsExtrapolation( x1 , y1 , z1 , v1, &
                                      x2 , y2 , z2 , v2, &
                                      x3 , y3 , z3 , v3, &
                                      xe , ye , ze         )


      real (kind = rdf) :: x1, y1, z1, v1
      real (kind = rdf) :: x2, y2, z2, v2
      real (kind = rdf) :: x3, y3, z3, v3
      real (kind = rdf) :: xe, ye, ze
      real (kind = rdf) :: ThreePointsExtrapolation

      ! local variables
      real (kind = rdf) :: xints, yints, zints, vints
      real (kind = rdf), dimension(3) :: Pints
      real (kind = rdf) :: m ! interpolation slope

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! This function returns the extrapolated value at a point Pe
      ! in a quadrilateral, using the information of the three other 
      ! vertices. The corresponding vertices are sorted as follows:
      ! 
      ! 
      ! Nodes are numbered using right-hand rule
      !  
      !    P3, v3         P2, v2
      !       *----------*
      !      .             . 
      !     .      * Pints   .
      !    .                   .
      !   x---------------------*
      ! Pe,ve                    P1, v1
      !
      ! 
      ! Diagonals: 
      ! P1 - P3 line
      ! P2 - Pe line
      ! Pints : point where the diagonals intersect
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      Pints = TwoLinesIntersection(  x1 , y1 , z1 , &
                                     x3 , y3 , z3 , &
                                     x2 , y2 , z2 , &
                                     xe , ye , ze      )
    
      xints = Pints(1)
      yints = Pints(2)
      zints = Pints(3)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! FIRST STEP
      ! The value of the variable at the intersection point
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      m =    sqrt( ( xints - x1 )**2 + ( yints - y1 )**2 + ( zints - z1 )**2 ) &
           / sqrt( ( x3    - x1 )**2 + ( y3    - y1 )**2 + ( z3    - z1 )**2 )

      vints = v1 + m * ( v3 - v1 )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! SECOND STEP
      ! The value of the variable at the point Pe = (xe,ye,ze)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      m =    sqrt( ( xe    - x2 )**2 + ( ye    - y2 )**2 + ( ze    - z2 )**2 ) &
           / sqrt( ( xints - x2 )**2 + ( yints - y2 )**2 + ( zints - z2 )**2 )

      ThreePointsExtrapolation = v2 + m * ( vints - v2 ) 

   end function ThreePointsExtrapolation


   function TwoLinesIntersection( x1 , y1 , z1 , &
                                  x2 , y2 , z2 , &
                                  x3 , y3 , z3 , &
                                  x4 , y4 , z4      )
   
      ! input
      real (kind = rdf) :: x1, y1, z1
      real (kind = rdf) :: x2, y2, z2
      real (kind = rdf) :: x3, y3, z3
      real (kind = rdf) :: x4, y4, z4

      ! output
      real (kind = rdf), dimension(3) :: TwoLinesIntersection

      ! local variables
      real (kind = rdf), dimension(3) :: V1, V2, P1, P2, P3, P4
      real (kind = rdf) :: t , Numerator , Denominator

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! This routines return the position vector of the intersection point
      ! between lines P1-P2 and P3-P4
      ! It was retrieved from: https://math.stackexchange.com/questions/
      ! 28503/how-to-find-intersection-of-two-lines-in-3d
      !
      !  \ P1
      !   *      P4
      !    \    o
      !     \  /
      !      \/ Pints  
      !      /\ 
      !     /  \
      !    /    \
      !   o      * P2
      !  P3       \
      !
      ! Line 1 formed by P1 = (x1,y1,z1) and P2 = (x2,y2,z2)
      ! Line 2 formed by P3 = (x3,y3,z3) and P4 = (x4,y4,z4)
      ! Pints is the intersection point
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      ! P1,P2,P3 and P4 in vector form
      P1(1) = x1; P1(2) = y1 ; P1(3) = z1
      P2(1) = x2; P2(2) = y2 ; P2(3) = z2
      P3(1) = x3; P3(2) = y3 ; P3(3) = z3
      P4(1) = x4; P4(2) = y4 ; P4(3) = z4


      ! Direction vector 1

      V1 = P2 - P1
      V2 = P4 - P3

      Numerator   = ( ( P3(1) - P1(1) ) * ( -V2(2) ) - ( P3(2) - P1(2) )*( -V2(1)) )
      Denominator = ( ( V1(1) ) * ( -V2(2) ) - ( V1(2) ) * ( -V2(1) ) )

      if ( abs(Denominator) < eps_sims ) then
         Numerator   = ( ( P3(2) - P1(2) ) * ( -V2(3) ) - ( P3(3) - P1(3) )*( -V2(2)) )
         Denominator = ( ( V1(2) ) * ( -V2(3) ) - ( V1(3) ) * ( -V2(2) ) )
      end if

      if ( abs(Denominator) < eps_sims ) then
         Numerator   = ( ( P3(3) - P1(3) ) * ( -V2(1) ) - ( P3(1) - P1(1) )*( -V2(3)) )
         Denominator = ( ( V1(3) ) * ( -V2(1) ) - ( V1(1) ) * ( -V2(3) ) )
      end if

      t = Numerator / Denominator 

      TwoLinesIntersection = P1 + t * V1

   end function TwoLinesIntersection


   function GetSignPhi0 (phi0 , epslsm)

      real(kind = rdf) :: phi0, epslsm
      real(kind = rdf) :: GetSignPhi0

      if ( abs(phi0) > epslsm) then       
         GetSignPhi0 = sign( one, phi0 )
      else
         GetSignPhi0 =  phi0 / epslsm - one/pi * sin( ( pi * phi0 ) / epslsm )           
      end if   

   end function GetSignPhi0


  ! minmod function

   function minmod(a, b)
      
      real(kind = rdf) , intent(in) :: a, b
      real(kind = rdf) :: minmod

      minmod = zero

      if (a > zero .and. b > zero) minmod = min(a, b)
      if (a < zero .and. b < zero) minmod = max(a, b) 

   end function minmod


   ! Van Leer limiter

   function VanLeerLimiterDerivative(f0,f1,f2,BiasDirection) result(f_prime)


      real (kind = rdf) , intent(in) :: f0, f1, f2
      integer , intent(in) :: BiasDirection
      real (kind = rdf) :: f_prime, r
      
      ! Van Leer limiter constants
      real (kind = rdf) , parameter :: eps = 1.e-10
      real (kind = rdf) :: beta, phi
      

      if (BiasDirection == 1) then ! left boundary

         ! compute right and left differences
         r = (f2 - f1)
         f_prime = (f1 - f0)
         
         ! apply Van Leer limiter
         beta = (r + eps) / (f_prime + eps)
         phi = (beta + abs(beta)) / (one + beta)
         
         ! compute third-order finite difference approximation with limiter
         f_prime = f_prime * phi * (two - phi) + ((two - phi) / two) * (f2 - two * f1 + f0)
      
      end if


      if (BiasDirection == -1) then ! right boundary

         ! compute right and left differences
         r = (f0 - f1)
         f_prime = (f0 - f1)
         
         ! apply Van Leer limiter
         beta = (r + eps) / (f_prime + eps)
         phi = (beta + abs(beta)) / (one + beta)
         
         ! compute third-order finite difference approximation with limiter
         f_prime = f_prime * phi * (one - phi) - ((one - phi) / two) * (f0 - two * f1 + f2)
      
      end if

   end function VanLeerLimiterDerivative



   function ENO2Boundary(f0, f1, f2, BiasDirection) result(f_prime)
      
      implicit none
      
      ! input
      real ( kind = rdf ) , intent(in) :: f0, f1, f2
      integer, intent(in) :: BiasDirection

      ! output and local variables
      real ( kind = rdf ) :: f_prime, D1, D2, D3, as
      
      as = real(BiasDirection, kind = rdf)
      
      ! Initialize D1, D2, D3
      D1 = as * (f1 - f0) 
      D2 = as * (f2 - f1) 
      D3 = as * (f2 - f0) / two
      
      ! Compute the weights for the stencil

           if ( abs(D1 - D2)   < eps_sims )  then

        ! first-order derivative
        f_prime = D1

      else if ( abs(D1)        < abs(D2)  )  then

        ! second-order derivative
        f_prime = (two*D1 + D2) / three

      else if ( abs(D2)        < abs(D1)  )  then
        
        ! second-order derivative
        f_prime = (D1 + two*D2) / three

      else if ( abs(D1)       == abs(D2)  )  then

        f_prime = (four*D1 + D2) / five

      else if ( abs(D1) < abs(D3)  .and. abs(D2) <  abs(D3) ) then

        f_prime = (three*D1 + six*D2 - D3) / seven

      else if ( abs(D1) >= abs(D3) .and. abs(D2) >= abs(D3) ) then

        f_prime = D3

      else

        f_prime = D1

      endif
     
   end function ENO2Boundary
   

   subroutine WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR    , &
                                           fLL, fL , fC, fR, fRR    , &
                                           dc                       , &
                                           BiasDirection            , & 
                                           exsign                   , & 
                                           df_dcsi                    &
                                         )


      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !  Description: 
      ! 
      !  This subroutine computes a derivative with different formulations (1st order, 
      !  2nd order; centred differnce, one-sided difference ) depending on the phase of the 
      !  nodes of the stencil (given by rLL, rL, rC, rR, rRR )
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Input variables

      real ( kind = rdf ) , intent (in) :: rLL , rL , rC , rR , rRR
      real ( kind = rdf ) , intent (in) :: fLL , fL , fC , fR , fRR
      real ( kind = rdf ) , intent (in) :: dc
      integer             , intent (in) :: BiasDirection

      ! Output variables

      real ( kind = rdf ) , intent (out) :: exsign
      real ( kind = rdf ) , intent (out) :: df_dcsi


      ! Local variables

      integer :: StencilBitArray
      integer :: i
      real ( kind = rdf ) :: dc2 

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
      dc2 = one_half * dc

      exsign  = zero
      df_dcsi = zero

      ! The central node is not within the water phase
      if ( rC < one_half ) return

      StencilBitArray = 0

      if( rLL > one_half ) StencilBitArray = ibset ( StencilBitArray , 0 )
      if( rL  > one_half ) StencilBitArray = ibset ( StencilBitArray , 1 )
      if( rC  > one_half ) StencilBitArray = ibset ( StencilBitArray , 2 )
      if( rR  > one_half ) StencilBitArray = ibset ( StencilBitArray , 3 )
      if( rRR > one_half ) StencilBitArray = ibset ( StencilBitArray , 4 )

      if ( StencilBitArray == 4 ) return ! A----A----W----A----A : only rC is water

      select case ( BiasDirection ) 

         case ( 0 ) ! No Bias

            select case ( StencilBitArray )

                
               case ( 14 , 15 , 30 , 31 ) 
                  
                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !                                  |
                  ! W----W----W----W----W            |            31
                  ! LL   L    C    R    RR           |
                  !                                  |
                  ! X----W----W----W----X            |         14, 15, 30  
                  ! LL   L    C    R    RR           |
                  !                                  |
                  ! X : Any                          |
                  !---------------------------------------------------------
                  
                  exsign  = one
                  df_dcsi = dc2 * ( fR - fL )
                  return

               case ( 28 , 29 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !           /
                  !          /                       |
                  ! A----A--•-W----W----W            |            28
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !
                  !           /
                  !          /                       |
                  ! W----A--•-W----W----W            |            29
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !---------------------------------------------------------

                  exsign = one
                  df_dcsi = dc2 * ( - three * fC + four * fR - fRR )
                  return

               case ( 12 , 13 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !           /
                  !          /                       |
                  ! A----A--•-W----W----A            |            12
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !
                  !           /
                  !          /                       |
                  ! W----A--•-W----W----A            |            13
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !---------------------------------------------------------

                  exsign = one
                  ! First order approximation
                  df_dcsi = dc * ( fR - fC )
                  return


               case ( 7 , 23 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !                /
                  !               /                  |
                  ! W----W----W--•-A----A            |            7
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !
                  !                /
                  !               /                  |
                  ! W----W----W--•-A----W            |            23
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !---------------------------------------------------------


                  exsign = one
                  df_dcsi = dc2 * ( three * fC - four * fL + fLL )
                  return


               case ( 6 , 22 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !                /
                  !               /                  |
                  ! A----W----W--•-A----A            |            6
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !
                  !                /
                  !               /                  |
                  ! A----W----W--•-A----W            |            22
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !---------------------------------------------------------


                  exsign = one
                  ! First order approximation
                  df_dcsi = dc * ( fC - fL )
                  return

            end select 



         case ( -1 ) ! Left biased derivatives

            ! Set R and RR bits to zero
            StencilBitArray = ibclr ( StencilBitArray , 3 )
            StencilBitArray = ibclr ( StencilBitArray , 4 )
            
            select case ( StencilBitArray )

               case ( 7 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !                /
                  !               /                  |
                  ! W----W----W--•-D----D            |            7
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !
                  ! D : Dummy
                  !---------------------------------------------------------


                  exsign = one
                  df_dcsi = dc2 * ( three * fC - four * fL + fLL )
                  return

               case ( 6 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !                /
                  !               /                  |
                  ! A----W----W--•-D----D            |            6
                  ! LL   L    C /  R    RR           |
                  !            /                     |
                  !
                  ! D : Dummy
                  !---------------------------------------------------------

                  exsign = one
                  ! First order approximation
                  df_dcsi = dc * ( fC - fL )
                  return


            end select

         case (  1 ) ! Right biased derivatives


            ! Set R and RR bits to zero
            StencilBitArray = ibclr ( StencilBitArray , 3 )
            StencilBitArray = ibclr ( StencilBitArray , 4 )

            select case ( StencilBitArray )


               case ( 28 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !           /
                  !          /                       |
                  ! D----D--•-W----W----W            |            28
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !
                  ! D : Dummy
                  !---------------------------------------------------------


                  exsign = one
                  df_dcsi = dc2 * ( - three * fC + four * fR - fRR )
                  return

               case ( 12 )   

                  !---------------------------------------------------------               
                  ! Phase stencil                    |    Resulting integer
                  ! --------------------------------------------------------
                  !           /
                  !          /                       |
                  ! D----D--•-W----W----A            |            12
                  ! LL   L /  C    R    RR           |
                  !       /                          |
                  !
                  ! D : Dummy
                  !---------------------------------------------------------


                  exsign = one
                  ! First order approximation
                  df_dcsi = dc * ( fR - fC )
                  return

            end select

      end select

   end subroutine WaterPhase_FirstDerivative

      subroutine NearSurface_P_FirstDerivative ( phiLL , phiL , phiC , phiR , phiRR , &
                                                 pLL   , pL   , pC   , pR   , pRR   , &
                                                 dc                                 , &
                                                 BiasDirection                      , & 
                                                 exsign                             , & 
                                                 dp_dcsi                              &
                                               )

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      !  Description: 
      ! 
      !  This subroutine computes a derivative with different formulations (1st order, 
      !  2nd order; centred differnce, one-sided difference ) depending on the phase of the 
      !  nodes of the stencil (given by phiLL, phiL, phiC, phiR, phiRR )
      !
      !  It asumes that the pressure at free surface is zero
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Input variables

      real ( kind = rdf ) , intent (in) :: phiLL , phiL , phiC , phiR , phiRR
      real ( kind = rdf ) , intent (in) :: pLL , pL , pC , pR , pRR
      real ( kind = rdf ) , intent (in) :: dc
      integer             , intent (in) :: BiasDirection

      ! Output variables
      real ( kind = rdf ) , intent (out) :: exsign
      real ( kind = rdf ) , intent (out) :: dp_dcsi

      ! Local variables
      real ( kind = rdf ) :: dc2 
      real ( kind = rdf ) :: dcsi , dcsiL , dcsiR 
      real ( kind = rdf ) :: InterpCoeff
      real ( kind = rdf ) :: a , b , c
      real ( kind = rdf ) :: tol, told
      real ( kind = rdf ) :: d1, d2, lambda, denom

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
      ! remember: dc = 1/ Δξ --> dc2 = 1/ ( 2Δξ ) 
      dc2 = one_half * dc

      tol  = 1.0E-10
      told = 1.0E-03

      exsign  = zero
      dp_dcsi = zero

      select case ( BiasDirection ) 

         case ( 0 ) ! No Bias

            !! case 1
            if ( phiL > eps_sims .and. &
                 phiC > eps_sims .and. &
                 phiR > eps_sims            ) then

            !if ( phiL > tol .and. &
            !     phiC > tol .and. &
            !     phiR > tol            ) then

               exsign  = one
               dp_dcsi = dc2 * ( pR - pL )

               if ( abs( dp_dcsi ) > two ) then

                  print *, 'case 1 '
                  print *, 'dp_dcsi = ', dp_dcsi
                  stop

               end if

               return

            !! case 2
            else if ( phiL >  eps_sims    .and. &
                      phiC >  eps_sims    .and. &
                      phiR < -eps_sims               ) then

            ! case 2
            !else if ( phiL >  tol    .and. &
            !          phiC >  tol    .and. &
            !          phiR < -tol               ) then

               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               !  RIGHT NODE WITHIN THE AIR PHASE                     
               !                       
               !                water            ϕ = 0      air
               !                (ϕ>0)               \      (ϕ<0)
               !                                     \   
               !     ----o----------------o-----------x-------o---  --> ξ
               !        i-1               i         fs \     i+1
               !        (L)              (C)            \    (R)
               !
               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
               ! Interpolation coefficient:
               ! Free surface approaches node C ==>  InterpCoeff --> 0
               ! Free surface approaches node R ==>  InterpCoeff --> 1
               
               InterpCoeff = abs( phiC / ( phiC - phiR  ) )
               
               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiL = one / dc ! ΔξL = Δξ
               dcsiR = InterpCoeff * ( one / dc ) ! ΔξR = InterpCoeff * Δξ 

               d1     = dcsiL * dcsiR * ( dcsiL + dcsiR )
               d2     = dcsiL + dcsiR 
               lambda = ( one - sign( one , dcsiR - told * dcsiL  ) ) / two ! 1 or 0

               ! if dcsiR is too small, I keep d2 (1st order acc). Otherwise, I keep d1 (2nd order acc)
               denom = d1 + lambda * ( d2 - d1 )

               a = -dcsiR**2
               a = a + lambda * (-one - a)
               a = a / denom 

               b = ( dcsiR**2 - dcsiL**2 )
               b = b + lambda * (zero - b)
               b = b/denom

               c = dcsiL**2
               c = c + lambda * (zero - c)
               c = c/denom 

               exsign  = one
               dp_dcsi = a * pL + b * pC

               if ( abs( dp_dcsi ) > two ) then

                  stop
   
               end if

               return


            ! case 3
            else if ( phiL < -eps_sims    .and. &
                      phiC >  eps_sims    .and. &
                      phiR >  eps_sims               ) then

            !else if ( phiL < -tol    .and. &
            !          phiC >  tol    .and. &
            !          phiR >  tol               ) then

               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               !  LEFT NODE WITHIN THE AIR PHASE                     
               !                       
               !        air       ϕ = 0    water
               !       (ϕ<0)       /       (ϕ>0)
               !                  /
               !     ----o-------x--------o------------------o---  --> ξ
               !        i-1     /fs       i                 i+1
               !        (L)    /         (C)                (R)
               !
               ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
               ! Interpolation coefficient:
               ! Free surface approaches node C ==>  InterpCoeff --> 0
               ! Free surface approaches node L ==>  InterpCoeff --> 1
         
               InterpCoeff = abs( phiC / ( phiC - phiR  ) )
               
               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiL = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 
               dcsiR = one / dc ! ΔξR = Δξ

               d1     = dcsiL * dcsiR * ( dcsiL + dcsiR )
               d2     = dcsiL + dcsiR 
               lambda = ( one - sign( one , dcsiL - told * dcsiR  ) ) / two ! 1 or 0

               ! if dcsiR is too small, I keep d2 (1st order acc). Otherwise, I keep d1 (2nd order acc)
               denom = d1 + lambda * ( d2 - d1 )

               a = -dcsiR**2
               a = a + lambda * (zero - a)
               a = a / denom 

               b = ( dcsiR**2 - dcsiL**2 )
               b = b + lambda * (zero - b)
               b = b/denom

               c = dcsiL**2
               c = c + lambda * (one - c)
               c = c/denom 

               exsign  = one
               dp_dcsi = b * pC + c * pR

               if ( abs( dp_dcsi ) > two ) then

                  stop
   
               end if

               return


            ! case 4
            else if ( phiL < -eps_sims    .and. &
                      phiC >  eps_sims    .and. &
                      phiR < -eps_sims               ) then

            !else if ( phiL < -tol    .and. &
            !          phiC >  tol    .and. &
            !          phiR < -tol               ) then

               InterpCoeff = abs( phiC / ( phiC - phiL ) )

               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiL = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 

               InterpCoeff = abs( phiC / ( phiC - phiR  ) )

               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiR = InterpCoeff * ( one / dc ) ! ΔξR = InterpCoeff * Δξ 

               a =   -dcsiR**2              / ( ( dcsiR + tol ) * ( dcsiL + tol ) * ( dcsiR + dcsiL + tol ) )  
               b = (  dcsiR**2 - dcsiL**2 ) / ( ( dcsiR + tol ) * ( dcsiL + tol ) * ( dcsiR + dcsiL + tol ) )  
               c =    dcsiL**2              / ( ( dcsiR + tol ) * ( dcsiL + tol ) * ( dcsiR + dcsiL + tol ) )  

               exsign  = one
               dp_dcsi = a* pL + b * pC + c * pR

   
               print *, 'case 4 '
               print *, 'phiLL, phiL, phiC, phiR, phiRR =', phiLL, phiL, phiC, phiR, phiRR
               print *, 'pLL, pL, pC, pR, pRR =', pLL, pL, pC, pR, pRR
               print *, 'dcsiL, dcsiR = ', dcsiL, dcsiR
               print *, 'InterpCoeff = ', InterpCoeff
               print *, 'a,b,c = ', a,b,c
               print *, 'dp_dcsi = ', dp_dcsi
               print *, ' '
                     
               if ( abs( dp_dcsi ) > two ) then
                     
                  stop
      
               end if

               return

            ! case 5
            else if (       phiL   >  eps_sims    .and. &
                            phiC   >  eps_sims    .and. &
                      abs ( phiR ) <  eps_sims               ) then

            !else if (       phiL   >  tol    .and. &
            !                phiC   >  tol    .and. &
            !          abs ( phiR ) <  tol               ) then

               exsign  = one
               dp_dcsi = dc2 * ( -pL )

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 5 '
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 6
            else if ( phiL   >  eps_sims    .and. &
                      phiC   >  zero        .and. &
                      phiC   <  eps_sims    .and. &
                      phiR   < -eps_sims               ) then

            !else if ( phiL   >  tol         .and. &
            !          phiC   >  zero        .and. &
            !          phiC   <  tol         .and. &
            !          phiR   < -tol               ) then

               exsign  = one
               dp_dcsi = dc * ( -pL )

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 6'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 7
            else if ( phiL   < -eps_sims    .and. &
                      phiC   >  zero        .and. &
                      phiC   <  eps_sims    .and. &
                      phiR   >  eps_sims               ) then

            !else if ( phiL   < -tol    .and. &
            !          phiC   >  zero   .and. &
            !          phiC   <  tol    .and. &
            !          phiR   >  tol               ) then

               exsign  = one
               dp_dcsi = dc * ( pR )

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 7'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 8
            else if ( abs ( phiL ) <  eps_sims    .and. &
                            phiC   >  eps_sims    .and. &
                            phiR   >  eps_sims               ) then

            !else if ( abs ( phiL ) <  tol    .and. &
            !                phiC   >  tol    .and. &
            !                phiR   >  tol               ) then

               exsign  = one
               dp_dcsi = dc2 * ( pR )

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 8'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return


            ! case 9
            else if ( ( abs ( phiL ) <  eps_sims    .and. &
                        abs ( phiC ) <  eps_sims  ) .or.  &
                      ( abs ( phiC ) <  eps_sims    .and. &
                        abs ( phiR ) <  eps_sims  ) .or.  &
                      ( abs ( phiL ) <  eps_sims    .and. &
                        abs ( phiC ) <  eps_sims    .and. &
                        abs ( phiR ) <  eps_sims  )          ) then

            ! case 9
            !else if ( ( abs ( phiL ) <  tol    .and. &
            !            abs ( phiC ) <  tol  ) .or.  &
            !          ( abs ( phiC ) <  tol    .and. &
            !            abs ( phiR ) <  tol  ) .or.  &
            !          ( abs ( phiL ) <  tol    .and. &
            !            abs ( phiC ) <  tol    .and. &
            !            abs ( phiR ) <  tol  )          ) then

               print *, 'case 9'
               exsign  = one
               dp_dcsi = zero
               return

            end if

         case ( -1 ) ! Left biased derivatives

            ! case 10
            if ( phiLL > eps_sims .and. &
                 phiL  > eps_sims .and. &
                 phiC  > eps_sims            ) then


            !if ( phiLL > tol .and. &
            !     phiL  > tol .and. &
            !     phiC  > tol            ) then

               exsign  = one 
               dp_dcsi = dc * BiasedDerivative( pC , pL , pLL , zero , 2 , BiasDirection ) 

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 10'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 11
            elseif ( phiLL < -eps_sims .and. &
                     phiL  >  eps_sims .and. &
                     phiC  >  eps_sims            ) then


            !elseif ( phiLL < -tol .and. &
            !         phiL  >  tol .and. &
            !         phiC  >  tol            ) then

               InterpCoeff = abs( phiL / ( phiL - phiLL ) )

               ! If there are errors calculating InterpCoeff, I just return
               !if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               !dcsiL = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 
               !dcsi  = one / dc

               !exsign  = one
               !dp_dcsi = - ( pL * ( dcsi + dcsiL )**two - pC * ( ( dcsi + dcsiL )**two - dcsi**two ) )/ &
               !            ( ( dcsiL + tol ) * ( dcsi + tol ) * ( dcsi + dcsiL + tol) )


               
               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiL = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 
               dcsiR = one / dc ! ΔξR = Δξ

               d1     = dcsiL * dcsiR * ( dcsiL + dcsiR )
               d2     = dcsiL + dcsiR 
               lambda = ( one - sign( one , dcsiL - told * dcsiR  ) ) / two ! 1 or 0

               ! if dcsiR is too small, I keep d2 (1st order acc). Otherwise, I keep d1 (2nd order acc)
               denom = d1 + lambda * ( d2 - d1 )

               a = -dcsiR**2
               a = a + lambda * (zero - a)
               a = a / denom 

               b = ( dcsiR**2 - dcsiL**2 )
               b = b + lambda * (zero - b)
               b = b/denom

               c = dcsiL**2
               c = c + lambda * (one - c)
               c = c/denom 

               exsign  = one
               dp_dcsi = b * pL + c * pC


               if ( abs( dp_dcsi ) > two ) then
   
                  print *, 'case 11' 
                  print *, 'dp_dcsi = ', dp_dcsi
                  stop
   
               end if

               return

            ! case 12
            elseif ( phiLL < -eps_sims .and. &
                     phiL  < -eps_sims .and. &
                     phiC  >  eps_sims            ) then

            !elseif ( phiLL < -tol .and. &
            !         phiL  < -tol .and. &
            !         phiC  >  tol            ) then

               InterpCoeff = abs( phiC / ( phiC - phiL ) )

               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return

               dcsiL = InterpCoeff * ( one / dc ) 

               exsign  = one
               dp_dcsi = pC / dcsiL

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 12'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 13
            elseif (      phiLL   < -eps_sims .and. &
                     abs( phiL  ) <  eps_sims .and. &
                     abs( phiC  ) <  eps_sims            ) then

            !elseif (      phiLL   < -tol .and. &
            !         abs( phiL  ) <  tol .and. &
            !         abs( phiC  ) <  tol            ) then

               print *, 'case 13'
               exsign  = one
               dp_dcsi = zero
               return

            end if

         case (  1 ) ! Right biased derivatives


            ! case 14
            if ( phiC  > eps_sims .and. &
                 phiR  > eps_sims .and. &
                 phiRR > eps_sims            ) then

            !if ( phiC  > tol .and. &
            !     phiR  > tol .and. &
            !     phiRR > tol            ) then

               exsign  = one
               dp_dcsi = dc2 * ( -three * pC + four * pR - one * pRR ) 

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 14'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 15
            elseif ( phiC  >  eps_sims .and. &
                     phiR  >  eps_sims .and. &
                     phiRR < -eps_sims            ) then

            !elseif ( phiC  >  tol .and. &
            !         phiR  >  tol .and. &
            !         phiRR < -tol            ) then

               InterpCoeff = abs( phiR / ( phiR - phiRR ) )

               ! If there are errors calculating InterpCoeff, I just return
               !if ( InterpCoeff < zero .or. InterpCoeff > one ) return
               
               !dcsiR = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 
               !dcsi  = one / dc

               !exsign  = one
               !dp_dcsi =  ( pR * ( dcsi + dcsiR )**two - pC * ( ( dcsi + dcsiR )**two - dcsi**two ) )/ &
               !           ( ( dcsiR + tol ) * ( dcsi + tol ) * ( dcsi + dcsiR + tol ) )

               dcsiL = one / dc ! ΔξL = Δξ
               dcsiR = InterpCoeff * ( one / dc ) ! ΔξR = InterpCoeff * Δξ 

               d1     = dcsiL * dcsiR * ( dcsiL + dcsiR )
               d2     = dcsiL + dcsiR 
               lambda = ( one - sign( one , dcsiR - told * dcsiL  ) ) / two ! 1 or 0

               ! if dcsiR is too small, I keep d2 (1st order acc). Otherwise, I keep d1 (2nd order acc)
               denom = d1 + lambda * ( d2 - d1 )

               a = -dcsiR**2
               a = a + lambda * (-one - a)
               a = a / denom 

               b = ( dcsiR**2 - dcsiL**2 )
               b = b + lambda * (zero - b)
               b = b/denom

               c = dcsiL**2
               c = c + lambda * (zero - c)
               c = c/denom 

               exsign  = one
               dp_dcsi = a * pC + b * pR


               if ( abs( dp_dcsi ) > two ) then
   
                  print *, 'case 15'
                  print *, 'dp_dcsi = ', dp_dcsi
                  stop
   
               end if

               return

            ! case 16
            elseif ( phiC  >  eps_sims .and. &
                     phiR  < -eps_sims .and. &
                     phiRR < -eps_sims            ) then

            !elseif ( phiC  >  tol .and. &
            !         phiR  < -tol .and. &
            !         phiRR < -tol            ) then

               InterpCoeff = abs( phiC / ( phiC - phiR ) )

               ! If there are errors calculating InterpCoeff, I just return
               if ( InterpCoeff < zero .or. InterpCoeff > one ) return
               
               dcsiR = InterpCoeff * ( one / dc ) ! ΔξL = InterpCoeff * Δξ 

               exsign  = one
               dp_dcsi = ( -pC ) / ( dcsiR + tol )

                  if ( abs( dp_dcsi ) > two ) then
   
                     print *, 'case 16'
                     print *, 'dp_dcsi = ', dp_dcsi
                     stop
   
                  end if

               return

            ! case 17
            elseif ( abs( phiC  ) <  eps_sims .and. &
                     abs( phiR  ) <  eps_sims .and. &
                          phiRR   < -eps_sims            ) then

            !elseif ( abs( phiC  ) <  tol .and. &
            !         abs( phiR  ) <  tol .and. &
            !              phiRR   < -tol            ) then

               print *, 'case 17'
               exsign  = one
               dp_dcsi = zero
               return

            end if

      end select

   end subroutine NearSurface_P_FirstDerivative

   function LinearExtrapolation(x1,y1,z1,x2,y2,z2,f1,f2,xe,ye,ze) result( fe )

      real(kind = rdf) :: x1,y1,z1,x2,y2,z2,f1,f2,xe,ye,ze,fe
      real(kind = rdf) :: deltai, deltae
      real(kind = rdf) :: alpha1 , alpha2


      ! Distance between interior points
      deltai = sqrt( (x2-x1)**two + (y2-y1)**two + (z2-z1)**two )

      ! Distance between extrpolated point and point 1
      deltae = sqrt( (xe-x1)**two + (ye-y1)**two + (ze-z1)**two )

      alpha1 = ( deltai + deltae ) / deltai
      alpha2 =   deltae / deltai 

      fe = alpha1 * f1 - alpha2 * f2

      print *, 'a1 =', alpha1 , ' , a2', alpha2

   end function LinearExtrapolation


   subroutine GetLinearExtrapolationCoefficients( x1 , y1 , z1 ,     &
                                                  x2 , y2 , z2 ,     &
                                                  xe , ye , ze ,     &
                                                  alpha1 , alpha2    &
                                                )

      real(kind = rdf) , intent(in)  :: x1,y1,z1,x2,y2,z2,xe,ye,ze
      real(kind = rdf) , intent(out) :: alpha1 , alpha2

      real(kind = rdf) :: deltai, deltae

      ! Distance between interior points
      deltai = sqrt( (x2-x1)**two + (y2-y1)**two + (z2-z1)**two )

      ! Distance between extrpolated point and point 1
      deltae = sqrt( (xe-x1)**two + (ye-y1)**two + (ze-z1)**two )

      alpha1 = ( deltai + deltae ) / deltai
      alpha2 =  -deltae / deltai 

   end subroutine GetLinearExtrapolationCoefficients

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! ########  ######## ########  ##     ##  ######    ######   #### ##    ##  ######   
! ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ###   ## ##    ##  
! ##     ## ##       ##     ## ##     ## ##        ##         ##  ####  ## ##        
! ##     ## ######   ########  ##     ## ##   #### ##   ####  ##  ## ## ## ##   #### 
! ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##  #### ##    ##  
! ##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##   ### ##    ##  
! ########  ######## ########   #######   ######    ######   #### ##    ##  ######   
! 
! ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######      
! ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##     
! ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##           
! ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######      
! ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##     
! ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##     
! ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######     
! 
! 
! Some functions to debug benchmark cases
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


function AnalyticalPhiSphere( xrequest, yrequest, zrequest                              , &
                              xCenterOfRotation , yCenterOfRotation , zCenterOfRotation , & 
                              Radius, x0 , y0 , z0 , tstep , DeltaTime , Tperiod             ) &
                              result( PhiSphere )

   real ( kind = rdf ), intent(in) :: xrequest , yrequest, zrequest ! points where I want phi
   real ( kind = rdf ), intent(in) :: xCenterOfRotation, yCenterOfRotation, zCenterOfRotation
   real ( kind = rdf ), intent(in) :: Radius           ! Sphere's radius
   real ( kind = rdf ), intent(in) :: x0 , y0 , z0     ! Initial location of the sphere's center
   real ( kind = rdf ), intent(in) :: tstep, DeltaTime ! current time step (from 1 to Ntimesteps) and Δt
   real ( kind = rdf ), intent(in) :: Tperiod          ! Rotation period of the sphere
   real ( kind = rdf) :: PhiSphere

   ! local variables
   real ( kind = rdf ) :: xaux , yaux , zaux
   real ( kind = rdf ) :: xcs , ycs , zcs
   real ( kind = rdf ) :: ttime

   ! tstep = 0 is the initial solution where the sphere centre is located at x0, y0, z0
   
   ttime = tstep * DeltaTime
   
   ! Translation  
   xaux = x0 - xCenterOfRotation
   yaux = y0 - yCenterOfRotation
   zaux = z0 - zCenterOfRotation

   ! Apply the rotation
   xcs = xaux * cos(2*pi/Tperiod * ttime) - yaux * sin(2*pi/Tperiod * ttime)
   ycs = xaux * sin(2*pi/Tperiod * ttime) + yaux * cos(2*pi/Tperiod * ttime)
   zcs = zaux
   
   ! Re-translate the rotated sphere to the original reference frame

   xcs = xcs + xCenterOfRotation
   ycs = ycs + yCenterOfRotation
   zcs = zcs + zCenterOfRotation

   PhiSphere = Radius - sqrt( (xrequest - xcs)**two + (yrequest - ycs)**two + (zrequest - zcs)**two )

end function AnalyticalPhiSphere


function AnalyticalPhiSphereError( xrequest, yrequest, zrequest, PhiRequest ,       &
                                   xCenterOfRotation , yCenterOfRotation , zCenterOfRotation , & 
                                   Radius, x0 , y0 , z0 , tstep , DeltaTime , Tperiod ) result( PhiSphereError )

   
   real ( kind = rdf ), intent(in) :: xrequest , yrequest, zrequest, PhiRequest ! points where I want phi error
   real ( kind = rdf ), intent(in) :: xCenterOfRotation, yCenterOfRotation, zCenterOfRotation
   real ( kind = rdf ), intent(in) :: Radius           ! Sphere's radius
   real ( kind = rdf ), intent(in) :: x0 , y0 , z0     ! Initial location of the sphere's center
   real ( kind = rdf ), intent(in) :: tstep, DeltaTime ! current time step (from 1 to Ntimesteps) and Δt
   real ( kind = rdf ), intent(in) :: Tperiod          ! Rotation period of the sphere
   real ( kind = rdf) :: PhiSphereError
   
   PhiSphereError = PhiRequest - AnalyticalPhiSphere( xrequest, yrequest, zrequest, &
                                                      xCenterOfRotation , yCenterOfRotation , zCenterOfRotation , & 
                                                      Radius, x0 , y0 , z0 , tstep , DeltaTime , Tperiod ) 


end function AnalyticalPhiSphereError


subroutine PhiBiasFakeAdvectionAndDistorsion ( phi , il, iu, jl, ju, kl, ku, distortion )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   !
   ! This subroutine simulates the advection of the free surface and introduces a distorsion. 
   ! The idea of this subroutine is to generate a water domain to apply the volume correction method
   ! to see where it corrects more and where it does less.
   !
   !        .---.                    ---.          .---  
   !       /     \                       \        /
   !      /       \       =====>          \      /  
   !   __/         \__                     \____/
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   implicit none

   integer, intent(in) :: il, iu, jl, ju, kl, ku
   real ( kind = rdf ), dimension( il:iu , jl:ju , kl:ku ) , intent(inout) :: phi
   real ( kind = rdf ) :: distortion

   ! local variables
   real ( kind = rdf ), dimension( il:iu , jl:ju , kl:ku ) :: phi_aux
   integer :: k

   ! subroutine body
   
   phi_aux = phi

   do k=kl,ku
      phi(:,:,k) = -phi_aux(:,:,ku-k+1)  
   end do

   phi(2,:,7) = -distortion/3
   phi(4,:,3) = -distortion/3
   phi(6,:,3) = -distortion/3
   phi(8,:,7) = -distortion/3


end subroutine PhiBiasFakeAdvectionAndDistorsion


end module AdvectionMethods

























