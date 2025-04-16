subroutine trilinear_interpolation(nodal_var_values,coeff,var)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! This routine computes an interpolated value of nodal defined variables using Lagrange
! polynomials described in 
! 
! Delandmeter, P., & Sebille, E. V. (2019). The Parcels v2. 0 Lagrangian framework: new 
! field interpolation schemes. Geoscientific Model Development, 12(8), 3571-3584.
!
! It uses the interpolation coefficients computed in get_trilinear_intrp_coefficients.F90 
!
! This subroutine returns the interpolated variable var
!
! Jorge Sandoval, UoE/PUC. Edinburgh, October 18th, 2021.
! j.sandoval@ed.ac.uk / jcsandov@uc.cl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   implicit none

   real (kind = rdf), intent(in)  :: nodal_var_values(8), coeff(3)   
   real (kind = rdf), intent(out) :: var

   ! local variables
   real (kind = rdf)  :: lagrange_polynomial(8)   
   integer :: v
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! NOTES:
   ! Nodes nomenclature for interpolation subroutines
   ! v1 = v000, v2 = v100, v3 = v110, v4 = v010
   ! v5 = v001, v6 = v101, v7 = v111, v8 = v011
   !
   ! where v101 means v at (i+1,j,k+1)
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! Polynomials defined in Delandmeter et al., (2019)

   lagrange_polynomial(1) = (one-coeff(1))  *  (one-coeff(2))   *  (one-coeff(3))
   lagrange_polynomial(2) = coeff(1)        *  (one-coeff(2))   *  (one-coeff(3))
   lagrange_polynomial(3) = coeff(1)        *  coeff(2)         *  (one-coeff(3))
   lagrange_polynomial(4) = (one-coeff(1))  *  coeff(2)         *  (one-coeff(3))
   lagrange_polynomial(5) = (one-coeff(1))  *  (one-coeff(2))   *  coeff(3)
   lagrange_polynomial(6) = coeff(1)        *  (one-coeff(2))   *  coeff(3)
   lagrange_polynomial(7) = coeff(1)        *  coeff(2)         *  coeff(3)
   lagrange_polynomial(8) = (one-coeff(1))  *  coeff(2)         *  coeff(3)
   
   ! interpolated var initialisation
   var = zero

   do v = 1,8
      var = var + lagrange_polynomial(v)*nodal_var_values(v)
   end do
   
end subroutine trilinear_interpolation