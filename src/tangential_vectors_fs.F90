
subroutine tangential_vectors_fs( i , j , k , tvec, svec)
   
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! Computation of tangential vector system
   ! 
   ! Tangential vectors are computed using the main curvature directions as suggested  
   ! by Dr. Watanabe. The main curvature direction are the eigenvectors of phi Hessian 
   ! matrix in curvilinear coordinates.
   !
   ! The expressions for maximum and minimum curvatures and tangential vectores are
   ! detailed in
   !
   ! Albin, E., Knikker, R., Xin, S., Paschereit, C. O., & d’Angelo, Y. (2016, June). 
   ! Computational assessment of curvatures and principal directions of implicit 
   ! surfaces from 3D scalar data. In International Conference on Mathematical Methods 
   ! for Curves and Surfaces (pp. 1-22). Springer, Cham.
   !
   ! available in: https://www.archives-ouvertes.fr/hal-01486547/document
   !
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! input - output variables
   !
   ! values at (0,0,0) means values at (i,j,k), so values at (-1,0,1) imply values at
   ! (i-1,j,k+1). This structure is to avoid big arrays communication between program
   ! units when it is no needed.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   integer, intent(in) :: i,j,k   
   real (kind = rdf), dimension(3), intent(out) :: tvec, svec

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! local variables

   real (kind = rdf) ::   phi_local(-1:1,-1:1,-1:1)           , &
                          phi_grad_local(1:3,-1:1,-1:1,-1:1)  , &
                          nvec(1:3)                           , &
                          csi_local(1:3,-1:1,-1:1,-1:1)       , &
                          eta_local(1:3,-1:1,-1:1,-1:1)       , &
                          zet_local(1:3,-1:1,-1:1,-1:1)       , &
                          aj_local(-1:1,-1:1,-1:1)


   real (kind = rdf) :: phi_hessian(1:3,1:3)
   real (kind = rdf) :: tvec_guess(1:3), svec_guess(1:3)
   real (kind = rdf) :: tol, aux(1:3), aux1, aux2
   real (kind = rdf) :: phi_n, phi_uu, phi_vv, phi_uv
   real (kind = rdf) :: kh, kk, kmin, kmax, k1, k2, zet_sign
   real (kind = rdf) :: tvec_guess_mod
   real (kind = rdf) :: tvec_mod, svec_mod
   real (kind = rdf) :: t1(1:3), t2(1:3)
   integer :: idx
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   nvec = -phi_gradient(1:3,i,j,k) / norm2( phi_gradient(1:3,i,j,k) )

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   phi_local(-1:1,-1:1,-1:1)      = phi(i-1:i+1,j-1:j+1,k-1:k+1)     
   csi_local(1:3,-1:1,-1:1,-1:1)  = csi(1:3,i-1:i+1,j-1:j+1,k-1:k+1) 
   eta_local(1:3,-1:1,-1:1,-1:1)  = eta(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
   zet_local(1:3,-1:1,-1:1,-1:1)  = zet(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
   aj_local(-1:1,-1:1,-1:1)       = aj(i-1:i+1,j-1:j+1,k-1:k+1)

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! inital guesses for tangent vectors
   
   tol            = ten**(-14.0_rdf) ! tol = 10^(-14)
   tvec_guess     = zero 
   svec_guess     = zero
   t1             = zero
   t2             = zero
   tvec           = zero 
   svec           = zero
   tvec_guess_mod = zero
   phi_hessian    = zero
   aux            = zero
   
   phi_n          = zero
   phi_uu         = zero
   phi_uv         = zero
   phi_vv         = zero

   tvec_guess(1) =  nvec(3)
   tvec_guess(2) =  nvec(3)
   tvec_guess(3) = -nvec(1)-nvec(2) 
   
   tvec_guess_mod = sqrt(tvec_guess(1)**2+tvec_guess(2)**2+tvec_guess(3)**2)
   
   if (tvec_guess_mod.le.tol) then
      tvec_guess(1) = -nvec(2)-nvec(3)
      tvec_guess(2) = nvec(1)
      tvec_guess(3) = nvec(1) 
   
      tvec_guess_mod = sqrt(tvec_guess(1)**2+tvec_guess(2)**2+tvec_guess(3)**2)
   
      if (tvec_guess_mod.le.tol) then
         tvec_guess(1) = nvec(2)
         tvec_guess(2) = -nvec(3)-nvec(1)
         tvec_guess(3) = nvec(2) 
      end if
   
   end if
   
   ! svec_guess = nvec x tvec_guess
   
   svec_guess(1)=tvec_guess(2)*nvec(3)-tvec_guess(3)*nvec(2);
   svec_guess(2)=tvec_guess(3)*nvec(1)-tvec_guess(1)*nvec(3);
   svec_guess(3)=tvec_guess(1)*nvec(2)-tvec_guess(2)*nvec(1);
   
   call hessian_matrix_curvilinear( phi_local     , & 
                                    csi_local     , & 
                                    eta_local     , & 
                                    zet_local     , & 
                                    dc, de, dz    , & 
                                    aj_local      , & 
                                    phi_hessian     &
                                  )
   
   ! principal curvature directions
   
   do idx = 1,3
      phi_n = phi_n + phi_gradient(idx,i,j,k) * nvec(idx)
   end do
   
   aux = matmul(phi_hessian, svec_guess)
   
   do idx = 1,3
      phi_uv = phi_uv + tvec_guess(idx) * aux(idx)
   end do
   
   aux = matmul(phi_hessian, tvec_guess)
   
   do idx = 1,3
      phi_uu = phi_uu + tvec_guess(idx)*aux(idx)
   end do
   
   aux = matmul(phi_hessian, svec_guess)
   
   do idx = 1,3
      phi_vv = phi_vv + svec_guess(idx)*aux(idx)
   end do
   
   ! curvature computation
   
   ! Mean curvature
   
   kh = (phi_uu+phi_vv)/(two*abs(phi_n))
   
   ! Gaussian - Curvature
   
   kk = (phi_uu*phi_vv-phi_uv**two)*(phi_n**(-two))
   
   ! min-max curvatures
   
   kmin = kh-sqrt(abs(kh**two-kk))
   kmax = kh+sqrt(abs(kh**two-kk))
   
   ! zet_sign-criterior
   
   if(abs(kmin*phi_n-phi_uu).ge.abs(kmin*phi_n-phi_vv)) then
      zet_sign = one
   else
      zet_sign = -one
   end if

   ! Principal curvatures
   
   k1 = kh-sqrt(abs(kh**two-kk))*zet_sign
   k2 = kh+sqrt(abs(kh**two-kk))*zet_sign
   
   ! Principal directions
   
   t1(1) = (k1*phi_n-phi_uu)*svec_guess(1) + phi_uv*tvec_guess(1)
   t1(2) = (k1*phi_n-phi_uu)*svec_guess(2) + phi_uv*tvec_guess(2)
   t1(3) = (k1*phi_n-phi_uu)*svec_guess(3) + phi_uv*tvec_guess(3)
   
   t2(1) = (k2*phi_n-phi_vv)*tvec_guess(1) + phi_uv*svec_guess(1)
   t2(2) = (k2*phi_n-phi_vv)*tvec_guess(2) + phi_uv*svec_guess(2)
   t2(3) = (k2*phi_n-phi_vv)*tvec_guess(3) + phi_uv*svec_guess(3)
   
   ! tmin tmax
   
   aux1 = (one+zet_sign)/two ; aux2 = (one-zet_sign)/two

   do idx = 1,3
      tvec(idx) = t1(idx) * aux1 + t2(idx) * aux2
      svec(idx) = t1(idx) * aux2 - t2(idx) * aux1
   end do

   tvec_mod = sqrt(tvec(1)**2+tvec(2)**2+tvec(3)**2)
   svec_mod = sqrt(svec(1)**2+svec(2)**2+svec(3)**2)

   ! if the tangential vectors are near zero, we use the guess vector
   ! as a tangential vectors system. It could happen, for instance, when
   ! phi gradient is constant (and therefore its Hessian matrix is a zero
   ! matrix)

   if (            norm2(tvec) < tol .or.            norm2(svec) < tol .or. &
        dot_product(nvec,tvec) > tol .or. dot_product(nvec,svec) > tol ) then

      tvec(1) = tvec_guess(1)
      tvec(2) = tvec_guess(2)
      tvec(3) = tvec_guess(3)
      
      svec(1) = svec_guess(1)
      svec(2) = svec_guess(2)
      svec(3) = svec_guess(3)

      ! vectors normalisation

      tvec = tvec / norm2(tvec)
      svec = svec / norm2(svec)

   else ! tangential vector comes from principal directions

      ! vectors normalisation

      tvec = tvec / norm2(tvec)
      svec = svec / norm2(svec)

   end if

end subroutine tangential_vectors_fs