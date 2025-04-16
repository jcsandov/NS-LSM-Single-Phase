subroutine nts_vectors(phi, phi_grad, csi, eta, zet, dc, de, dz, aj, nvec, tvec, svec)
   
   implicit none

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! input - output variables
   !
   ! values at (0,0,0) mean values at (i,j,k), so values at (-1,0,1) imply values at
   ! (i-1,j,k+1). This structure is to avoid big arrays communication between program
   ! units when it is no needed.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   

   real (kind = rdf), intent(in) ::    phi(-1:1,-1:1,-1:1)           , &
                                    &  phi_grad(1:3,-1:1,-1:1,-1:1)  , &
                                    &  csi(1:3,-1:1,-1:1,-1:1)       , &
                                    &  eta(1:3,-1:1,-1:1,-1:1)       , &
                                    &  zet(1:3,-1:1,-1:1,-1:1)       , &
                                    &  dc, de, dz                    , &
                                    &  aj(-1:1,-1:1,-1:1)

   real (kind = rdf), dimension(3), intent(out) :: nvec, tvec, svec

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! local variables

   real (kind = rdf) :: phi_gradient_norm

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


   ! normal vector
   
   phi_gradient_norm = sqrt( phi_grad(1,0,0,0)**two + &
                             phi_grad(2,0,0,0)**two + &
                             phi_grad(3,0,0,0)**two      )
   
   nvec = zero

   nvec(1) = -phi_grad(1,0,0,0)/phi_gradient_norm
   nvec(2) = -phi_grad(2,0,0,0)/phi_gradient_norm
   nvec(3) = -phi_grad(3,0,0,0)/phi_gradient_norm
   
   ! tangential vectors
   
   call tangential_vectors(phi, phi_grad, nvec, csi, eta, zet, dc, de, dz, aj, tvec, svec )

end subroutine nts_vectors