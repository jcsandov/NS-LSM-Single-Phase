subroutine hessian_matrix_curvilinear(var, csi, eta, zet, dc, de, dz, aj, var_hessian)

   ! subroutine for computing local Hessian-Matrix of a var in curvilinear-coordinates
   !
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! input - output variables
   !
   ! values at (0,0,0) mean values at (i,j,k), so values at (-1,0,1) imply values at
   ! (i-1,j,k+1). This structure is to avoid big arrays communication between program
   ! units when it is no needed.
   ! 

   real (kind = rdf), intent(in) ::    var(-1:1,-1:1,-1:1)       , &
                                    &  csi(1:3,-1:1,-1:1,-1:1)   , &
                                    &  eta(1:3,-1:1,-1:1,-1:1)   , &
                                    &  zet(1:3,-1:1,-1:1,-1:1)   , &
                                    &  dc, de, dz                , &
                                    &  aj(-1:1,-1:1,-1:1)

   real (kind = rdf), intent(out) ::  var_hessian(1:3,1:3)
   

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! local variables

   real (kind = rdf) :: metric_tensor(1:3,1:3,-1:1,-1:1,-1:1)   
   real (kind = rdf) :: metric_tensor_im(1:3,1:3)   
   real (kind = rdf) :: metric_tensor_ip(1:3,1:3)   
   real (kind = rdf) :: metric_tensor_jm(1:3,1:3)   
   real (kind = rdf) :: metric_tensor_jp(1:3,1:3)   
   real (kind = rdf) :: metric_tensor_km(1:3,1:3)   
   real (kind = rdf) :: metric_tensor_kp(1:3,1:3)   

   real (kind = rdf) ::      aux_dcsi, aux_deta, aux_dzet   
   real (kind = rdf) ::      aj_im, aj_ip, aj_jm, aj_jp, aj_km, aj_kp   , &
                           & dvar_dcsi_im, dvar_dcsi_ip                 , & 
                           & dvar_deta_im, dvar_deta_ip                 , &
                           & dvar_dzet_im, dvar_dzet_ip                 , &
                           & dvar_dcsi_jm, dvar_dcsi_jp                 , &
                           & dvar_deta_jm, dvar_deta_jp                 , &
                           & dvar_dzet_jm, dvar_dzet_jp                 , &
                           & dvar_dcsi_km, dvar_dcsi_kp                 , &
                           & dvar_deta_km, dvar_deta_kp                 , &
                           & dvar_dzet_km, dvar_dzet_kp                 , &
                           & metric_1_xp_im, metric_1_xp_ip             , &
                           & metric_1_xq_im, metric_1_xq_ip             , &
                           & metric_2_xp_jm, metric_2_xp_jp             , &
                           & metric_2_xq_jm, metric_2_xq_jp             , &
                           & metric_3_xp_km, metric_3_xp_kp             , &
                           & metric_3_xq_km, metric_3_xq_kp             , &
                           & dc4, de4, dz4

   integer :: ph,qh

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! local variables initialisation

   metric_tensor = zero
   var_hessian = zero

   aj_im = zero
   aj_ip = zero
   aj_jm = zero
   aj_jp = zero
   aj_km = zero
   aj_kp = zero

   dvar_dcsi_im = zero
   dvar_dcsi_ip = zero
   dvar_deta_im = zero
   dvar_deta_ip = zero
   dvar_dzet_im = zero
   dvar_dzet_ip = zero

   dvar_dcsi_jm = zero
   dvar_dcsi_jp = zero
   dvar_deta_jm = zero
   dvar_deta_jp = zero
   dvar_dzet_jm = zero
   dvar_dzet_jp = zero

   dvar_dcsi_km = zero
   dvar_dcsi_kp = zero
   dvar_deta_km = zero
   dvar_deta_kp = zero
   dvar_dzet_km = zero
   dvar_dzet_kp = zero

   metric_tensor_im = zero
   metric_tensor_ip = zero
   metric_tensor_jm = zero
   metric_tensor_jp = zero
   metric_tensor_km = zero
   metric_tensor_kp = zero

   dc4 = zero; de4 = zero; dz4 = zero;

   ! metric tensor at i,j,k
   metric_tensor (1,1,0,0,0) = csi(1,0,0,0) ! ∂ξ/∂x
   metric_tensor (1,2,0,0,0) = csi(2,0,0,0) ! ∂ξ/∂y
   metric_tensor (1,3,0,0,0) = csi(3,0,0,0) ! ∂ξ/∂z

   metric_tensor (2,1,0,0,0) = eta(1,0,0,0) ! ∂η/∂x
   metric_tensor (2,2,0,0,0) = eta(2,0,0,0) ! ∂η/∂y
   metric_tensor (2,3,0,0,0) = eta(3,0,0,0) ! ∂η/∂z

   metric_tensor (3,1,0,0,0) = zet(1,0,0,0) ! ∂ζ/∂x
   metric_tensor (3,2,0,0,0) = zet(2,0,0,0) ! ∂ζ/∂y
   metric_tensor (3,3,0,0,0) = zet(3,0,0,0) ! ∂ζ/∂z
   
   !i-1,j,k
   metric_tensor (1,1,-1,0,0) = csi(1,-1,0,0) 
   metric_tensor (1,2,-1,0,0) = csi(2,-1,0,0) 
   metric_tensor (1,3,-1,0,0) = csi(3,-1,0,0) 
   metric_tensor (2,1,-1,0,0) = eta(1,-1,0,0) 
   metric_tensor (2,2,-1,0,0) = eta(2,-1,0,0) 
   metric_tensor (2,3,-1,0,0) = eta(3,-1,0,0) 
   metric_tensor (3,1,-1,0,0) = zet(1,-1,0,0) 
   metric_tensor (3,2,-1,0,0) = zet(2,-1,0,0) 
   metric_tensor (3,3,-1,0,0) = zet(3,-1,0,0)
   
   !i+1,j,k
   metric_tensor (1,1,1,0,0) = csi(1,1,0,0) 
   metric_tensor (1,2,1,0,0) = csi(2,1,0,0) 
   metric_tensor (1,3,1,0,0) = csi(3,1,0,0) 
   metric_tensor (2,1,1,0,0) = eta(1,1,0,0) 
   metric_tensor (2,2,1,0,0) = eta(2,1,0,0) 
   metric_tensor (2,3,1,0,0) = eta(3,1,0,0) 
   metric_tensor (3,1,1,0,0) = zet(1,1,0,0) 
   metric_tensor (3,2,1,0,0) = zet(2,1,0,0) 
   metric_tensor (3,3,1,0,0) = zet(3,1,0,0)
   
   !i,j-1,k
   metric_tensor (1,1,0,-1,0) = csi(1,0,-1,0) 
   metric_tensor (1,2,0,-1,0) = csi(2,0,-1,0) 
   metric_tensor (1,3,0,-1,0) = csi(3,0,-1,0) 
   metric_tensor (2,1,0,-1,0) = eta(1,0,-1,0) 
   metric_tensor (2,2,0,-1,0) = eta(2,0,-1,0) 
   metric_tensor (2,3,0,-1,0) = eta(3,0,-1,0) 
   metric_tensor (3,1,0,-1,0) = zet(1,0,-1,0) 
   metric_tensor (3,2,0,-1,0) = zet(2,0,-1,0) 
   metric_tensor (3,3,0,-1,0) = zet(3,0,-1,0)
   
   !i,j+1,k
   metric_tensor (1,1,0,1,0) = csi(1,0,1,0) 
   metric_tensor (1,2,0,1,0) = csi(2,0,1,0) 
   metric_tensor (1,3,0,1,0) = csi(3,0,1,0) 
   metric_tensor (2,1,0,1,0) = eta(1,0,1,0) 
   metric_tensor (2,2,0,1,0) = eta(2,0,1,0) 
   metric_tensor (2,3,0,1,0) = eta(3,0,1,0) 
   metric_tensor (3,1,0,1,0) = zet(1,0,1,0) 
   metric_tensor (3,2,0,1,0) = zet(2,0,1,0) 
   metric_tensor (3,3,0,1,0) = zet(3,0,1,0)
   
   !i,j,k-1
   metric_tensor (1,1,0,0,-1) = csi(1,0,0,-1) 
   metric_tensor (1,2,0,0,-1) = csi(2,0,0,-1) 
   metric_tensor (1,3,0,0,-1) = csi(3,0,0,-1) 
   metric_tensor (2,1,0,0,-1) = eta(1,0,0,-1) 
   metric_tensor (2,2,0,0,-1) = eta(2,0,0,-1) 
   metric_tensor (2,3,0,0,-1) = eta(3,0,0,-1) 
   metric_tensor (3,1,0,0,-1) = zet(1,0,0,-1) 
   metric_tensor (3,2,0,0,-1) = zet(2,0,0,-1) 
   metric_tensor (3,3,0,0,-1) = zet(3,0,0,-1)
   
   !i,j,k+1
   metric_tensor (1,1,0,0,1) = csi(1,0,0,1) 
   metric_tensor (1,2,0,0,1) = csi(2,0,0,1) 
   metric_tensor (1,3,0,0,1) = csi(3,0,0,1) 
   metric_tensor (2,1,0,0,1) = eta(1,0,0,1) 
   metric_tensor (2,2,0,0,1) = eta(2,0,0,1) 
   metric_tensor (2,3,0,0,1) = eta(3,0,0,1) 
   metric_tensor (3,1,0,0,1) = zet(1,0,0,1) 
   metric_tensor (3,2,0,0,1) = zet(2,0,0,1) 
   metric_tensor (3,3,0,0,1) = zet(3,0,0,1)
   
   
   ! We compute the derivatives of the variable at half-node points
   ! using central difference in curvilinear coordinates. The detail of
   ! the numerical implementation is described in:
   ! TO DO: generate a small report about Hessian and principal curvature 
   ! tangential space.
   ! 
   ! the nomenclature is as follows:
   ! 

   ! im = i-1/2, j, k
   ! ip = i+1/2, j, k
   ! jm = i, j-1/2, k
   ! jp = i, j+1/2, k
   ! km = i, j, k-1/2
   ! kp = i, j, k+1/2

   dc4 = one_fourth * dc
   de4 = one_fourth * de
   dz4 = one_fourth * dz

   ! i - direction
   dvar_dcsi_im = dc  * (var(0,0,0) - var(-1,0,0))
   dvar_dcsi_ip = dc  * (var(1,0,0) - var(0,0,0))
   dvar_deta_im = de4 * (var(0,1,0) - var(0,-1,0) + var(-1,1,0) - var(-1,-1,0))
   dvar_deta_ip = de4 * (var(1,1,0) - var(1,-1,0) + var(0,1,0)  - var(0,-1,0))
   dvar_dzet_im = dz4 * (var(0,0,1) - var(0,0,-1) + var(-1,0,1) - var(-1,0,-1))
   dvar_dzet_ip = dz4 * (var(1,0,1) - var(1,0,-1) + var(0,0,1)  - var(0,0,-1))

   ! j - direction
   dvar_dcsi_jm = dc4  * (var(1,0,0) - var(-1,0,0) + var(1,-1,0) - var(-1,-1,0))
   dvar_dcsi_jp = dc4  * (var(1,1,0) - var(-1,1,0) + var(1,0,0)  - var(-1,0,0))
   dvar_deta_jm = de   * (var(0,0,0) - var(0,-1,0))
   dvar_deta_jp = de   * (var(0,1,0) - var(0,0,0))
   dvar_dzet_jm = dz4  * (var(0,0,1) - var(0,0,-1) + var(0,-1,1) - var(0,-1,-1))
   dvar_dzet_jp = dz4  * (var(0,1,1) - var(0,1,-1) + var(0,0,1)  - var(0,0,-1))

   ! k - direction
   dvar_dcsi_km = dc4  * (var(1,0,0) - var(-1,0,0) + var(1,0,-1) - var(-1,0,-1))
   dvar_dcsi_kp = dc4  * (var(1,0,1) - var(-1,0,1) + var(1,0,0)  - var(-1,0,0))
   dvar_deta_km = de4  * (var(0,1,0) - var(0,-1,0) + var(0,1,-1) - var(0,-1,-1))
   dvar_deta_kp = de4  * (var(0,1,1) - var(0,-1,1) + var(0,1,0)  - var(0,-1,0))
   dvar_dzet_km = dz   * (var(0,0,0) - var(0,0,-1))
   dvar_dzet_kp = dz   * (var(0,0,1) - var(0,0,0))

   ! Here we compute the jacobian and the metric tensor at half-node points

   aj_im = one_half * (aj(-1,0,0) + aj(0,0,0))
   aj_ip = one_half * (aj(0,0,0)  + aj(1,0,0))
   aj_jm = one_half * (aj(0,-1,0) + aj(0,0,0))
   aj_jp = one_half * (aj(0,0,0)  + aj(0,1,0))
   aj_km = one_half * (aj(0,0,-1) + aj(0,0,0))
   aj_kp = one_half * (aj(0,0,0)  + aj(0,0,1))

   metric_tensor_im(1:3,1:3) = one_half * (metric_tensor(1:3,1:3,-1,0,0) + &
                                           metric_tensor(1:3,1:3, 0,0,0))

   metric_tensor_ip(1:3,1:3) = one_half * (metric_tensor(1:3,1:3, 0,0,0) + &
                                           metric_tensor(1:3,1:3, 1,0,0))

   metric_tensor_jm(1:3,1:3) = one_half * (metric_tensor(1:3,1:3,0,-1,0) + &
                                           metric_tensor(1:3,1:3,0, 0,0))

   metric_tensor_jp(1:3,1:3) = one_half * (metric_tensor(1:3,1:3,0, 0,0) + &
                                           metric_tensor(1:3,1:3,0, 1,0))

   metric_tensor_km(1:3,1:3) = one_half * (metric_tensor(1:3,1:3,0,0,-1) + &
                                           metric_tensor(1:3,1:3,0,0, 0))

   metric_tensor_kp(1:3,1:3) = one_half * (metric_tensor(1:3,1:3,0,0, 0) + &
                                           metric_tensor(1:3,1:3,0,0, 1))


   ! The Hessian matrix in curvilinear coordinates for a variable var 
   ! (consider Einstein sumation rule):
   ! H(p,q) = J * ∂/∂ξ^m(1/J * ∂ξ^m/∂x_p * ∂ξ^l/∂x_q * ∂var/∂ξ_l)
   
   do ph = 1,3 ! ph =  i-hessian (row hessian index)
      do qh = 1,3 ! qh =  j-hessian (column hessian index)

         aux_dcsi = zero; aux_deta = zero; aux_dzet = zero;

         aux_dcsi = dc * ((one/aj_ip * metric_tensor_ip(1,ph) * metric_tensor_ip(1,qh) * dvar_dcsi_ip  - &
                           one/aj_im * metric_tensor_im(1,ph) * metric_tensor_im(1,qh) * dvar_dcsi_im) + &
                          (one/aj_ip * metric_tensor_ip(1,ph) * metric_tensor_ip(2,qh) * dvar_deta_ip  - &
                           one/aj_im * metric_tensor_im(1,ph) * metric_tensor_im(2,qh) * dvar_deta_im) + &
                          (one/aj_ip * metric_tensor_ip(1,ph) * metric_tensor_ip(3,qh) * dvar_dzet_ip  - &
                           one/aj_im * metric_tensor_im(1,ph) * metric_tensor_im(3,qh) * dvar_dzet_im))

         aux_deta = de * ((one/aj_jp * metric_tensor_jp(2,ph) * metric_tensor_jp(1,qh) * dvar_dcsi_jp  - &
                           one/aj_jm * metric_tensor_jm(2,ph) * metric_tensor_jm(1,qh) * dvar_dcsi_jm) + &
                          (one/aj_jp * metric_tensor_jp(2,ph) * metric_tensor_jp(2,qh) * dvar_deta_jp  - &
                           one/aj_jm * metric_tensor_jm(2,ph) * metric_tensor_jm(2,qh) * dvar_deta_jm) + &
                          (one/aj_jp * metric_tensor_jp(2,ph) * metric_tensor_jp(3,qh) * dvar_dzet_jp  - &
                           one/aj_jm * metric_tensor_jm(2,ph) * metric_tensor_jm(3,qh) * dvar_dzet_jm))

         aux_dzet = dz * ((one/aj_kp * metric_tensor_kp(3,ph) * metric_tensor_kp(1,qh) * dvar_dcsi_kp  - &
                           one/aj_km * metric_tensor_km(3,ph) * metric_tensor_km(1,qh) * dvar_dcsi_km) + &
                          (one/aj_kp * metric_tensor_kp(3,ph) * metric_tensor_kp(2,qh) * dvar_deta_kp  - &
                           one/aj_km * metric_tensor_km(3,ph) * metric_tensor_km(2,qh) * dvar_deta_km) + &
                          (one/aj_kp * metric_tensor_kp(3,ph) * metric_tensor_kp(3,qh) * dvar_dzet_kp  - &
                           one/aj_km * metric_tensor_km(3,ph) * metric_tensor_km(3,qh) * dvar_dzet_km))

         var_hessian(ph,qh) = aj(0,0,0)*(aux_dcsi + aux_deta + aux_dzet)

      end do
   end do

end subroutine hessian_matrix_curvilinear
