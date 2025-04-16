module global_obstacle
!variables globales para obstaculo

use precision
implicit none

!real(kind = rdf), dimension(:), allocatable :: phi_zero, phi_static
!real(kind = rdf), dimension(:), allocatable :: phi,phi_n

!real(kind = rdf) :: epslsm !parametro malla epsilon

integer , dimension(8) :: i_la_o, i_lb_o, j_la_o, j_lb_o, k_la_o, k_lb_o !indices que recorren el borde del obstaculo.
integer :: is_obstacle !1 si el procesador tiene parte de la malla con obstaculo, 0 sino
integer ,dimension(8) :: obstacle_zone  !1 si es que el procesador tiene obstaculo, 0 si no
integer, dimension(4,8) :: i_idx_o, j_idx_o  !indices para exprapolacion lineal
integer, dimension(:,:,:), allocatable :: csi_obs, eta_obs  !indices para convec_quick
integer, dimension(:,:,:,:), allocatable :: obstacle_lsm_ad !identifica zona donde se debe recalcular RH en rutina de advection LSM
integer, dimension(:,:,:,:), allocatable :: boundary_obstacle !1 si estoy en el borde 1, 2,3,4 (excel)
integer, dimension(:,:,:), allocatable :: act_obstacle_ad !identifica nodos donde se aplica actualizacion advection ( valor == 1)
integer, dimension(:,:,:), allocatable :: act_obstacle_rn !idem pero rutinas de reinicializacion

integer :: li_obs_ia, li_obs_ib, li_obs_ja, li_obs_jb, li_obs_ka, li_obs_kb! Indexes to save the local set of indexes that are inside
! the obstacle to set their velocities to zero in solver_daf (J. Sandoval, 2020)

end module global_obstacle
