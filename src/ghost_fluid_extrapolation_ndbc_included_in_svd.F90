subroutine ghost_fluid_extrapolation( il,iu            , &
                                      jl,ju            , &
                                      kl,ku            , &
                                      igp, jgp, kgp    , &
                                      dc, de, dz       , &
                                      q                , &
                                      rh               , &
                                      csi              , &
                                      eta              , &
                                      zet              , &
                                      aj               , &
                                      phi              , & 
                                      rsign            , &
                                      x,y,z            , &
                                      iteraciontiempo    &      
                                    )

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute the velocity on the nodes adjacent to the free surface by the extrapolation 
! procedure described in Watanabe, Saruwatari & Ingram (JCP, 2008). A Least-Square 
! Method is used to compute a linear correction-function to extrapolate velocities based 
! on normal and tangential (zero-shear) dynamic boundary conditions.
! 
! The method is adapted here for a generalised curvilinear-system and collocated grid 
! configuration.
!
! For details, refer to the following report: 
! https://drive.google.com/file/d/1LM7Q2jkNWmI8R3Apj76ihLYL5RxRlrSX/view?usp=sharing 
!
! The original method is detailed in: Watanabe, Y., Saruwatari, A., & Ingram, D. M.(2008). 
! Free-surface flows under impacting droplets. Journal of Computational Physics, 227(4), 
! 2344-2365.
!
! Jorge Sandoval, UoE/PUC. Edinburgh, September 21th, 2021.
! j.sandoval@ed.ac.uk / jcsandov@uc.cl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

use global_param
use global_app
use global_mpi
use global_lsm, only: mms, nms, sweep_lsqm, radius_lsqm, phi_outputiter
use global_debug

implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! input-output arguments
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, intent(in) :: il,iu,jl,ju,kl,ku ! external nodes
integer, intent(in) :: igp, jgp, kgp ! # of ghostpoint at local processor 
real (kind = rdf), intent(in) :: dc,de,dz ! dcsi, deta, dzet
real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi , eta, zet ! metrics
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: aj ! |J|: jacobian
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: x,y,z ! grid coordinates
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi, rsign

! In-out arguments. q vector is modified in its velocity components
real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: q

! The first update to rhs are the pressure gradient terms
real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(out) :: rh

!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
integer , intent(in) :: iteraciontiempo
character (len = 256) :: debugname
!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! local  variables
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: i,j,k ! main-loop indexes
integer :: i_mysta, &
           j_mysta, &
           k_mysta, &
           i_myend, &
           j_myend, &
           k_myend

integer :: ilsearch, jlsearch, klsearch, iusearch, jusearch, kusearch ! serching for near layer

real (kind = rdf), dimension(:,:,:,:), allocatable :: phi_gradient ! ∂phi/∂x_j
real (kind = rdf), dimension(:,:,:)  , allocatable :: rsign_aux 
real (kind = rdf) :: sum_aux

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! normal-tangential vector system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
real (kind = rdf), dimension(1:3,-1:1,-1:1,-1:1) :: csi_work, eta_work, zet_work
real (kind = rdf), dimension(-1:1,-1:1,-1:1) :: aj_work
real (kind = rdf) :: dc2, de2, dz2
!real (kind = rdf) :: dphi_dcsi, dphi_deta, dphi_dzet
real (kind = rdf), dimension(-1:1,-1:1,-1:1) :: phi_work
real (kind = rdf), dimension(1:3,-1:1,-1:1,-1:1) :: phi_grad_work
real (kind = rdf), dimension(1:3) :: nvec, tvec, svec ! local normal and tangential vectors
real (kind = rdf), dimension(1:3,1:3) :: nj_ti, nj_si, nj_ni ! normal-tangential products for matrix
                                                      ! system
real (kind = rdf), dimension(1:3) :: alpha_vec ! position vectorsto the fs
real (kind = rdf), dimension(1:3) :: alpha_vec_extp ! position vectors from the extrapolation 
                                                    ! node to the fs
real (kind = rdf) :: xs, ys, zs ! free-surface location coordinates
integer :: m,n,l ! work arrays indexes

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix - lsqm system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real (kind = rdf), dimension(1:nms) :: ttvec, tsvec , tnvec ! Tt and Ts vectors of the SVD system
real (kind = rdf), dimension(0:3) :: alpha_local    ! flagged local position vector to the fs 
real (kind = rdf), dimension(1:3) :: rdiff          ! local position vector to the fs
real (kind = rdf) :: rdiff_norm ! rdiff vector norm                     
real (kind = rdf) :: dx ! local cell size
real (kind = rdf) :: exsign, aux_least_dis ! flag varibles for lsqm
real (kind = rdf) :: least_dis ! least distance from ii,jj,kk nodes to the free-surface
real (kind = rdf), dimension(1:3,1:3) :: velocity_curv_gradient !∂ui/∂xi^j
integer :: ii,jj,kk ! indexes for lsqm 
integer :: inearest,jnearest,knearest ! indexes identify local neasrest node to the free-surface 
integer :: ivs,jvs,kvs, p, coeff_loop ! indexes for Tt, Ts, Bt, Bs, A - vector-system
integer :: ims,jms,kms ! indexes for T, B, A - matrix-system
integer :: cont ! auxiliar variable
integer :: col, row ! auxiliar variable
real (kind = rdf), dimension(1:3,1:3) :: du_dx_fs_lsqm ! velocity gradient at the
                                                       ! free-surface
real (kind = rdf), dimension(1:3) :: u_fs_lsqm ! velocity at the
                                               ! free-surface

! Big arrays will be allocatables (I allocate memory when I use them only)

real (kind = rdf), dimension(:,:), allocatable   :: t_matrix_system
real (kind = rdf), dimension(:)  , allocatable   :: a_coeff_vector
real (kind = rdf), dimension(:)  , allocatable   :: b_matrix_system

! continuity correction

real (kind = rdf) :: coc ! velocity divergence
integer :: icoc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ghost fluid method
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! least_dis_extp is an array that contains the distance from the free-surface to a node that is
! candidate to be extrapolated. The array is updated if that distance is less than previous ones   
real (kind = rdf), dimension(:,:,:), allocatable :: least_dis_extp

real (kind = rdf), dimension(1:3,1:3) :: metric_tensor_vg ! local metric tensor
real (kind = rdf), dimension(1:3,1:3) :: velocity_gradient ! local velocity gradients
real (kind = rdf), dimension(1:3,1:3) :: velocity_gradient_extp ! local vel grad used for extp

! indexes for local velocity gradient computation
integer :: ivg, jvg, msum

! error arrays
real (kind = rdf), dimension(:,:,:), allocatable :: error_tdir, error_sdir, error_ndir

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Normal Dynamic Boundary Condition
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real (kind = rdf), dimension(:,:,:), allocatable :: pJ

real (kind = rdf), dimension(1:3,1:3) :: GradVelL, GradVelC, GradVelR ! velocity gradient tensor
real (kind = rdf), dimension(1:3)     :: PressureFluxAux ! ∂/∂ξ (p/J * ∂ξ/∂xj)

! ----------------------------------------------------------------------------------
! Singular - Value Decomposition (SVD) parameters
! ----------------------------------------------------------------------------------

! SGELSD
integer :: info, rank, nlvl, lwork, liwork
integer, parameter :: smlsiz=25
integer, parameter :: nrhs=1
real(kind = rdf), dimension(:), allocatable :: iwork,work
integer, parameter:: lda=max(1,mms)
integer, parameter:: ldb=max(1,max(mms,nms))
real(kind = rdf):: s(min(mms,nms)),ersvdm !!,sm,res,resf
real(kind = rdf), parameter:: rcond=0.00001_rdf

external:: dgelsd, sgelsd,ilaenv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! These variables are defined as the dgelsd documentation suggests
! https://extras.csc.fi/math/nag/mark21/pdf/F08/f08kcf.pdf
! (don't modify if not necessary)

nlvl = max(0,int( log(real(min(mms,nms))/real(smlsiz+1))/log(2.) )+1)
if(mms.ge.nms) then      
  lwork  = 12*nms+2*nms*smlsiz+8*nms*nlvl+ nms*nrhs+(smlsiz+1)*(smlsiz+1)
else
  lwork  = 12*mms+2*mms*smlsiz+8*mms*nlvl+ mms*nrhs+(smlsiz+1)*(smlsiz+1)
end if

liwork = max(1,3*min(mms,nms)*nlvl+11*min(mms,nms))

allocate(work(max(1,lwork)) , iwork(max(1,liwork)))

! ----------------------------------------------------------------------------------

i_mysta = il + igp
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp
j_myend = ju - jgp
k_myend = ku - kgp
         
! processes on the domain boundaries

if (myback  == mpi_proc_null)  i_mysta = il + igp + 1
if (myleft  == mpi_proc_null)  j_mysta = jl + jgp + 1
if (mydown  == mpi_proc_null)  k_mysta = kl + kgp + 1

if (myfront == mpi_proc_null)  i_myend = iu - igp - 1
if (myright == mpi_proc_null)  j_myend = ju - jgp - 1
if (myup    == mpi_proc_null)  k_myend = ku - kgp - 1

! memory allocation

allocate ( phi_gradient   ( 1:3, il:iu , jl:ju , kl:ku )  , &
           least_dis_extp (      il:iu , jl:ju , kl:ku )  , &
           error_tdir     (      il:iu , jl:ju , kl:ku )  , &
           error_sdir     (      il:iu , jl:ju , kl:ku )  , &
           error_ndir     (      il:iu , jl:ju , kl:ku )  , &
           pJ             (      il:iu , jl:ju , kl:ku )    & 
         )


! Error arrays initialisation

error_tdir = zero
error_sdir = zero
error_ndir = zero

! phi gradient computation
! phi_gradient = (dϕ/dx, dϕ/dx, dϕ/dz)
! TO DO: gradient and hessian should be contained in different external files
! TO DO: compute the gradient and the hessian just in the area near the free surface

least_dis_extp = 40.0_rdf
phi_gradient   = zero
dc2            = one_half * dc
de2            = one_half * de
dz2            = one_half * dz

! we get the phi gradient over the whole domain (of the current processor)
call calc_phi_gradient()


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                          
!
!  #    # ###### #       ####  ##### #  ####  ##### #   #                              
!  #    # #      #      #    #   #   # #    #   #    # #                               
!  #    # #####  #      #    #   #   # #        #     #                                
!  #    # #      #      #    #   #   # #        #     #                                
!   #  #  #      #      #    #   #   # #    #   #     #                                
!    ##   ###### ######  ####    #   #  ####    #     #                                                        
!                                                                                                                                                                           
!  ###### #    # ##### #####    ##   #####   ####  #        ##   ##### #  ####  #    # 
!  #       #  #    #   #    #  #  #  #    # #    # #       #  #    #   # #    # ##   # 
!  #####    ##     #   #    # #    # #    # #    # #      #    #   #   # #    # # #  # 
!  #        ##     #   #####  ###### #####  #    # #      ######   #   # #    # #  # # 
!  #       #  #    #   #   #  #    # #      #    # #      #    #   #   # #    # #   ## 
!  ###### #    #   #   #    # #    # #       ####  ###### #    #   #   #  ####  #    #                                                                                      
!                                                                                      
!  #    # ######   ##   #####     #        ##   #   # ###### #####                     
!  ##   # #       #  #  #    #    #       #  #   # #  #      #    #                    
!  # #  # #####  #    # #    #    #      #    #   #   #####  #    #                    
!  #  # # #      ###### #####     #      ######   #   #      #####                     
!  #   ## #      #    # #   #     #      #    #   #   #      #   #                     
!  #    # ###### #    # #    #    ###### #    #   #   ###### #    #                    
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                          

allocate ( rsign_aux ( il:iu , jl:ju , kl:ku ) )

rsign_aux = rsign

do k = il , iu
do j = jl , ju
do i = kl , ku
         
   ! ----------------------------------------------------------------------------
   ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
   ! ----------------------------------------------------------------------------
   ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
   ! phase (rsign = 1), then this if is true. This condition identifies air nodes
   ! next to the free-surface
   ! ----------------------------------------------------------------------------

   ilsearch = max( il , i-1 )
   iusearch = min( iu , i+1 )

   jlsearch = max( jl , j-1 )
   jusearch = min( ju , j+1 )

   klsearch = max( kl , k-1 )
   kusearch = min( ku , k+1 )

   ! The first if is to reduce as much as I can the area of the domain where I'm gonna
   ! the adjacency to the free-surface. This is because the search over the 
   ! rsign(i-1:i+1,j-1:j+1,k-1:k+1) neighbourhood is computationally costly.

   if ( rsign(i,j,k) > one_half .and. phi(i,j,k) < ten  ) then

      ! I check the i±1,j±1,k±1 directions. If any of them lie within the air phase,
      ! then AirNodeAround < 1/2 == .true. and the interpolation is performed

      if (   rsign( ilsearch , j , k   ) * rsign( iusearch , j , k ) * &
             rsign( i , jlsearch , k   ) * rsign( i , jusearch , k ) * &
             rsign( i , j   , klsearch ) * rsign( i , j , kusearch )  < one_half ) then

         rsign_aux(i,j,k) = zero

      end if

   end if

end do
end do
end do


do k = k_mysta, k_myend
do j = j_mysta, j_myend
do i = i_mysta, i_myend

   ! Only the layer near the free-surface
   if( ( one - rsign_aux(i,j,k) ) * rsign(i,j,k) > one_half ) then

      exsign = 0
      velocity_curv_gradient = zero ! variable initialisation before update

      ! velocity_curv_gradient = ∂u_i/∂ξ^m tensor                  
      call velocity_curv_gradient_tensor_near_layer(i,j,k-1, velocity_curv_gradient, exsign)

      metric_tensor_vg = zero ! vg means velocity gradient

      ! metric_tensor_vg = ∂ξ^m/∂x_j

      metric_tensor_vg(1,1) = csi(1,i,j,k-1)
      metric_tensor_vg(1,2) = csi(2,i,j,k-1)
      metric_tensor_vg(1,3) = csi(3,i,j,k-1)
      metric_tensor_vg(2,1) = eta(1,i,j,k-1)
      metric_tensor_vg(2,2) = eta(2,i,j,k-1)
      metric_tensor_vg(2,3) = eta(3,i,j,k-1)
      metric_tensor_vg(3,1) = zet(1,i,j,k-1)
      metric_tensor_vg(3,2) = zet(2,i,j,k-1)
      metric_tensor_vg(3,3) = zet(3,i,j,k-1)

      velocity_gradient = zero

      ! velocity_gradient = ∂u_i/∂ξ^m * ∂ξ^m/∂x_j (Einstein's summation)

      do jvg = 1,3
      do ivg = 1,3
      do msum = 1, 3

         velocity_gradient(ivg,jvg) =  velocity_gradient(ivg,jvg) + &
                                       velocity_curv_gradient(ivg,msum) * &
                                       metric_tensor_vg(msum,jvg)                                 
      end do
      end do
      end do

      ! velocity update

      q(2,i,j,k) = q(2,i,j,k-1) + exsign * velocity_gradient(1,3) * (z(i,j,k) - z(i,j,k-1))
      q(3,i,j,k) = q(3,i,j,k-1) + exsign * velocity_gradient(2,3) * (z(i,j,k) - z(i,j,k-1))
      q(4,i,j,k) = q(4,i,j,k-1) + exsign * velocity_gradient(3,3) * (z(i,j,k) - z(i,j,k-1))


   end if

end do
end do
end do


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                          
!
!    #    # ###### #       ####   ####  # ##### #   #    
!    #    # #      #      #    # #    # #   #    # #     
!    #    # #####  #      #    # #      #   #     #      
!    #    # #      #      #    # #      #   #     #      
!     #  #  #      #      #    # #    # #   #     #      
!      ##   ###### ######  ####   ####  #   #     #      
!                                                     
!   ###### #    # ##### #####    ##   #####   ####  #        ##   ##### #  ####  #    # 
!   #       #  #    #   #    #  #  #  #    # #    # #       #  #    #   # #    # ##   # 
!   #####    ##     #   #    # #    # #    # #    # #      #    #   #   # #    # # #  # 
!   #        ##     #   #####  ###### #####  #    # #      ######   #   # #    # #  # # 
!   #       #  #    #   #   #  #    # #      #    # #      #    #   #   # #    # #   ## 
!   ###### #    #   #   #    # #    # #       ####  ###### #    #   #   #  ####  #    # 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                          


! loop over the whole local domain looking for nodes to be interpolated

do k = k_mysta, k_myend
do j = j_mysta, j_myend
do i = i_mysta, i_myend
         
   ! ----------------------------------------------------------------------------
   ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
   ! ----------------------------------------------------------------------------
   ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
   ! phase (rsign = 1), then this if is true. This condition identifies air nodes
   ! next to the free-surface
   ! ----------------------------------------------------------------------------

   if ( rsign(i,j,k) < one_half ) then ! Air-phase

      sum_aux = zero
   
      do kk = k-1 , k+1
      do jj = j-1 , j+1
      do ii = i-1 , i+1
   
         sum_aux = sum_aux + rsign(ii,jj,kk)
   
      end do
      end do
      end do

   if( sum_aux > one_half ) then
            
      ! Allocating matrix arrays for SVD solver

      allocate(  t_matrix_system ( 1:mms , 1:nms  ) , &
                 a_coeff_vector  ( 1:nms )          , &
                 b_matrix_system ( 1:max(mms,nms) )   &
               )

      t_matrix_system = zero
      a_coeff_vector  = zero
      b_matrix_system = zero

      ! ----------------------------------------------------------------------
      ! 2. Normal and tangential vectors computation 
      !
      ! * Unit normal vector computation: n = grad(phi)/gradient_norm(phi)
      ! * tangential vectors computation:
      ! 
      ! tangential vectors are computed using the main 
      ! curvature directions as suggested by Dr. Watanabe. 
      ! The main curvature direction are the eigenvectors 
      ! of phi Hessian matrix in curvilinear coordinates
      ! ----------------------------------------------------------------------

      nvec = zero; tvec = zero; svec = zero;

      ! I get the part of the arrays I will need for extrapolation
      ! for enhancing the memory access management
      ! TO DO: maybe consider a better memory managemente for this call

      csi_work = zero; eta_work = zero; zet_work = zero; aj_work = zero; 
      phi_work = zero; phi_grad_work = zero;

      csi_work     (1:3,-1:1,-1:1,-1:1)  = csi(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
      eta_work     (1:3,-1:1,-1:1,-1:1)  = eta(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
      zet_work     (1:3,-1:1,-1:1,-1:1)  = zet(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
      aj_work          (-1:1,-1:1,-1:1)  = aj(i-1:i+1,j-1:j+1,k-1:k+1)
      phi_work         (-1:1,-1:1,-1:1)  = phi(i-1:i+1,j-1:j+1,k-1:k+1)
      phi_grad_work(1:3,-1:1,-1:1,-1:1)  = phi_gradient(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
      

      ! nj_ti tensor = nj*ti + ni*tj (it's a symetric tensor)
      ! where n = normal vector, t = tangential vector 1

      call nts_vectors(  phi_work, phi_grad_work,csi_work, eta_work, zet_work, &
                       & dc, de, dz, aj_work, nvec, tvec, svec)

      ! n_j * t_i tensor = nj*ti + ni*tj
      ! where n = normal vector, t = tangential vector 
      !(it's a symetric tensor)
      
      nj_ti(1,1) = nvec(1) * tvec(1) + nvec(1) * tvec(1)
      nj_ti(1,2) = nvec(2) * tvec(1) + nvec(1) * tvec(2)
      nj_ti(1,3) = nvec(3) * tvec(1) + nvec(1) * tvec(3)
      nj_ti(2,1) = nvec(1) * tvec(2) + nvec(2) * tvec(1)
      nj_ti(2,2) = nvec(2) * tvec(2) + nvec(2) * tvec(2)
      nj_ti(2,3) = nvec(3) * tvec(2) + nvec(2) * tvec(3)
      nj_ti(3,1) = nvec(1) * tvec(3) + nvec(3) * tvec(1)
      nj_ti(3,2) = nvec(2) * tvec(3) + nvec(3) * tvec(2)
      nj_ti(3,3) = nvec(3) * tvec(3) + nvec(3) * tvec(3)

      ! nj_si tensor = nj*si + ni*sj (it's a symetric tensor)
      ! where n = normal vector, s = tangential vector 2

      nj_si(1,1) = nvec(1) * svec(1) + nvec(1) * svec(1)
      nj_si(1,2) = nvec(2) * svec(1) + nvec(1) * svec(2)
      nj_si(1,3) = nvec(3) * svec(1) + nvec(1) * svec(3)
      nj_si(2,1) = nvec(1) * svec(2) + nvec(2) * svec(1)
      nj_si(2,2) = nvec(2) * svec(2) + nvec(2) * svec(2)
      nj_si(2,3) = nvec(3) * svec(2) + nvec(2) * svec(3)
      nj_si(3,1) = nvec(1) * svec(3) + nvec(3) * svec(1)
      nj_si(3,2) = nvec(2) * svec(3) + nvec(3) * svec(2)
      nj_si(3,3) = nvec(3) * svec(3) + nvec(3) * svec(3)


      nj_ni(1,1) = nvec(1) * nvec(1) + nvec(1) * nvec(1)
      nj_ni(1,2) = nvec(2) * nvec(1) + nvec(1) * nvec(2)
      nj_ni(1,3) = nvec(3) * nvec(1) + nvec(1) * nvec(3)
      nj_ni(2,1) = nvec(1) * nvec(2) + nvec(2) * nvec(1)
      nj_ni(2,2) = nvec(2) * nvec(2) + nvec(2) * nvec(2)
      nj_ni(2,3) = nvec(3) * nvec(2) + nvec(2) * nvec(3)
      nj_ni(3,1) = nvec(1) * nvec(3) + nvec(3) * nvec(1)
      nj_ni(3,2) = nvec(2) * nvec(3) + nvec(3) * nvec(2)
      nj_ni(3,3) = nvec(3) * nvec(3) + nvec(3) * nvec(3)


      ! --------------------------------------------------------------
      ! 3. Determination of the free-surface location
      ! 
      ! we compute the alpha vector, which is the position 
      ! vector from the free surface to the grid point (i,j,k) 
      !
      ! alpha = -ϕ * n 
      ! 
      ! --------------------------------------------------------------

      alpha_vec(1) = -phi(i,j,k) * nvec(1) 
      alpha_vec(2) = -phi(i,j,k) * nvec(2) 
      alpha_vec(3) = -phi(i,j,k) * nvec(3) 

      ! x_ijk = xs + α --> xs = x_ijk - α
      ! (xs,ys,zs) : position vector at the free-surface

      xs = x(i,j,k) - alpha_vec(1)
      ys = y(i,j,k) - alpha_vec(2)
      zs = z(i,j,k) - alpha_vec(3)

      ! ---------------------------------------------------------
      ! 4. Computation of T and B matrices using Least-Square  
      ! method (lsqm).
      !
      ! The following loop sweep the search area defined by +/-  
      ! sweep_lsqm nodes in every direction. For every evaluated   
      ! node, Tt, Ts and Bt, Bs vectors are evaluated
      !
      ! sweep_lsqm: # of nodes used for computing Tt, Ts and Bt, 
      ! Bs by least-square method (lsqm). It was set originally to
      ! 3. Now it is given as an input parameter in control.dat
      !
      ! For the vector system formed by the dynamic tangential
      ! boundary condition:
      ! 
      ! Tt vector: ttvec
      ! Ts vector: tsvec
      ! Bt vector: btvec
      ! Bs vector: bsvec
      !
      ! For the final matrix system: T*A = B to be solved using
      ! singular-value decomposition (SVD):
      !
      ! T matrix : t_matrix_system
      ! B vector : b_matrix_system
      ! A vector : a_coeff_vector
      !----------------------------------------------------------
      
      ! loop for lsqm computation
      ! ii, jj, kk : indexes for i,j,k neighbourhood swept. Only nodes in the
      ! water phase are considered for the lsqm (exsign flag variable controls
      ! that)

      ! dx = J^(-1/3) = (dV_ijk)^(1/3): a local length scale 
      ! for cell size

      dx = aj(i,j,k)**( -one_third )  

      ! (inearest, jnearest, knearest) are the nearest location to the free-surface
      ! in the water phase, swept for lsqm. We initialise them at the i,j,k nodes

      inearest = i
      jnearest = j     
      knearest = k

      least_dis      = 40.0_rdf ! 20.0: a huge number for initialisation
      aux_least_dis  = zero             ! flag_variable
      
      velocity_gradient_extp = zero 
      alpha_vec_extp         = zero

      do kk = max( k_mysta , k-sweep_lsqm ), min( k_myend , k+sweep_lsqm )
      do jj = max( j_mysta , j-sweep_lsqm ), min( j_myend , j+sweep_lsqm )
      do ii = max( i_mysta , i-sweep_lsqm ), min( i_myend , i+sweep_lsqm )

        ! rdiff: position vector from node ii,jj,kk to the free
        ! surface  

        rdiff(1) = x(ii,jj,kk) - xs  
        rdiff(2) = y(ii,jj,kk) - ys  
        rdiff(3) = z(ii,jj,kk) - zs  

        rdiff_norm = norm2( rdiff ) !sqrt(rdiff(1)**2+rdiff(2)**2+rdiff(3)**2)

        ! exsign is a flag variable whose value can be 1 or 0
        ! 
        ! 1: node ii,jj,kk is within the water phase and its
        ! position lies inside the defined radius. Originally,
        ! radius_lsqm = 3.5. It is now given as an input parameter
        ! in control.dat
        ! 
        ! 0: otherwise
        ! 
        ! when exsign = 0 (air-phase/outbound range), the local T and B terms
        ! added to the global summation are both zero (not considered for lsqm)

        exsign = rsign(ii,jj,kk) * ( sign( one , radius_lsqm * dx - rdiff_norm ) + one )/two

        ! alpha_local is the position vector between xs and the point ii,jj,kk  

        ! alpha_vec = (alpha, beta, gamma)^T
        ! a, b, c, d : linear correction-function coefficients

        alpha_local(0) = exsign *  one       ! because it multiplies a_ij coefficient
        alpha_local(1) = exsign * (rdiff(1)) ! b_ij -> alpha_ii,jj,kk
        alpha_local(2) = exsign * (rdiff(2)) ! c_ij -> beta_ii,jj,kk
        alpha_local(3) = exsign * (rdiff(3)) ! d_ij -> gamma_ii,jj,kk

        !-----------------------------------------------------------------
        ! These two loops (jvs and ivs) construct the Tt and Ts vectors
        ! at ii,jj,kk point (ivs = i-vector system)
        !-----------------------------------------------------------------

        ! counter for Tt and Ts construction. cont=1..36
        cont  = 1 
        
        ! Tt and Ts vectors initialisation at node ii,jj,kk
        ttvec = zero
        tsvec = zero
        tnvec = zero

        do jvs = 1,3 ! jvs : j-vector system
           do ivs = 1,3 ! ivs : i-vector system
              do coeff_loop = 0,3 ! a, b, c, d components loop for Tt and Ts vecs
                    
                 ! sizes: ttvec(1:36), tsvec(1:36) 
                 
                 ttvec(cont) = ttvec(cont) +  nj_ti(ivs,jvs) &
                                         &  * alpha_local(coeff_loop)
                    
                 tsvec(cont) = tsvec(cont) +  nj_si(ivs,jvs) &
                                         &  * alpha_local(coeff_loop)

                 tnvec(cont) = tnvec(cont) +  nj_ni(ivs,jvs) &
                                         &  * alpha_local(coeff_loop)

                 ! counter update
                 cont=cont+1 

              end do
           end do 
        end do

        !----------------------------------------------------------------------
        ! T-matrix construction. It is separated in two loops to enhance
        ! memory-access locality 
        !
        ! Tmatrix_{IJ} = sum_N(Tt(I)*Tt(J))...
        !                .
        !                .
        !                .
        !                sum_N(Ts(I)*Ts(J))) (two vertical blocks)
        !
        ! mms = m-matrix sytem. # of rows of the LSQM matrix system = # eqs of
        ! the system (overdetermined system).
        ! nms = n-matrix sytem. # of columns of the LSQM matrix system = # of
        ! A vector coefficients (by default = 36)
        !----------------------------------------------------------------------

        ! vertical block 1 (Tt vector)
        
        do col = 1,nms
           do row = 1,nms
              t_matrix_system(row,col) = t_matrix_system(row,col) + &
                                         & ttvec(row)*ttvec(col) 
           end do
        end do

        !  vertical block 2 (Ts vector)

        do col = 1,nms
           do row = 1,nms
              t_matrix_system(nms+row,col) = t_matrix_system(nms+row,col) + &
                                            & tsvec(row)*tsvec(col) 
           end do
        end do


        !  vertical block 3 (Tn vector)

        do col = 1,nms
           do row = 1,nms
              t_matrix_system(2*nms+row,col) = t_matrix_system(2*nms+row,col) + &
                                               tnvec(row)*tnvec(col) 
           end do
        end do

        !----------------------------------------------------------------------
        ! B-vector construction. It is separated in two loops to enhance
        ! memory-access locality
        !
        !         /  ∂ui               ∂uj            \
        ! Bt  =  (  ----- (xs + α) +  ----- (xs + α)   ) * nj * ti         
        !         \  ∂xj               ∂xi            /
        !       
        !         /  ∂ui               ∂uj            \
        ! Bs  =  (  ----- (xs + α) +  ----- (xs + α)   ) * nj * si         
        !         \  ∂xj               ∂xi            /
        !         
        ! ims, jms : i-matrix system, j-matrix system
        !----------------------------------------------------------------------

        ! compute the velocity "curvilinear-gradient" tensor at node ii,jj,kk
        ! (a tensor with the derivatives of cartesian velocities in curvilinear
        ! directions)
        ! the flag variable exsign is updated in the routine 
        ! velocity_curv_gradient_tensor to take into account if the gradient is
        ! computed by using only water-phase nodes (in that case exsign = 1). 
        ! If it is not possible, exsign = 0, the routines returns a zero 
        ! curvilinear velocity gradient, and therefore the node ii,jj,kk is not 
        ! considered as a valid node for the lsqm (it adds zero to the B vector
        ! computation at b_matrix_system loop).

        velocity_curv_gradient = zero ! variable initialisation before update

        ! velocity_curv_gradient = ∂u_i/∂ξ^m tensor                  
        call velocity_curv_gradient_tensor(ii,jj,kk, velocity_curv_gradient, exsign)

        metric_tensor_vg = zero ! vg means velocity gradient

        ! metric_tensor_vg = ∂ξ^m/∂x_j

        metric_tensor_vg(1,1) = csi(1,ii,jj,kk)
        metric_tensor_vg(1,2) = csi(2,ii,jj,kk)
        metric_tensor_vg(1,3) = csi(3,ii,jj,kk)
        metric_tensor_vg(2,1) = eta(1,ii,jj,kk)
        metric_tensor_vg(2,2) = eta(2,ii,jj,kk)
        metric_tensor_vg(2,3) = eta(3,ii,jj,kk)
        metric_tensor_vg(3,1) = zet(1,ii,jj,kk)
        metric_tensor_vg(3,2) = zet(2,ii,jj,kk)
        metric_tensor_vg(3,3) = zet(3,ii,jj,kk)

        velocity_gradient = zero

        ! velocity_gradient = ∂u_i/∂ξ^m * ∂ξ^m/∂x_j (Einstein's summation)

        do jvg = 1,3
           do ivg = 1,3
              do msum = 1, 3

                 velocity_gradient(ivg,jvg) =  velocity_gradient(ivg,jvg) + &
                                               velocity_curv_gradient(ivg,msum) * &
                                               metric_tensor_vg(msum,jvg)                                 
              end do
           end do
        end do

        ! aux_least_dis = 1 if exsign = 1 and the distance between xs and (ii,jj,kk)   
        !                   location is less than the distance between xs and previous  
        !                   (ii,jj,kk) location. If aux_least_dis  = 1, then inearest,
        !                   jnearest and knearest are updated to:
        !                   inearest = ii, jnearest = jj and knearest = kk
        !                   and the velocity gradient for extrapolation is updated to
        !                   the velocity gradient at ii,jj,kk
        !   
        ! aux_least_dis = 0 otherwise

        aux_least_dis =  exsign * (sign(one , least_dis - rdiff_norm) + one)/two

        ! least distance updating (if aux_least_dis = 1)
        least_dis = least_dis + aux_least_dis * (-least_dis + rdiff_norm)

        ! indexes for the nearest position to i,j,k updating (if aux_least_dis = 1)
        inearest = inearest + nint(aux_least_dis) * (-inearest + ii)
        jnearest = jnearest + nint(aux_least_dis) * (-jnearest + jj)
        knearest = knearest + nint(aux_least_dis) * (-knearest + kk)

        ! indexes for gradient velocity from nearest position updating (if aux_least_dis = 1)
        velocity_gradient_extp =   velocity_gradient_extp &
                                 + aux_least_dis * (-velocity_gradient_extp + velocity_gradient)

        ! alpha_vec_extp is the position vector from the fluid node from where the 
        ! extrapolation to the free-surface is carried out (inearest, jnearest, knearest)

        alpha_vec_extp = alpha_vec_extp + aux_least_dis * (-alpha_vec_extp + alpha_local(1:3))

        ! vertical block 1 (Bt vector)

        do row = 1,nms
           do jms = 1,3
              do ims = 1,3
                 b_matrix_system(row) =    b_matrix_system(row) +           &
                                           (velocity_gradient(ims,jms) +    &
                                            velocity_gradient(jms,ims))*    &
                                            nvec(jms)*tvec(ims)*ttvec(row)
              end do
           end do
        end do

        !  vertical block 2 (Bs vector)

        do row = 1,nms
           do jms = 1,3
              do ims = 1,3
                 b_matrix_system(nms+row) =     b_matrix_system(nms+row) +      &
                                               (velocity_gradient(ims,jms) +    &
                                                velocity_gradient(jms,ims))*    &
                                                nvec(jms)*svec(ims)*tsvec(row)
              end do
           end do
        end do

        !  vertical block 3 (Bn vector)

        do row = 1,nms
           do jms = 1,3
              do ims = 1,3
                 b_matrix_system(2*nms+row) =     b_matrix_system(2*nms+row) +      &
                                                 (velocity_gradient(ims,jms) +    &
                                                  velocity_gradient(jms,ims))*    &
                                                  nvec(jms)*nvec(ims)*tnvec(row)
              end do
           end do
        end do


        ! continuity correction coc = ∂u_i/∂x_i
        ! it is something that is implemented in the Hokkaido code, but for our tests
        ! it didn't make a big difference. We keep it commented in the code in case we
        ! need to activate in the future (2022/03/14)   

        coc = velocity_gradient(1,1) + velocity_gradient(2,2) + velocity_gradient(3,3) 
   
        b_matrix_system(2*nms+1) = b_matrix_system(2*nms+1) + coc * alpha_local(1)
        b_matrix_system(2*nms+2) = b_matrix_system(2*nms+2) + coc * alpha_local(2)
        b_matrix_system(2*nms+3) = b_matrix_system(2*nms+3) + coc * alpha_local(3)


        do icoc = 1,3

           t_matrix_system(3*nms+icoc,2) =    t_matrix_system(3*nms+icoc,2) &
                                            + alpha_local(icoc) * alpha_local(1)
           t_matrix_system(3*nms+icoc,3) =    t_matrix_system(3*nms+icoc,3) &
                                            + alpha_local(icoc) * alpha_local(2)
           t_matrix_system(3*nms+icoc,4) =    t_matrix_system(3*nms+icoc,4) &
                                            + alpha_local(icoc) * alpha_local(3)

           t_matrix_system(3*nms+icoc,18) =   t_matrix_system(3*nms+icoc,18) &
                                            + alpha_local(icoc) * alpha_local(1)
           t_matrix_system(3*nms+icoc,19) =   t_matrix_system(3*nms+icoc,19) &
                                            + alpha_local(icoc) * alpha_local(2)
           t_matrix_system(3*nms+icoc,20) =   t_matrix_system(3*nms+icoc,20) &
                                            + alpha_local(icoc) * alpha_local(3)

           t_matrix_system(3*nms+icoc,34) =   t_matrix_system(3*nms+icoc,34) &
                                            + alpha_local(icoc) * alpha_local(1)
           t_matrix_system(3*nms+icoc,35) =   t_matrix_system(3*nms+icoc,35) &
                                            + alpha_local(icoc) * alpha_local(2)
           t_matrix_system(3*nms+icoc,36) =   t_matrix_system(3*nms+icoc,36) &
                                            + alpha_local(icoc) * alpha_local(3)

        end do

      end do ! ii
      end do ! jj
      end do ! kk , end of lsqm loop over i,j,k neighbourhood (ii,jj,kk)

      !----------------------------------------------------------------------
      ! 5. Singular Value Decomposition (SVD) step
      !
      ! For technical specification of the input/output of the sgelsd routine
      ! the reader is referred to: https://www.ibm.com/docs/en/essl/6.2?topic
      ! =llss-sgelsd-dgelsd-cgelsd-zgelsd-linear-least-squares-solution-
      ! general-matrix-using-singular-value-decomposition
      ! 
      !----------------------------------------------------------------------
      ! The sgelsd solver compute the linear least squares solution of TA=B
      ! for a general matrix T using the singular value decomposition. 
      ! the s at the beginning of sgelsd indicates single precision. In case
      ! of double precision we must use dgelsd.
      !
      ! In this case:
      ! T : t_matrix_system
      ! B : b_matrix_system. B is overwritten in the method by the n-size  
      ! solution vector A
      ! A : a_coeff_vector
      !----------------------------------------------------------------------

      s  = zero ! sgelsd output parameter

      !call sgelsd ( mms, nms, nrhs, t_matrix_system,           &
      !              lda, b_matrix_system, ldb, s, rcond, rank, &
      !              work, lwork, iwork, info )

      call dgelsd ( mms, nms, nrhs, t_matrix_system,           &
                    lda, b_matrix_system, ldb, s, rcond, rank, &
                    work, lwork, iwork, info )

      ! We assign the A vector (a_ij coeffs) with the solution of SVD 
      ! system overwritten to b_matrix_system
      
      a_coeff_vector = zero

      do ims = 1, nms
         a_coeff_vector(ims) = b_matrix_system(ims)
      end do

      !----------------------------------------------------------------------
      ! 6. Ghost-fluid Method (gfm) step
      !----------------------------------------------------------------------
      ! Here, we use the coefficients obtained by using lsqm to extrapolate
      ! the velocity gradient to the free-surface and by means of this,   
      ! extrapolate the velocities to the free-surface and air-phase nodes
      ! 
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      ! Velocity-gradient extrapolation to the free-surface
      !----------------------------------------------------------------------

      ! du_dx_fs_lsqm : velocity gradient at the free-surface computed
      ! using the coefficients obtained by using lsqm. 
      !
      ! It goes into velocity_gradient_extrapolation_free_surface_lsqm subroutine
      ! as an inout argument

      du_dx_fs_lsqm = zero

      call velocity_gradient_extrapolation_free_surface_lsqm ( a_coeff_vector          , &
                                                               alpha_vec_extp          , & 
                                                               velocity_gradient_extp  , &
                                                               du_dx_fs_lsqm             &
                                                            )

      !---------------------------------------------------
      ! Velocity extrapolation to the free-surface
      !---------------------------------------------------

      u_fs_lsqm = zero
         
      call velocity_extrapolation_free_surface_lsqm ( alpha_vec_extp  , &
                                                      inearest        , & 
                                                      jnearest        , & 
                                                      knearest        , &
                                                      du_dx_fs_lsqm   , & 
                                                      u_fs_lsqm         &
                                                    )

      !---------------------------------------------------
      ! ghost-nodes velocity extrapolation
      !---------------------------------------------------

      call ghost_nodes_velocity_extrapolation( i  , j  , k   , &
                                               xs , ys , zs  , &
                                               u_fs_lsqm     , &
                                               du_dx_fs_lsqm   &
                                             )

      !---------------------------------------------------
      ! zero - shear condition error computation
      !---------------------------------------------------

      call error_gfm(i, j, k, nvec, tvec, svec, du_dx_fs_lsqm)

      ! deallocate big arrays
      deallocate( t_matrix_system , a_coeff_vector , b_matrix_system )

   end if ! if the node i,j,k is extrapolated
   end if ! rsign(i,j,k) == 0
end do ! i
end do ! j
end do ! k

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
! #####  #####  ######  ####   ####  #    # #####  ######     ####  #####    ##   #####  # ###### #    # #####    
! #    # #    # #      #      #      #    # #    # #         #    # #    #  #  #  #    # # #      ##   #   #      
! #    # #    # #####   ####   ####  #    # #    # #####     #      #    # #    # #    # # #####  # #  #   #      
! #####  #####  #           #      # #    # #####  #         #  ### #####  ###### #    # # #      #  # #   #      
! #      #   #  #      #    # #    # #    # #   #  #         #    # #   #  #    # #    # # #      #   ##   #      
! #      #    # ######  ####   ####   ####  #    # ######     ####  #    # #    # #####  # ###### #    #   #      
!                                                                                                                                                                                                                                                       
!
! ###### #      #    # #    #    ##### ###### #####  #    # 
! #      #      #    #  #  #       #   #      #    # ##  ## 
! #####  #      #    #   ##        #   #####  #    # # ## # 
! #      #      #    #   ##        #   #      #####  #    # 
! #      #      #    #  #  #       #   #      #   #  #    # 
! #      ######  ####  #    #      #   ###### #    # #    # 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! 
! In this code section, the linear flux term related to the pressure is calculated and incorporated into the 
! rh vector
!
!  rh = rh + PressureFlux
!  
! where
!
!                       _                                                             _
!                      |                            0                                  |
!                      | ∂/∂ξ (p/J * ∂ξ/∂x) + ∂/∂η (p/J * ∂η/∂x) + ∂/∂ζ (p/J * ∂ζ/∂x)  |
!     PressureFlux =   | ∂/∂ξ (p/J * ∂ξ/∂y) + ∂/∂η (p/J * ∂η/∂y) + ∂/∂ζ (p/J * ∂ζ/∂y)  |             
!                      | ∂/∂ξ (p/J * ∂ξ/∂z) + ∂/∂η (p/J * ∂η/∂z) + ∂/∂ζ (p/J * ∂ζ/∂z)  |
!                       -                                                             - 
!
! I'm adding the terms of every component of the vector in cumulative steps. First ∂/∂ξ terms, then ∂/∂η 
! terms and finally ∂/∂ζ terms.
! 
! If i,j,k node is next to the free surface then I need to apply the dynamic normal boundary condition at  
! the free surface to obtain the pressure at the free-surface and then the flux at i,j,k. Normal vector and  
! velocity gradients at the free surface are obtained by linear interpolation.             
!
!
!         2      ∂ui            |    
! pfs  + ---- * ----- * ni * nj |    = 0 
!         Re     ∂xj            |fs
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! pJ = p/J
pJ = zero 

! I calculate p/J within the water phase
where ( rsign > one_half ) pJ = q(1,:,:,:) / aj

do k = k_mysta, k_myend
do j = j_mysta, j_myend
do i = i_mysta, i_myend

   ! I just compute this term within the water phase
   if ( rsign(i,j,k) > one_half ) then

      PressureFluxAux = zero

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! ξ - direction
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! i-1 and i+1 nodes are within the water phase
      if ( rsign(i-1,j,k) > one_half .and. rsign(i+1,j,k) > one_half ) then

         rh(2,i,j,k) = dc2 * ( pJ(i+1,j,k) * csi(1,i+1,j,k) - pJ(i-1,j,k) * csi(1,i-1,j,k) )
         rh(3,i,j,k) = dc2 * ( pJ(i+1,j,k) * csi(2,i+1,j,k) - pJ(i-1,j,k) * csi(2,i-1,j,k) )
         rh(4,i,j,k) = dc2 * ( pJ(i+1,j,k) * csi(3,i+1,j,k) - pJ(i-1,j,k) * csi(3,i-1,j,k) )

      else

         ! Velocity gradient tensor initialisation. They are updated only if the nodes is within
         ! the air phase as it is needed to compute pfs at the free surface from the NDBC

         GradVelL = zero ; GradVelC = zero ; GradVelR = zero  
 
         ! Velocity gradient tensor at i,j,k
         call VelocityGradientTensor( i , j , k , GradVelC )

         ! Velocity gradient tensor at i-1,j,k if the air phase is in that direction
         if( rsign(i-1,j,k) < one_half ) call VelocityGradientTensor( i-1 , j , k , GradVelL ) 
         ! Velocity gradient tensor at i+1,j,k if the air phase is in that direction
         if( rsign(i+1,j,k) < one_half ) call VelocityGradientTensor( i+1 , j , k , GradVelR ) 

         ! To compute ∂/∂ξ (p/J * ∂ξ/∂xj) I use the NDBC to calculate the pressure at the free
         ! surface
         PressureFluxAux = PressureLinearFlux_NDBC( q(1,i-1,j,k)            , & ! pL
                                                    q(1,i  ,j,k)            , & ! pC
                                                    q(1,i+1,j,k)            , & ! pR
                                                    aj(i-1,j,k)             , & ! JL
                                                    aj(i  ,j,k)             , & ! JC
                                                    aj(i+1,j,k)             , & ! JR
                                                    csi(:,i-1,j,k)          , & ! ∂ξ/∂xj_L 
                                                    csi(:,i  ,j,k)          , & ! ∂ξ/∂xj_C
                                                    csi(:,i+1,j,k)          , & ! ∂ξ/∂xj_R
                                                    dc                      , & ! Δξ
                                                    phi(i-1,j,k)            , & ! ϕL
                                                    phi(i  ,j,k)            , & ! ϕC
                                                    phi(i+1,j,k)            , & ! ϕR
                                                    phi_gradient(:,i-1,j,k) , & ! ∇ϕL                                  
                                                    phi_gradient(:,i  ,j,k) , & ! ∇ϕC
                                                    phi_gradient(:,i+1,j,k) , & ! ∇ϕR
                                                    GradVelL                , & ! ∇uL
                                                    GradVelC                , & ! ∇uC
                                                    GradVelR                  & ! ∇uR  
                                                  )

         rh(2,i,j,k) = PressureFluxAux(1)
         rh(3,i,j,k) = PressureFluxAux(2)
         rh(4,i,j,k) = PressureFluxAux(3)


      end if
      
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! η - direction
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! j-1 and j+1 nodes are within the water phase
      if ( rsign(i,j-1,k) > one_half .and. rsign(i,j+1,k) > one_half ) then

         rh(2,i,j,k) = rh(2,i,j,k) + de2 * ( pJ(i,j+1,k) * eta(1,i,j+1,k) - pJ(i,j-1,k) * eta(1,i,j-1,k) )
         rh(3,i,j,k) = rh(3,i,j,k) + de2 * ( pJ(i,j+1,k) * eta(2,i,j+1,k) - pJ(i,j-1,k) * eta(2,i,j-1,k) )
         rh(4,i,j,k) = rh(4,i,j,k) + de2 * ( pJ(i,j+1,k) * eta(3,i,j+1,k) - pJ(i,j-1,k) * eta(3,i,j-1,k) )

      else

         ! Velocity gradient tensor initialisation. They are updated only if the nodes is within
         ! the air phase as it is needed to compute pfs at the free surface from the NDBC

         GradVelL = zero ; GradVelC = zero ; GradVelR = zero  
 
         ! Velocity gradient tensor at i,j,k
         call VelocityGradientTensor( i , j , k , GradVelC )

         ! Velocity gradient tensor at i,j-1,k if the air phase is in that direction
         if( rsign(i,j-1,k) < one_half ) call VelocityGradientTensor( i , j-1 , k , GradVelL ) 
         ! Velocity gradient tensor at i,j+1,k if the air phase is in that direction
         if( rsign(i,j+1,k) < one_half ) call VelocityGradientTensor( i , j+1 , k , GradVelR ) 

         ! To compute ∂/∂η (p/J * ∂η/∂xj) I use the NDBC to calculate the pressure at the free
         ! surface
         PressureFluxAux = PressureLinearFlux_NDBC( q(1,i,j-1,k)            , & ! pL
                                                    q(1,i,j  ,k)            , & ! pC
                                                    q(1,i,j+1,k)            , & ! pR
                                                    aj(i,j-1,k)             , & ! JL
                                                    aj(i,j  ,k)             , & ! JC
                                                    aj(i,j+1,k)             , & ! JR
                                                    eta(:,i,j-1,k)          , & ! ∂η/∂xj_L 
                                                    eta(:,i,j  ,k)          , & ! ∂η/∂xj_C
                                                    eta(:,i,j+1,k)          , & ! ∂η/∂xj_R
                                                    de                      , & ! Δη
                                                    phi(i,j-1,k)            , & ! ϕL
                                                    phi(i,j  ,k)            , & ! ϕC
                                                    phi(i,j+1,k)            , & ! ϕR
                                                    phi_gradient(:,i,j-1,k) , & ! ∇ϕL                                  
                                                    phi_gradient(:,i,j  ,k) , & ! ∇ϕC
                                                    phi_gradient(:,i,j+1,k) , & ! ∇ϕR
                                                    GradVelL                , & ! ∇uL
                                                    GradVelC                , & ! ∇uC
                                                    GradVelR                  & ! ∇uR  
                                                  )

         rh(2,i,j,k) = rh(2,i,j,k) + PressureFluxAux(1)
         rh(3,i,j,k) = rh(3,i,j,k) + PressureFluxAux(2)
         rh(4,i,j,k) = rh(4,i,j,k) + PressureFluxAux(3)

      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! ζ - direction
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! k-1 and k+1 nodes are within the water phase
      if ( rsign(i,j,k-1) > one_half .and. rsign(i,j,k+1) > one_half ) then

         rh(2,i,j,k) = rh(2,i,j,k) + dz2 * ( pJ(i,j,k+1)*zet(1,i,j,k+1) - pJ(i,j,k-1)*zet(1,i,j,k-1) )
         rh(3,i,j,k) = rh(3,i,j,k) + dz2 * ( pJ(i,j,k+1)*zet(2,i,j,k+1) - pJ(i,j,k-1)*zet(2,i,j,k-1) )
         rh(4,i,j,k) = rh(4,i,j,k) + dz2 * ( pJ(i,j,k+1)*zet(3,i,j,k+1) - pJ(i,j,k-1)*zet(3,i,j,k-1) )

      else

         ! Velocity gradient tensor initialisation. They are updated only if the nodes is within
         ! the air phase as it is needed to compute pfs at the free surface from the NDBC

         GradVelL = zero ; GradVelC = zero ; GradVelR = zero  
 
         ! Velocity gradient tensor at i,j,k
         call VelocityGradientTensor( i , j , k , GradVelC )

         ! Velocity gradient tensor at i,j,k-1 if the air phase is in that direction
         if( rsign(i,j,k-1) < one_half ) call VelocityGradientTensor( i , j , k-1 , GradVelL ) 
         ! Velocity gradient tensor at i,j+1,k if the air phase is in that direction
         if( rsign(i,j,k+1) < one_half ) call VelocityGradientTensor( i , j , k+1 , GradVelR ) 

         ! To compute ∂/∂ζ (p/J * ∂ζ/∂xj) I use the NDBC to calculate the pressure at the free
         ! surface
         PressureFluxAux = PressureLinearFlux_NDBC( q(1,i,j,k-1)            , & ! pL
                                                    q(1,i,j,k  )            , & ! pC
                                                    q(1,i,j,k+1)            , & ! pR
                                                    aj(i,j,k-1)             , & ! JL
                                                    aj(i,j,k  )             , & ! JC
                                                    aj(i,j,k+1)             , & ! JR
                                                    zet(:,i,j,k-1)          , & ! ∂ζ/∂xj_L 
                                                    zet(:,i,j,k  )          , & ! ∂ζ/∂xj_C
                                                    zet(:,i,j,k+1)          , & ! ∂ζ/∂xj_R
                                                    dz                      , & ! Δζ
                                                    phi(i,j,k-1)            , & ! ϕL
                                                    phi(i,j,k  )            , & ! ϕC
                                                    phi(i,j,k+1)            , & ! ϕR
                                                    phi_gradient(:,i,j,k-1) , & ! ∇ϕL                                  
                                                    phi_gradient(:,i,j,k  ) , & ! ∇ϕC
                                                    phi_gradient(:,i,j,k+1) , & ! ∇ϕR
                                                    GradVelL                , & ! ∇uL
                                                    GradVelC                , & ! ∇uC
                                                    GradVelR                  & ! ∇uR
                                                  )

         rh(2,i,j,k) = rh(2,i,j,k) + PressureFluxAux(1)
         rh(3,i,j,k) = rh(3,i,j,k) + PressureFluxAux(2)
         rh(4,i,j,k) = rh(4,i,j,k) + PressureFluxAux(3)

      end if

   end if

end do
end do
end do


print *, ' '
print *, 'max error tangential BC gfm : ' , max(maxval(abs(error_tdir)),maxval(abs(error_sdir)))
print *, 'max error normal BC gfm     : ' , maxval(abs(error_ndir))
print *, ' '

write(debugname, fmt ='(a,i6.6)') 'error_gfm_tdir',iteraciontiempo
call outputD3_real( error_tdir , debugname ) 

write(debugname, fmt ='(a,i6.6)') 'error_gfm_sdir',iteraciontiempo
call outputD3_real( error_sdir , debugname ) 

write(debugname, fmt ='(a,i6.6)') 'error_gfm_ndir',iteraciontiempo
call outputD3_real( error_ndir , debugname ) 


! Deallocate big arrays
deallocate( phi_gradient , least_dis_extp , error_tdir, error_sdir, error_ndir,  pJ , rsign_aux )

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

contains

include 'calc_phi_gradient.F90'
include 'VelocityGradientTensor.F90'
include 'PhiGradientVector.F90'


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
   
   call tangential_vectors(phi, phi_grad, nvec, csi, eta, zet, dc, de, dz, aj, tvec, svec)

end subroutine nts_vectors



subroutine tangential_vectors(phi, phi_grad, nvec, csi, eta, zet, dc, de, dz, aj, tvec, svec)
   
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


   real (kind = rdf), intent(in) ::    phi(-1:1,-1:1,-1:1)           , &
                                    &  phi_grad(1:3,-1:1,-1:1,-1:1)  , &
                                    &  nvec(1:3)                     , &
                                    &  csi(1:3,-1:1,-1:1,-1:1)       , &
                                    &  eta(1:3,-1:1,-1:1,-1:1)       , &
                                    &  zet(1:3,-1:1,-1:1,-1:1)       , &
                                    &  dc, de, dz                    , &
                                    &  aj(-1:1,-1:1,-1:1)
   
   real (kind = rdf), dimension(3), intent(out) :: tvec, svec

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! local variables

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

   

   ! inital guesses for tangent vectors
   
   tol = ten**(-14.0_rdf) ! tol = 10^(-14)
   tvec_guess = zero; svec_guess = zero;
   t1 = zero; t2 = zero;
   tvec = zero; svec = zero;
   tvec_guess_mod = zero;
   phi_hessian = zero;

   tvec_guess(1) = nvec(3)
   tvec_guess(2) = nvec(3)
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
   
   call hessian_matrix_curvilinear(phi, csi, eta, zet, dc, de, dz, aj, phi_hessian)
   
   ! principal curvature directions
   
   do idx = 1,3
      phi_n = phi_n + phi_grad(idx,0,0,0)*nvec(idx)
   end do
   
   aux = matmul(phi_hessian, svec_guess)
   
   do idx = 1,3
      phi_uv = phi_uv + tvec_guess(idx)*aux(idx)
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
      tvec(idx) = t1(idx)*aux1 + t2(idx)*aux2
      svec(idx) = t1(idx)*aux2 - t2(idx)*aux1
   end do

   tvec_mod = sqrt(tvec(1)**2+tvec(2)**2+tvec(3)**2)
   svec_mod = sqrt(svec(1)**2+svec(2)**2+svec(3)**2)

   ! if the tangential vectors are near zero, we use the guess vector
   ! as a tangential vectors system. It could happen, for instance, when
   ! phi gradient is constant (and therefore its Hessian matrix is a zero
   ! matrix)

   if (tvec_mod.le.tol.or.svec_mod.le.tol) then

      tvec(1) = tvec_guess(1)
      tvec(2) = tvec_guess(2)
      tvec(3) = tvec_guess(3)
      
      svec(1) = svec_guess(1)
      svec(2) = svec_guess(2)
      svec(3) = svec_guess(3)

      ! vectors normalisation

      tvec = tvec/sqrt(tvec(1)**2+tvec(2)**2+tvec(3)**2)
      svec = svec/sqrt(svec(1)**2+svec(2)**2+svec(3)**2)

   else ! tangential vector comes from principal directions

   ! vectors normalisation

      tvec = tvec/tvec_mod
      svec = svec/svec_mod

   end if

end subroutine tangential_vectors




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




subroutine velocity_curv_gradient_tensor(i, j, k, velocity_curv_gradient, exsign)
   
   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
   
   ! local variables

   ! rsign_im, rsign_ic and rsign_ip are flag variables to decide which stencil to
   ! use depending on whether the node is and their neighbors. 
   ! For the variables:   m = minus/upwind   (i-2,i-1,i or i-1,i stencils)
   !                      c = central        (i-1,i+1 stencil)
   !                      p = plus/downwind  (i,i+1,i+2 or i,i+1 stencils).
   !
   ! Only one of {rsign_im, rsign_ic,rsign_ip} variables can be 1.0 for a node 
   ! i,j,k at the same time

   real (kind = rdf) ::   rsign_im, rsign_ic, rsign_ip, &
                        & rsign_jm, rsign_jc, rsign_jp, &
                        & rsign_km, rsign_kc, rsign_kp

   ! velocity derivatives for the different used stencils

   real (kind = rdf) ::   du_dcsi_m, du_dcsi_c, du_dcsi_p, &
                        & du_deta_m, du_deta_c, du_deta_p, &
                        & du_dzet_m, du_dzet_c, du_dzet_p, &
                        & dv_dcsi_m, dv_dcsi_c, dv_dcsi_p, &
                        & dv_deta_m, dv_deta_c, dv_deta_p, &
                        & dv_dzet_m, dv_dzet_c, dv_dzet_p, &
                        & dw_dcsi_m, dw_dcsi_c, dw_dcsi_p, &
                        & dw_deta_m, dw_deta_c, dw_deta_p, &
                        & dw_dzet_m, dw_dzet_c, dw_dzet_p

   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   real (kind = rdf) :: aux_sum ! sum over all the stencil flag variables

   du_dcsi_m = zero; du_dcsi_c = zero; du_dcsi_p = zero;
   du_deta_m = zero; du_deta_c = zero; du_deta_p = zero;
   du_dzet_m = zero; du_dzet_c = zero; du_dzet_p = zero;
   dv_dcsi_m = zero; dv_dcsi_c = zero; dv_dcsi_p = zero;
   dv_deta_m = zero; dv_deta_c = zero; dv_deta_p = zero;
   dv_dzet_m = zero; dv_dzet_c = zero; dv_dzet_p = zero;
   dw_dcsi_m = zero; dw_dcsi_c = zero; dw_dcsi_p = zero;
   dw_deta_m = zero; dw_deta_c = zero; dw_deta_p = zero;
   dw_dzet_m = zero; dw_dzet_c = zero; dw_dzet_p = zero;

   ! Only one of {rsign_im, rsign_ic, rsign_ip} can be one for a node i,j,k
   ! at the same time

   ! rsign_im = 1 --> rsign_ic = rsign_ip = 0 --> upwind scheme is used
   ! rsign_ic = 1 --> rsign_im = rsign_ip = 0 --> centred scheme is used
   ! rsign_ip = 1 --> rsign_im = rsign_ic = 0 --> downwind scheme is used

   ! if i-1,i,i+1 = water --> rsign_ic = 1. Otherwise, rsign_ic = 0
   rsign_ic = rsign(i-1,j,k)   * rsign(i,j,k) * rsign(i+1,j,k)
   ! if i = water, i-1 = water and rsign_ic = 0 --> rsign_im = 1. Otherwise, rsign_im = 0
   rsign_im = (one - rsign_ic) * rsign(i,j,k) * rsign(i-1,j,k)
   ! if i = water, i+1 = water and rsign_ic = 0 --> rsign_ip = 1. Otherwise, rsign_ip = 0
   rsign_ip = (one - rsign_ic) * rsign(i,j,k) * rsign(i+1,j,k)

   rsign_jc = rsign(i,j-1,k)   * rsign(i,j,k) * rsign(i,j+1,k)
   rsign_jm = (one - rsign_jc) * rsign(i,j,k) * rsign(i,j-1,k)
   rsign_jp = (one - rsign_jc) * rsign(i,j,k) * rsign(i,j+1,k)

   rsign_kc = rsign(i,j,k-1)   * rsign(i,j,k) * rsign(i,j,k+1)
   rsign_km = (one - rsign_kc) * rsign(i,j,k) * rsign(i,j,k-1)
   rsign_kp = (one - rsign_kc) * rsign(i,j,k) * rsign(i,j,k+1)

   ! if there is any direction where the flag variable is zero, the gradient
   ! shouldn't be considered for lsqm (the gradient is set to zero and the
   ! flag variable, exsign, too)

   aux_sum =  rsign_im + rsign_ic + rsign_ip &
            + rsign_jm + rsign_jc + rsign_jp &
            + rsign_km + rsign_kc + rsign_kp

   ! if aux_sum<3, exsign = 0

   exsign = one - (sign( one , 2.5_rdf - aux_sum) + one) / two

   if (exsign>0) then ! if exsign = 0, I don't compute any derivative
   
      ! - - - - - - - - - - - - - - 
      ! i - direction
      ! - - - - - - - - - - - - - - 
   
      du_dcsi_c = dc2 * (q(2,i+1,j,k)-q(2,i-1,j,k))
      dv_dcsi_c = dc2 * (q(3,i+1,j,k)-q(3,i-1,j,k))
      dw_dcsi_c = dc2 * (q(4,i+1,j,k)-q(4,i-1,j,k))
   
      if (i == i_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_dcsi_m = dc * (q(2,i,j,k) - q(2,i-1,j,k))
         dv_dcsi_m = dc * (q(3,i,j,k) - q(3,i-1,j,k))
         dw_dcsi_m = dc * (q(4,i,j,k) - q(4,i-1,j,k))
   
         ! downwind scheme
         du_dcsi_p = dc2 * (  - one   * q(2,i+2,j,k) &
                            & + four  * q(2,i+1,j,k) &
                            & - three * q(2,i  ,j,k) )
   
         dv_dcsi_p = dc2 * (  - one   * q(3,i+2,j,k) &
                            & + four  * q(3,i+1,j,k) &
                            & - three * q(3,i  ,j,k) )
   
         dw_dcsi_p = dc2 * (  - one   * q(4,i+2,j,k) &
                            & + four  * q(4,i+1,j,k) &
                            & - three * q(4,i  ,j,k) )
   
         ! if i+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dcsi_p = du_dcsi_p + (one - rsign(i+2,j,k)) * (- du_dcsi_p + dc * (q(2,i+1,j,k)-q(2,i,j,k)))
         dv_dcsi_p = dv_dcsi_p + (one - rsign(i+2,j,k)) * (- dv_dcsi_p + dc * (q(3,i+1,j,k)-q(3,i,j,k)))
         dw_dcsi_p = dw_dcsi_p + (one - rsign(i+2,j,k)) * (- dw_dcsi_p + dc * (q(4,i+1,j,k)-q(4,i,j,k)))
   
      else if (i == i_myend) then
   
         ! upwind scheme
         du_dcsi_m = dc2 * (    one   * q(2,i-2,j,k) &
                            & - four  * q(2,i-1,j,k) &
                            & + three * q(2,i  ,j,k) )
   
         dv_dcsi_m = dc2 * (    one   * q(3,i-2,j,k) &
                            & - four  * q(3,i-1,j,k) &
                            & + three * q(3,i  ,j,k) )
   
         dw_dcsi_m = dc2 * (    one   * q(4,i-2,j,k) &
                            & - four  * q(4,i-1,j,k) &
                            & + three * q(4,i  ,j,k) )
   
         ! if i-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dcsi_m = du_dcsi_m + (one - rsign(i-2,j,k)) * (- du_dcsi_m + dc * (q(2,i,j,k)-q(2,i-1,j,k)))
         dv_dcsi_m = dv_dcsi_m + (one - rsign(i-2,j,k)) * (- dv_dcsi_m + dc * (q(3,i,j,k)-q(3,i-1,j,k)))
         dw_dcsi_m = dw_dcsi_m + (one - rsign(i-2,j,k)) * (- dw_dcsi_m + dc * (q(4,i,j,k)-q(4,i-1,j,k)))
   
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_dcsi_p = dc * (q(2,i+1,j,k) - q(2,i,j,k))
         dv_dcsi_p = dc * (q(3,i+1,j,k) - q(3,i,j,k))
         dw_dcsi_p = dc * (q(4,i+1,j,k) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_dcsi_m = dc2 * (    one   * q(2,i-2,j,k) &
                            & - four  * q(2,i-1,j,k) &
                            & + three * q(2,i  ,j,k) )
   
         dv_dcsi_m = dc2 * (    one   * q(3,i-2,j,k) &
                            & - four  * q(3,i-1,j,k) &
                            & + three * q(3,i  ,j,k) )
   
         dw_dcsi_m = dc2 * (    one   * q(4,i-2,j,k) &
                            & - four  * q(4,i-1,j,k) &
                            & + three * q(4,i  ,j,k) )
   
         ! if i-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dcsi_m = du_dcsi_m + (one - rsign(i-2,j,k)) * (- du_dcsi_m + dc * (q(2,i,j,k)-q(2,i-1,j,k)))
         dv_dcsi_m = dv_dcsi_m + (one - rsign(i-2,j,k)) * (- dv_dcsi_m + dc * (q(3,i,j,k)-q(3,i-1,j,k)))
         dw_dcsi_m = dw_dcsi_m + (one - rsign(i-2,j,k)) * (- dw_dcsi_m + dc * (q(4,i,j,k)-q(4,i-1,j,k)))
   
   
         ! downwind derivatives
         du_dcsi_p = dc2 * (  - one   * q(2,i+2,j,k) &
                            & + four  * q(2,i+1,j,k) &
                            & - three * q(2,i  ,j,k) )
   
         dv_dcsi_p = dc2 * (  - one   * q(3,i+2,j,k) &
                            & + four  * q(3,i+1,j,k) &
                            & - three * q(3,i  ,j,k) )
   
         dw_dcsi_p = dc2 * (  - one   * q(4,i+2,j,k) &
                            & + four  * q(4,i+1,j,k) &
                            & - three * q(4,i  ,j,k) )
   
         ! if i+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dcsi_p = du_dcsi_p + (one - rsign(i+2,j,k)) * (- du_dcsi_p + dc * (q(2,i+1,j,k)-q(2,i,j,k)))
         dv_dcsi_p = dv_dcsi_p + (one - rsign(i+2,j,k)) * (- dv_dcsi_p + dc * (q(3,i+1,j,k)-q(3,i,j,k)))
         dw_dcsi_p = dw_dcsi_p + (one - rsign(i+2,j,k)) * (- dw_dcsi_p + dc * (q(4,i+1,j,k)-q(4,i,j,k)))
   
      end if
     
      ! - - - - - - - - - - - - - - 
      ! j - direction
      ! - - - - - - - - - - - - - - 
   
      du_deta_c = de2 * (q(2,i,j+1,k)-q(2,i,j-1,k))
      dv_deta_c = de2 * (q(3,i,j+1,k)-q(3,i,j-1,k))
      dw_deta_c = de2 * (q(4,i,j+1,k)-q(4,i,j-1,k))
   
      if (j == j_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_deta_m = de * (q(2,i,j,k) - q(2,i,j-1,k))
         dv_deta_m = de * (q(3,i,j,k) - q(3,i,j-1,k))
         dw_deta_m = de * (q(4,i,j,k) - q(4,i,j-1,k))
   
         ! downwind scheme
         du_deta_p = de2 * (  - one   * q(2,i,j+2,k) &
                            & + four  * q(2,i,j+1,k) &
                            & - three * q(2,i,j  ,k) )
   
         dv_deta_p = de2 * (  - one   * q(3,i,j+2,k) &
                            & + four  * q(3,i,j+1,k) &
                            & - three * q(3,i,j  ,k) )
   
         dw_deta_p = de2 * (  - one   * q(4,i,j+2,k) &
                            & + four  * q(4,i,j+1,k) &
                            & - three * q(4,i,j  ,k) )
   
         ! if j+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_deta_p = du_deta_p + (one - rsign(i,j+2,k)) * (- du_deta_p + de * (q(2,i,j+1,k)-q(2,i,j,k)))
         dv_deta_p = dv_deta_p + (one - rsign(i,j+2,k)) * (- dv_deta_p + de * (q(3,i,j+1,k)-q(3,i,j,k)))
         dw_deta_p = dw_deta_p + (one - rsign(i,j+2,k)) * (- dw_deta_p + de * (q(4,i,j+1,k)-q(4,i,j,k)))
   
   
      else if (j == j_myend) then
   
         ! upwind scheme
         du_deta_m = dc2 * (    one   * q(2,i,j-2,k) &
                            & - four  * q(2,i,j-1,k) &
                            & + three * q(2,i,j  ,k) )
   
         dv_deta_m = dc2 * (    one   * q(3,i,j-2,k) &
                            & - four  * q(3,i,j-1,k) &
                            & + three * q(3,i,j  ,k) )
   
         dw_deta_m = dc2 * (    one   * q(4,i,j-2,k) &
                            & - four  * q(4,i,j-1,k) &
                            & + three * q(4,i,j  ,k) )
   
         ! if j-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_deta_m = du_deta_m + (one - rsign(i,j-2,k)) * (- du_deta_m + de * (q(2,i,j,k)-q(2,i,j-1,k)))
         dv_deta_m = dv_deta_m + (one - rsign(i,j-2,k)) * (- dv_deta_m + de * (q(3,i,j,k)-q(3,i,j-1,k)))
         dw_deta_m = dw_deta_m + (one - rsign(i,j-2,k)) * (- dw_deta_m + de * (q(4,i,j,k)-q(4,i,j-1,k)))
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_deta_p = de * (q(2,i,j+1,k) - q(2,i,j,k))
         dv_deta_p = de * (q(3,i,j+1,k) - q(3,i,j,k))
         dw_deta_p = de * (q(4,i,j+1,k) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_deta_m = dc2 * (    one   * q(2,i,j-2,k) &
                            & - four  * q(2,i,j-1,k) &
                            & + three * q(2,i,j  ,k) )
   
         dv_deta_m = dc2 * (    one   * q(3,i,j-2,k) &
                            & - four  * q(3,i,j-1,k) &
                            & + three * q(3,i,j  ,k) )
   
         dw_deta_m = dc2 * (    one   * q(4,i,j-2,k) &
                            & - four  * q(4,i,j-1,k) &
                            & + three * q(4,i,j  ,k) )
   
         ! if j-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_deta_m = du_deta_m + (one - rsign(i,j-2,k)) * (- du_deta_m + de * (q(2,i,j,k)-q(2,i,j-1,k)))
         dv_deta_m = dv_deta_m + (one - rsign(i,j-2,k)) * (- dv_deta_m + de * (q(3,i,j,k)-q(3,i,j-1,k)))
         dw_deta_m = dw_deta_m + (one - rsign(i,j-2,k)) * (- dw_deta_m + de * (q(4,i,j,k)-q(4,i,j-1,k)))
   
   
         ! downwind derivatives
         du_deta_p = de2 * (  - one   * q(2,i,j+2,k) &
                            & + four  * q(2,i,j+1,k) &
                            & - three * q(2,i,j  ,k) )
   
         dv_deta_p = de2 * (  - one   * q(3,i,j+2,k) &
                            & + four  * q(3,i,j+1,k) &
                            & - three * q(3,i,j  ,k) )
   
         dw_deta_p = de2 * (  - one   * q(4,i,j+2,k) &
                            & + four  * q(4,i,j+1,k) &
                            & - three * q(4,i,j  ,k) )
   
         ! if j+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_deta_p = du_deta_p + (one - rsign(i,j+2,k)) * (- du_deta_p + de * (q(2,i,j+1,k)-q(2,i,j,k)))
         dv_deta_p = dv_deta_p + (one - rsign(i,j+2,k)) * (- dv_deta_p + de * (q(3,i,j+1,k)-q(3,i,j,k)))
         dw_deta_p = dw_deta_p + (one - rsign(i,j+2,k)) * (- dw_deta_p + de * (q(4,i,j+1,k)-q(4,i,j,k)))
   
      end if
     
      ! - - - - - - - - - - - - - - 
      ! k - direction
      ! - - - - - - - - - - - - - - 
   
      du_dzet_c = de2 * (q(2,i,j,k+1)-q(2,i,j,k-1))
      dv_dzet_c = de2 * (q(3,i,j,k+1)-q(3,i,j,k-1))
      dw_dzet_c = de2 * (q(4,i,j,k+1)-q(4,i,j,k-1))
   
      if (k == k_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_dzet_m = dz * (q(2,i,j,k) - q(2,i,j,k-1))
         dv_dzet_m = dz * (q(3,i,j,k) - q(3,i,j,k-1))
         dw_dzet_m = dz * (q(4,i,j,k) - q(4,i,j,k-1))
   
         ! downwind
         du_dzet_p = dz2 * (  - one   * q(2,i,j,k+2) &
                            & + four  * q(2,i,j,k+1) &
                            & - three * q(2,i,j,k  ) )
   
         dv_dzet_p = dz2 * (  - one   * q(3,i,j,k+2) &
                            & + four  * q(3,i,j,k+1) &
                            & - three * q(3,i,j,k  ) )
   
         dw_dzet_p = dz2 * (  - one   * q(4,i,j,k+2) &
                            & + four  * q(4,i,j,k+1) &
                            & - three * q(4,i,j,k  ) )
   
         ! if k+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dzet_p = du_dzet_p + (one - rsign(i,j,k+2)) * (- du_dzet_p + dz * (q(2,i,j,k+1)-q(2,i,j,k)))
         dv_dzet_p = dv_dzet_p + (one - rsign(i,j,k+2)) * (- dv_dzet_p + dz * (q(3,i,j,k+1)-q(3,i,j,k)))
         dw_dzet_p = dw_dzet_p + (one - rsign(i,j,k+2)) * (- dw_dzet_p + dz * (q(4,i,j,k+1)-q(4,i,j,k)))
   
   
      else if (k == k_myend) then
   
         ! upwind scheme
         du_dzet_m = dz2 * (    one   * q(2,i,j,k-2) &
                            & - four  * q(2,i,j,k-1) &
                            & + three * q(2,i,j,k  ) )
   
         dv_dzet_m = dz2 * (    one   * q(3,i,j,k-2) &
                            & - four  * q(3,i,j,k-1) &
                            & + three * q(3,i,j,k  ) )
   
         dw_dzet_m = dz2 * (    one   * q(4,i,j,k-2) &
                            & - four  * q(4,i,j,k-1) &
                            & + three * q(4,i,j,k  ) )
   
         ! if k-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dzet_m = du_dzet_m + (one - rsign(i,j,k-2)) * (- du_dzet_m + dz * (q(2,i,j,k)-q(2,i,j,k-1)))
         dv_dzet_m = dv_dzet_m + (one - rsign(i,j,k-2)) * (- dv_dzet_m + dz * (q(3,i,j,k)-q(3,i,j,k-1)))
         dw_dzet_m = dw_dzet_m + (one - rsign(i,j,k-2)) * (- dw_dzet_m + dz * (q(4,i,j,k)-q(4,i,j,k-1)))
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_dzet_p = dz * (q(2,i,j,k+1) - q(2,i,j,k))
         dv_dzet_p = dz * (q(3,i,j,k+1) - q(3,i,j,k))
         dw_dzet_p = dz * (q(4,i,j,k+1) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_dzet_m = dz2 * (    one   * q(2,i,j,k-2) &
                            & - four  * q(2,i,j,k-1) &
                            & + three * q(2,i,j,k  ) )
   
         dv_dzet_m = dz2 * (    one   * q(3,i,j,k-2) &
                            & - four  * q(3,i,j,k-1) &
                            & + three * q(3,i,j,k  ) )
   
         dw_dzet_m = dz2 * (    one   * q(4,i,j,k-2) &
                            & - four  * q(4,i,j,k-1) &
                            & + three * q(4,i,j,k  ) )
   
         ! if k-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dzet_m = du_dzet_m + (one - rsign(i,j,k-2)) * (- du_dzet_m + dz * (q(2,i,j,k)-q(2,i,j,k-1)))
         dv_dzet_m = dv_dzet_m + (one - rsign(i,j,k-2)) * (- dv_dzet_m + dz * (q(3,i,j,k)-q(3,i,j,k-1)))
         dw_dzet_m = dw_dzet_m + (one - rsign(i,j,k-2)) * (- dw_dzet_m + dz * (q(4,i,j,k)-q(4,i,j,k-1)))
   
         ! downwind derivatives
         du_dzet_p = dz2 * (  - one   * q(2,i,j,k+2) &
                            & + four  * q(2,i,j,k+1) &
                            & - three * q(2,i,j,k  ) )
   
         dv_dzet_p = dz2 * (  - one   * q(3,i,j,k+2) &
                            & + four  * q(3,i,j,k+1) &
                            & - three * q(3,i,j,k  ) )
   
         dw_dzet_p = dz2 * (  - one   * q(4,i,j,k+2) &
                            & + four  * q(4,i,j,k+1) &
                            & - three * q(4,i,j,k  ) )
   
         ! if k+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dzet_p = du_dzet_p + (one - rsign(i,j,k+2)) * (- du_dzet_p + dz * (q(2,i,j,k+1)-q(2,i,j,k)))
         dv_dzet_p = dv_dzet_p + (one - rsign(i,j,k+2)) * (- dv_dzet_p + dz * (q(3,i,j,k+1)-q(3,i,j,k)))
         dw_dzet_p = dw_dzet_p + (one - rsign(i,j,k+2)) * (- dw_dzet_p + dz * (q(4,i,j,k+1)-q(4,i,j,k)))
   
      end if

      ! - - - - - - - - - - - - - - 
      ! GRADIENT COMPUTATION
      ! - - - - - - - - - - - - - - 
   
      ! - - - - - - - - - - - - - - 
      ! u - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂u/∂ξ
      velocity_curv_gradient(1,1) = rsign_im * du_dcsi_m + &
                                    rsign_ic * du_dcsi_c + &
                                    rsign_ip * du_dcsi_p
   
      ! ∂u/∂η
      velocity_curv_gradient(1,2) = rsign_jm * du_deta_m + &
                                    rsign_jc * du_deta_c + &
                                    rsign_jp * du_deta_p
   
      ! ∂u/∂ζ
      velocity_curv_gradient(1,3) = rsign_km * du_dzet_m + &
                                    rsign_kc * du_dzet_c + &
                                    rsign_kp * du_dzet_p
      ! - - - - - - - - - - - - - - 
      ! v - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂v/∂ξ
      velocity_curv_gradient(2,1) = rsign_im * dv_dcsi_m + &
                                    rsign_ic * dv_dcsi_c + &
                                    rsign_ip * dv_dcsi_p
   
      ! ∂v/∂η
      velocity_curv_gradient(2,2) = rsign_jm * dv_deta_m + &
                                    rsign_jc * dv_deta_c + &
                                    rsign_jp * dv_deta_p
   
      ! ∂v/∂ζ
      velocity_curv_gradient(2,3) = rsign_km * dv_dzet_m + &
                                    rsign_kc * dv_dzet_c + &
                                    rsign_kp * dv_dzet_p
      ! - - - - - - - - - - - - - - 
      ! w - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂w/∂ξ
      velocity_curv_gradient(3,1) = rsign_im * dw_dcsi_m + &
                                    rsign_ic * dw_dcsi_c + &
                                    rsign_ip * dw_dcsi_p
   
      ! ∂w/∂η
      velocity_curv_gradient(3,2) = rsign_jm * dw_deta_m + &
                                    rsign_jc * dw_deta_c + &
                                    rsign_jp * dw_deta_p
   
      ! ∂w/∂ζ
      velocity_curv_gradient(3,3) = rsign_km * dw_dzet_m + &
                                    rsign_kc * dw_dzet_c + &
                                    rsign_kp * dw_dzet_p
   
   end if ! (if exsign > 0)

   ! if there is a direction where the flag variable is zero, the gradient
   ! is not considered for lsqm (the gradient is set to zero by setting exsign to zero)

   velocity_curv_gradient = exsign * velocity_curv_gradient

end subroutine velocity_curv_gradient_tensor 

subroutine velocity_curv_gradient_tensor_near_layer(i, j, k, velocity_curv_gradient, exsign)
   
   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
   
   ! local variables

   ! rsign_im, rsign_ic and rsign_ip are flag variables to decide which stencil to
   ! use depending on whether the node is and their neighbors. 
   ! For the variables:   m = minus/upwind   (i-2,i-1,i or i-1,i stencils)
   !                      c = central        (i-1,i+1 stencil)
   !                      p = plus/downwind  (i,i+1,i+2 or i,i+1 stencils).
   !
   ! Only one of {rsign_im, rsign_ic,rsign_ip} variables can be 1.0 for a node 
   ! i,j,k at the same time

   real (kind = rdf) ::   rsign_im, rsign_ic, rsign_ip, &
                        & rsign_jm, rsign_jc, rsign_jp, &
                        & rsign_km, rsign_kc, rsign_kp

   ! velocity derivatives for the different used stencils

   real (kind = rdf) ::   du_dcsi_m, du_dcsi_c, du_dcsi_p, &
                        & du_deta_m, du_deta_c, du_deta_p, &
                        & du_dzet_m, du_dzet_c, du_dzet_p, &
                        & dv_dcsi_m, dv_dcsi_c, dv_dcsi_p, &
                        & dv_deta_m, dv_deta_c, dv_deta_p, &
                        & dv_dzet_m, dv_dzet_c, dv_dzet_p, &
                        & dw_dcsi_m, dw_dcsi_c, dw_dcsi_p, &
                        & dw_deta_m, dw_deta_c, dw_deta_p, &
                        & dw_dzet_m, dw_dzet_c, dw_dzet_p

   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   real (kind = rdf) :: aux_sum ! sum over all the stencil flag variables

   du_dcsi_m = zero; du_dcsi_c = zero; du_dcsi_p = zero;
   du_deta_m = zero; du_deta_c = zero; du_deta_p = zero;
   du_dzet_m = zero; du_dzet_c = zero; du_dzet_p = zero;
   dv_dcsi_m = zero; dv_dcsi_c = zero; dv_dcsi_p = zero;
   dv_deta_m = zero; dv_deta_c = zero; dv_deta_p = zero;
   dv_dzet_m = zero; dv_dzet_c = zero; dv_dzet_p = zero;
   dw_dcsi_m = zero; dw_dcsi_c = zero; dw_dcsi_p = zero;
   dw_deta_m = zero; dw_deta_c = zero; dw_deta_p = zero;
   dw_dzet_m = zero; dw_dzet_c = zero; dw_dzet_p = zero;

   ! Only one of {rsign_im, rsign_ic, rsign_ip} can be one for a node i,j,k
   ! at the same time

   ! rsign_im = 1 --> rsign_ic = rsign_ip = 0 --> upwind scheme is used
   ! rsign_ic = 1 --> rsign_im = rsign_ip = 0 --> centred scheme is used
   ! rsign_ip = 1 --> rsign_im = rsign_ic = 0 --> downwind scheme is used

   ! if i-1,i,i+1 = water --> rsign_ic = 1. Otherwise, rsign_ic = 0
   rsign_ic = rsign_aux(i-1,j,k)   * rsign_aux(i,j,k) * rsign_aux(i+1,j,k)
   ! if i = water, i-1 = water and rsign_ic = 0 --> rsign_im = 1. Otherwise, rsign_im = 0
   rsign_im = (one - rsign_ic) * rsign_aux(i,j,k) * rsign_aux(i-1,j,k)
   ! if i = water, i+1 = water and rsign_ic = 0 --> rsign_ip = 1. Otherwise, rsign_ip = 0
   rsign_ip = (one - rsign_ic) * rsign_aux(i,j,k) * rsign_aux(i+1,j,k)

   rsign_jc = rsign_aux(i,j-1,k)   * rsign_aux(i,j,k) * rsign_aux(i,j+1,k)
   rsign_jm = (one - rsign_jc) * rsign_aux(i,j,k) * rsign_aux(i,j-1,k)
   rsign_jp = (one - rsign_jc) * rsign_aux(i,j,k) * rsign_aux(i,j+1,k)

   rsign_kc = rsign_aux(i,j,k-1)   * rsign_aux(i,j,k) * rsign_aux(i,j,k+1)
   rsign_km = (one - rsign_kc) * rsign_aux(i,j,k) * rsign_aux(i,j,k-1)
   rsign_kp = (one - rsign_kc) * rsign_aux(i,j,k) * rsign_aux(i,j,k+1)

   ! if there is any direction where the flag variable is zero, the gradient
   ! shouldn't be considered for lsqm (the gradient is set to zero and the
   ! flag variable, exsign, too)

   aux_sum =  rsign_im + rsign_ic + rsign_ip &
            + rsign_jm + rsign_jc + rsign_jp &
            + rsign_km + rsign_kc + rsign_kp

   ! if aux_sum<3, exsign = 0

   exsign = one - (sign( one , 2.5_rdf - aux_sum) + one) / two

   if (exsign>0) then ! if exsign = 0, I don't compute any derivative
   
      ! - - - - - - - - - - - - - - 
      ! i - direction
      ! - - - - - - - - - - - - - - 
   
      du_dcsi_c = dc2 * (q(2,i+1,j,k)-q(2,i-1,j,k))
      dv_dcsi_c = dc2 * (q(3,i+1,j,k)-q(3,i-1,j,k))
      dw_dcsi_c = dc2 * (q(4,i+1,j,k)-q(4,i-1,j,k))
   
      if (i == i_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_dcsi_m = dc * (q(2,i,j,k) - q(2,i-1,j,k))
         dv_dcsi_m = dc * (q(3,i,j,k) - q(3,i-1,j,k))
         dw_dcsi_m = dc * (q(4,i,j,k) - q(4,i-1,j,k))
   
         ! downwind scheme
         du_dcsi_p = dc2 * (  - one   * q(2,i+2,j,k) &
                            & + four  * q(2,i+1,j,k) &
                            & - three * q(2,i  ,j,k) )
   
         dv_dcsi_p = dc2 * (  - one   * q(3,i+2,j,k) &
                            & + four  * q(3,i+1,j,k) &
                            & - three * q(3,i  ,j,k) )
   
         dw_dcsi_p = dc2 * (  - one   * q(4,i+2,j,k) &
                            & + four  * q(4,i+1,j,k) &
                            & - three * q(4,i  ,j,k) )
   
         ! if i+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dcsi_p = du_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- du_dcsi_p + dc * (q(2,i+1,j,k)-q(2,i,j,k)))
         dv_dcsi_p = dv_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- dv_dcsi_p + dc * (q(3,i+1,j,k)-q(3,i,j,k)))
         dw_dcsi_p = dw_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- dw_dcsi_p + dc * (q(4,i+1,j,k)-q(4,i,j,k)))
   
      else if (i == i_myend) then
   
         ! upwind scheme
         du_dcsi_m = dc2 * (    one   * q(2,i-2,j,k) &
                            & - four  * q(2,i-1,j,k) &
                            & + three * q(2,i  ,j,k) )
   
         dv_dcsi_m = dc2 * (    one   * q(3,i-2,j,k) &
                            & - four  * q(3,i-1,j,k) &
                            & + three * q(3,i  ,j,k) )
   
         dw_dcsi_m = dc2 * (    one   * q(4,i-2,j,k) &
                            & - four  * q(4,i-1,j,k) &
                            & + three * q(4,i  ,j,k) )
   
         ! if i-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dcsi_m = du_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- du_dcsi_m + dc * (q(2,i,j,k)-q(2,i-1,j,k)))
         dv_dcsi_m = dv_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- dv_dcsi_m + dc * (q(3,i,j,k)-q(3,i-1,j,k)))
         dw_dcsi_m = dw_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- dw_dcsi_m + dc * (q(4,i,j,k)-q(4,i-1,j,k)))
   
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_dcsi_p = dc * (q(2,i+1,j,k) - q(2,i,j,k))
         dv_dcsi_p = dc * (q(3,i+1,j,k) - q(3,i,j,k))
         dw_dcsi_p = dc * (q(4,i+1,j,k) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_dcsi_m = dc2 * (    one   * q(2,i-2,j,k) &
                            & - four  * q(2,i-1,j,k) &
                            & + three * q(2,i  ,j,k) )
   
         dv_dcsi_m = dc2 * (    one   * q(3,i-2,j,k) &
                            & - four  * q(3,i-1,j,k) &
                            & + three * q(3,i  ,j,k) )
   
         dw_dcsi_m = dc2 * (    one   * q(4,i-2,j,k) &
                            & - four  * q(4,i-1,j,k) &
                            & + three * q(4,i  ,j,k) )
   
         ! if i-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dcsi_m = du_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- du_dcsi_m + dc * (q(2,i,j,k)-q(2,i-1,j,k)))
         dv_dcsi_m = dv_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- dv_dcsi_m + dc * (q(3,i,j,k)-q(3,i-1,j,k)))
         dw_dcsi_m = dw_dcsi_m + (one - rsign_aux(i-2,j,k)) * (- dw_dcsi_m + dc * (q(4,i,j,k)-q(4,i-1,j,k)))
   
   
         ! downwind derivatives
         du_dcsi_p = dc2 * (  - one   * q(2,i+2,j,k) &
                            & + four  * q(2,i+1,j,k) &
                            & - three * q(2,i  ,j,k) )
   
         dv_dcsi_p = dc2 * (  - one   * q(3,i+2,j,k) &
                            & + four  * q(3,i+1,j,k) &
                            & - three * q(3,i  ,j,k) )
   
         dw_dcsi_p = dc2 * (  - one   * q(4,i+2,j,k) &
                            & + four  * q(4,i+1,j,k) &
                            & - three * q(4,i  ,j,k) )
   
         ! if i+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dcsi_p = du_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- du_dcsi_p + dc * (q(2,i+1,j,k)-q(2,i,j,k)))
         dv_dcsi_p = dv_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- dv_dcsi_p + dc * (q(3,i+1,j,k)-q(3,i,j,k)))
         dw_dcsi_p = dw_dcsi_p + (one - rsign_aux(i+2,j,k)) * (- dw_dcsi_p + dc * (q(4,i+1,j,k)-q(4,i,j,k)))
   
      end if
     
      ! - - - - - - - - - - - - - - 
      ! j - direction
      ! - - - - - - - - - - - - - - 
   
      du_deta_c = de2 * (q(2,i,j+1,k)-q(2,i,j-1,k))
      dv_deta_c = de2 * (q(3,i,j+1,k)-q(3,i,j-1,k))
      dw_deta_c = de2 * (q(4,i,j+1,k)-q(4,i,j-1,k))
   
      if (j == j_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_deta_m = de * (q(2,i,j,k) - q(2,i,j-1,k))
         dv_deta_m = de * (q(3,i,j,k) - q(3,i,j-1,k))
         dw_deta_m = de * (q(4,i,j,k) - q(4,i,j-1,k))
   
         ! downwind scheme
         du_deta_p = de2 * (  - one   * q(2,i,j+2,k) &
                            & + four  * q(2,i,j+1,k) &
                            & - three * q(2,i,j  ,k) )
   
         dv_deta_p = de2 * (  - one   * q(3,i,j+2,k) &
                            & + four  * q(3,i,j+1,k) &
                            & - three * q(3,i,j  ,k) )
   
         dw_deta_p = de2 * (  - one   * q(4,i,j+2,k) &
                            & + four  * q(4,i,j+1,k) &
                            & - three * q(4,i,j  ,k) )
   
         ! if j+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_deta_p = du_deta_p + (one - rsign_aux(i,j+2,k)) * (- du_deta_p + de * (q(2,i,j+1,k)-q(2,i,j,k)))
         dv_deta_p = dv_deta_p + (one - rsign_aux(i,j+2,k)) * (- dv_deta_p + de * (q(3,i,j+1,k)-q(3,i,j,k)))
         dw_deta_p = dw_deta_p + (one - rsign_aux(i,j+2,k)) * (- dw_deta_p + de * (q(4,i,j+1,k)-q(4,i,j,k)))
   
   
      else if (j == j_myend) then
   
         ! upwind scheme
         du_deta_m = dc2 * (    one   * q(2,i,j-2,k) &
                            & - four  * q(2,i,j-1,k) &
                            & + three * q(2,i,j  ,k) )
   
         dv_deta_m = dc2 * (    one   * q(3,i,j-2,k) &
                            & - four  * q(3,i,j-1,k) &
                            & + three * q(3,i,j  ,k) )
   
         dw_deta_m = dc2 * (    one   * q(4,i,j-2,k) &
                            & - four  * q(4,i,j-1,k) &
                            & + three * q(4,i,j  ,k) )
   
         ! if j-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_deta_m = du_deta_m + (one - rsign_aux(i,j-2,k)) * (- du_deta_m + de * (q(2,i,j,k)-q(2,i,j-1,k)))
         dv_deta_m = dv_deta_m + (one - rsign_aux(i,j-2,k)) * (- dv_deta_m + de * (q(3,i,j,k)-q(3,i,j-1,k)))
         dw_deta_m = dw_deta_m + (one - rsign_aux(i,j-2,k)) * (- dw_deta_m + de * (q(4,i,j,k)-q(4,i,j-1,k)))
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_deta_p = de * (q(2,i,j+1,k) - q(2,i,j,k))
         dv_deta_p = de * (q(3,i,j+1,k) - q(3,i,j,k))
         dw_deta_p = de * (q(4,i,j+1,k) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_deta_m = dc2 * (    one   * q(2,i,j-2,k) &
                            & - four  * q(2,i,j-1,k) &
                            & + three * q(2,i,j  ,k) )
   
         dv_deta_m = dc2 * (    one   * q(3,i,j-2,k) &
                            & - four  * q(3,i,j-1,k) &
                            & + three * q(3,i,j  ,k) )
   
         dw_deta_m = dc2 * (    one   * q(4,i,j-2,k) &
                            & - four  * q(4,i,j-1,k) &
                            & + three * q(4,i,j  ,k) )
   
         ! if j-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_deta_m = du_deta_m + (one - rsign_aux(i,j-2,k)) * (- du_deta_m + de * (q(2,i,j,k)-q(2,i,j-1,k)))
         dv_deta_m = dv_deta_m + (one - rsign_aux(i,j-2,k)) * (- dv_deta_m + de * (q(3,i,j,k)-q(3,i,j-1,k)))
         dw_deta_m = dw_deta_m + (one - rsign_aux(i,j-2,k)) * (- dw_deta_m + de * (q(4,i,j,k)-q(4,i,j-1,k)))
   
   
         ! downwind derivatives
         du_deta_p = de2 * (  - one   * q(2,i,j+2,k) &
                            & + four  * q(2,i,j+1,k) &
                            & - three * q(2,i,j  ,k) )
   
         dv_deta_p = de2 * (  - one   * q(3,i,j+2,k) &
                            & + four  * q(3,i,j+1,k) &
                            & - three * q(3,i,j  ,k) )
   
         dw_deta_p = de2 * (  - one   * q(4,i,j+2,k) &
                            & + four  * q(4,i,j+1,k) &
                            & - three * q(4,i,j  ,k) )
   
         ! if j+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_deta_p = du_deta_p + (one - rsign_aux(i,j+2,k)) * (- du_deta_p + de * (q(2,i,j+1,k)-q(2,i,j,k)))
         dv_deta_p = dv_deta_p + (one - rsign_aux(i,j+2,k)) * (- dv_deta_p + de * (q(3,i,j+1,k)-q(3,i,j,k)))
         dw_deta_p = dw_deta_p + (one - rsign_aux(i,j+2,k)) * (- dw_deta_p + de * (q(4,i,j+1,k)-q(4,i,j,k)))
   
      end if
     
      ! - - - - - - - - - - - - - - 
      ! k - direction
      ! - - - - - - - - - - - - - - 
   
      du_dzet_c = de2 * (q(2,i,j,k+1)-q(2,i,j,k-1))
      dv_dzet_c = de2 * (q(3,i,j,k+1)-q(3,i,j,k-1))
      dw_dzet_c = de2 * (q(4,i,j,k+1)-q(4,i,j,k-1))
   
      if (k == k_mysta)  then
   
         ! if the node is next to the border, the upwind
         ! scheme is just first - order accurate
   
         du_dzet_m = dz * (q(2,i,j,k) - q(2,i,j,k-1))
         dv_dzet_m = dz * (q(3,i,j,k) - q(3,i,j,k-1))
         dw_dzet_m = dz * (q(4,i,j,k) - q(4,i,j,k-1))
   
         ! downwind
         du_dzet_p = dz2 * (  - one   * q(2,i,j,k+2) &
                            & + four  * q(2,i,j,k+1) &
                            & - three * q(2,i,j,k  ) )
   
         dv_dzet_p = dz2 * (  - one   * q(3,i,j,k+2) &
                            & + four  * q(3,i,j,k+1) &
                            & - three * q(3,i,j,k  ) )
   
         dw_dzet_p = dz2 * (  - one   * q(4,i,j,k+2) &
                            & + four  * q(4,i,j,k+1) &
                            & - three * q(4,i,j,k  ) )
   
         ! if k+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dzet_p = du_dzet_p + (one - rsign_aux(i,j,k+2)) * (- du_dzet_p + dz * (q(2,i,j,k+1)-q(2,i,j,k)))
         dv_dzet_p = dv_dzet_p + (one - rsign_aux(i,j,k+2)) * (- dv_dzet_p + dz * (q(3,i,j,k+1)-q(3,i,j,k)))
         dw_dzet_p = dw_dzet_p + (one - rsign_aux(i,j,k+2)) * (- dw_dzet_p + dz * (q(4,i,j,k+1)-q(4,i,j,k)))
   
   
      else if (k == k_myend) then
   
         ! upwind scheme
         du_dzet_m = dz2 * (    one   * q(2,i,j,k-2) &
                            & - four  * q(2,i,j,k-1) &
                            & + three * q(2,i,j,k  ) )
   
         dv_dzet_m = dz2 * (    one   * q(3,i,j,k-2) &
                            & - four  * q(3,i,j,k-1) &
                            & + three * q(3,i,j,k  ) )
   
         dw_dzet_m = dz2 * (    one   * q(4,i,j,k-2) &
                            & - four  * q(4,i,j,k-1) &
                            & + three * q(4,i,j,k  ) )
   
         ! if k-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dzet_m = du_dzet_m + (one - rsign_aux(i,j,k-2)) * (- du_dzet_m + dz * (q(2,i,j,k)-q(2,i,j,k-1)))
         dv_dzet_m = dv_dzet_m + (one - rsign_aux(i,j,k-2)) * (- dv_dzet_m + dz * (q(3,i,j,k)-q(3,i,j,k-1)))
         dw_dzet_m = dw_dzet_m + (one - rsign_aux(i,j,k-2)) * (- dw_dzet_m + dz * (q(4,i,j,k)-q(4,i,j,k-1)))
   
         ! if the node is next to the border, the downwind
         ! scheme is just first - order accurate
   
         du_dzet_p = dz * (q(2,i,j,k+1) - q(2,i,j,k))
         dv_dzet_p = dz * (q(3,i,j,k+1) - q(3,i,j,k))
         dw_dzet_p = dz * (q(4,i,j,k+1) - q(4,i,j,k))
   
      else
   
         ! upwind derivatives
         du_dzet_m = dz2 * (    one   * q(2,i,j,k-2) &
                            & - four  * q(2,i,j,k-1) &
                            & + three * q(2,i,j,k  ) )
   
         dv_dzet_m = dz2 * (    one   * q(3,i,j,k-2) &
                            & - four  * q(3,i,j,k-1) &
                            & + three * q(3,i,j,k  ) )
   
         dw_dzet_m = dz2 * (    one   * q(4,i,j,k-2) &
                            & - four  * q(4,i,j,k-1) &
                            & + three * q(4,i,j,k  ) )
   
         ! if k-2 node is within the air-phase, the upwind scheme is 
         ! replaced by a first-order one
         
         du_dzet_m = du_dzet_m + (one - rsign_aux(i,j,k-2)) * (- du_dzet_m + dz * (q(2,i,j,k)-q(2,i,j,k-1)))
         dv_dzet_m = dv_dzet_m + (one - rsign_aux(i,j,k-2)) * (- dv_dzet_m + dz * (q(3,i,j,k)-q(3,i,j,k-1)))
         dw_dzet_m = dw_dzet_m + (one - rsign_aux(i,j,k-2)) * (- dw_dzet_m + dz * (q(4,i,j,k)-q(4,i,j,k-1)))
   
         ! downwind derivatives
         du_dzet_p = dz2 * (  - one   * q(2,i,j,k+2) &
                            & + four  * q(2,i,j,k+1) &
                            & - three * q(2,i,j,k  ) )
   
         dv_dzet_p = dz2 * (  - one   * q(3,i,j,k+2) &
                            & + four  * q(3,i,j,k+1) &
                            & - three * q(3,i,j,k  ) )
   
         dw_dzet_p = dz2 * (  - one   * q(4,i,j,k+2) &
                            & + four  * q(4,i,j,k+1) &
                            & - three * q(4,i,j,k  ) )
   
         ! if k+2 node is within the air-phase, the downdind scheme is 
         ! replaced by a first-order one
   
         du_dzet_p = du_dzet_p + (one - rsign_aux(i,j,k+2)) * (- du_dzet_p + dz * (q(2,i,j,k+1)-q(2,i,j,k)))
         dv_dzet_p = dv_dzet_p + (one - rsign_aux(i,j,k+2)) * (- dv_dzet_p + dz * (q(3,i,j,k+1)-q(3,i,j,k)))
         dw_dzet_p = dw_dzet_p + (one - rsign_aux(i,j,k+2)) * (- dw_dzet_p + dz * (q(4,i,j,k+1)-q(4,i,j,k)))
   
      end if

      ! - - - - - - - - - - - - - - 
      ! GRADIENT COMPUTATION
      ! - - - - - - - - - - - - - - 
   
      ! - - - - - - - - - - - - - - 
      ! u - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂u/∂ξ
      velocity_curv_gradient(1,1) = rsign_im * du_dcsi_m + &
                                    rsign_ic * du_dcsi_c + &
                                    rsign_ip * du_dcsi_p
   
      ! ∂u/∂η
      velocity_curv_gradient(1,2) = rsign_jm * du_deta_m + &
                                    rsign_jc * du_deta_c + &
                                    rsign_jp * du_deta_p
   
      ! ∂u/∂ζ
      velocity_curv_gradient(1,3) = rsign_km * du_dzet_m + &
                                    rsign_kc * du_dzet_c + &
                                    rsign_kp * du_dzet_p
      ! - - - - - - - - - - - - - - 
      ! v - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂v/∂ξ
      velocity_curv_gradient(2,1) = rsign_im * dv_dcsi_m + &
                                    rsign_ic * dv_dcsi_c + &
                                    rsign_ip * dv_dcsi_p
   
      ! ∂v/∂η
      velocity_curv_gradient(2,2) = rsign_jm * dv_deta_m + &
                                    rsign_jc * dv_deta_c + &
                                    rsign_jp * dv_deta_p
   
      ! ∂v/∂ζ
      velocity_curv_gradient(2,3) = rsign_km * dv_dzet_m + &
                                    rsign_kc * dv_dzet_c + &
                                    rsign_kp * dv_dzet_p
      ! - - - - - - - - - - - - - - 
      ! w - derivatives
      ! - - - - - - - - - - - - - - 
   
      ! ∂w/∂ξ
      velocity_curv_gradient(3,1) = rsign_im * dw_dcsi_m + &
                                    rsign_ic * dw_dcsi_c + &
                                    rsign_ip * dw_dcsi_p
   
      ! ∂w/∂η
      velocity_curv_gradient(3,2) = rsign_jm * dw_deta_m + &
                                    rsign_jc * dw_deta_c + &
                                    rsign_jp * dw_deta_p
   
      ! ∂w/∂ζ
      velocity_curv_gradient(3,3) = rsign_km * dw_dzet_m + &
                                    rsign_kc * dw_dzet_c + &
                                    rsign_kp * dw_dzet_p
   
   end if ! (if exsign > 0)

   ! if there is a direction where the flag variable is zero, the gradient
   ! is not considered for lsqm (the gradient is set to zero by setting exsign to zero)

   velocity_curv_gradient = exsign * velocity_curv_gradient

end subroutine velocity_curv_gradient_tensor_near_layer


subroutine velocity_gradient_extrapolation_free_surface_lsqm(a_coeff_vector, alpha_vec, & 
                                                              velocity_gradient, du_dx_fs_lsqm)

   implicit none

   real (kind = rdf), dimension(0:3,1:9), intent(in) :: a_coeff_vector ! in matrix form for an 
                                                                       ! easier handling 

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! xneasrest
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: velocity_gradient ! velocity_gradient at
                                                                          ! inearest, jnearest, knearest

   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: du_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface

   ! local variables

   integer :: igfs, jgfs
   real (kind = rdf), dimension(1:3,1:3) :: fmatrix ! f_ij matrix for linear function approximation 

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                     __                                   __
   !                    |  a11 a21 a31 a12 a22 a32 a13 a23 a33  |
   !                    |  b11 b21 b31 b12 b22 b32 b13 b23 b33  |
   !  a_coeff_vector =  |  c11 c21 c31 c12 c22 c32 c13 c23 c33  |
   !                    |  d11 d21 d31 d12 d22 d32 d13 d23 d33  |
   !                     --                                   --
   !
   !
   ! aij is the coefficient for extrapolating ∂u_i/∂x_j
   !
   !
   ! du_dx_fs_lsqm : cartesian velocity gradient at the free-surface computed
   ! using the coefficients obtained by lsqm.
   ! 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! f_ij = a_ij + b_ij * α + c_ij * β + d_ij * γ

   fmatrix(1,1) = a_coeff_vector(0,1) + a_coeff_vector(1,1)*alpha_vec(1) &
                                    & + a_coeff_vector(2,1)*alpha_vec(2) &
                                    & + a_coeff_vector(3,1)*alpha_vec(3)

   fmatrix(2,1) = a_coeff_vector(0,2) + a_coeff_vector(1,2)*alpha_vec(1) &
                                    & + a_coeff_vector(2,2)*alpha_vec(2) &
                                    & + a_coeff_vector(3,2)*alpha_vec(3)

   fmatrix(3,1) = a_coeff_vector(0,3) + a_coeff_vector(1,3)*alpha_vec(1) &
                                    & + a_coeff_vector(2,3)*alpha_vec(2) &
                                    & + a_coeff_vector(3,3)*alpha_vec(3)

   fmatrix(1,2) = a_coeff_vector(0,4) + a_coeff_vector(1,4)*alpha_vec(1) &
                                    & + a_coeff_vector(2,4)*alpha_vec(2) &
                                    & + a_coeff_vector(3,4)*alpha_vec(3)

   fmatrix(2,2) = a_coeff_vector(0,5) + a_coeff_vector(1,5)*alpha_vec(1) &
                                    & + a_coeff_vector(2,5)*alpha_vec(2) &
                                    & + a_coeff_vector(3,5)*alpha_vec(3)

   fmatrix(3,2) = a_coeff_vector(0,6) + a_coeff_vector(1,6)*alpha_vec(1) &
                                    & + a_coeff_vector(2,6)*alpha_vec(2) &
                                    & + a_coeff_vector(3,6)*alpha_vec(3)

   fmatrix(1,3) = a_coeff_vector(0,7) + a_coeff_vector(1,7)*alpha_vec(1) &
                                    & + a_coeff_vector(2,7)*alpha_vec(2) &
                                    & + a_coeff_vector(3,7)*alpha_vec(3)

   fmatrix(2,3) = a_coeff_vector(0,8) + a_coeff_vector(1,8)*alpha_vec(1) &
                                    & + a_coeff_vector(2,8)*alpha_vec(2) &
                                    & + a_coeff_vector(3,8)*alpha_vec(3)

   fmatrix(3,3) = a_coeff_vector(0,9) + a_coeff_vector(1,9)*alpha_vec(1) &
                                    & + a_coeff_vector(2,9)*alpha_vec(2) &
                                    & + a_coeff_vector(3,9)*alpha_vec(3)


   ! xs : free-surface location
   !
   ! ∂u_i/∂x^j (xs) = ∂u_i/∂x_j (xs + α) - f_ij
   ! Where ∂u_i/∂x_j (xs + α) is the velocity gradient at the nearest water-phase node
   ! with a valid computed gradient (using the rules defined in the velocity_curv_gradient_tensor
   ! subroutine)

   do jgfs = 1,3 ! j-index for the velocity gradient at the free-surface loop
      do igfs= 1,3 ! i-index for for the velocity gradient the free-surface loop

            du_dx_fs_lsqm(igfs,jgfs) = velocity_gradient(igfs,jgfs) - fmatrix(igfs,jgfs)

      end do
   end do

end subroutine velocity_gradient_extrapolation_free_surface_lsqm




subroutine velocity_extrapolation_free_surface_lsqm ( alpha_vec, i, j, k, du_dx_fs_lsqm, u_fs_lsqm)

   implicit none

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! x(i,j,k) (xs and x as vectors)
   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface
   
   ! I declared this variable as inout to see it and initialise it in the main code before calling                                                                         
   real (kind = rdf), dimension(1:3), intent(inout) :: u_fs_lsqm ! velocity at the free-surface

   ! local variables

   integer :: iufs, tsum

   do iufs = 1,3 ! index for velocity free-surface computation (ufs). iufs = 1 -> u, iufs = 2 -> v... 
      
      ! u_fs_lsqm(iufs) = q(iufs+1) because q(2,i,j,k) = u(i,j,k), q(3,i,j,k) = v(i,j,k) ...         

      u_fs_lsqm(iufs) = q(iufs+1, i, j, k) 

      do tsum = 1,3 ! summation index from Taylor expansion: u(xs) = u(xs+α) - ∂u/∂x*α - ∂u/∂y*β - ∂u/∂z*γ
         
         u_fs_lsqm(iufs) = u_fs_lsqm(iufs) - du_dx_fs_lsqm(iufs,tsum)*alpha_vec(tsum)

      end do
   end do

end subroutine velocity_extrapolation_free_surface_lsqm




subroutine ghost_nodes_velocity_extrapolation(i,j,k, xs, ys, zs, u_fs_lsqm, du_dx_fs_lsqm)

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(in) :: u_fs_lsqm
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm
   real (kind = rdf), intent(in) :: xs, ys, zs ! free-surface position closest to (i,j,k) node

   ! local variables
   real (kind = rdf), dimension(1:3) :: alpha_local
   real (kind = rdf) :: uextp, vextp, wextp ! extrapolated velocity using Taylor expansion
   real (kind = rdf) :: rdiff_norm, exsign  

   integer :: ii, jj, kk, extp_swept

   ! extp_swept = # of nodes around i,j,k to be candidates for extrapolation
   extp_swept = 3 ! TO DO: it could be set from control.dat

   ! max(k_mysta-1,k-extp_swept),min(k_myend+1,k+extp_swept) means that it even extrapolates
   ! the velocity to the air nodes at the boundaries of the domain. This is gonna be useful to
   ! compute velocity gradients in the air phase for the normal boundary condition.

   !do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   !do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   !do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

   do kk = max( kl , k-extp_swept ) , min( ku , k+extp_swept ) 
   do jj = max( jl , j-extp_swept ) , min( ju , j+extp_swept )
   do ii = max( il , i-extp_swept ) , min( iu , i+extp_swept )

      ! if the (ii,jj,kk) is in the air-phase, it's extrapolated using
      ! a Taylor expansion
      
      if ( rsign(ii,jj,kk) < one_half ) then ! (air-phase)

         ! position vector from the free-surface to the extrapolated node 

         alpha_local = zero

         alpha_local(1) = x(ii,jj,kk) - xs
         alpha_local(2) = y(ii,jj,kk) - ys
         alpha_local(3) = z(ii,jj,kk) - zs

         ! distance between the extrapolation candidate (ii,jj,kk) node and the
         ! free - surface
         rdiff_norm = norm2( alpha_local )

         ! exsign: flag variable
         !
         ! exsign = 1 ,if current distance between extrapolation candidate and free 
         ! surface is less than a previous stored one (previous calls of the 
         ! ghost_nodes_velocity_extrapolation subroutine).
         !
         ! exsign = 0 otherwise

         exsign =  (sign(one , least_dis_extp(ii,jj,kk) - rdiff_norm) + one)/two
         
         ! if exsign = 1, the new least distance is stored in least_dis_extp(ii,jj,kk)
         least_dis_extp(ii,jj,kk) =            least_dis_extp(ii,jj,kk) + & 
                                    exsign * (-least_dis_extp(ii,jj,kk) + rdiff_norm)

         ! extrapolated velocities using Taylor expansion
         uextp = zero; vextp = zero; wextp = zero;

         uextp = u_fs_lsqm(1) +  du_dx_fs_lsqm(1,1)*alpha_local(1) + & 
                                 du_dx_fs_lsqm(1,2)*alpha_local(2) + & 
                                 du_dx_fs_lsqm(1,3)*alpha_local(3)

         vextp = u_fs_lsqm(2) +  du_dx_fs_lsqm(2,1)*alpha_local(1) + & 
                                 du_dx_fs_lsqm(2,2)*alpha_local(2) + & 
                                 du_dx_fs_lsqm(2,3)*alpha_local(3)

         wextp = u_fs_lsqm(3) +  du_dx_fs_lsqm(3,1)*alpha_local(1) + & 
                                 du_dx_fs_lsqm(3,2)*alpha_local(2) + & 
                                 du_dx_fs_lsqm(3,3)*alpha_local(3)

         ! if exsign = 1 (least distanced node), the velocity at ii,jj,kk
         ! is replaced by the extrapolated one

         if ( exsign > one_half ) then
      
            q(2,ii,jj,kk) = q(2,ii,jj,kk) + exsign * (-q(2,ii,jj,kk) + uextp)
            q(3,ii,jj,kk) = q(3,ii,jj,kk) + exsign * (-q(3,ii,jj,kk) + vextp)
            q(4,ii,jj,kk) = q(4,ii,jj,kk) + exsign * (-q(4,ii,jj,kk) + wextp)

         end if
      
      end if

   end do
   end do
   end do

end subroutine ghost_nodes_velocity_extrapolation


subroutine error_gfm(i, j, k, nvec, tvec, svec, du_dx_fs_lsqm)

   implicit none

   integer, intent(in) :: i,j,k 
   real (kind = rdf), dimension(1:3), intent(in) :: nvec, tvec, svec  ! vector system
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm ! velocity gradient at the
                                                                      ! free-surface

   ! local variables
   integer :: isum, jsum ! indexes for Einstein sumation

   error_tdir(i,j,k) = zero
   error_sdir(i,j,k) = zero
   error_ndir(i,j,k) = zero
   
   do isum = 1,3
      do jsum = 1,3

         error_tdir(i,j,k) = error_tdir(i,j,k) + (  du_dx_fs_lsqm(isum,jsum)     &
                                                  + du_dx_fs_lsqm(jsum,isum) ) * &
                                                    nvec(isum)*tvec(jsum)   

         error_sdir(i,j,k) = error_sdir(i,j,k) + (  du_dx_fs_lsqm(isum,jsum)     &
                                                  + du_dx_fs_lsqm(jsum,isum) ) * &
                                                    nvec(isum)*svec(jsum)   

         error_ndir(i,j,k) = error_sdir(i,j,k) + (  du_dx_fs_lsqm(isum,jsum)     &
                                                  + du_dx_fs_lsqm(jsum,isum) ) * &
                                                    nvec(isum)*nvec(jsum)   

      end do
   end do


end subroutine error_gfm



function PressureLinearFlux_NDBC(  pL       , pC       , pR        , &
                                   JL       , JC       , JR        , &
                                   metricsL , metricsC , metricsR  , &
                                   dcsi                            , & 
                                   phiL     , phiC     , phiR      , & 
                                   GradPhiL , GradPhiC , GradPhiR  , &  
                                   GradVelL , GradVelC , GradVelR  , &
                                   PrintResults                      &
                                 )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! This function computes the flux related to the pressure in the linear flux
   ! term of Navier-Stokes equations for nodes next to free surface. It uses the 
   ! normal dynamic boundary condition to compute the pressure at the free surface.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   !                                                          these derivatives, may need to be adjusted  
   !                 _                   _     _        _     depending on the phase of the L,R sides.                                             
   !                | ∂/∂ξ (p/J * ∂ξ/∂x)  |   | ∂f(1)/∂ξ |    To compute them we use the "irregular stars"
   ! PressureFlux = | ∂/∂ξ (p/J * ∂ξ/∂y)  | = | ∂f(2)/∂ξ | -> algorithm described in "A Computer Study of                                
   !                | ∂/∂ξ (p/J * ∂ξ/∂z)  |   | ∂f(3)/∂ξ |    Finite-Amplitude Water Waves" by Chan &  
   !                •-                   -•   •-        -•    Street (JCP, 1970).                    
   !
   !
   !                    ϕ = 0 
   !                     /  
   !        air         /      water
   !       (ϕ<0)       /       (ϕ>0)
   !                  /
   !     ----o-------x--------o------------------o---  --> ξ
   !        i-1     /fs       i                 i+1
   !        (L)    /         (C)                (R)
   !         |    /           |                  |
   !         |   /            |                  | 
   !         |  /             |                  |
   !         | |              |                  |
   !         | |<--- ΔξL  --->|<-----  ΔξR ----->|
   !         
   !         
   !  If it is necessary to modify the stencil, the general expression for the derivative becomes
   !  ( Ferziger, Computational Methods for Fluid Dynamics, 4th edition, section 3.4 )  
   !
   !  ∂f/∂ξ = a * fL + b * fC + c * fR 
   !  
   !  Where 
   !  
   !  a = (       - ΔξR^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !  b = ( ΔξR^2 - ΔξL^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !  c = (         ΔξL^2 ) / (ΔξR * ΔξL * ( ΔξL + ΔξR ) )
   !
   ! As fL or fR may lie on the free surface, the pressure needs to be computed using the normal 
   ! dynamic boundary condition (NDBC) in non-dimensional form
   ! 
   !
   !
   !         2      ∂ui            |    
   ! pfs  + ---- * ----- * ni * nj |    = 0 
   !         Re     ∂xj            |fs
   !
   !
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Input variables
   real ( kind = rdf ), intent(in) :: pL, pC, pR ! pressure
   real ( kind = rdf ), intent(in) :: JL, JC, JR ! Jacobian
   real ( kind = rdf ), dimension(3)   , intent(in) :: metricsL, metricsC, metricsR ! ∂ξ/∂xj
   real ( kind = rdf ), intent(in) :: dcsi ! Δξ
   real ( kind = rdf ), intent(in) :: phiL , phiC, phiR ! ϕ
   real ( kind = rdf ), dimension(3)   , intent(in) :: GradPhiL , GradPhiC , GradPhiR ! ∇ϕ
   real ( kind = rdf ), dimension(3,3) , intent(in) :: GradVelL , GradVelC , GradVelR ! ∇u
   logical, optional :: PrintResults

   ! Local variables

   integer :: i,j
   real ( kind = rdf ) :: dcsi_dx_L ,  dcsi_dy_L , dcsi_dz_L ! ∂ξ/∂xj_L 
   real ( kind = rdf ) :: dcsi_dx_C ,  dcsi_dy_C , dcsi_dz_C ! ∂ξ/∂xj_C
   real ( kind = rdf ) :: dcsi_dx_R ,  dcsi_dy_R , dcsi_dz_R ! ∂ξ/∂xj_R
   real ( kind = rdf ), dimension(3) :: fL, fC, fR ! p/J * ∂ξ/∂xj
   real ( kind = rdf ) :: InterpCoeff ! Interpolation coefficient (0 < InterpCoeff < 1)
   real ( kind = rdf ) :: dcsiL , dcsiR
   real ( kind = rdf ) :: tol

   ! Variables at the free surface

   real ( kind = rdf ) :: J_fs , dcsi_dx_fs , dcsi_dy_fs , dcsi_dz_fs ! Jacobian and metrics at the fs
   real ( kind = rdf ), dimension(3)   :: GradPhi_fs , nvec_fs ! ∇ϕ_fs & normal vector at the fs
   real ( kind = rdf ), dimension(3,3) :: GradVel_fs ! ∇u_fs
   real ( kind = rdf ) :: dui_dxj_ni_nj ! ∂ui/∂xj * ni * nj
   real ( kind = rdf ) :: pfs ! Pressure at the free surface obtained by NDBC
   real ( kind = rdf ) :: a , b, c ! Weights for the stencil of general derivative

   ! Output variable
   real ( kind = rdf ), dimension(3) :: PressureLinearFlux_NDBC ! ∂/∂ξ (p/J * ∂ξ/∂xj)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Set the tolerance

   tol = ten**(-eight)

   ! Metrics ∂ξ/∂x, ∂ξ/∂y & ∂ξ/∂z at L,C & R nodes

   dcsi_dx_L = metricsL(1)
   dcsi_dy_L = metricsL(2)
   dcsi_dz_L = metricsL(3)

   dcsi_dx_C = metricsC(1)
   dcsi_dy_C = metricsC(2)
   dcsi_dz_C = metricsC(3)

   dcsi_dx_R = metricsR(1)
   dcsi_dy_R = metricsR(2)
   dcsi_dz_R = metricsR(3)

   ! Initialise ΔξL and ΔξR as the mesh grid size
   dcsiL = dcsi
   dcsiR = dcsi

   ! Initialise fL, fc and fR without considering phase change
   
   ! p/J * ∂ξ/∂x
   fL(1) = pL / JL * dcsi_dx_L 
   fC(1) = pC / JC * dcsi_dx_C 
   fR(1) = pR / JR * dcsi_dx_R 

   ! p/J * ∂ξ/∂y
   fL(2) = pL / JL * dcsi_dy_L 
   fC(2) = pC / JC * dcsi_dy_C 
   fR(2) = pR / JR * dcsi_dy_R 

   ! p/J * ∂ξ/∂z
   fL(3) = pL / JL * dcsi_dz_L 
   fC(3) = pC / JC * dcsi_dz_C 
   fR(3) = pR / JR * dcsi_dz_R 

!   if( phiL < -eps_sims ) then
   if( phiL < zero ) then
   
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

      InterpCoeff = abs( phiC / ( phiC - phiL + eps_sims ) )
      !InterpCoeff = abs( phiC / abs( phiC - phiL ) )

      ! I modify dcsiL interpolating Δξ to the free surface
      dcsiL = InterpCoeff * dcsi 
   
      ! I compute J, ∂ξ/∂xj, ∇ϕ (normal vector), and ∇u at the free surface as an weighted combination
      ! of the values at nodes L and C
   
      J_fs        = InterpCoeff * JL        + ( one - InterpCoeff ) * JC       
      
      dcsi_dx_fs  = InterpCoeff * dcsi_dx_L + ( one - InterpCoeff ) * dcsi_dx_C  
      dcsi_dy_fs  = InterpCoeff * dcsi_dy_L + ( one - InterpCoeff ) * dcsi_dy_C  
      dcsi_dz_fs  = InterpCoeff * dcsi_dz_L + ( one - InterpCoeff ) * dcsi_dz_C  
      
      GradPhi_fs  = InterpCoeff * GradPhiL  + ( one - InterpCoeff ) * GradPhiC 
      GradVel_fs  = InterpCoeff * GradVelL  + ( one - InterpCoeff ) * GradVelC 

      ! unitary normal vector at the free surface n = -∇ϕ / |∇ϕ|
      !nvec_fs = -GradPhi_fs / ( norm2( GradPhi_fs ) + eps_sims )
      nvec_fs = -GradPhi_fs / ( norm2( GradPhi_fs ) )

      ! Computation of ∂ui/∂xj * ni * nj

      dui_dxj_ni_nj = zero

      do j = 1,3
         do i = 1,3
            
            dui_dxj_ni_nj = dui_dxj_ni_nj + GradVel_fs(i,j) * nvec_fs(i) * nvec_fs(j) 

         end do
      end do

      ! Normal Dynamic Boundary Condition p_fs = -2/Re * ( ∂ui / ∂xj * ni * nj )_fs
      !pfs = -two / ren * dui_dxj_ni_nj
      pfs = zero

      ! Update of fL using the values at the free surface
      fL(1) = pfs / J_fs * dcsi_dx_fs 
      fL(2) = pfs / J_fs * dcsi_dy_fs 
      fL(3) = pfs / J_fs * dcsi_dz_fs 

   end if


!   if ( phiR < -eps_sims ) then
   if ( phiR < zero ) then

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

      InterpCoeff = abs( phiC / ( phiC - phiR + eps_sims ) )
      !InterpCoeff = abs( phiC / ( phiC - phiR ) )

      ! I modify dcsiR interpolating Δξ to the free surface
      dcsiR = InterpCoeff * dcsi 
   
      ! I compute J, ∂ξ/∂xj, ∇ϕ (normal vector), and ∇u at the free surface as an weighted combination
      ! of the values at nodes L and C
   
      J_fs        = InterpCoeff * JR        + ( one - InterpCoeff ) * JC

      dcsi_dx_fs  = InterpCoeff * dcsi_dx_R + ( one - InterpCoeff ) * dcsi_dx_C  
      dcsi_dy_fs  = InterpCoeff * dcsi_dy_R + ( one - InterpCoeff ) * dcsi_dy_C  
      dcsi_dz_fs  = InterpCoeff * dcsi_dz_R + ( one - InterpCoeff ) * dcsi_dz_C  
      
      GradPhi_fs  = InterpCoeff * GradPhiR  + ( one - InterpCoeff ) * GradPhiC 
      GradVel_fs  = InterpCoeff * GradVelR  + ( one - InterpCoeff ) * GradVelC 

      ! unitary normal vector at the free surface n = -∇ϕ / |∇ϕ|
      nvec_fs = -GradPhi_fs / norm2( GradPhi_fs )

      ! Computation of ∂ui/∂xj * ni * nj

      dui_dxj_ni_nj = zero

      do j = 1,3
         do i = 1,3
            
            dui_dxj_ni_nj = dui_dxj_ni_nj + GradVel_fs(i,j) * nvec_fs(i) * nvec_fs(j) 

         end do
      end do

      ! Normal Dynamic Boundary Condition p_fs = -2/Re * ( ∂ui / ∂xj * ni * nj )_fs
      !pfs = -two / ren * dui_dxj_ni_nj
      pfs = zero

      ! Update of fR using the values at the free surface
      fR(1) = pfs / J_fs * dcsi_dx_fs 
      fR(2) = pfs / J_fs * dcsi_dy_fs 
      fR(3) = pfs / J_fs * dcsi_dz_fs 
   
   end if

   ! Stencil coefficients for ∂f/∂ξ = a * fL + b * fC + c * fR (Ferziger)
   
   ! I think there's a typo error with this coefficients, as a should go with
   ! dcsiL instead of dcsiR
   !a =   -dcsiR**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
   !b = (  dcsiR**2 - dcsiL**2 ) / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
   !c =    dcsiL**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  


   if ( tol * dcsiR > dcsiL ) then

      PressureLinearFlux_NDBC = (fR - fC) / dcsi

   else if ( tol * dcsiL > dcsiR ) then

      PressureLinearFlux_NDBC = (fC - fL) / dcsi

   else

      a =   -dcsiR**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
      b = (  dcsiR**2 - dcsiL**2 ) / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  
      c =    dcsiL**2              / ( dcsiR * dcsiL * ( dcsiR + dcsiL ) )  

      ! Finally, I compute PressureLinearFlux_NDBC = ∂f/∂ξ 
      PressureLinearFlux_NDBC = (a * fL) + (b * fC) + (c * fR)

   end if

end function PressureLinearFlux_NDBC



end subroutine ghost_fluid_extrapolation