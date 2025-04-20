subroutine ghost_fluid_extrapolation( il,iu          , &
                                      jl,ju          , &
                                      kl,ku          , &
                                      igp, jgp, kgp  , &
                                      dc, de, dz     , &
                                      q              , &
                                      xnut           , &
                                      csi            , &
                                      eta            , &
                                      zet            , &
                                      aj             , &
                                      phi            , & 
                                      phi_n          , & 
                                      phi_gradient   , & 
                                      x,y,z            &                                    
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
use global_lsm, only: mms, nms, sweep_lsqm, radius_lsqm, phi_outputiter, FrLSM, BigPhi, &
                      limit_ghost_velocities, zero_pressure_fs
use global_debug

implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! input-output arguments
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, intent(in) :: il,iu,jl,ju,kl,ku ! external nodes
integer, intent(in) :: igp, jgp, kgp     ! # of ghostpoint at local processor 

real (kind = rdf), intent(in) :: dc,de,dz ! dcsi, deta, dzet
real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) , intent(in)    :: csi , eta, zet ! metrics
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: aj ! |J|: jacobian
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: x,y,z ! grid coordinates

real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(inout) :: xnut
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: phi
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: phi_n
!real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: h , hn
!real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(inout) :: rsign
real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) , intent(in)    :: phi_gradient ! ∂phi/∂x_j

! In-out arguments. q vector is modified in its velocity components
real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: q

!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
!
!real (kind = rdf) , dimension(1:3,il:iu,jl:ju,kl:ku) :: nvecfield
!real (kind = rdf) , dimension(1:3,il:iu,jl:ju,kl:ku) :: tvecfield
!real (kind = rdf) , dimension(1:3,il:iu,jl:ju,kl:ku) :: svecfield
!
!real (kind = rdf) , dimension(il:iu,jl:ju,kl:ku) :: du_dx_tvec
!real (kind = rdf) , dimension(il:iu,jl:ju,kl:ku) :: du_dx_svec
!
!integer, intent(in) :: itc
!
!real (kind = rdf) :: f11,f21,f31, &
!                     f12,f22,f32, &
!                     f13,f23,f33
!
!character(len = 256) :: debugname

!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

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

integer :: ista, iend, jsta, jend, ksta, kend
integer :: im,ip,jm,jp,km,kp
real (kind = rdf) :: nWaterNeighbours
real (kind = rdf) :: dmax
logical :: PointWithinCell

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
real (kind = rdf), dimension(1:3,1:3) :: nj_ti, nj_si ! normal-tangential products for matrix
                                                      ! system
real (kind = rdf), dimension(1:3) :: alpha_vec ! position vectorsto the fs
real (kind = rdf), dimension(1:3) :: alpha_vec_extp ! position vectors from the extrapolation 
                                                    ! node to the fs
real (kind = rdf) :: xs, ys, zs ! free-surface location coordinates
integer :: m,n,l ! work arrays indexes

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix - lsqm system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real (kind = rdf), dimension(1:36) :: ttvec, tsvec     ! Tt and Ts vectors of the SVD system
real (kind = rdf), dimension(1:12) :: ttvec_p, tsvec_p ! Tt and Ts vectors of the SVD system
real (kind = rdf), dimension(0:3) :: alpha_local       ! flagged local position vector to the fs 
real (kind = rdf), dimension(1:3) :: rdiff             ! local position vector to the fs
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

real (kind = rdf), dimension(1:3) :: dp_dx_fs_lsqm ! pressure gradient at the
                                                   ! free-surface

real (kind = rdf), dimension(1:3) :: u_fs_lsqm ! velocity at the
                                               ! free-surface

! Big arrays will be allocatables (I allocate memory when I use them only)

real (kind = rdf), dimension(:,:), allocatable   :: t_matrix_system
real (kind = rdf), dimension(:)  , allocatable   :: a_coeff_vector
real (kind = rdf), dimension(:)  , allocatable   :: b_matrix_system

real (kind = rdf) :: weight
real (kind = rdf) :: xnut_fs


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

real (kind = rdf) :: max_vel_norm , max_vel_norm_global ! 
real (kind = rdf) :: max_vel_norm_neighbour ! 


! indexes for local velocity gradient computation
integer :: ivg, jvg, msum


real (kind = rdf) :: pfs !  pressure at the free-surface

logical :: BlankingFlag

! error arrays ( commented when we're not debbuging )
! real (kind = rdf), dimension(:,:,:), allocatable :: error_tdir, error_sdir

! ----------------------------------------------------------------------------------
! Singular - Value Decomposition (SVD) parameters
! ----------------------------------------------------------------------------------

!integer, parameter :: mmsu = 75
!integer, parameter :: nmsu = 36

integer, parameter :: mmsu = 57 ! 2 * 27 + 3 * continuity = 57
integer, parameter :: nmsu = 27

integer, parameter :: mmsp = 32
integer, parameter :: nmsp = 12

! SGELSD
integer :: info, rank
integer :: nlvlu, lworku, liworku


integer, parameter :: smlsiz=25
integer, parameter :: nrhs=1

real(kind = rdf), dimension(:), allocatable :: iworku , worku
real(kind = rdf), dimension(:), allocatable :: iworkp , workp

integer, parameter:: ldau=max(1,mmsu)

integer, parameter:: ldbu=max(1,max(mmsu,nmsu))

real(kind = rdf) :: ersvdm !!,sm,res,resf
real(kind = rdf) :: su(min(mmsu,nmsu)) 

real(kind = rdf), parameter:: rcond = 0.00001_rdf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!list variables

integer, save :: gfmnodes
integer, dimension(:,:), allocatable , save :: gfmnodes_list
integer :: g

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

external:: dgelsd , sgelsd , ilaenv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! These variables are defined as the dgelsd documentation suggests
! https://extras.csc.fi/math/nag/mark21/pdf/F08/f08kcf.pdf
! (don't modify if not necessary)

nlvlu = max(0,int( log(real(min(mmsu,nmsu))/real(smlsiz+1))/log(2.) )+1)
if(mmsu.ge.nmsu) then      
  lworku  = 12*nmsu+2*nmsu*smlsiz+8*nmsu*nlvlu+ nmsu*nrhs+(smlsiz+1)*(smlsiz+1)
else
  lworku  = 12*mmsu+2*mmsu*smlsiz+8*mmsu*nlvlu+ mmsu*nrhs+(smlsiz+1)*(smlsiz+1)
end if

liworku = max(1,3*min(mmsu,nmsu)*nlvlu+11*min(mmsu,nmsu))

allocate(worku(max(1,lworku)) , iworku(max(1,liworku)))


! ----------------------------------------------------------------------------------

i_mysta = il + igp
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp
j_myend = ju - jgp
k_myend = ku - kgp
         
! processes on the domain boundaries

if ( myback  == mpi_proc_null )  i_mysta = il + igp + 1
if ( myleft  == mpi_proc_null )  j_mysta = jl + jgp + 1
if ( mydown  == mpi_proc_null )  k_mysta = kl + kgp + 1

if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1
if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1
if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1

! Physical boundaries

ista = il ; jsta = jl ; ksta = kl 
iend = iu ; jend = ju ; kend = ku 

if ( myback  == mpi_proc_null )  ista = il + igp 
if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

if ( myfront == mpi_proc_null )  iend = iu - igp
if ( myright == mpi_proc_null )  jend = ju - jgp
if ( myup    == mpi_proc_null )  kend = ku - kgp


! memory allocation
allocate ( least_dis_extp ( il:iu , jl:ju , kl:ku ) )

! Error arrays initialisation
! error_tdir = zero
! error_sdir = zero

! Variables initialisation
least_dis_extp = 200.0_rdf
dc2            = one_half * dc
de2            = one_half * de
dz2            = one_half * dz


! q correction performs a linear extrapolation to the nodes that weren't resolve
! the last time step, but now are part of the water phase. It is only needed the
! first pseudo-time iteration as their flow variables are calculated from then on
if ( itc == 1 ) then

  !call p_correction_post_advection()

  ! Now q correction includes the pressure (that's why I commented p_correction_post_advection)
  call q_correction_post_advection()

  ! update ghost nodes q and xnut
  call rhs_exchng3_4d ( q )
  call rhs_exchng3_3d ( xnut )

end if

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

! Proc-local maximum velocity magnitude 
!max_vel_norm = zero

!do k = k_mysta , k_myend
!do j = j_mysta , j_myend
!do i = i_mysta , i_myend
!
!  if ( rsign(i,j,k) > one_half ) then
!    max_vel_norm = max( max_vel_norm , norm2( (/q(2,i,j,k) , q(3,i,j,k) , q(4,i,j,k)/) ) )
!  end if
!
!end do 
!end do
!end do

! Global maximum velocity magnitude
!call mpi_allreduce( max_vel_norm , max_vel_norm_global , 1 , MPI_REAL_TYPE , &
!                    mpi_max , mpi_comm_world , ierr )

! loop over the whole local domain looking for nodes to be interpolated
if ( itc == 1 ) call get_gfm_list()

! Loop over the ghost fluid method nodes
do g = 1, gfmnodes

  i = gfmnodes_list(g,1)
  j = gfmnodes_list(g,2)
  k = gfmnodes_list(g,3)

  ! Allocating matrix arrays for SVD solver

  allocate(  t_matrix_system ( 1:mmsu , 1:nmsu  ) , &
             a_coeff_vector  ( 1:nmsu )           , &
             b_matrix_system ( 1:max(mmsu,nmsu) )   &
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

  csi_work     (1:3,-1:1,-1:1,-1:1) = csi(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
  eta_work     (1:3,-1:1,-1:1,-1:1) = eta(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
  zet_work     (1:3,-1:1,-1:1,-1:1) = zet(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
  aj_work          (-1:1,-1:1,-1:1) = aj(i-1:i+1,j-1:j+1,k-1:k+1)
  phi_work         (-1:1,-1:1,-1:1) = phi(i-1:i+1,j-1:j+1,k-1:k+1)
  phi_grad_work(1:3,-1:1,-1:1,-1:1) = phi_gradient(1:3,i-1:i+1,j-1:j+1,k-1:k+1)
      
  ! nj_ti tensor = nj*ti + ni*tj (it's a symetric tensor)
  ! where n = normal vector, t = tangential vector 1

  ! x_ijk = xs + α --> xs = x_ijk - α
  ! (xs,ys,zs) : position vector at the free-surface
  ! as phi(i,j,k) is negative, and phi_gradient points towards from the free-surface
  ! so alpha_vec points from the free-surface to (i,j,k) in the air phase

  alpha_vec(1) = phi(i,j,k) * ( phi_gradient(1,i,j,k) / norm2( phi_gradient(1:3,i,j,k) ) ) 
  alpha_vec(2) = phi(i,j,k) * ( phi_gradient(2,i,j,k) / norm2( phi_gradient(1:3,i,j,k) ) ) 
  alpha_vec(3) = phi(i,j,k) * ( phi_gradient(3,i,j,k) / norm2( phi_gradient(1:3,i,j,k) ) ) 

  xs = x(i,j,k) - alpha_vec(1)
  ys = y(i,j,k) - alpha_vec(2)
  zs = z(i,j,k) - alpha_vec(3)

  PointWithinCell = .false.
  
  ! Get the normal and tangential vectors (and xnut) at the free surface
  call nts_vectors_fs( i    , j    , k    , &
                       xs   , ys   , zs   , &
                       nvec , tvec , svec , &
                       xnut_fs            , &
                       PointWithinCell      &
                      ) 

  ! if the interpolation to the free-surface failed, I just use the local nts 
  ! as in the original formulation of Watanabe et al.
  if ( .not. PointWithinCell ) then
  
    call nts_vectors(  phi_work, phi_grad_work,csi_work, eta_work, zet_work, &
                       dc, de, dz, aj_work, nvec, tvec, svec )

    xnut_fs = xnut(i,j,k)
  
  end if

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


  ! ---------------------------------------------------------
  ! 3. Computation of T and B matrices using Least-Square  
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

  least_dis      = 200.0_rdf ! 20.0: a huge number for initialisation
  aux_least_dis  = zero             ! flag_variable
      
  velocity_gradient_extp = zero 
  alpha_vec_extp         = zero

  ! Neighbourhood indexes definition
  im = max( ista , i-sweep_lsqm ) ; ip = min( iend , i+sweep_lsqm )
  jm = max( jsta , j-sweep_lsqm ) ; jp = min( jend , j+sweep_lsqm )
  km = max( ksta , k-sweep_lsqm ) ; kp = min( kend , k+sweep_lsqm )

  ! Max distance 
  ! dmax = hundred
  ! call get_max_distance_neighbourhood(xs,ys,zs,im,ip,jm,jp,km,kp,dmax)
  
  max_vel_norm_neighbour = zero
  
  do kk = km , kp
  do jj = jm , jp
  do ii = im , ip

    BlankingFlag = .false. 

    if ( nblke/=0 ) then
      
      do nb = 1,nblke

        if ( ii > le_blk_ia(1,nb) .and. ii < le_blk_ib(1,nb) .and. & 
             jj > le_blk_ja(1,nb) .and. jj < le_blk_jb(1,nb) .and. &
             kk > le_blk_ka(1,nb) .and. kk < le_blk_kb(1,nb) ) then

          BlankingFlag = .true.

        end if
      
      end do

    end if

    ! skip the air nodes
    if ( rsign( phi(ii,jj,kk) ) < one_half .or. BlankingFlag ) cycle

    max_vel_norm_neighbour = max( max_vel_norm_neighbour , &
                                  norm2( (/q(2,ii,jj,kk) , q(3,ii,jj,kk) , q(4,ii,jj,kk)/) ) )

    ! Don't use not resolved nodes for extrapolation
    if ( itc == 1 .and. phi(ii,jj,kk) * phi_n(ii,jj,kk) < zero ) cycle

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

    exsign = rsign( phi(ii,jj,kk) ) ! * ( sign( one , radius_lsqm * dx - rdiff_norm ) + one )/two
    
    ! I'll skip nodes where phi_gradient is large because it's probably taking
    ! nodes from the phi-blanked region (>3 nodes away from the free surface) 
    if ( exsign < one_half .or. &
         norm2( phi_gradient(1:3,ii,jj,kk) ) > two ) cycle ! water phase + searching area

    velocity_curv_gradient = zero ! variable initialisation before update

    ! velocity_curv_gradient = ∂u_i/∂ξ^m tensor                  
    call velocity_curv_gradient_tensor(ii,jj,kk,velocity_curv_gradient,exsign)

    if ( exsign < one_half ) cycle ! water phase + searching area

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

    ! counter for Tt and Ts construction. cont=1..27
    cont  = 1 
        
    ! Tt and Ts vectors initialisation at node ii,jj,kk
    ttvec = zero
    tsvec = zero

    do jvs = 1,3 ! jvs : j-vector system
    do ivs = 1,3 ! ivs : i-vector system
              
      do coeff_loop = 1,3 ! b, c, d components loop for Tt and Ts vecs
                
        ! sizes: ttvec(1:27), tsvec(1:27) 
             
        ttvec(cont) = ttvec(cont) +  nj_ti(ivs,jvs) * alpha_local(coeff_loop)
        tsvec(cont) = tsvec(cont) +  nj_si(ivs,jvs) * alpha_local(coeff_loop)
        
        ! counter update
        cont = cont + 1 

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
    
    ! weight of the current node values
    call get_weight_LSM_GFM( xs , ys , zs , nvec , dmax , ii , jj , kk , weight )

    do col = 1, nmsu
    do row = 1, nmsu
      t_matrix_system(row,col)     =    t_matrix_system ( row , col ) &
                                      + weight * ttvec(row) * ttvec(col) 
    end do
    end do

    !  vertical block 2 (Ts vector)

    do col = 1,nmsu
    do row = 1,nmsu
      t_matrix_system(nmsu+row,col) =   t_matrix_system ( nmsu + row , col ) &
                                      + weight * tsvec(row) * tsvec(col) 
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
    
    velocity_gradient = zero

    ! ∂u      ∂u   ∂ξ     ∂u   ∂η     ∂u   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂x      ∂ξ   ∂x     ∂η   ∂x     ∂ζ   ∂x 
    velocity_gradient(1,1) = velocity_curv_gradient(1,1) * csi(1,ii,jj,kk) + &
                             velocity_curv_gradient(1,2) * eta(1,ii,jj,kk) + &
                             velocity_curv_gradient(1,3) * zet(1,ii,jj,kk)
    
    ! ∂u      ∂u   ∂ξ     ∂u   ∂η     ∂u   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂y      ∂ξ   ∂y     ∂η   ∂y     ∂ζ   ∂y 
    velocity_gradient(1,2) = velocity_curv_gradient(1,1) * csi(2,ii,jj,kk) + &
                             velocity_curv_gradient(1,2) * eta(2,ii,jj,kk) + &
                             velocity_curv_gradient(1,3) * zet(2,ii,jj,kk)
    
    ! ∂u      ∂u   ∂ξ     ∂u   ∂η     ∂u   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂z      ∂ξ   ∂z     ∂η   ∂z     ∂ζ   ∂z 
    velocity_gradient(1,3) = velocity_curv_gradient(1,1) * csi(3,ii,jj,kk) + &
                             velocity_curv_gradient(1,2) * eta(3,ii,jj,kk) + &
                             velocity_curv_gradient(1,3) * zet(3,ii,jj,kk)
    
    ! ∂v      ∂v   ∂ξ     ∂v   ∂η     ∂v   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂x      ∂ξ   ∂x     ∂η   ∂x     ∂ζ   ∂x 
    velocity_gradient(2,1) = velocity_curv_gradient(2,1) * csi(1,ii,jj,kk) + &
                             velocity_curv_gradient(2,2) * eta(1,ii,jj,kk) + &
                             velocity_curv_gradient(2,3) * zet(1,ii,jj,kk)
    
    ! ∂v      ∂v   ∂ξ     ∂v   ∂η     ∂v   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂y      ∂ξ   ∂y     ∂η   ∂y     ∂ζ   ∂y 
    velocity_gradient(2,2) = velocity_curv_gradient(2,1) * csi(2,ii,jj,kk) + &
                             velocity_curv_gradient(2,2) * eta(2,ii,jj,kk) + &
                             velocity_curv_gradient(2,3) * zet(2,ii,jj,kk)
    
    ! ∂v      ∂v   ∂ξ     ∂v   ∂η     ∂v   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂z      ∂ξ   ∂z     ∂η   ∂z     ∂ζ   ∂z 
    velocity_gradient(2,3) = velocity_curv_gradient(2,1) * csi(3,ii,jj,kk) + &
                             velocity_curv_gradient(2,2) * eta(3,ii,jj,kk) + &
                             velocity_curv_gradient(2,3) * zet(3,ii,jj,kk)
    
    ! ∂w      ∂w   ∂ξ     ∂w   ∂η     ∂w   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂x      ∂ξ   ∂x     ∂η   ∂x     ∂ζ   ∂x 
    velocity_gradient(3,1) = velocity_curv_gradient(3,1) * csi(1,ii,jj,kk) + &
                             velocity_curv_gradient(3,2) * eta(1,ii,jj,kk) + &
                             velocity_curv_gradient(3,3) * zet(1,ii,jj,kk)
    
    ! ∂w      ∂w   ∂ξ     ∂w   ∂η     ∂w   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂y      ∂ξ   ∂y     ∂η   ∂y     ∂ζ   ∂y 
    velocity_gradient(3,2) = velocity_curv_gradient(3,1) * csi(2,ii,jj,kk) + &
                             velocity_curv_gradient(3,2) * eta(2,ii,jj,kk) + &
                             velocity_curv_gradient(3,3) * zet(2,ii,jj,kk)
    
    ! ∂w      ∂w   ∂ξ     ∂w   ∂η     ∂w   ∂ζ 
    !---- =  ---- ---- + ---- ---- + ---- ----   
    ! ∂z      ∂ξ   ∂z     ∂η   ∂z     ∂ζ   ∂z 
    velocity_gradient(3,3) = velocity_curv_gradient(3,1) * csi(3,ii,jj,kk) + &
                             velocity_curv_gradient(3,2) * eta(3,ii,jj,kk) + &
                             velocity_curv_gradient(3,3) * zet(3,ii,jj,kk)


    ! This routine calculates the gradient with the Jacobian inside                         
    !call velocity_curv_gradient_tensor2(ii,jj,kk,velocity_gradient,exsign)

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
    inearest = inearest + nint( aux_least_dis ) * ( -inearest + ii )
    jnearest = jnearest + nint( aux_least_dis ) * ( -jnearest + jj )
    knearest = knearest + nint( aux_least_dis ) * ( -knearest + kk )


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

        b_matrix_system(row) =    b_matrix_system(row) +  weight * &
                                  (velocity_gradient(ims,jms) +    &
                                   velocity_gradient(jms,ims))*    &
                                   nvec(jms)*tvec(ims)*ttvec(row)
      end do
      end do
        
    end do

    !  vertical block 2 (Bs vector)
    do row = 1,nmsu
      
      do jms = 1,3
      do ims = 1,3
        
        b_matrix_system(nmsu+row) =     b_matrix_system(nmsu+row) +  weight * &
                                        (velocity_gradient(ims,jms) +    &
                                         velocity_gradient(jms,ims))*    &
                                         nvec(jms)*svec(ims)*tsvec(row)
      end do
      end do
    
    end do

    ! continuity correction coc = ∂u_i/∂x_i
    ! it is something that is implemented in the Hokkaido code, but for our tests
    ! it didn't make a big difference. We keep it commented in the code in case we
    ! need to activate in the future (2022/03/14)   

    coc = velocity_gradient(1,1) + velocity_gradient(2,2) + velocity_gradient(3,3) 
   
    b_matrix_system(2*nmsu+1) = b_matrix_system(2*nmsu+1) + weight * coc * alpha_local(1)
    b_matrix_system(2*nmsu+2) = b_matrix_system(2*nmsu+2) + weight * coc * alpha_local(2)
    b_matrix_system(2*nmsu+3) = b_matrix_system(2*nmsu+3) + weight * coc * alpha_local(3)

    do icoc = 1,3

      t_matrix_system(2*nmsu+icoc,1) =   t_matrix_system(2*nmsu+icoc,1) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,2) =   t_matrix_system(2*nmsu+icoc,2) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,3) =   t_matrix_system(2*nmsu+icoc,3) &
                                       + weight * alpha_local(icoc) * alpha_local(3)

      t_matrix_system(2*nmsu+icoc,13) =  t_matrix_system(2*nmsu+icoc,13) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,14) =  t_matrix_system(2*nmsu+icoc,14) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,15) =  t_matrix_system(2*nmsu+icoc,15) &
                                       + weight * alpha_local(icoc) * alpha_local(3)

      t_matrix_system(2*nmsu+icoc,25) =  t_matrix_system(2*nmsu+icoc,25) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,26) =  t_matrix_system(2*nmsu+icoc,26) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,27) =  t_matrix_system(2*nmsu+icoc,27) &
                                       + weight * alpha_local(icoc) * alpha_local(3)
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

  su  = zero ! sgelsd output parameter

  call dgelsd ( mmsu, nmsu, nrhs, t_matrix_system,           &
                ldau, b_matrix_system, ldbu, su, rcond, rank, &
                worku, lworku, iworku, info )

  ! We assign the A vector (a_ij coeffs) with the solution of SVD 
  ! system overwritten to b_matrix_system
      
  a_coeff_vector = zero

  do ims = 1, nmsu
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

  ! With the velocity gradient at the free surface, I apply the NDBC
  call normal_dynamic_bc( du_dx_fs_lsqm , xnut_fs , nvec, pfs )

  ! With the pressure at the fs from the NDBC, I compute the pressure gradient at the 
  ! free surface using a Weighted Least Square Method
  dp_dx_fs_lsqm = zero
  call free_surface_pressure_gradient( i , j , k , xs, ys, zs, pfs , nvec , dp_dx_fs_lsqm )

  
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

  !call ghost_nodes_velocity_extrapolation( i  , j  , k   , &
  !                                         xs , ys , zs  , &
  !                                         u_fs_lsqm     , &
  !                                         du_dx_fs_lsqm   &
  !                                       )


   call ghost_nodes_extrapolation( i  , j  , k    , &
                                   xs , ys , zs   , &
                                   u_fs_lsqm      , &
                                   pfs            , &
                                   du_dx_fs_lsqm  , &
                                   dp_dx_fs_lsqm    &
                                 )

  !---------------------------------------------------
  ! zero - shear condition error computation
  !---------------------------------------------------

  ! commented to avoid unnecesary calculations
  !call error_gfm(i, j, k, nvec, tvec, svec, du_dx_fs_lsqm)

  ! deallocate big arrays
  deallocate( t_matrix_system   , a_coeff_vector   , b_matrix_system   )

end do ! gfm nodes

!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
! I'll use the previously implemented pressure extrapolation to
! see if it makes any difference for the solitary wave test case
! call pressure_extrapolation()
!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


! Deallocate big arrays
!deallocate( least_dis_extp , error_tdir, error_sdir, pJ , rsign_aux )
deallocate( least_dis_extp )

! deaollocate svd arrays
deallocate( worku , iworku )


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

contains

include 'get_gfm_list.F90'
include 'nts_vectors.F90'
include 'nts_vectors_fs.F90'
include 'tangential_vectors.F90'
include 'tangential_vectors_fs.F90'
include 'hessian_matrix_curvilinear.F90'

!include 'pressure_extrapolation.F90'
!include 'p_correction_post_advection.F90'
include 'q_correction_post_advection.F90'
include 'VelocityGradientTensor.F90'
include 'PhiGradientVector.F90'
include 'rhs_exchng3_3d.F90'
include 'rhs_exchng3_4d.F90'
include 'get_weight_LSM_GFM.F90'
include 'velocity_curv_gradient_tensor.F90'
!include 'get_max_distance_neighbourhood.F90'

include 'velocity_gradient_extrapolation_free_surface_lsqm.F90'
include 'normal_dynamic_bc.F90'
include 'free_surface_pressure_gradient.F90'
include 'velocity_extrapolation_free_surface_lsqm.F90'
include 'ghost_nodes_velocity_extrapolation.F90'
include 'ghost_nodes_extrapolation.F90'


!Uncomment for debugging
!subroutine error_gfm(i, j, k, nvec, tvec, svec, du_dx_fs_lsqm)
!
!   implicit none
!
!   integer, intent(in) :: i,j,k 
!   real (kind = rdf), dimension(1:3), intent(in) :: nvec, tvec, svec  ! vector system
!   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm ! velocity gradient at the
!                                                                      ! free-surface
!
!   ! local variables
!   integer :: isum, jsum ! indexes for Einstein sumation
!
!   error_tdir(i,j,k) = zero
!   error_sdir(i,j,k) = zero
!   
!   do isum = 1,3
!      do jsum = 1,3
!
!         error_tdir(i,j,k) = error_tdir(i,j,k) + (  du_dx_fs_lsqm(isum,jsum)     &
!                                                  + du_dx_fs_lsqm(jsum,isum) ) * &
!                                                    nvec(isum)*tvec(jsum)   
!
!         error_sdir(i,j,k) = error_sdir(i,j,k) + (  du_dx_fs_lsqm(isum,jsum)     &
!                                                  + du_dx_fs_lsqm(jsum,isum) ) * &
!                                                    nvec(isum)*svec(jsum)   
!      end do
!   end do
!
!
!end subroutine error_gfm


end subroutine ghost_fluid_extrapolation