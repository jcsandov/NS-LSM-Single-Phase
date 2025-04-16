subroutine ghost_fluid_extrapolation( il,iu          , &
                                      jl,ju          , &
                                      kl,ku          , &
                                      igp, jgp, kgp  , &
                                      dc, de, dz     , &
                                      q              , &
                                      csi            , &
                                      eta            , &
                                      zet            , &
                                      aj             , &
                                      phi            , & 
                                      phi_n          , & 
                                      phi_gradient   , & 
                                      h              , & 
                                      hn             , & 
                                      rsign          , &
                                      x,y,z          , &                                    
                                      iteraciontiempo  &
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
use global_lsm, only: mms, nms, sweep_lsqm, radius_lsqm, phi_outputiter, FrLSM
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

real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: phi
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: phi_n
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(in)    :: h , hn
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)     , intent(inout) :: rsign
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
integer :: iteraciontiempo
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

real (kind = rdf), dimension(:,:), allocatable   :: t_matrix_system_p
real (kind = rdf), dimension(:)  , allocatable   :: a_coeff_vector_p
real (kind = rdf), dimension(:)  , allocatable   :: b_matrix_system_p
real (kind = rdf) :: weight


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
real (kind = rdf), dimension(1:3)     :: pressure_gradient ! local pressure gradients

real (kind = rdf), dimension(1:3,1:3) :: velocity_gradient_extp ! local vel grad used for extp
real (kind = rdf), dimension(1:3)     :: pressure_gradient_extp ! local p grad used for extp

! indexes for local velocity gradient computation
integer :: ivg, jvg, msum

! error arrays ( commented when we're not debbuging )
! real (kind = rdf), dimension(:,:,:), allocatable :: error_tdir, error_sdir

! ----------------------------------------------------------------------------------
! Singular - Value Decomposition (SVD) parameters
! ----------------------------------------------------------------------------------

integer, parameter :: mmsu = 75
integer, parameter :: nmsu = 36

integer, parameter :: mmsp = 32
integer, parameter :: nmsp = 12

! SGELSD
integer :: info, rank
integer :: nlvlu, lworku, liworku
integer :: nlvlp, lworkp, liworkp


integer, parameter :: smlsiz=25
integer, parameter :: nrhs=1

real(kind = rdf), dimension(:), allocatable :: iworku , worku
real(kind = rdf), dimension(:), allocatable :: iworkp , workp

integer, parameter:: ldau=max(1,mmsu)
integer, parameter:: ldap=max(1,mmsp)

integer, parameter:: ldbu=max(1,max(mmsu,nmsu))
integer, parameter:: ldbp=max(1,max(mmsp,nmsp))

real(kind = rdf) :: ersvdm !!,sm,res,resf
real(kind = rdf) :: su(min(mmsu,nmsu)) 
real(kind = rdf) :: sp(min(mmsp,nmsp)) 

real(kind = rdf), parameter:: rcond = 0.00001_rdf

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


nlvlp = max(0,int( log(real(min(mmsp,nmsp))/real(smlsiz+1))/log(2.) )+1)
if(mmsp.ge.nmsp) then      
  lworkp  = 12*nmsp+2*nmsp*smlsiz+8*nmsp*nlvlp+ nmsp*nrhs+(smlsiz+1)*(smlsiz+1)
else
  lworkp  = 12*mmsp+2*mmsp*smlsiz+8*mmsp*nlvlp+ mmsp*nrhs+(smlsiz+1)*(smlsiz+1)
end if

liworkp = max(1,3*min(mmsp,nmsp)*nlvlp+11*min(mmsp,nmsp))

allocate(workp(max(1,lworkp)) , iworkp(max(1,liworkp)))

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


! correct pressure after level-set advection for nodes who kept
! negative pressures after being crossed by the free surface (it does
! it only the first pseudo-time iteration)
! if ( itc  == 1 ) then

call p_correction_post_advection()

  ! return the pressure to dynamic pressure
  !if ( hydraulic_mode ) q(1,:,:,:) = q(1,:,:,:) - ( one / FrLSM**2 ) * hn(:,:,:)   

! end if

! Extrapolate the pressure to the first layer of ghost nodes
call pressure_extrapolation()

! update ghost nodes for pressure
call rhs_exchng3_3d ( q(1,:,:,:) )


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
do k = k_mysta , k_myend
do j = j_mysta , j_myend
do i = i_mysta , i_myend
         
  ! ----------------------------------------------------------------------------
  ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
  ! ----------------------------------------------------------------------------
  ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
  ! phase (rsign = 1), then this if is true. This condition identifies air nodes
  ! next to the free-surface
  ! ----------------------------------------------------------------------------

  if ( rsign(i,j,k) > one_half ) cycle ! water-phase

  nWaterNeighbours = zero

  search_water_neighbour : &   
  do kk = k-1 , k+1
  do jj = j-1 , j+1
  do ii = i-1 , i+1
      
    !nWaterNeighbours = nWaterNeighbours + rsign(ii,jj,kk)
         
    if ( rsign(ii,jj,kk) > one_half ) then

      nWaterNeighbours = one 
      exit search_water_neighbour

    end if

  end do
  end do
  end do search_water_neighbour

  if ( nWaterNeighbours < one_half ) cycle  

  ! Allocating matrix arrays for SVD solver

  allocate(  t_matrix_system ( 1:mmsu , 1:nmsu  ) , &
             a_coeff_vector  ( 1:nmsu )           , &
             b_matrix_system ( 1:max(mmsu,nmsu) )   &
           )    

  t_matrix_system = zero
  a_coeff_vector  = zero
  b_matrix_system = zero


  allocate(  t_matrix_system_p ( 1:mmsp , 1:nmsp  ) , &
             a_coeff_vector_p  ( 1:nmsp )           , &
             b_matrix_system_p ( 1:max(mmsp,nmsp) )   &
           )    

  t_matrix_system_p = zero
  a_coeff_vector_p  = zero
  b_matrix_system_p = zero


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
  
  call nts_vectors_fs( i , j , k , xs , ys , zs, nvec , tvec , svec , PointWithinCell ) 

  ! if the interpolation to the free-surface failed, I just use the local nts 
  ! as in the original formulation
  if ( .not. PointWithinCell ) then
  
    call nts_vectors(  phi_work, phi_grad_work,csi_work, eta_work, zet_work, &
                       dc, de, dz, aj_work, nvec, tvec, svec )
  end if

  !nvecfield(1:3,i,j,k) = nvec
  !tvecfield(1:3,i,j,k) = tvec
  !svecfield(1:3,i,j,k) = svec


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
  pressure_gradient_extp = zero 
  alpha_vec_extp         = zero

  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
  !if ( i == 220 .and. j < 41 ) then
  !    
  !  print *, '============================================================================='
  !  print *, ' '
  !  print *, '                   i,j,k = ', i,j,k
  !  print *, ' '
  !  print *, '============================================================================='
  !
  !end if
  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

  im = max( ista , i-sweep_lsqm ) ; ip = min( iend , i+sweep_lsqm )
  jm = max( jsta , j-sweep_lsqm ) ; jp = min( jend , j+sweep_lsqm )
  km = max( ksta , k-sweep_lsqm ) ; kp = min( kend , k+sweep_lsqm )

  ! Max distance 
  dmax = hundred
  call get_max_distance_neighbourhood(xs,ys,zs,im,ip,jm,jp,km,kp,dmax)
  
  do kk = km , kp
  do jj = jm , jp
  do ii = im , ip

    ! skip the air nodes
    if ( rsign(ii,jj,kk) < one_half ) cycle

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

    exsign = rsign(ii,jj,kk) ! * ( sign( one , radius_lsqm * dx - rdiff_norm ) + one )/two
    if ( exsign < one_half ) cycle ! water phase + searching area

    velocity_curv_gradient = zero ! variable initialisation before update

    ! velocity_curv_gradient = ∂u_i/∂ξ^m tensor                  
    !call velocity_curv_gradient_tensor(ii,jj,kk,velocity_curv_gradient,exsign)
    call velocity_curv_gradient_tensor3(ii,jj,kk,velocity_curv_gradient,exsign)

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

    ! counter for Tt and Ts construction. cont=1..36
    cont  = 1 
        
    ! Tt and Ts vectors initialisation at node ii,jj,kk
    ttvec = zero
    tsvec = zero

    do jvs = 1,3 ! jvs : j-vector system
    do ivs = 1,3 ! ivs : i-vector system
              
      do coeff_loop = 0,3 ! a, b, c, d components loop for Tt and Ts vecs
                
        ! sizes: ttvec(1:36), tsvec(1:36) 
             
        ttvec(cont) = ttvec(cont) +  nj_ti(ivs,jvs) * alpha_local(coeff_loop)
        tsvec(cont) = tsvec(cont) +  nj_si(ivs,jvs) * alpha_local(coeff_loop)
        
        ! counter update
        cont = cont + 1 

      end do
        
    end do 
    end do

    ! Pressure extrapolation matrices
    ! counter for Tt and Ts construction. cont=1..36
    cont  = 1 
        
    ! Tt and Ts vectors initialisation at node ii,jj,kk
    ttvec_p = zero
    tsvec_p = zero

    do ivs = 1,3 ! ivs : i-vector system
      do coeff_loop = 0,3 ! a, b, c, d components loop for Tt and Ts vecs
                 
        ! sizes: ttvec(1:36), tsvec(1:36) 
              
        ttvec_p(cont) = ttvec_p(cont) +  tvec(ivs) * alpha_local(coeff_loop)
        tsvec_p(cont) = tsvec_p(cont) +  svec(ivs) * alpha_local(coeff_loop)
                 
        ! counter update
        cont = cont + 1 

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
    
    !weight = abs( dot_product( nvec , rdiff / rdiff_norm ) )
    call get_weight_LSQM( xs , ys , zs , nvec , dmax , ii , jj , kk , weight )

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

    ! vertical block 1 (Tt vector)
    
    do col = 1,nmsp
    do row = 1,nmsp
      t_matrix_system_p(row,col) =   t_matrix_system_p(row,col)    &
                                   + weight * ttvec_p(row) * ttvec_p(col) 
    end do
    end do

    !  vertical block 2 (Ts vector)

    do col = 1,nmsp
    do row = 1,nmsp
      t_matrix_system_p(nmsp+row,col) =   t_matrix_system_p(nmsp+row,col)  &
                                        + weight * tsvec_p(row) * tsvec_p(col) 
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

    !if ( exsign < one_half ) cycle

    ! if exsign returns as zero, then velocity_curv_gradient = 0
    velocity_curv_gradient = exsign * velocity_curv_gradient        
    
    velocity_gradient = zero

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
   
    ! velocity_gradient = ∂u_i/∂ξ^m * ∂ξ^m/∂x_j (Einstein's summation)

    !do jvg = 1,3
    !do ivg = 1,3
    !        
    !  do msum = 1, 3
    !
    !    velocity_gradient(ivg,jvg) =  velocity_gradient(ivg,jvg)       + &
    !                                  velocity_curv_gradient(ivg,msum) * &
    !                                  metric_tensor_vg(msum,jvg)                                 
    !  end do
    !      
    !end do
    !end do

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

      t_matrix_system(2*nmsu+icoc,2) =    t_matrix_system(2*nmsu+icoc,2) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,3) =    t_matrix_system(2*nmsu+icoc,3) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,4) =    t_matrix_system(2*nmsu+icoc,4) &
                                       + weight * alpha_local(icoc) * alpha_local(3)

      t_matrix_system(2*nmsu+icoc,18) =   t_matrix_system(2*nmsu+icoc,18) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,19) =   t_matrix_system(2*nmsu+icoc,19) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,20) =   t_matrix_system(2*nmsu+icoc,20) &
                                       + weight * alpha_local(icoc) * alpha_local(3)

      t_matrix_system(2*nmsu+icoc,34) =   t_matrix_system(2*nmsu+icoc,34) &
                                       + weight * alpha_local(icoc) * alpha_local(1)
      t_matrix_system(2*nmsu+icoc,35) =   t_matrix_system(2*nmsu+icoc,35) &
                                       + weight * alpha_local(icoc) * alpha_local(2)
      t_matrix_system(2*nmsu+icoc,36) =   t_matrix_system(2*nmsu+icoc,36) &
                                       + weight * alpha_local(icoc) * alpha_local(3)

    end do

    ! - - - - - - - - - - - - - -- - - - - - -- - - - - - -- - - - - - -- - - - - - -
    !   Pressure LSM: Tp * avec = bvec
    ! - - - - - - - - - - - - - -- - - - - - -- - - - - - -- - - - - - -- - - - - - -

    !print *,' '
    !print *, 'ii,jj,kk = ', ii,jj,kk
    !print *,' '

    call calc_pressure_gradient( ii , jj , kk , pressure_gradient , exsign)
    !call calc_pressure_gradient2( ii , jj , kk , pressure_gradient , exsign)
    
    if ( exsign < one_half ) cycle
    
    do row = 1,nmsp
      do ims = 1,3
     
        b_matrix_system_p( row ) =   b_matrix_system_p( row )                                   &
                                   + weight * pressure_gradient(ims) * tvec(ims) * ttvec_p (row)
    
      end do
    end do
       
    !  vertical block 2 (Bs vector)
   
    do row = 1,nmsp
      do ims = 1,3
     
        b_matrix_system_p( nmsp + row ) =   b_matrix_system_p( nmsp + row )                     &
                                          + weight * pressure_gradient(ims) * svec(ims) * tsvec_p (row)
    
      end do
    end do

    !! neares5 pressure gradient to use for extrapolation to the free surface
    pressure_gradient_extp =   pressure_gradient_extp &
                             + aux_least_dis * (-pressure_gradient_extp + pressure_gradient)

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

  ! SVD solver for the pressure system
  sp  = zero ! sgelsd output parameter

  call dgelsd ( mmsp, nmsp, nrhs, t_matrix_system_p,          &
                ldap, b_matrix_system_p, ldbp, sp, rcond, rank, &
                workp, lworkp, iworkp, info )

  
  a_coeff_vector_p = zero
  
  do ims = 1, nmsp
     a_coeff_vector_p(ims) = b_matrix_system_p(ims)
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

  dp_dx_fs_lsqm = zero
  
  call pressure_gradient_extrapolation_free_surface_lsqm ( a_coeff_vector_p        , &
                                                           alpha_vec_extp          , & 
                                                           pressure_gradient_extp  , &
                                                           dp_dx_fs_lsqm             &
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

  !call ghost_nodes_velocity_extrapolation( i  , j  , k   , &
  !                                         xs , ys , zs  , &
  !                                         u_fs_lsqm     , &
  !                                         du_dx_fs_lsqm   &
  !                                       )


   call ghost_nodes_extrapolation( i  , j  , k    , &
                                   xs , ys , zs   , &
                                   u_fs_lsqm      , &
                                   du_dx_fs_lsqm  , &
                                   dp_dx_fs_lsqm   &
                                 )

  !---------------------------------------------------
  ! zero - shear condition error computation
  !---------------------------------------------------

  ! commented to avoid unnecesary calculations
  !call error_gfm(i, j, k, nvec, tvec, svec, du_dx_fs_lsqm)

  ! deallocate big arrays
  deallocate( t_matrix_system   , a_coeff_vector   , b_matrix_system   )
  deallocate( t_matrix_system_p , a_coeff_vector_p , b_matrix_system_p )

end do ! i
end do ! j
end do ! k

! Deallocate big arrays
!deallocate( least_dis_extp , error_tdir, error_sdir, pJ , rsign_aux )
deallocate( least_dis_extp )

! deaollocate svd arrays
deallocate( worku , iworku )


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

contains

include 'pressure_extrapolation.F90'
include 'p_correction_post_advection.F90'
include 'VelocityGradientTensor.F90'
include 'PhiGradientVector.F90'
include 'rhs_exchng3_3d.F90'

subroutine get_max_distance_neighbourhood(xs,ys,zs,im,ip,jm,jp,km,kp,dmax)

  !------------------------------------------------------------------------
  ! Description:
  ! -------------
  ! get_max_distance_neighbourhood returns the largest distance
  ! among the water nodes in the region defined by the nodes
  ! im,ip,jm,jp,km,kp
  !------------------------------------------------------------------------

  implicit none
  ! input/output variables
  real(kind = rdf) , intent(in) :: xs,ys,zs
  integer , intent(in) :: im,ip,jm,jp,km,kp
  real(kind = rdf) :: dmax

  !local variables
  real(kind=rdf),dimension(1:3) :: delta_r
  integer :: ii,jj,kk
  real(kind=rdf) :: aux_sign

  dmax = zero

  do kk = km , kp
  do jj = jm , jp
  do ii = im , ip

    aux_sign = one

    ! If I check an air node, the Δr vector is set to zero
    if( rsign(ii,jj,kk) < one_half ) aux_sign = zero

    delta_r(1) = aux_sign * ( xs - x(ii,jj,kk) )
    delta_r(2) = aux_sign * ( ys - y(ii,jj,kk) )
    delta_r(3) = aux_sign * ( zs - z(ii,jj,kk) )

    dmax = max( dmax , norm2(delta_r) )

  end do
  end do
  end do

end subroutine get_max_distance_neighbourhood


subroutine get_weight_LSQM( xs , ys , zs , nvec , dmax , ii , jj , kk , weight )

  implicit none

  real(kind = rdf) , intent(in) :: xs,ys,zs
  real(kind=rdf),dimension(1:3) , intent(in) :: nvec
  integer , intent(in) :: ii,jj,kk
  real(kind = rdf) , intent(in) :: dmax
  real(kind = rdf) :: weight

  !local variables
  real(kind=rdf),dimension(1:3) :: delta_r
  real(kind=rdf) :: num , denom
  real(kind=rdf), parameter :: eps_denom = 1.0E-6

  if ( rsign(ii,jj,kk) < one_half ) return

  !---------------------------------------------------
  !   
  !          | r̂ • n̂ |
  !   w  = -------------- 
  !            d̄ + ε   
  !
  ! where:  n̂ : normal vector at the free surface
  !         r̂ = r / ‖ r ‖
  !         d̄ = d / dmax
  !         ε : parameter to avoid division by zero  
  !--------------------------------------------------

  delta_r(1) = x(ii,jj,kk) - xs
  delta_r(2) = y(ii,jj,kk) - ys
  delta_r(3) = z(ii,jj,kk) - zs

  num   = abs( dot_product( delta_r / norm2(delta_r)  , nvec ) )
  denom = ( norm2(delta_r) / dmax ) + eps_denom

  weight = num / denom

end subroutine get_weight_LSQM

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


subroutine nts_vectors_fs( i,j,k, xs,ys,zs,nvec, tvec, svec, PointWithinCell)
   
   use InterpolationMethods

   implicit none

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! input - output variables
   !
   ! values at (0,0,0) mean values at (i,j,k), so values at (-1,0,1) imply values at
   ! (i-1,j,k+1). This structure is to avoid big arrays communication between program
   ! units when it is no needed.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   integer, intent(in) :: i,j,k
   real (kind = rdf) , intent(in) :: xs,ys,zs
   real (kind = rdf), dimension(3), intent(out) :: nvec, tvec, svec
   logical :: PointWithinCell 

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   ! local variables

   real (kind = rdf) :: phi_gradient_norm

   ! local variables interpolation

   integer                                      :: iv1 , jv1 , kv1
   integer                                      :: iv  , jv  , kv 
   real ( kind = rdf ) , dimension(8,3)         :: CellVerticesCoordinates
   real ( kind = rdf ) , dimension(6)           :: DistanceToFaces

   integer :: CellVertexLoop

   ! Phi gradietnt
   real ( kind = rdf ) , dimension(8)           :: dphidx_CellArray , dphidy_CellArray , dphidz_CellArray
   ! Tangential vectors
   real ( kind = rdf ) , dimension(8)           :: t1_CellArray , t2_CellArray , t3_CellArray
   real ( kind = rdf ) , dimension(8)           :: s1_CellArray , s2_CellArray , s3_CellArray

   integer , parameter                          :: nvars = 9
   real ( kind = rdf ) , dimension( nvars , 8 ) :: VarsToInterpolate_CellArray
   real ( kind = rdf ) , dimension( nvars )     :: nts_interpolated

   integer, parameter :: iOffset(1:8) = (/0,1,1,0,0,1,1,0/) ! i - index offset
   integer, parameter :: jOffset(1:8) = (/0,0,1,1,0,0,1,1/) ! j - index offset
   integer, parameter :: kOffset(1:8) = (/0,0,0,0,1,1,1,1/) ! k - index offset
  

   PointWithinCell = .false.

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! normal vector interpolation
   
   ! I look over the eight octants aroun i,j,k
   searchloop_extp : &
   do iv1 = i-1 , i 
   do jv1 = j-1 , j 
   do kv1 = k-1 , k 

      ! 
      !     i,j+1,k+1 ----i+1,j+1,k+1    
      !      /|(8)           /|(7)                            
      !     / |             / |                                  
      !  i,j,k+1-------i+1,j,k+1                   
      !    |(5)           |(6)|                                   
      !    |  |           |   |                                
      !    |  |           |   |                                
      !    |  i,j+1,k-----|-i+1,j+1,k           
      !    | /(4)         |  /(3)                 
      !    |/             | /
      !  i,j,k-----------i+1,j,k
      !   (1)               (2)
      ! 
      
      ! Coordinates of the eight vertices of the cell
      CellVerticesCoordinates(1,:) = (/ x( iv1   , jv1   , kv1   ) , &
                                        y( iv1   , jv1   , kv1   ) , &
                                        z( iv1   , jv1   , kv1   ) /)

      CellVerticesCoordinates(2,:) = (/ x( iv1+1 , jv1   , kv1   ) , &
                                        y( iv1+1 , jv1   , kv1   ) , &
                                        z( iv1+1 , jv1   , kv1   ) /)

      CellVerticesCoordinates(3,:) = (/ x( iv1+1 , jv1+1 , kv1   ) , &
                                        y( iv1+1 , jv1+1 , kv1   ) , &
                                        z( iv1+1 , jv1+1 , kv1   ) /)

      CellVerticesCoordinates(4,:) = (/ x( iv1   , jv1+1 , kv1   ) , &
                                        y( iv1   , jv1+1 , kv1   ) , &
                                        z( iv1   , jv1+1 , kv1   ) /)
   
      CellVerticesCoordinates(5,:) = (/ x( iv1   , jv1   , kv1+1 ) , &
                                        y( iv1   , jv1   , kv1+1 ) , &
                                        z( iv1   , jv1   , kv1+1 ) /)

      CellVerticesCoordinates(6,:) = (/ x( iv1+1 , jv1   , kv1+1 ) , &
                                        y( iv1+1 , jv1   , kv1+1 ) , &
                                        z( iv1+1 , jv1   , kv1+1 ) /)

      CellVerticesCoordinates(7,:) = (/ x( iv1+1 , jv1+1 , kv1+1 ) , &
                                        y( iv1+1 , jv1+1 , kv1+1 ) , &
                                        z( iv1+1 , jv1+1 , kv1+1 ) /)

      CellVerticesCoordinates(8,:) = (/ x( iv1   , jv1+1 , kv1+1 ) , &
                                        y( iv1   , jv1+1 , kv1+1 ) , &
                                        z( iv1   , jv1+1 , kv1+1 ) /)

      ! verification if the fs lies on the current cell

      call PointWithinCellCheck ( CellVerticesCoordinates , &
                                  (/xs,ys,zs/)            , &
                                  PointWithinCell         , &
                                  DistanceToFaces           &
                                )

      if ( PointWithinCell ) exit searchloop_extp

      
   end do
   end do
   end do searchloop_extp

   if ( .not. PointWithinCell ) return

   ! I came out this searchloop with iv1, jv1 and kv1 which tells me the v1 node
   ! of the cell where the free surface is located. 

   ! ---------------------------------------------------------------------------
   ! ϕ GRADIENT FOR NORMAL VECTOR CALCULATION
   ! ---------------------------------------------------------------------------


    do CellVertexLoop = 1,8
             
      ! Current vertex
      iv = iv1 + iOffset( CellVertexLoop )
      jv = jv1 + jOffset( CellVertexLoop )
      kv = kv1 + kOffset( CellVertexLoop )

      dphidx_CellArray( CellVertexLoop ) = phi_gradient( 1 , iv , jv , kv )
      dphidy_CellArray( CellVertexLoop ) = phi_gradient( 2 , iv , jv , kv )
      dphidz_CellArray( CellVertexLoop ) = phi_gradient( 3 , iv , jv , kv )

      ! I calculate the tangential vectors away from the walls 

      if( iv /= ista .and. iv /= iend  .and. &
          jv /= jsta .and. jv /= jend  .and. &
          kv /= ksta .and. kv /= kend ) call tangential_vectors2( iv , jv , kv , tvec , svec )


      t1_CellArray( CellVertexLoop ) = tvec(1)
      t2_CellArray( CellVertexLoop ) = tvec(2)
      t3_CellArray( CellVertexLoop ) = tvec(3)

      s1_CellArray( CellVertexLoop ) = svec(1)
      s2_CellArray( CellVertexLoop ) = svec(2)
      s3_CellArray( CellVertexLoop ) = svec(3)

    end do

    if ( iv == ista ) then
    
      t1_CellArray(1) = t1_CellArray(2) ; t2_CellArray(1) = t2_CellArray(2) ; t3_CellArray(1) = t3_CellArray(2)
      t1_CellArray(4) = t1_CellArray(3) ; t2_CellArray(4) = t2_CellArray(3) ; t3_CellArray(4) = t3_CellArray(3)
      t1_CellArray(5) = t1_CellArray(6) ; t2_CellArray(5) = t2_CellArray(6) ; t3_CellArray(5) = t3_CellArray(6)
      t1_CellArray(8) = t1_CellArray(7) ; t2_CellArray(8) = t2_CellArray(7) ; t3_CellArray(8) = t3_CellArray(7)

      s1_CellArray(1) = s1_CellArray(2) ; s2_CellArray(1) = s2_CellArray(2) ; s3_CellArray(1) = s3_CellArray(2)
      s1_CellArray(4) = s1_CellArray(3) ; s2_CellArray(4) = s2_CellArray(3) ; s3_CellArray(4) = s3_CellArray(3)
      s1_CellArray(5) = s1_CellArray(6) ; s2_CellArray(5) = s2_CellArray(6) ; s3_CellArray(5) = s3_CellArray(6)
      s1_CellArray(8) = s1_CellArray(7) ; s2_CellArray(8) = s2_CellArray(7) ; s3_CellArray(8) = s3_CellArray(7)

    end if

    if ( iv == iend ) then
    
      t1_CellArray(2) = t1_CellArray(1) ; t2_CellArray(2) = t2_CellArray(1) ; t3_CellArray(2) = t3_CellArray(1)
      t1_CellArray(3) = t1_CellArray(4) ; t2_CellArray(3) = t2_CellArray(4) ; t3_CellArray(3) = t3_CellArray(4)
      t1_CellArray(6) = t1_CellArray(5) ; t2_CellArray(6) = t2_CellArray(5) ; t3_CellArray(6) = t3_CellArray(5)
      t1_CellArray(7) = t1_CellArray(8) ; t2_CellArray(7) = t2_CellArray(8) ; t3_CellArray(7) = t3_CellArray(8)

      s1_CellArray(2) = s1_CellArray(1) ; s2_CellArray(2) = s2_CellArray(1) ; s3_CellArray(2) = s3_CellArray(1)
      s1_CellArray(3) = s1_CellArray(4) ; s2_CellArray(3) = s2_CellArray(4) ; s3_CellArray(3) = s3_CellArray(4)
      s1_CellArray(6) = s1_CellArray(5) ; s2_CellArray(6) = s2_CellArray(5) ; s3_CellArray(6) = s3_CellArray(5)
      s1_CellArray(7) = s1_CellArray(8) ; s2_CellArray(7) = s2_CellArray(8) ; s3_CellArray(7) = s3_CellArray(8)

    end if

    if ( jv == jsta ) then
    
      t1_CellArray(1) = t1_CellArray(4) ; t2_CellArray(1) = t2_CellArray(4) ; t3_CellArray(1) = t3_CellArray(4)
      t1_CellArray(2) = t1_CellArray(3) ; t2_CellArray(2) = t2_CellArray(3) ; t3_CellArray(2) = t3_CellArray(3)
      t1_CellArray(5) = t1_CellArray(8) ; t2_CellArray(5) = t2_CellArray(8) ; t3_CellArray(5) = t3_CellArray(8)
      t1_CellArray(6) = t1_CellArray(7) ; t2_CellArray(6) = t2_CellArray(7) ; t3_CellArray(6) = t3_CellArray(7)

      s1_CellArray(1) = s1_CellArray(4) ; s2_CellArray(1) = s2_CellArray(4) ; s3_CellArray(1) = s3_CellArray(4)
      s1_CellArray(2) = s1_CellArray(3) ; s2_CellArray(2) = s2_CellArray(3) ; s3_CellArray(2) = s3_CellArray(3)
      s1_CellArray(5) = s1_CellArray(8) ; s2_CellArray(5) = s2_CellArray(8) ; s3_CellArray(5) = s3_CellArray(8)
      s1_CellArray(6) = s1_CellArray(7) ; s2_CellArray(6) = s2_CellArray(7) ; s3_CellArray(6) = s3_CellArray(7)

    end if

    if ( jv == jend ) then

      t1_CellArray(4) = t1_CellArray(1) ; t2_CellArray(4) = t2_CellArray(1) ; t3_CellArray(4) = t3_CellArray(1)
      t1_CellArray(3) = t1_CellArray(2) ; t2_CellArray(3) = t2_CellArray(2) ; t3_CellArray(3) = t3_CellArray(2)
      t1_CellArray(8) = t1_CellArray(5) ; t2_CellArray(8) = t2_CellArray(5) ; t3_CellArray(8) = t3_CellArray(5)
      t1_CellArray(7) = t1_CellArray(6) ; t2_CellArray(7) = t2_CellArray(6) ; t3_CellArray(7) = t3_CellArray(6)

      s1_CellArray(4) = s1_CellArray(1) ; s2_CellArray(4) = s2_CellArray(1) ; s3_CellArray(4) = s3_CellArray(1)
      s1_CellArray(3) = s1_CellArray(2) ; s2_CellArray(3) = s2_CellArray(2) ; s3_CellArray(3) = s3_CellArray(2)
      s1_CellArray(8) = s1_CellArray(5) ; s2_CellArray(8) = s2_CellArray(5) ; s3_CellArray(8) = s3_CellArray(5)
      s1_CellArray(7) = s1_CellArray(6) ; s2_CellArray(7) = s2_CellArray(6) ; s3_CellArray(7) = s3_CellArray(6)    

    end if

    if ( kv == ksta ) then
    
      t1_CellArray(1) = t1_CellArray(5) ; t2_CellArray(1) = t2_CellArray(5) ; t3_CellArray(1) = t3_CellArray(5)
      t1_CellArray(2) = t1_CellArray(6) ; t2_CellArray(2) = t2_CellArray(6) ; t3_CellArray(2) = t3_CellArray(6)
      t1_CellArray(3) = t1_CellArray(7) ; t2_CellArray(3) = t2_CellArray(7) ; t3_CellArray(3) = t3_CellArray(7)
      t1_CellArray(4) = t1_CellArray(8) ; t2_CellArray(4) = t2_CellArray(8) ; t3_CellArray(4) = t3_CellArray(8)

      s1_CellArray(1) = s1_CellArray(5) ; s2_CellArray(1) = s2_CellArray(5) ; s3_CellArray(1) = s3_CellArray(5)
      s1_CellArray(2) = s1_CellArray(6) ; s2_CellArray(2) = s2_CellArray(6) ; s3_CellArray(2) = s3_CellArray(6)
      s1_CellArray(3) = s1_CellArray(7) ; s2_CellArray(3) = s2_CellArray(7) ; s3_CellArray(3) = s3_CellArray(7)
      s1_CellArray(4) = s1_CellArray(8) ; s2_CellArray(4) = s2_CellArray(8) ; s3_CellArray(4) = s3_CellArray(8)

    end if

    if ( kv == kend ) then
    
      t1_CellArray(5) = t1_CellArray(1) ; t2_CellArray(5) = t2_CellArray(1) ; t3_CellArray(5) = t3_CellArray(1)
      t1_CellArray(6) = t1_CellArray(2) ; t2_CellArray(6) = t2_CellArray(2) ; t3_CellArray(6) = t3_CellArray(2)
      t1_CellArray(7) = t1_CellArray(3) ; t2_CellArray(7) = t2_CellArray(3) ; t3_CellArray(7) = t3_CellArray(3)
      t1_CellArray(8) = t1_CellArray(4) ; t2_CellArray(8) = t2_CellArray(4) ; t3_CellArray(8) = t3_CellArray(4)

      s1_CellArray(5) = s1_CellArray(1) ; s2_CellArray(5) = s2_CellArray(1) ; s3_CellArray(5) = s3_CellArray(1)
      s1_CellArray(6) = s1_CellArray(2) ; s2_CellArray(6) = s2_CellArray(2) ; s3_CellArray(6) = s3_CellArray(2)
      s1_CellArray(7) = s1_CellArray(3) ; s2_CellArray(7) = s2_CellArray(3) ; s3_CellArray(7) = s3_CellArray(3)
      s1_CellArray(8) = s1_CellArray(4) ; s2_CellArray(8) = s2_CellArray(4) ; s3_CellArray(8) = s3_CellArray(4)

    end if

   ! minus sign to align it with the normal vector
   VarsToInterpolate_CellArray( 1  , 1:8 ) =  - dphidx_CellArray 
   VarsToInterpolate_CellArray( 2  , 1:8 ) =  - dphidy_CellArray
   VarsToInterpolate_CellArray( 3  , 1:8 ) =  - dphidz_CellArray

   VarsToInterpolate_CellArray( 4  , 1:8 ) =  t1_CellArray 
   VarsToInterpolate_CellArray( 5  , 1:8 ) =  t2_CellArray
   VarsToInterpolate_CellArray( 6  , 1:8 ) =  t3_CellArray

   VarsToInterpolate_CellArray( 7  , 1:8 ) =  s1_CellArray 
   VarsToInterpolate_CellArray( 8  , 1:8 ) =  s2_CellArray
   VarsToInterpolate_CellArray( 9  , 1:8 ) =  s3_CellArray

   call TrilinearInterpolation( CellVerticesCoordinates     , &
                                (/xs,ys,zs/)                , &
                                DistanceToFaces             , &
                                nvars                       , &
                                VarsToInterpolate_CellArray , &
                                nts_interpolated              &
                              )   

   nvec = nts_interpolated(1:3) / norm2 ( nts_interpolated (1:3) )
   tvec = nts_interpolated(4:6) / norm2 ( nts_interpolated (4:6) )
   svec = nts_interpolated(7:9) / norm2 ( nts_interpolated (7:9) )

end subroutine nts_vectors_fs


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

!   if ( tvec_mod < tol .or. svec_mod < tol .or. &
!        dot_product(nvec,tvec) > tol .or. dot_product(nvec,svec) > tol ) then

      tvec(1) = tvec_guess(1)
      tvec(2) = tvec_guess(2)
      tvec(3) = tvec_guess(3)
      
      svec(1) = svec_guess(1)
      svec(2) = svec_guess(2)
      svec(3) = svec_guess(3)

      ! vectors normalisation

      tvec = tvec / norm2( tvec ) 
      svec = svec / norm2( svec )

!   else ! tangential vector comes from principal directions

   ! vectors normalisation

!      tvec = tvec / tvec_mod
!      svec = svec / svec_mod

!   end if

end subroutine tangential_vectors


subroutine tangential_vectors2( i , j , k , tvec, svec)
   
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

end subroutine tangential_vectors2


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
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: uLL , uL , uC , uR , uRR
   real ( kind = rdf ) :: vLL , vL , vC , vR , vRR
   real ( kind = rdf ) :: wLL , wL , wC , wR , wRR

   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
      
   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = zero ! dummy

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( i == iend ) then
         
      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   end if

   exsign = zero

   ! ∂u/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,1)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,1)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,1)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j   , k )  
      rR  = rsign ( i , j+1 , k )  
      rRR = rsign ( i , j+2 , k )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;
   
   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;

   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = zero ! dummy

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( j == jend ) then
         
      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;

   end if

   exsign = zero

   ! 

   ! ∂u/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,2)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,2)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,2)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then
 
      bias_zet = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j , k   )  
      rR  = rsign ( i , j , k+1 )  
      rRR = rsign ( i , j , k+2 )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;
      
   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;

   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = zero ! dummy

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;

   end if

   exsign = zero

   ! 

   ! ∂u/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,3)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,3)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,3)    &
                                    )

   if ( exsign < one_half ) return                                    


end subroutine velocity_curv_gradient_tensor

subroutine velocity_curv_gradient_tensor2(i, j, k, velocity_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: u1LL , u1L , u1C , u1R , u1RR
   real ( kind = rdf ) :: u2LL , u2L , u2C , u2R , u2RR
   real ( kind = rdf ) :: u3LL , u3L , u3C , u3R , u3RR
   real ( kind = rdf ) :: u4LL , u4L , u4C , u4R , u4RR
   real ( kind = rdf ) :: u5LL , u5L , u5C , u5R , u5RR
   real ( kind = rdf ) :: u6LL , u6L , u6C , u6R , u6RR
   real ( kind = rdf ) :: u7LL , u7L , u7C , u7R , u7RR
   real ( kind = rdf ) :: u8LL , u8L , u8C , u8R , u8RR
   real ( kind = rdf ) :: u9LL , u9L , u9C , u9R , u9RR

   real ( kind = rdf ) :: v1LL , v1L , v1C , v1R , v1RR
   real ( kind = rdf ) :: v2LL , v2L , v2C , v2R , v2RR
   real ( kind = rdf ) :: v3LL , v3L , v3C , v3R , v3RR
   real ( kind = rdf ) :: v4LL , v4L , v4C , v4R , v4RR
   real ( kind = rdf ) :: v5LL , v5L , v5C , v5R , v5RR
   real ( kind = rdf ) :: v6LL , v6L , v6C , v6R , v6RR
   real ( kind = rdf ) :: v7LL , v7L , v7C , v7R , v7RR
   real ( kind = rdf ) :: v8LL , v8L , v8C , v8R , v8RR
   real ( kind = rdf ) :: v9LL , v9L , v9C , v9R , v9RR

   real ( kind = rdf ) :: w1LL , w1L , w1C , w1R , w1RR
   real ( kind = rdf ) :: w2LL , w2L , w2C , w2R , w2RR
   real ( kind = rdf ) :: w3LL , w3L , w3C , w3R , w3RR
   real ( kind = rdf ) :: w4LL , w4L , w4C , w4R , w4RR
   real ( kind = rdf ) :: w5LL , w5L , w5C , w5R , w5RR
   real ( kind = rdf ) :: w6LL , w6L , w6C , w6R , w6RR
   real ( kind = rdf ) :: w7LL , w7L , w7C , w7R , w7RR
   real ( kind = rdf ) :: w8LL , w8L , w8C , w8R , w8RR
   real ( kind = rdf ) :: w9LL , w9L , w9C , w9R , w9RR

   real ( kind = rdf ) :: du_dx , du1_dcsi , du4_deta , du7_dzet  
   real ( kind = rdf ) :: du_dy , du2_dcsi , du5_deta , du8_dzet  
   real ( kind = rdf ) :: du_dz , du3_dcsi , du6_deta , du9_dzet  

   real ( kind = rdf ) :: dv_dx , dv1_dcsi , dv4_deta , dv7_dzet  
   real ( kind = rdf ) :: dv_dy , dv2_dcsi , dv5_deta , dv8_dzet  
   real ( kind = rdf ) :: dv_dz , dv3_dcsi , dv6_deta , dv9_dzet  

   real ( kind = rdf ) :: dw_dx , dw1_dcsi , dw4_deta , dw7_dzet  
   real ( kind = rdf ) :: dw_dy , dw2_dcsi , dw5_deta , dw8_dzet  
   real ( kind = rdf ) :: dw_dz , dw3_dcsi , dw6_deta , dw9_dzet  


   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   velocity_gradient = zero
   
   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      u1LL = zero                                    
      u1L  = zero                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = zero                                    
      u2L  = zero                                    
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = zero                                    
      u3L  = zero                                    
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = zero                                    
      v1L  = zero                                    
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = zero                                    
      v2L  = zero                                    
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = zero                                    
      v3L  = zero                                    
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = zero                                    
      w1L  = zero                                    
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = zero                                    
      w2L  = zero                                    
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = zero                                    
      w3L  = zero                                    
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 


   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      u1LL = zero                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = zero                                    
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = zero                                    
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = zero                                    
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = zero                                    
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = zero                                    
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = zero                                    
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = zero                                    
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = zero                                    
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 
   
   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = zero ! dummy

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = zero

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = zero 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = zero 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = zero 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = zero

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = zero 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = zero 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = zero 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = zero 


   else if ( i == iend ) then
         
      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = zero 
      u1RR = zero 

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = zero 
      u2RR = zero 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = zero 
      u3RR = zero 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = zero 
      v1RR = zero 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = zero 
      v2RR = zero 

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = zero 
      v3RR = zero 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = zero 
      w1RR = zero 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = zero 
      w2RR = zero 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = zero 
      w3RR = zero 

   else

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k)
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 


   end if

   exsign = zero

   ! ∂u/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u1LL, u1L , u1C, u1R, u1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u2LL, u2L , u2C, u2R, u2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u3LL, u3L , u3C, u3R, u3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du3_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v1LL, v1L , v1C, v1R, v1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v2LL, v2L , v2C, v2R, v2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v3LL, v3L , v3C, v3R, v3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv3_dcsi                       &
                                    )

   if ( exsign < one_half ) return   

   ! ∂w/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w1LL, w1L , w1C, w1R, w1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w2LL, w2L , w2C, w2R, w2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w3LL, w3L , w3C, w3R, w3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw3_dcsi                       &
                                    )

   if ( exsign < one_half ) return   

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j   , k )  
      rR  = rsign ( i , j+1 , k )  
      rRR = rsign ( i , j+2 , k )  

      u4LL = zero                                    
      u4L  = zero                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = zero 
      u5L  = zero 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = zero 
      u6L  = zero 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = zero 
      v4L  = zero 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = zero 
      v5L  = zero 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = zero 
      w4L  = zero 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = zero 
      w5L  = zero 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = zero 
      w6L  = zero 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 


   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      u4LL = zero                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = zero 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = zero 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = zero 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = zero 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = zero 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = zero 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = zero 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = zero 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

   
   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = zero ! dummy

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = zero

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = zero

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = zero 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = zero

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = zero

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = zero

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = zero

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = zero

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = zero

   else if ( j == jend ) then
         
      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = zero 
      u4RR = zero 

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = zero 
      u5RR = zero 

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = zero 
      u6RR = zero 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = zero 
      v4RR = zero 

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = zero 
      v5RR = zero 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = zero 
      v6RR = zero 

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = zero 
      w4RR = zero 

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = zero 
      w5RR = zero 

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = zero 
      w6RR = zero 

   else

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 


   end if

   exsign = zero

   ! ∂u/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u4LL, u4L , u4C, u4R, u4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u5LL, u5L , u5C, u5R, u5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u6LL, u6L , u6C, u6R, u6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du6_deta                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v4LL, v4L , v4C, v4R, v4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v5LL, v5L , v5C, v5R, v5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v6LL, v6L , v6C, v6R, v6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv6_deta                       &
                                    )

   if ( exsign < one_half ) return     

   ! ∂w/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w4LL, w4L , w4C, w4R, w4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w5LL, w5L , w5C, w5R, w5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w6LL, w6L , w6C, w6R, w6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw6_deta                       &
                                    )

   if ( exsign < one_half ) return                                    


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then

      bias_zet = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j , k   )  
      rR  = rsign ( i , j , k+1 )  
      rRR = rsign ( i , j , k+2 )  

      u7LL = zero                                    
      u7L  = zero                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = zero 
      u8L  = zero 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = zero 
      u9L  = zero 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = zero 
      v7L  = zero 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = zero 
      v8L  = zero 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = zero 
      v9L  = zero 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = zero 
      w7L  = zero 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = zero 
      w8L  = zero 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = zero 
      w9L  = zero 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 



   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      u7LL = zero                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = zero 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = zero 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = zero 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = zero 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = zero 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = zero 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = zero 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = zero 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

   
   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = zero ! dummy

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = zero 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = zero 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = zero 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = zero 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = zero 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = zero

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = zero 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = zero 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = zero 

   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = zero 
      u7RR = zero 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = zero 
      u8RR = zero 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = zero 
      u9RR = zero 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = zero 
      v7RR = zero 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = zero 
      v8RR = zero 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = zero 
      v9RR = zero 

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = zero 
      w7RR = zero 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = zero 
      w8RR = zero 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = zero 
      w9RR = zero 

   else

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 


   end if

   exsign = zero

   ! ∂u/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u7LL, u7L , u7C, u7R, u7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u8LL, u8L , u8C, u8R, u8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u9LL, u9L , u9C, u9R, u9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du9_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v7LL, v7L , v7C, v7R, v7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v8LL, v8L , v8C, v8R, v8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v9LL, v9L , v9C, v9R, v9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv9_dzet                       &
                                    )

   if ( exsign < one_half ) return        

   ! ∂w/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w7LL, w7L , w7C, w7R, w7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w8LL, w8L , w8C, w8R, w8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w9LL, w9L , w9C, w9R, w9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw9_dzet                       &
                                    )

   if ( exsign < one_half ) return                                                                  



   du_dx = aj(i,j,k) * ( du1_dcsi + du4_deta + du7_dzet ) 
   du_dy = aj(i,j,k) * ( du2_dcsi + du5_deta + du8_dzet ) 
   du_dz = aj(i,j,k) * ( du3_dcsi + du6_deta + du9_dzet ) 

   dv_dx = aj(i,j,k) * ( dv1_dcsi + dv4_deta + dv7_dzet ) 
   dv_dy = aj(i,j,k) * ( dv2_dcsi + dv5_deta + dv8_dzet ) 
   dv_dz = aj(i,j,k) * ( dv3_dcsi + dv6_deta + dv9_dzet ) 

   dw_dx = aj(i,j,k) * ( dw1_dcsi + dw4_deta + dw7_dzet ) 
   dw_dy = aj(i,j,k) * ( dw2_dcsi + dw5_deta + dw8_dzet ) 
   dw_dz = aj(i,j,k) * ( dw3_dcsi + dw6_deta + dw9_dzet ) 

   velocity_gradient(1,1) = du_dx
   velocity_gradient(1,2) = du_dy
   velocity_gradient(1,3) = du_dz

   velocity_gradient(2,1) = dv_dx
   velocity_gradient(2,2) = dv_dy
   velocity_gradient(2,3) = dv_dz

   velocity_gradient(3,1) = dw_dx
   velocity_gradient(3,2) = dw_dy
   velocity_gradient(3,3) = dw_dz

end subroutine velocity_curv_gradient_tensor2


subroutine velocity_curv_gradient_tensor3(i, j, k, velocity_curv_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: uLL , uL , uC , uR , uRR
   real ( kind = rdf ) :: vLL , vL , vC , vR , vRR
   real ( kind = rdf ) :: wLL , wL , wC , wR , wRR

   real ( kind = rdf ) :: ucon_1 , ucon_2 , ucon_3
   real ( kind = rdf ) :: up1 , um1 , up2 , um2 , up3 , um3

   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0


   !  Contravariant velocity
   ucon_1 = csi(1,i,j,k) * q(2,i,j,k) + &
            csi(2,i,j,k) * q(3,i,j,k) + &
            csi(3,i,j,k) * q(4,i,j,k)    

   ucon_2 = eta(1,i,j,k) * q(2,i,j,k) + &
            eta(2,i,j,k) * q(3,i,j,k) + &
            eta(3,i,j,k) * q(4,i,j,k)    

   ucon_3 = zet(1,i,j,k) * q(2,i,j,k) + &
            zet(2,i,j,k) * q(3,i,j,k) + &
            zet(3,i,j,k) * q(4,i,j,k)    


   up1 = zero
   um1 = zero

   up2 = zero
   um2 = zero

   up3 = zero
   um3 = zero

   if ( ucon_1 > 0 ) then
      up1 = one
   else
      um1 = one
   end if

   if ( ucon_2 > 0 ) then
      up2 = one
   else
      um2 = one
   end if

   if ( ucon_3 > 0 ) then
      up3 = one
   else
      um3 = one
   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
 
      ! ∂u/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,1)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,1)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,1)    &
                                       )
   
      if ( exsign < one_half ) return                                    

   else if (i == iend) then

      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

      ! ∂u/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,1)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,1)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,1)    &
                                       )
   
      if ( exsign < one_half ) return     

   else

      rL = ( rsign ( i-1 , j , k ) + abs(rsign ( i-1 , j , k )) ) / two
      rC = ( rsign ( i   , j , k ) + abs(rsign ( i   , j , k )) ) / two
      rR = ( rsign ( i+1 , j , k ) + abs(rsign ( i+1 , j , k )) ) / two

      exsign = up1 * ( rC * rL ) + &
               um1 * ( rR * rC )
      
      if ( exsign < one_half ) return

      velocity_curv_gradient(1,1) = up1 * dc * ( q(2,i,j,k)   - q(2,i-1,j,k)  ) + &
                                    um1 * dc * ( q(2,i+1,j,k) - q(2,i,j,k)  )

      velocity_curv_gradient(2,1) = up1 * dc * ( q(3,i,j,k)   - q(3,i-1,j,k)  ) + &
                                    um1 * dc * ( q(3,i+1,j,k) - q(3,i,j,k)  )

      velocity_curv_gradient(3,1) = up1 * dc * ( q(4,i,j,k)   - q(4,i-1,j,k)  ) + &
                                    um1 * dc * ( q(4,i+1,j,k) - q(4,i,j,k)  )
   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i, j    , k )  
      rR  = rsign ( i, j+1  , k )  
      rRR = rsign ( i, j+2  , k )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;

      ! ∂u/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,2)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,2)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dc                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,2)    &
                                       )
   
      if ( exsign < one_half ) return      
   
   else if ( j == jend ) then

      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

      ! ∂u/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,2)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,2)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dc                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,2)    &
                                       )
   
      if ( exsign < one_half ) return   

   else

      rL = ( rsign ( i , j-1 , k ) + abs(rsign ( i , j-1 , k )) ) / two
      rC = ( rsign ( i , j   , k ) + abs(rsign ( i , j   , k )) ) / two
      rR = ( rsign ( i , j+1 , k ) + abs(rsign ( i , j+1 , k )) ) / two

      exsign = up2 * ( rC * rL ) + um2 * ( rC * rR )
      
      if ( exsign < one_half ) return

      velocity_curv_gradient(1,2) = up2 * de * ( q(2,i,j,k)   - q(2,i,j-1,k)  ) + &
                                    um2 * de * ( q(2,i,j+1,k) - q(2,i,j,k)  )

      velocity_curv_gradient(2,2) = up2 * de * ( q(3,i,j,k)   - q(3,i,j-1,k)  ) + &
                                    um2 * de * ( q(3,i,j+1,k) - q(3,i,j,k)  )

      velocity_curv_gradient(3,2) = up2 * de * ( q(4,i,j,k)   - q(4,i,j-1,k)  ) + &
                                    um2 * de * ( q(4,i,j+1,k) - q(4,i,j,k)  )
   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then

      bias_zet = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j , k   )  
      rR  = rsign ( i , j , k+1 )  
      rRR = rsign ( i , j , k+2 )  

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;
 
      ! ∂u/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,3)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,3)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,3)    &
                                       )
   
      if ( exsign < one_half ) return   

   else if ( k == kend ) then

      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

      ! ∂u/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        uLL, uL , uC, uR, uRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(1,3)    &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
      ! ∂v/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        vLL, vL , vC, vR, vRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(2,3)    &
                                       )
   
      if ( exsign < one_half ) return                                    
   
   
      ! ∂w/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        wLL, wL , wC, wR, wRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        velocity_curv_gradient(3,3)    &
                                       )
   
      if ( exsign < one_half ) return  

   else

      rL = ( rsign ( i , j , k-1 ) + abs(rsign ( i , j , k-1 )) ) / two
      rC = ( rsign ( i , j , k   ) + abs(rsign ( i , j , k   )) ) / two
      rR = ( rsign ( i , j , k+1 ) + abs(rsign ( i , j , k+1 )) ) / two

      exsign = up3 * ( rC * rL ) + um3 * ( rC * rR )
      
      if ( exsign < one_half ) return

      velocity_curv_gradient(1,3) = up3 * dz * ( q(2,i,j,k)   - q(2,i,j,k-1)  ) + &
                                    um3 * dz * ( q(2,i,j,k+1) - q(2,i,j,k)  )

      velocity_curv_gradient(2,3) = up3 * dz * ( q(3,i,j,k)   - q(3,i,j,k-1)  ) + &
                                    um3 * dz * ( q(3,i,j,k+1) - q(3,i,j,k)  )

      velocity_curv_gradient(3,3) = up3 * dz * ( q(4,i,j,k)   - q(4,i,j,k-1)  ) + &
                                    um3 * dz * ( q(4,i,j,k+1) - q(4,i,j,k)  )
   end if


end subroutine velocity_curv_gradient_tensor3


subroutine calc_pressure_gradient(i, j, k, PressureGradient, exsign ) 
   
  ! Input and output variables
  integer, intent(in) :: i,j,k
  real ( kind = rdf ), intent(inout) :: exsign
  real ( kind = rdf ), dimension(3), intent(out) :: PressureGradient

   ! Local variables
   real ( kind = rdf ), dimension(3)   :: PressureCurvilinearGradient
   real ( kind = rdf ), dimension(3,3) :: MetricsTensor
   real ( kind = rdf ) :: dc2, de2, dz2
   integer :: m,l,p
   integer :: ista , jsta , ksta , iend , jend , kend


   ! Physical boundaries

   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp

   PressureGradient            = zero
   PressureCurvilinearGradient = zero
   MetricsTensor               = zero

   dc2 = one_half * dc
   de2 = one_half * de
   dz2 = one_half * dz 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                                      _         _                                          
   !                                     |   ∂p/∂ξ   |  
   !       PressureCurvilinearGradient = |   ∂p/∂η   |                               
   !                                     |   ∂p/∂ζ   |  
   !                                      •-        -• 
   !       We calculate the derivatives using a second order centred
   !       difference scheme. If the stencil is bounded, then a second
   !       order biased derivative is applied.
   ! 
   !       It assumes the routine that calls it has ista, jsta, ksta and
   !       iend, jend, kend, defined.
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ξ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( myback == mpi_proc_null .and. i == ista ) then

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1)   =  dc2 * (  - one   * q(1,i+2,j,k) &
                                                   + four  * q(1,i+1,j,k) &
                                                   - three * q(1,i  ,j,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i+2,j,k) * rsign(i+1,j,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myfront == mpi_proc_null .and. i == iend) then

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1)   =  dc2 * (    one   * q(1,i-2,j,k) &
                                                   - four  * q(1,i-1,j,k) &
                                                   + three * q(1,i  ,j,k) )

      exsign = abs( rsign(i-2,j,k) * rsign(i-1,j,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1) = dc2 * ( q(1,i+1,j,k) - q(1,i-1,j,k) )

      exsign = abs( rsign(i+1,j,k) * rsign(i-1,j,k) )

      if ( exsign < one_half ) return                                    

   end if
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            η - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( myleft == mpi_proc_null .and. j == jsta ) then

      ! ∂p/∂η
      PressureCurvilinearGradient(2)   =  de2 * (  - one   * q(1,i,j+2,k) &
                                                   + four  * q(1,i,j+1,k) &
                                                   - three * q(1,i,j  ,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j+2,k) * rsign(i,j+1,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myright == mpi_proc_null .and. j == jend) then

      ! ∂p/∂η
      PressureCurvilinearGradient(2)   =  de2 * (    one   * q(1,i,j-2,k) &
                                                   - four  * q(1,i,j-1,k) &
                                                   + three * q(1,i,j  ,k) )
      
      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j-2,k) * rsign(i,j-1,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂η
      PressureCurvilinearGradient(2) = de2 * ( q(2,i,j+1,k) - q(2,i,j-1,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j+1,k) * rsign(i,j-1,k)  )

      if ( exsign < one_half ) return                                    

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ζ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( mydown == mpi_proc_null .and. k == ksta ) then

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3)   =  dz2 * (  - one   * q(1,i,j,k+2) &
                                                   + four  * q(1,i,j,k+1) &
                                                   - three * q(1,i,j,k  ) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k+2) * rsign(i,j,k+1) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myup == mpi_proc_null .and. k == kend) then

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3)   =  dz2 * (    one   * q(1,i,j,k-2) &
                                                   - four  * q(1,i,j,k-1) &
                                                   + three * q(1,i,j,k  ) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k-2) * rsign(i,j,k-1) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3) = dz2 * ( q(1,i,j,k+1) - q(1,i,j,k-1) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k+1) * rsign(i,j,k-1) )

      if ( exsign < one_half ) return                                    

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                                _                    _                                          
   !                               | ∂ξ/∂x  ∂ξ/∂y  ∂ξ/∂z  |  
   !               MetricsTensor = | ∂η/∂x  ∂η/∂y  ∂η/∂z  |                               
   !                               | ∂ζ/∂x  ∂ζ/∂y  ∂ζ/∂z  |  
   !                               •-                    -•  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   MetricsTensor(1,1) = csi(1,i,j,k)
   MetricsTensor(1,2) = csi(2,i,j,k)
   MetricsTensor(1,3) = csi(3,i,j,k)

   MetricsTensor(2,1) = eta(1,i,j,k)
   MetricsTensor(2,2) = eta(2,i,j,k)
   MetricsTensor(2,3) = eta(3,i,j,k)

   MetricsTensor(3,1) = zet(1,i,j,k)
   MetricsTensor(3,2) = zet(2,i,j,k)
   MetricsTensor(3,3) = zet(3,i,j,k)

   ! PressureGradient(m) = ∂p/∂xm = ( ∂p/∂ξp ) * ( ∂ξp/∂xm )

   do m = 1,3
    do p = 1,3
       
       PressureGradient(m) =   PressureGradient(m) &
                             + PressureCurvilinearGradient(p) * MetricsTensor(p,m)
    end do
   end do

end subroutine calc_pressure_gradient


subroutine calc_pressure_gradient2(i, j, k, pressure_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(out) :: pressure_gradient 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   real (kind = rdf), dimension(1:3) :: pressure_curv_gradient    
   real ( kind = rdf ), dimension(3,3) :: MetricsTensor
   integer :: m,p 

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: pLL , pL , pC , pR , pRR

   real ( kind = rdf ) :: ucon_1 , ucon_2 , ucon_3
   real ( kind = rdf ) :: up1 , um1 , up2 , um2 , up3 , um3

   integer :: bias_csi , bias_eta , bias_zet


   pressure_gradient = zero
   
   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0


   !  Contravariant velocity
   ucon_1 = csi(1,i,j,k) * q(2,i,j,k) + &
            csi(2,i,j,k) * q(3,i,j,k) + &
            csi(3,i,j,k) * q(4,i,j,k)    

   ucon_2 = eta(1,i,j,k) * q(2,i,j,k) + &
            eta(2,i,j,k) * q(3,i,j,k) + &
            eta(3,i,j,k) * q(4,i,j,k)    

   ucon_3 = zet(1,i,j,k) * q(2,i,j,k) + &
            zet(2,i,j,k) * q(3,i,j,k) + &
            zet(3,i,j,k) * q(4,i,j,k)    


   up1 = zero
   um1 = zero

   up2 = zero
   um2 = zero

   up3 = zero
   um3 = zero

   if ( ucon_1 > 0 ) then
      up1 = one
   else
      um1 = one
   end if

   if ( ucon_2 > 0 ) then
      up2 = one
   else
      um2 = one
   end if

   if ( ucon_3 > 0 ) then
      up3 = one
   else
      um3 = one
   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      pLL = zero                
      pL  = zero                
      pC  = q ( 1, i   , j , k )
      pR  = q ( 1, i+1 , j , k )
      pRR = q ( 1, i+2 , j , k )
 
      ! ∂p/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(1)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                                                      

   else if (i == iend) then

      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i-2 , j , k ) 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂u/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(1)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   else

      rL = ( rsign ( i-1 , j , k ) + abs(rsign ( i-1 , j , k )) ) / two
      rC = ( rsign ( i   , j , k ) + abs(rsign ( i   , j , k )) ) / two
      rR = ( rsign ( i+1 , j , k ) + abs(rsign ( i+1 , j , k )) ) / two

      exsign = up1 * ( rC * rL ) + um1 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(1) = up1 * dc * ( q(1,i,j,k)   - q(1,i-1,j,k)  ) + &
                                  um1 * dc * ( q(1,i+1,j,k) - q(1,i,j,k)    )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i, j    , k )  
      rR  = rsign ( i, j+1  , k )  
      rRR = rsign ( i, j+2  , k )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = q ( 1, i , j+2 , k ) 

      ! ∂p/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(2)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
   else if ( j == jend ) then

      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j-2 , k ) 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂p/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(2)      &
                                       )
   
      if ( exsign < one_half ) return                                    


   else

      rL = ( rsign ( i , j-1 , k ) + abs(rsign ( i , j-1 , k )) ) / two
      rC = ( rsign ( i , j   , k ) + abs(rsign ( i , j   , k )) ) / two
      rR = ( rsign ( i , j+1 , k ) + abs(rsign ( i , j+1 , k )) ) / two

      exsign = up2 * ( rC * rL ) + um2 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(2) = up2 * de * ( q(1,i,j,k)   - q(1,i,j-1,k)  ) + &
                                  um2 * de * ( q(1,i,j+1,k) - q(1,i,j,k)    )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then

      bias_zet = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j , k   )  
      rR  = rsign ( i , j , k+1 )  
      rRR = rsign ( i , j , k+2 )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j , k   ) 
      pR  = q ( 1, i , j , k+1 ) 
      pRR = q ( 1, i , j , k+2 ) 
 
      ! ∂p/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(3)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   else if ( k == kend ) then

      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j , k-2 ) 
      pL  = q ( 1, i , j , k-1 ) 
      pC  = q ( 1, i , j , k   ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂p/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(3)      &
                                       )
   
      if ( exsign < one_half ) return                                    

   else

      rL = ( rsign ( i , j , k-1 ) + abs(rsign ( i , j , k-1 )) ) / two
      rC = ( rsign ( i , j , k   ) + abs(rsign ( i , j , k   )) ) / two
      rR = ( rsign ( i , j , k+1 ) + abs(rsign ( i , j , k+1 )) ) / two

      exsign = up3 * ( rC * rL ) + um3 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(3) = up3 * dz * ( q(1,i,j,k)   - q(1,i,j,k-1)  ) + &
                                  um3 * dz * ( q(1,i,j,k+1) - q(1,i,j,k)  )

   end if

   MetricsTensor(1,1) = csi(1,i,j,k)
   MetricsTensor(1,2) = csi(2,i,j,k)
   MetricsTensor(1,3) = csi(3,i,j,k)

   MetricsTensor(2,1) = eta(1,i,j,k)
   MetricsTensor(2,2) = eta(2,i,j,k)
   MetricsTensor(2,3) = eta(3,i,j,k)

   MetricsTensor(3,1) = zet(1,i,j,k)
   MetricsTensor(3,2) = zet(2,i,j,k)
   MetricsTensor(3,3) = zet(3,i,j,k)

   ! PressureGradient(m) = ∂p/∂xm = ( ∂p/∂ξp ) * ( ∂ξp/∂xm )

   do m = 1,3
   do p = 1,3
       
      pressure_gradient(m) =   pressure_gradient(m) &
                             + pressure_curv_gradient(p) * MetricsTensor(p,m)
   end do
   end do



end subroutine calc_pressure_gradient2

subroutine velocity_gradient_extrapolation_free_surface_lsqm( a_coeff_vector    , alpha_vec      , & 
                                                              velocity_gradient , du_dx_fs_lsqm      )

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


subroutine pressure_gradient_extrapolation_free_surface_lsqm( a_coeff_vector_p, alpha_vec, & 
                                                              pressure_gradient, dp_dx_fs_lsqm)

   implicit none

   real (kind = rdf), dimension(1:12), intent(in) :: a_coeff_vector_p ! in matrix form for an 
                                                                         ! easier handling 

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! xneasrest
   real (kind = rdf), dimension(1:3), intent(in) :: pressure_gradient ! velocity_gradient at
                                                                          ! inearest, jnearest, knearest

   real (kind = rdf), dimension(1:3), intent(inout) :: dp_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface

   ! local variables

   integer :: igfs, jgfs

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                       __        __
   !                      |  a1 a2 a3  |
   !                      |  b1 b2 b3  |
   !  a_coeff_vector_p =  |  c1 c2 c3  |
   !                      |  d1 d2 d3  |
   !                       --        --
   !
   !
   ! aij is the coefficient for extrapolating ∂p_i/∂x_j
   !
   !
   ! dp_dx_fs_lsqm : cartesian pressure gradient at the free-surface computed
   ! using the coefficients obtained by lsqm.
   ! 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! xs : free-surface location
   !
   ! ∂p/∂x_i (xs) = ∂p/∂x_i (xs + α) - f_i
   ! Where ∂u_i/∂x_j (xs + α) is the velocity gradient at the nearest water-phase node
   ! with a valid computed gradient (using the rules defined in the velocity_curv_gradient_tensor
   ! subroutine)


   dp_dx_fs_lsqm(1) = pressure_gradient(1) - ( a_coeff_vector_p(1) + a_coeff_vector_p(2)  * alpha_vec(1) &
                                                                   + a_coeff_vector_p(3)  * alpha_vec(2) &
                                                                   + a_coeff_vector_p(4)  * alpha_vec(3) )


   dp_dx_fs_lsqm(2) = pressure_gradient(2) - ( a_coeff_vector_p(5) + a_coeff_vector_p(6)  * alpha_vec(1) &
                                                                   + a_coeff_vector_p(7)  * alpha_vec(2) &
                                                                   + a_coeff_vector_p(8)  * alpha_vec(3) )

   dp_dx_fs_lsqm(3) = pressure_gradient(3) - ( a_coeff_vector_p(9) + a_coeff_vector_p(10) * alpha_vec(1) &
                                                                   + a_coeff_vector_p(11) * alpha_vec(2) &
                                                                   + a_coeff_vector_p(12) * alpha_vec(3) )


end subroutine pressure_gradient_extrapolation_free_surface_lsqm


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

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

   ! Interior nodes including domain boundaries
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! extp_swept = # of nodes around i,j,k to be candidates for extrapolation
   extp_swept = 3 ! TO DO: it could be set from control.dat

   ! max(k_mysta-1,k-extp_swept),min(k_myend+1,k+extp_swept) means that it even extrapolates
   ! the velocity to the air nodes at the boundaries of the domain. This is gonna be useful to
   ! compute velocity gradients in the air phase for the normal boundary condition.

   do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

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


subroutine ghost_nodes_extrapolation(i,j,k, xs, ys, zs, u_fs_lsqm, du_dx_fs_lsqm, dp_dx_fs_lsqm)

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(in) :: u_fs_lsqm
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm
   real (kind = rdf), dimension(1:3), intent(in) :: dp_dx_fs_lsqm
   real (kind = rdf), intent(in) :: xs, ys, zs ! free-surface position closest to (i,j,k) node

   ! local variables
   real (kind = rdf), dimension(1:3) :: alpha_local
   real (kind = rdf) :: uextp, vextp, wextp ! extrapolated velocity using Taylor expansion
   real (kind = rdf) :: pextp
   real (kind = rdf) :: rdiff_norm, exsign  

   integer :: ii, jj, kk, extp_swept

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

  !Interior nodes including domain boundaries
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! extp_swept = # of nodes around i,j,k to be candidates for extrapolation
   extp_swept = 3 ! TO DO: it could be set from control.dat

   ! max(k_mysta-1,k-extp_swept),min(k_myend+1,k+extp_swept) means that it even extrapolates
   ! the velocity to the air nodes at the boundaries of the domain. This is gonna be useful to
   ! compute velocity gradients in the air phase for the normal boundary condition.

   !do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   !do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   !do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

   !do kk = max( kl , k-extp_swept ) , min( ku , k+extp_swept ) 
   !do jj = max( jl , j-extp_swept ) , min( ju , j+extp_swept )
   !do ii = max( il , i-extp_swept ) , min( iu , i+extp_swept )

   do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

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

         ! pressure extrapolation assuming pfs = 0
         pextp = zero

         pextp = dp_dx_fs_lsqm(1) * alpha_local(1) + &
                 dp_dx_fs_lsqm(2) * alpha_local(2) + &
                 dp_dx_fs_lsqm(3) * alpha_local(3)

         ! extrapolated velocities using Taylor expansion
         uextp = zero; vextp = zero; wextp = zero;

         uextp = u_fs_lsqm(1) +  du_dx_fs_lsqm(1,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(1,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(1,3) * alpha_local(3)

         vextp = u_fs_lsqm(2) +  du_dx_fs_lsqm(2,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(2,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(2,3) * alpha_local(3)

         wextp = u_fs_lsqm(3) +  du_dx_fs_lsqm(3,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(3,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(3,3) * alpha_local(3)

         ! if exsign = 1 (least distanced node), the velocity at ii,jj,kk
         ! is replaced by the extrapolated one

         if ( exsign > one_half ) then
      
            ! if rsign = -1, then I keep the geometrically extrapolated value
            ! if rsign =  0, then I use the LSQM extrapolation
            q(1,ii,jj,kk) =   q(1,ii,jj,kk)  &
                            + ( one + rsign(ii,jj,kk) ) * exsign * ( - q(1,ii,jj,kk) + pextp )

            ! Different from the one above, I will keep the LSQM extrapolation even for the 
            ! nondes next to the free surface
            !q(1,ii,jj,kk) = q(1,ii,jj,kk) + exsign * ( - q(1,ii,jj,kk) + pextp )
            q(2,ii,jj,kk) = q(2,ii,jj,kk) + exsign * ( - q(2,ii,jj,kk) + uextp )
            q(3,ii,jj,kk) = q(3,ii,jj,kk) + exsign * ( - q(3,ii,jj,kk) + vextp )
            q(4,ii,jj,kk) = q(4,ii,jj,kk) + exsign * ( - q(4,ii,jj,kk) + wextp )

         end if
      
      end if

   end do
   end do
   end do

end subroutine ghost_nodes_extrapolation


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