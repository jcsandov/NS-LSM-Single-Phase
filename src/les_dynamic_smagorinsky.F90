
subroutine les_dynamic_smagorinsky ( il, iu, jl, ju, kl, ku, &   
                                     igp, jgp, kgp,          & 
                                     dc, de, dz,             & 
                                     x,y,z,                  & 
                                     q,                      & 
                                     csi,                    & 
                                     eta,                    & 
                                     zet,                    & 
                                     aj,                     &        
                                     xnut)                                  

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Computes the turbulent viscosity (xnut) using the dynamic Smagorinsky  
  ! method. 
  ! 
  ! It receives the resolved velocity field, and returns xnut
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use global_param
  use global_app
  use global_mpi
  use AdvectionMethods

  implicit none

  !Inputs
  integer :: il,iu,jl,ju,kl,ku
  integer :: igp,jgp,kgp
  real (kind = rdf) :: dc,de,dz
  real (kind = rdf), dimension( 1:5 , il:iu , jl:ju , kl:ku ) :: q
  real (kind = rdf), dimension( 1:3 , il:iu , jl:ju , kl:ku ) :: csi , eta, zet
  real (kind = rdf), dimension(       il:iu , jl:ju , kl:ku ) :: aj, xnut
  real (kind = rdf), dimension(       il:iu , jl:ju , kl:ku ) :: x, y, z

  ! local variables
  !index
  integer :: i_mysta , i_myend
  integer :: j_mysta , j_myend
  integer :: k_mysta , k_myend
  integer :: ista , iend , jsta , jend , ksta , kend
  integer :: i, j, k
  integer :: a, b

  !Constants
  real (kind = rdf), parameter :: cul  = 0.05_rdf
  
  real (kind = rdf) :: s11, s22, s33
  real (kind = rdf) :: s12, s23, s13

  real (kind = rdf) :: du_dcsi , dv_dcsi , dw_dcsi
  real (kind = rdf) :: du_deta , dv_deta , dw_deta
  real (kind = rdf) :: du_dzet , dv_dzet , dw_dzet
  real (kind = rdf) :: du_dx   , du_dy   , du_dz
  real (kind = rdf) :: dv_dx   , dv_dy   , dv_dz
  real (kind = rdf) :: dw_dx   , dw_dy   , dw_dz
  
  real (kind = rdf) :: dc2, de2, dz2
  real (kind = rdf) :: suma, suma1, contador
  real (kind = rdf) :: grid_filter, test_filter

  real (kind = rdf) :: VelocityGradient(3,3)

  real (kind = rdf), dimension(:,:,:,:)  , allocatable :: S1ij, S2ij, S3ij, SS1ij, SS2ij, SS3ij
  real (kind = rdf), dimension(:,:,:,:,:), allocatable :: Sij_hat, SSij_hat
  real (kind = rdf), dimension(:,:,:,:,:), allocatable :: Lij, Mij
  real (kind = rdf), dimension(:,:,:)    , allocatable :: u_f, v_f, w_f, Sabs
  real (kind = rdf), dimension(:,:,:)    , allocatable :: uu, uv, uw, vv, vw, ww
  real (kind = rdf), dimension(:,:,:)    , allocatable :: uu_f, uv_f, uw_f, vv_f, vw_f, ww_f 
  real (kind = rdf), dimension(:,:,:)    , allocatable :: uuf_p, uvf_p, uwf_p, vvf_p, vwf_p, wwf_p
  real (kind = rdf), dimension(:,:,:)    , allocatable :: Sabs_hat
  real (kind = rdf), dimension(:,:,:)    , allocatable :: num, den, Cs
  real (kind = rdf), dimension(:,:,:)    , allocatable :: num_f, den_f

  !Allocate arrays
  allocate ( S1ij     (       1:3 , il:iu , jl:ju , kl:ku )  , &
             S2ij     (       1:3 , il:iu , jl:ju , kl:ku )  , &
             S3ij     (       1:3 , il:iu , jl:ju , kl:ku )  , &
             Sij_hat  ( 1:3 , 1:3 , il:iu , jl:ju , kl:ku )  , &
             SS1ij    (       1:3 , il:iu , jl:ju , kl:ku )  , &
             SS2ij    (       1:3 , il:iu , jl:ju , kl:ku )  , &
             SS3ij    (       1:3 , il:iu , jl:ju , kl:ku )  , &
             SSij_hat ( 1:3 , 1:3 , il:iu , jl:ju , kl:ku )  , &
             Lij      ( 1:3 , 1:3 , il:iu , jl:ju , kl:ku )  , &
             Mij      ( 1:3 , 1:3 , il:iu , jl:ju , kl:ku )    &
          )

  allocate (      u_f ( il:iu , jl:ju , kl:ku ) ,   v_f ( il:iu , jl:ju , kl:ku ) , &
                  w_f ( il:iu , jl:ju , kl:ku ) ,  Sabs ( il:iu , jl:ju , kl:ku ) )

  allocate (       uu ( il:iu , jl:ju , kl:ku ) ,    uv ( il:iu , jl:ju , kl:ku ) , &
                   uw ( il:iu , jl:ju , kl:ku ) ,    vv ( il:iu , jl:ju , kl:ku ) , &
                   vw ( il:iu , jl:ju , kl:ku ) ,    ww ( il:iu , jl:ju , kl:ku ) )

  allocate (     uu_f ( il:iu , jl:ju , kl:ku ) ,  uv_f ( il:iu , jl:ju , kl:ku ) , &
                 uw_f ( il:iu , jl:ju , kl:ku ) ,  vv_f ( il:iu , jl:ju , kl:ku ) , &
                 vw_f ( il:iu , jl:ju , kl:ku ) ,  ww_f ( il:iu , jl:ju , kl:ku ) )

  allocate (    uuf_p ( il:iu , jl:ju , kl:ku ) , uvf_p ( il:iu , jl:ju , kl:ku ) , &
                uwf_p ( il:iu , jl:ju , kl:ku ) , vvf_p ( il:iu , jl:ju , kl:ku ) , &
                vwf_p ( il:iu , jl:ju , kl:ku ) , wwf_p ( il:iu , jl:ju , kl:ku ) )

  allocate ( Sabs_hat ( il:iu , jl:ju , kl:ku ) ,    Cs ( il:iu , jl:ju , kl:ku ) , &
                  num ( il:iu , jl:ju , kl:ku ) , num_f ( il:iu , jl:ju , kl:ku ) , &
                  den ( il:iu , jl:ju , kl:ku ) , den_f ( il:iu , jl:ju , kl:ku ) )

  ! Variables initialisation

  u_f       = zero
  v_f       = zero
  w_f       = zero 
  Sabs      = zero 
  uu        = zero 
  uv        = zero 
  uw        = zero 
  vv        = zero 
  vw        = zero 
  ww        = zero 
  uu_f      = zero 
  uv_f      = zero 
  uw_f      = zero 
  vv_f      = zero 
  vw_f      = zero 
  ww_f      = zero 
  uuf_p     = zero 
  uvf_p     = zero 
  uwf_p     = zero 
  vvf_p     = zero 
  vwf_p     = zero 
  wwf_p     = zero 
  Sabs_hat  = zero
  Cs        = zero
  num       = zero
  den       = zero
  num_f     = zero
  den_f     = zero

  dc2 = pt5 * dc
  de2 = pt5 * de
  dz2 = pt5 * dz

  ! These interior nodes include the boundaries
  i_mysta = il + igp 
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp        
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! Physical boundaries

  ista = il ; jsta = jl ; ksta = kl 
  iend = iu ; jend = ju ; kend = ku 

  if ( myback  == mpi_proc_null )  ista = il + igp 
  if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
  if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

  if ( myfront == mpi_proc_null )  iend = iu - igp
  if ( myright == mpi_proc_null )  jend = ju - jgp
  if ( myup    == mpi_proc_null )  kend = ku - kgp


  !S tensor
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    ! Velocity gradient

    call VelocityGradientTensor( i , j , k , VelocityGradient )

    du_dx = VelocityGradient(1,1)
    du_dy = VelocityGradient(1,2)
    du_dz = VelocityGradient(1,3)

    dv_dx = VelocityGradient(2,1)
    dv_dy = VelocityGradient(2,2)
    dv_dz = VelocityGradient(2,3)

    dw_dx = VelocityGradient(3,1)
    dw_dy = VelocityGradient(3,2)
    dw_dz = VelocityGradient(3,3)

    !componentes del tensor S  en cada nodo
    !
    !          1   /  ∂ui     ∂uj  \
    !  Sij =  --- (  ----- + -----  )      
    !          2   \  ∂xj     ∂xi  /
    
    s11 = du_dx
    s22 = dv_dy
    s33 = dw_dz

    s12 = pt5 * ( du_dy + dv_dx )
    s13 = pt5 * ( du_dz + dw_dx )
    s23 = pt5 * ( dv_dz + dw_dy )

    !Tensor Sij
    S1ij(1,i,j,k) = s11 
    S1ij(2,i,j,k) = s12 
    S1ij(3,i,j,k) = s13 

    S2ij(1,i,j,k) = s12
    S2ij(2,i,j,k) = s22
    S2ij(3,i,j,k) = s23

    S3ij(1,i,j,k) = s13
    S3ij(2,i,j,k) = s23
    S3ij(3,i,j,k) = s33

    ! mean strain rate
    !          ___________
    ! | S | = √ 2 Sij Sij

    Sabs(i,j,k) = sqrt( two * (         ( s11 * s11 + s22 * s22 + s33 * s33 )  &
                                + two * ( s12 * s12 + s13 * s13 + s23 * s23 )    )  )


    !Tensor SSij = Sij|S|
    SS1ij(1,i,j,k) = s11 * Sabs(i,j,k)
    SS1ij(2,i,j,k) = s12 * Sabs(i,j,k)
    SS1ij(3,i,j,k) = s13 * Sabs(i,j,k)
    
    SS2ij(1,i,j,k) = s12 * Sabs(i,j,k)
    SS2ij(2,i,j,k) = s22 * Sabs(i,j,k)
    SS2ij(3,i,j,k) = s23 * Sabs(i,j,k)
    
    SS3ij(1,i,j,k) = s13 * Sabs(i,j,k)
    SS3ij(2,i,j,k) = s23 * Sabs(i,j,k)
    SS3ij(3,i,j,k) = s33 * Sabs(i,j,k)

  end do
  end do
  end do

  !Echange ghost point of Sabs, Sij and SSij in order to filter those variables
  
  call rhs_exchng3_3d ( Sabs )

  call rhs_exchng3_4d ( S1ij )
  call rhs_exchng3_4d ( S2ij )
  call rhs_exchng3_4d ( S3ij )

  call rhs_exchng3_4d ( SS1ij )
  call rhs_exchng3_4d ( SS2ij )
  call rhs_exchng3_4d ( SS3ij )

  ! ------------------------------------------------------------------------
  ! Mij tensor components
  !
  !       ___   ___   ___         _________
  ! Mij = ∆^2 * Sij * |S| - ∆^2 * Sij * |S|

  ! ___
  ! Sij
  Sij_hat(1,1,:,:,:) = testfilter_simpson_3d ( S1ij(1,:,:,:) , aj )
  Sij_hat(1,2,:,:,:) = testfilter_simpson_3d ( S1ij(2,:,:,:) , aj )
  Sij_hat(1,3,:,:,:) = testfilter_simpson_3d ( S1ij(3,:,:,:) , aj )

  Sij_hat(2,1,:,:,:) = testfilter_simpson_3d ( S2ij(1,:,:,:) , aj )
  Sij_hat(2,2,:,:,:) = testfilter_simpson_3d ( S2ij(2,:,:,:) , aj )
  Sij_hat(2,3,:,:,:) = testfilter_simpson_3d ( S2ij(3,:,:,:) , aj )
  
  Sij_hat(3,1,:,:,:) = testfilter_simpson_3d ( S3ij(1,:,:,:) , aj )
  Sij_hat(3,2,:,:,:) = testfilter_simpson_3d ( S3ij(2,:,:,:) , aj )
  Sij_hat(3,3,:,:,:) = testfilter_simpson_3d ( S3ij(3,:,:,:) , aj )

  ! ___
  ! |S|
  Sabs_hat = testfilter_simpson_3d( Sabs , aj )

  ! ____   _________
  ! SSij = Sij * |S|
  SSij_hat(1,1,:,:,:) = testfilter_simpson_3d ( SS1ij (1,:,:,:) , aj )
  SSij_hat(1,2,:,:,:) = testfilter_simpson_3d ( SS1ij (2,:,:,:) , aj )
  SSij_hat(1,3,:,:,:) = testfilter_simpson_3d ( SS1ij (3,:,:,:) , aj )

  SSij_hat(2,1,:,:,:) = testfilter_simpson_3d ( SS2ij (1,:,:,:) , aj )
  SSij_hat(2,2,:,:,:) = testfilter_simpson_3d ( SS2ij (2,:,:,:) , aj )
  SSij_hat(2,3,:,:,:) = testfilter_simpson_3d ( SS2ij (3,:,:,:) , aj )

  SSij_hat(3,1,:,:,:) = testfilter_simpson_3d ( SS3ij (1,:,:,:) , aj )
  SSij_hat(3,2,:,:,:) = testfilter_simpson_3d ( SS3ij (2,:,:,:) , aj )
  SSij_hat(3,3,:,:,:) = testfilter_simpson_3d ( SS3ij (3,:,:,:) , aj )

  ! ------------------------------------------------------------------------
  ! Lij tensor components
  !
  !       _______     __   __ 
  ! Lij = ui * uj  -  ui * uj

  ! _   _   _
  ! u , v , w
  u_f = testfilter_simpson_3d( q(2,:,:,:) , aj )
  v_f = testfilter_simpson_3d( q(3,:,:,:) , aj )
  w_f = testfilter_simpson_3d( q(4,:,:,:) , aj )

  ! ui * uj 
  uu(:,:,:) = q(2,:,:,:) * q(2,:,:,:)
  uv(:,:,:) = q(2,:,:,:) * q(3,:,:,:)
  uw(:,:,:) = q(2,:,:,:) * q(4,:,:,:) 
  vv(:,:,:) = q(3,:,:,:) * q(3,:,:,:)
  vw(:,:,:) = q(3,:,:,:) * q(4,:,:,:)
  ww(:,:,:) = q(4,:,:,:) * q(4,:,:,:)

  !  _______ 
  !  ui * uj   
  uu_f = testfilter_simpson_3d ( uu , aj )
  uv_f = testfilter_simpson_3d ( uv , aj )
  uw_f = testfilter_simpson_3d ( uw , aj )
  vv_f = testfilter_simpson_3d ( vv , aj )
  vw_f = testfilter_simpson_3d ( vw , aj )
  ww_f = testfilter_simpson_3d ( ww , aj )

  ! __   __
  ! ui * uj 
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    uuf_p(i,j,k) = u_f(i,j,k) * u_f(i,j,k)
    uvf_p(i,j,k) = u_f(i,j,k) * v_f(i,j,k)
    uwf_p(i,j,k) = u_f(i,j,k) * w_f(i,j,k)
    vvf_p(i,j,k) = v_f(i,j,k) * v_f(i,j,k)
    vwf_p(i,j,k) = v_f(i,j,k) * w_f(i,j,k)
    wwf_p(i,j,k) = w_f(i,j,k) * w_f(i,j,k)
    
  end do
  end do
  end do

  ! Leonard tensor 
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    Lij(1,1,i,j,k) = uu_f(i,j,k) - uuf_p(i,j,k)
    Lij(1,2,i,j,k) = uv_f(i,j,k) - uvf_p(i,j,k)
    Lij(1,3,i,j,k) = uw_f(i,j,k) - uwf_p(i,j,k)

    Lij(2,1,i,j,k) = uv_f(i,j,k) - uvf_p(i,j,k)
    Lij(2,2,i,j,k) = vv_f(i,j,k) - vvf_p(i,j,k)
    Lij(2,3,i,j,k) = vw_f(i,j,k) - vwf_p(i,j,k)
    
    Lij(3,1,i,j,k) = uw_f(i,j,k) - uwf_p(i,j,k)
    Lij(3,2,i,j,k) = vw_f(i,j,k) - vwf_p(i,j,k)
    Lij(3,3,i,j,k) = ww_f(i,j,k) - wwf_p(i,j,k)

  end do
  end do
  end do
 
  ! I need num and den at the boundaries, because I need to filter them
  ! afterwards
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    !Filters
    grid_filter = aj(i,j,k)**(-one_third)
    test_filter = two * grid_filter

    do a = 1 , 3
    do b = 1 , 3

      Mij(a,b,i,j,k) = two * ( ( grid_filter**two ) * SSij_hat(a,b,i,j,k)  -  &
                               ( test_filter**two ) * Sabs_hat(i,j,k) * Sij_hat(a,b,i,j,k) ) 
        
      num(i,j,k) = num(i,j,k) + Lij(a,b,i,j,k) * Mij(a,b,i,j,k)
      den(i,j,k) = den(i,j,k) + Mij(a,b,i,j,k) * Mij(a,b,i,j,k)

    end do
    end do

  end do
  end do
  end do

  !Ghost Points needed for the filtering
  call rhs_exchng3_3d ( num )
  call rhs_exchng3_3d ( den )

  !Promedio espacial para evitar singularidades
  !Promedio usando el filtro test
  num_f = testfilter_simpson_3d ( num , aj )
  den_f = testfilter_simpson_3d ( den , aj )

  ! Constante de Smagorinsky Dinámico
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    if ( abs( den_f(i,j,k) ) < eps_sims) then 
      
      Cs(i,j,k) = zero
    
    else

      Cs(i,j,k) = num_f(i,j,k) / den_f(i,j,k)

    end if

  end do
  end do
  end do

  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    if ( Cs(i,j,k) > 0.05_rdf ) then 
      Cs(i,j,k) = cul
    end if

    if ( Cs(i,j,k) < -0.05_rdf ) then 
      Cs(i,j,k) = -cul
    end if
  
  end do
  end do
  end do
  

  !Viscosidad turbulenta
  do k = k_mysta , k_myend
  do j = j_mysta , j_myend
  do i = i_mysta , i_myend

    grid_filter = aj(i,j,k)**(-one_third)

    ! νt = Cs * ∆^2 * |S|
    xnut(i,j,k) = Cs(i,j,k) * ( grid_filter**two ) * Sabs(i,j,k) 

  end do
  end do
  end do

  ! Update the ghost nodes
  call rhs_exchng3_3d ( xnut )  

  !Bordes
  call les_bcond()

  ! Update the ghost nodes of xnut
  call rhs_exchng3_3d ( xnut )


  !!Save Cs values
  !do k = k_mysta , k_myend
  !do j = j_mysta , j_myend
  !do i = i_mysta , i_myend
  !
  !  q(5,i,j,k) = Cs(i,j,k)
  !
  !end do
  !end do
  !end do

  ! deallocate variables
  deallocate ( S1ij, S2ij, S3ij, Sij_hat, SS1ij, SS2ij, SS3ij, SSij_hat, Lij, Mij )
  deallocate ( u_f, v_f, w_f, Sabs )
  deallocate ( uu, uv, uw, vv, vw, ww )
  deallocate ( uu_f, uv_f, uw_f, vv_f, vw_f, ww_f )
  deallocate ( uuf_p, uvf_p, uwf_p, vvf_p, vwf_p, wwf_p )
  deallocate ( Sabs_hat, Cs, num, num_f, den, den_f )
  
contains

  include 'testfilter_simpson_3d.F90'
  include 'les_bcond.F90'
  include 'rhs_exchng3_4d.F90'
  include 'rhs_exchng3_3d.F90'
  include 'VelocityGradientTensor.F90'

end subroutine les_dynamic_smagorinsky

