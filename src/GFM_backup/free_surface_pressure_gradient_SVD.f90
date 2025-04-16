
subroutine free_surface_pressure_gradient_SVD( i , j , k , xs, ys, zs, pfs , nvec , pgrad_fs )

   !--------------------------------------------------------------------------
   ! Description: 
   ! ------------
   ! 
   ! Calculates the pressure gradient at the free surface (xs,ys,zs) using a 
   ! Weighted Least Square Method (WLSM). 
   ! 
   ! the weights are a combination of distance and alignment with the normal 
   ! direction, such that nodes closer and/or more normal to the free surface
   ! will contribute more to the computation of the gradient
   !
   !--------------------------------------------------------------------------

   use InterpolationMethods 

   implicit none

   ! Input
   integer , intent(in) :: i,j,k ! Node in the air phase to have the pressure extrapolated
   real(kind = rdf) , intent(in) :: xs, ys , zs ! Free surface location
   real(kind = rdf) , intent(in) :: pfs ! Pressure at the free surface obtained from the NDBC
   real(kind = rdf), dimension(1:3), intent(in) :: nvec ! normal vector at the free surface

   ! Output
   real(kind = rdf), dimension(1:3), intent(out) :: pgrad_fs ! pressure gradient at the free surface


   ! Local variables

   ! SVD system
   real(kind = rdf), allocatable :: Rw_aux(:,:) ! M-matrix WLSM
   real(kind = rdf), allocatable :: bwvec_aux(:) ! rhs vector
   real(kind = rdf), allocatable :: Rw(:,:) ! M-matrix WLSM
   real(kind = rdf), allocatable :: bwvec(:) ! rhs vector
   real(kind = rdf), allocatable :: s(:) ! rhs vector
   real(kind = rdf), allocatable :: work(:) ! rhs vector
   real (kind=rdf) :: lw(1)
   integer :: m,n ! matrix dimensions
   integer, parameter :: nrhs = 1

   ! Local Scalars 
   real (Kind=rdf), parameter :: rcond = 0.00001_rdf

   ! Routine internal variables
   integer :: info, lda, ldb, liwork, lwork, rank
   
   ! More routine interanl variable that has to be set before calling
   integer, allocatable :: iwork(:)
   integer :: liw(1)
   integer :: smlsiz, nlvl


   integer :: ii,jj,kk
   real(kind = rdf) :: wlocal ! local weight
   real(kind = rdf) :: dx , dy , dz , dp 
   real(kind = rdf), dimension(1:3) :: dr ! Δr = rii - rfs
   real(kind = rdf) :: num , denom 
   real(kind = rdf) :: aux 
   real(kind=rdf), parameter :: eps_denom = 1.0E-6
   logical :: OK_FLAG

   real(kind=rdf) :: dsmin
   real(kind=rdf), parameter :: dsmin_fraction = 0.5

   integer :: counter


   pgrad_fs = zero

   ! Let's get the minimum ∆s around (i,j,k) (where ∆s can be
   ! either ∆x, ∆y or ∆z)

   ! Initialisation
   dsmin = ten * ten

   do kk = max( ksta , k-1 ), min( kend , k+1 )
   do jj = max( jsta , j-1 ), min( jend , j+1 )
   do ii = max( ista , i-1 ), min( iend , i+1 )

      if ( rsign(ii,jj,kk) < one_half ) cycle

      dx = x(i,j,k) - x(ii,jj,kk)
      dy = y(i,j,k) - y(ii,jj,kk)
      dz = z(i,j,k) - z(ii,jj,kk)

      dsmin = minval( (/ dsmin , abs(dx) , abs(dy) , abs(dz) /) )

   end do
   end do
   end do

   counter = 0

   allocate( Rw_aux(1:300,1:3) , bwvec_aux(1:300) )

   do kk = max( ksta , k-sweep_lsqm ), min( kend , k+sweep_lsqm )
   do jj = max( jsta , j-sweep_lsqm ), min( jend , j+sweep_lsqm )
   do ii = max( ista , i-sweep_lsqm ), min( iend , i+sweep_lsqm )

      if ( rsign(ii,jj,kk) < one_half ) cycle

      dx = x(ii,jj,kk) - xs
      dy = y(ii,jj,kk) - ys
      dz = z(ii,jj,kk) - zs

      dr = (/dx,dy,dz/)

      ! If the distance between the point ii,jj,kk and the fs is too
      ! small, then that gradient might be spurious
      if( norm2(dr) < dsmin_fraction * dsmin ) cycle

      ! Update the counter of the number of rows of Rw
      counter = counter + 1

      dp = q(1,ii,jj,kk) - pfs

      num   = abs( dot_product( dr/norm2(dr) , nvec ) )**(three/two) 
      denom = one !norm2(dr)**(three/two) + eps_denom

      wlocal = num / denom

      Rw_aux(counter,1) = wlocal * dx 
      Rw_aux(counter,2) = wlocal * dy 
      Rw_aux(counter,3) = wlocal * dz 

      bwvec_aux(counter) = wlocal * dp

   end do
   end do
   end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    ! All these variables are set following the guidelines from the routine documentation
    ! (in general, don't change)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    ! matrix dimensions
    m = counter
    n = 3

    smlsiz = 25
    nlvl   = max(0,int( log(real(min(m,n))/real(smlsiz+1))/log(2.) )+1)
    
    if(m.ge.n) then      
      lwork  = 12*n+2*n*smlsiz+8*n*nlvl+ n*nrhs+(smlsiz+1)*(smlsiz+1)
    else
      lwork  = 12*m+2*m*smlsiz+8*m*nlvl+ m*nrhs+(smlsiz+1)*(smlsiz+1)
    end if

    liwork = max(1,3*min(m,n)*nlvl+11*min(m,n))

    lda = max(1,m)
    ldb = max(1,max(m,n))

    allocate (work(lwork), iwork(liwork))
    allocate (Rw(lda,n), bwvec(max(m,n)), s(m))
    
    ! Data mapping
    Rw(1:m,1:3) = Rw_aux(1:m,1:3) 
    bwvec(1:m)  = bwvec_aux(1:m)

    deallocate( Rw_aux , bwvec_aux )

    call dgelsd(m, n, nrhs, Rw, lda, bwvec, ldb, s, rcond, rank, work, lwork, iwork, &
                info)

    pgrad_fs(1) = bwvec(1)
    pgrad_fs(2) = bwvec(2)
    pgrad_fs(3) = bwvec(3)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   

    deallocate(Rw)
    deallocate(bwvec)
    deallocate(s)


   !print *, ' '
   !print *, 'i,j,k: ',i,j,k,' , dp_dx theo  = ', aux , ' , pfs = ', pfs ,' , dp_dx WLSQM =' , pgrad_fs(1)    

end subroutine free_surface_pressure_gradient_SVD