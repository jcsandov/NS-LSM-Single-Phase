
subroutine mg_metrics

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! General curvilinear coordinates
   !
   ! for mp program
   !
   ! Calculate metrics of the geometric transformation for
   ! for the fine grid and each of the coarse grids
   !
   ! input
   !     x(nmxg)
   !     y(nmxg)
   !     z(nmxg)
   !
   ! output
   !     csi(1:3,nmxg) csi_x, csi_y and csi_z
   !     eta(1:3,nmxg) eta_x, eta_y and eta_z
   !     zet(1:3,nmxg) zet_x, zet_y and zet_z
   !      aj(nmxg)   Jacobian
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
   use global
   use global_param
   use global_app
   use global_mpi

   implicit none

   ! local variables
   real (kind = rdf) :: gj, dum, dum1
   real (kind = rdf) :: Sc, Se, Sz, SL

   integer :: n, i, j, k, l, p, m
   integer :: ic, ie, iz

   integer :: ka
   integer :: kb
   integer :: ja
   integer :: jb
   integer :: ia
   integer :: ib


   n = 1

   !------------------------------------------------------------------------
   !                          ξ - DIRECTION
   !------------------------------------------------------------------------

   dum1 = dc(n)
   dum  = one_half * dc(n)

   ia = le_ia(1)+1
   ib = le_ib(1)-1
   ja = le_ja(1)
   jb = le_jb(1)
   ka = le_ka(1)
   kb = le_kb(1)

   if ( myback  == mpi_proc_null ) ia = li_ia(1) + 1
   if ( myfront == mpi_proc_null ) ib = li_ib(1) - 1

   do k = ka , kb
   do j = ja , jb
   do i = ia , ib

      xc( le_idx(i,j,k,1) ) = dum * ( x( le_idx(i+1,j,k,1) ) - x( le_idx(i-1,j,k,1) ) )
      yc( le_idx(i,j,k,1) ) = dum * ( y( le_idx(i+1,j,k,1) ) - y( le_idx(i-1,j,k,1) ) )
      zc( le_idx(i,j,k,1) ) = dum * ( z( le_idx(i+1,j,k,1) ) - z( le_idx(i-1,j,k,1) ) ) 

   end do
   end do
   end do

   ! back boundary - 1st order derivative

   i = ia - 1

   do k = ka , kb
   do j = ja , jb

      xc( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i+1,j,k,1) ) - x( le_idx(i,j,k,1) ) )
      yc( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i+1,j,k,1) ) - y( le_idx(i,j,k,1) ) )
      zc( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i+1,j,k,1) ) - z( le_idx(i,j,k,1) ) ) 

   end do
   end do

   ! front boundary - 1st order derivative

   i = ib + 1

   do k = ka , kb
   do j = ja , jb

      xc( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i,j,k,1) ) - x( le_idx(i-1,j,k,1) ) )
      yc( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i,j,k,1) ) - y( le_idx(i-1,j,k,1) ) )
      zc( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i,j,k,1) ) - z( le_idx(i-1,j,k,1) ) ) 

   end do
   end do


   !------------------------------------------------------------------------
   !                          η - DIRECTION
   !------------------------------------------------------------------------

   dum1 = de(n)
   dum  = one_half * de(n)

   ia = le_ia(1)
   ib = le_ib(1)
   ja = le_ja(1)+1
   jb = le_jb(1)-1
   ka = le_ka(1)
   kb = le_kb(1)

   if ( myleft  == mpi_proc_null ) ja = li_ja(1) + 1
   if ( myright == mpi_proc_null ) jb = li_jb(1) - 1

   do k = ka , kb
   do j = ja , jb
   do i = ia , ib

      xe( le_idx(i,j,k,1) ) = dum * ( x( le_idx(i,j+1,k,1) ) - x( le_idx(i,j-1,k,1) ) )
      ye( le_idx(i,j,k,1) ) = dum * ( y( le_idx(i,j+1,k,1) ) - y( le_idx(i,j-1,k,1) ) )
      ze( le_idx(i,j,k,1) ) = dum * ( z( le_idx(i,j+1,k,1) ) - z( le_idx(i,j-1,k,1) ) ) 

   end do
   end do
   end do

   ! left boundary - 1st order derivative

   j = ja - 1

   do k = ka , kb
   do i = ia , ib

      xe( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i,j+1,k,1) ) - x( le_idx(i,j,k,1) ) )
      ye( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i,j+1,k,1) ) - y( le_idx(i,j,k,1) ) )
      ze( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i,j+1,k,1) ) - z( le_idx(i,j,k,1) ) ) 

   end do
   end do

   ! right boundary - 1st order derivative

   j = jb + 1

   do k = ka , kb
   do i = ia , ib

      xe( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i,j,k,1) ) - x( le_idx(i,j-1,k,1) ) )
      ye( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i,j,k,1) ) - y( le_idx(i,j-1,k,1) ) )
      ze( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i,j,k,1) ) - z( le_idx(i,j-1,k,1) ) ) 

   end do
   end do


   !------------------------------------------------------------------------
   !                          ζ - DIRECTION
   !------------------------------------------------------------------------

   dum1 = dz(n)
   dum  = one_half * dz(n)

   ia = le_ia(1)
   ib = le_ib(1)
   ja = le_ja(1)
   jb = le_jb(1)
   ka = le_ka(1)+1
   kb = le_kb(1)-1

   if ( mydown == mpi_proc_null ) ka = li_ka(1) + 1
   if ( myup   == mpi_proc_null ) kb = li_kb(1) - 1

   do k = ka , kb
   do j = ja , jb
   do i = ia , ib

      xz( le_idx(i,j,k,1) ) = dum * ( x( le_idx(i,j,k+1,1) ) - x( le_idx(i,j,k-1,1) ) )
      yz( le_idx(i,j,k,1) ) = dum * ( y( le_idx(i,j,k+1,1) ) - y( le_idx(i,j,k-1,1) ) )
      zz( le_idx(i,j,k,1) ) = dum * ( z( le_idx(i,j,k+1,1) ) - z( le_idx(i,j,k-1,1) ) ) 

   end do
   end do
   end do

   ! down boundary - 1st order derivative

   k = ka - 1

   do j = ja , jb
   do i = ia , ib

      xz( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i,j,k+1,1) ) - x( le_idx(i,j,k,1) ) )
      yz( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i,j,k+1,1) ) - y( le_idx(i,j,k,1) ) )
      zz( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i,j,k+1,1) ) - z( le_idx(i,j,k,1) ) ) 

   end do
   end do

   ! upper boundary - 1st order derivative

   k = kb + 1

   do j = ja , jb
   do i = ia , ib

      xz( le_idx(i,j,k,1) ) = dum1 * ( x( le_idx(i,j,k,1) ) - x( le_idx(i,j,k-1,1) ) )
      yz( le_idx(i,j,k,1) ) = dum1 * ( y( le_idx(i,j,k,1) ) - y( le_idx(i,j,k-1,1) ) )
      zz( le_idx(i,j,k,1) ) = dum1 * ( z( le_idx(i,j,k,1) ) - z( le_idx(i,j,k-1,1) ) ) 

   end do
   end do

     
   ! calculate Jacobian and the inverse transormation
   
   kb = le_kb(1)
   ka = le_ka(1)
   jb = le_jb(1)
   ja = le_ja(1)
   ib = le_ib(1)
   ia = le_ia(1)

   if ( mydown  == mpi_proc_null ) ka = li_ka(1)
   if ( myup    == mpi_proc_null ) kb = li_kb(1)
   if ( myleft  == mpi_proc_null ) ja = li_ja(1)
   if ( myright == mpi_proc_null ) jb = li_jb(1)
   if ( myback  == mpi_proc_null ) ia = li_ia(1)
   if ( myfront == mpi_proc_null ) ib = li_ib(1)

   do k = ka, kb
   do j = ja, jb
   do i = ia, ib
     
      l = le_idx(i,j,k,1)

      gj = xc(l) * ( ye(l) * zz(l) - yz(l) * ze(l) ) +  &
           xe(l) * ( yz(l) * zc(l) - yc(l) * zz(l) ) +  &
           xz(l) * ( yc(l) * ze(l) - ye(l) * zc(l) ) 
 
        if (gj == 0.0) print*, 'zero jacobian at ',n,i,j,k,gj

        aj(l) = one / gj

        csi(1,l) = aj(l) * ( ye(l) * zz(l) - yz(l) * ze(l) )
        csi(2,l) = aj(l) * ( xz(l) * ze(l) - xe(l) * zz(l) )
        csi(3,l) = aj(l) * ( xe(l) * yz(l) - xz(l) * ye(l) )

        eta(1,l) = aj(l) * ( yz(l) * zc(l) - yc(l) * zz(l) )
        eta(2,l) = aj(l) * ( xc(l) * zz(l) - xz(l) * zc(l) )
        eta(3,l) = aj(l) * ( xz(l) * yc(l) - xc(l) * yz(l) )

        zet(1,l) = aj(l) * ( yc(l) * ze(l) - ye(l) * zc(l) )
        zet(2,l) = aj(l) * ( xe(l) * zc(l) - xc(l) * ze(l) )
        zet(3,l) = aj(l) * ( xc(l) * ye(l) - xe(l) * yc(l) )
  
   end do
   end do
   end do
  

   n = 1

   ka = li_ka(ns)
   kb = li_kb(ns)
   ja = li_ja(ns)
   jb = li_jb(ns)
   ia = li_ia(ns)
   ib = li_ib(ns)

   ! calculate computational time step for des equaions
   if ( turbulence ) then

      do k = ka, kb
      do j = ja, jb
      do i = ia, ib
         l = le_idx(i,j,k,ns)

         Sc = sqrt(xc(l)*xc(l) + yc(l)*yc(l) + zc(l)*zc(l))
         Se = sqrt(xe(l)*xe(l) + ye(l)*ye(l) + ze(l)*ze(l))
         Sz = sqrt(xz(l)*xz(l) + yz(l)*yz(l) + zz(l)*zz(l))

         SL = min(Sc, Se, Sz)

         dtev(l) = SL * min(cfl2(myzone), vnn2(myzone)*ren*SL)

         ! Arnone et al. (1995)
         if ( unsteady ) dtev(l) = min(dtev(l),delti*two/three)

      end do
      end do
      end do

      ! blanking area
      if (nblk /= 0) then
         do nb = 1, nblk
            do k = li_blk_ka(ns,nb), li_blk_kb(ns,nb)
            do j = li_blk_ja(ns,nb), li_blk_jb(ns,nb)
            do i = li_blk_ia(ns,nb), li_blk_ib(ns,nb)
               l = le_idx(i,j,k,ns)
   
               dtev(l) = zero
               wd(l) = zero
            end do
            end do
            end do
         end do
      end if

   end if


  !==================================================
  !
  ! distance read in grid file not calculated because
  ! of limitations with mpi
  !
  !==================================================
  !
  ! Calculate variables determined only by the
  ! metrics of the grid
  !
  !if (daf) then

     ! calculate computational time step for daf only
     !
     !call met_daf_dtau

     ! calculate model matrices of Jacobian
     !
     !call met_model

  !end if

!contains

  !include 'met_model.F90'

  !include 'met_daf_dtau.F90'
  !include 'met_distance.F90'
  !include 'met_model_matrices.F90'

end subroutine mg_metrics

