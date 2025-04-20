
subroutine rhs_flux_sans_convec
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Orthogonal, curvlinear, cartesian coordinates

  ! Calculate pressure gradient and convective terms
  ! using second-order central differences.
  !
  ! CHANGE: For the level set version of the code, the pressure terms
  ! are calculated in the ghost fluid extrapolation routine (Jorge Sandoval,
  ! June 2023)

  ! Calculate divergence of the velocity for the continuity
  ! equation also using second-order central differences.

  ! Unlike previous codes we will add viscous, dissipation, and
  ! unsteady terms in the calling routine not in this one.

  ! input
  !     ucn_j(3,ijk)
  !     csi(3,ijk)
  !     aj(ijk)
  !     q(4,ijk)
  !     dc, de, dz

  ! output
  !     rh(4,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local arrays
  !
    real (kind = rdf), dimension(:,:,:), allocatable :: pJ !, pJh

  real (kind = rdf) :: dc2, de2, dz2
  
  ! local dummy reals
  !
  real (kind = rdf) :: dudx, dvdy, dwdz

  ! local dummy integers
  !
  integer :: i, j, k

  ! initialize
  ! 
  allocate (pJ(il:iu,jl:ju,kl:ku))
  !allocate (pJh(il:iu,jl:ju,kl:ku))

  pJ  = zero
  !pJh = zero

  ! dum --> dc2 or de2 or dz2
  dc2 = one_half * dc
  de2 = one_half * de
  dz2 = one_half * dz
  
  do k = k_mysta-1, k_myend+1
  do j = j_mysta-1, j_myend+1
  do i = i_mysta-1, i_myend+1
    pJ(i,j,k)  =   q(1,i,j,k) / aj(i,j,k) ! p/J
    !pJh(i,j,k) =  ( h(i,j,k)/(FrLSM**two) ) / aj(i,j,k) ! h / (J * Fr^2)
  end do
  end do
  end do

  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
  ! I'll deactivate this term for the moment:
  !pJh = zero
  !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

  ! TO DO: think about whole array operations
  ! csi direction
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
   
        rh(2,i,j,k) = dc2 * ( ( pJ(i+1,j,k) ) * csi(1,i+1,j,k)  - & 
                              ( pJ(i-1,j,k) ) * csi(1,i-1,j,k) )
        
        rh(3,i,j,k) = dc2 * ( ( pJ(i+1,j,k)  ) * csi(2,i+1,j,k)  - &
                              ( pJ(i-1,j,k)  ) * csi(2,i-1,j,k) )
        
        rh(4,i,j,k) = dc2 * (   pJ(i+1,j,k) * csi(3,i+1,j,k) - &
                                pJ(i-1,j,k) * csi(3,i-1,j,k) )
  
  end do
  end do
  end do
   
  ! eta direction
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
  
    rh(2,i,j,k) = rh(2,i,j,k) + de2 * ( ( pJ(i,j+1,k) ) * eta(1,i,j+1,k) - &
                                        ( pJ(i,j-1,k) ) * eta(1,i,j-1,k) )

    rh(3,i,j,k) = rh(3,i,j,k) + de2 * ( ( pJ(i,j+1,k) ) * eta(2,i,j+1,k) - &
                                        ( pJ(i,j-1,k) ) * eta(2,i,j-1,k) )

    rh(4,i,j,k) = rh(4,i,j,k) + de2 * (   pJ(i,j+1,k) * eta(3,i,j+1,k) - & 
                                          pJ(i,j-1,k) * eta(3,i,j-1,k) )
  
  end do
  end do
  end do
   
  ! zet direction
   
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

    rh(2,i,j,k) = rh(2,i,j,k) + dz2 * ( ( pJ(i,j,k+1) ) * zet(1,i,j,k+1)  - &
                                        ( pJ(i,j,k-1) ) * zet(1,i,j,k-1) )

    rh(3,i,j,k) = rh(3,i,j,k) + dz2 * ( ( pJ(i,j,k+1) ) * zet(2,i,j,k+1)  - &
                                        ( pJ(i,j,k-1) ) * zet(2,i,j,k-1) )
    
    rh(4,i,j,k) = rh(4,i,j,k) + dz2 * (   pJ(i,j,k+1) * zet(3,i,j,k+1) - &
                                          pJ(i,j,k-1) * zet(3,i,j,k-1) )
  
  end do
  end do
  end do

  ! local divergence

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     dudx = dc2 * ( ucn_j(1 , i+1 , j   , k  ) - ucn_j(1 , i-1 , j   , k  ) )
     dvdy = de2 * ( ucn_j(2 , i   , j+1 , k  ) - ucn_j(2 , i   , j-1 , k  ) )
     dwdz = dz2 * ( ucn_j(3 , i   , j   , k+1) - ucn_j(3 , i   , j   , k-1) )

     rh(1,i,j,k) = ( dudx + dvdy + dwdz )
     
  end do
  end do
  end do

  deallocate (pJ)
  !deallocate (pJh)

end subroutine rhs_flux_sans_convec

