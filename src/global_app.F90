module global_app

  use precision

  implicit none

  ! jet geometry
  !
  integer :: off_jet_num_pts
  integer, dimension(:,:), allocatable :: off_jet_index

  ! boundary type
  integer, dimension(:,:), allocatable :: btype
  integer, dimension(:,:), allocatable :: bdir

  ! boundary interpolation coefficients
  real (kind = rdf), dimension(:,:,:), allocatable :: rat

  ! diagonal coefficient matrices
  !
  real, dimension(:,:), allocatable :: sa, sb

  ! blanking in overset grids
  !
  integer, dimension(:,:,:), allocatable :: blanking
  integer, dimension(:,:,:), allocatable :: blktype
  integer, dimension(:),     allocatable :: nzblanking

  ! blanking in mpi+overset grids
  integer :: nb, nblk
  ! number of blanking regions within a proc, considering
  ! ghost nodes (e = exterior)
  integer :: nblke

  integer, dimension(:,:), allocatable :: li_blk_ia , li_blk_ib , &
                                          li_blk_ja , li_blk_jb , &
                                          li_blk_ka , li_blk_kb 

  ! The algorithms to search around the node need to be informed if
  ! the ghost nodes are part of a blank area or not
  integer, dimension(:,:), allocatable :: le_blk_ia , le_blk_ib , &
                                          le_blk_ja , le_blk_jb , &
                                          le_blk_ka , le_blk_kb 

  ! fixed pressure location
  ! 
  integer :: ifix
  integer :: jfix
  integer :: kfix
  integer :: nzfix

  integer :: local_kfix
  integer :: local_jfix
  integer :: local_ifix

  integer :: pfix_proc

end module global_app




