subroutine calc_h( il,iu          ,&
                   jl,ju          ,&
                   kl,ku          ,&
                   igp, jgp, kgp  ,&
                   z              ,&
                   phi            ,&
                   h               &
                  )

! Calculo de |grad(h0)| mediante ENO2, el cual es usado posteriormente para la correccion
! de Sussman. Tambien se calcula funcion S(h0).

use global_mpi
use global_app
use global_lsm, only : OrderReinitialisationBoundaries , ENOBCReinitialisation

implicit none


! Input variables

integer, intent(in) :: il,iu,jl,ju,kl,ku ! external nodes
integer, intent(in) :: igp, jgp, kgp ! # of ghostpoint at local processor 
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: z
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(in) :: phi

! Output variable
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) , intent(inout) :: h


!local

logical        :: FreeSurfaceWithinProc
integer        :: IDFreeSurfaceProc
real(kind=rdf) :: z_free_surface
real (kind = rdf) :: inter1, inter2

!index
integer  :: i_mysta,i_myend
integer  :: j_mysta,j_myend
integer  :: k_mysta,k_myend

integer  :: ista , iend
integer  :: jsta , jend
integer  :: ksta , kend

integer :: i , j , k 

! Switches (logicals) para forzar stencils bounded en los bordes del dominio

logical :: BackOutbounded   , LeftOutbounded   , DownOutbounded 
logical :: FrontOutbounded  , RightOutbounded  , UpOutbounded   

integer :: iBiasDirection, jBiasDirection, kBiasDirection

logical :: BlankingFlag

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!Nodos incluyendo el borde
i_mysta = il + igp
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp        
j_myend = ju - jgp
k_myend = ku - kgp

!Nodos sin incluir borde (el borde se actualiza en bcond_lsm.F90)

if (myback == mpi_proc_null)  i_mysta = il + igp + 1
if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

if (myfront == mpi_proc_null) i_myend = iu - igp - 1
if (myright == mpi_proc_null) j_myend = ju - jgp - 1
if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

! Physical boundaries of the processor
   
ista = il ; jsta = jl ; ksta = kl 
iend = iu ; jend = ju ; kend = ku 

if ( myback  == mpi_proc_null )  ista = il + igp 
if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

if ( myfront == mpi_proc_null )  iend = iu - igp
if ( myright == mpi_proc_null )  jend = ju - jgp
if ( myup    == mpi_proc_null )  kend = ku - kgp

inter1 =  four / three
inter2 =  -one / three


h = zero
z_free_surface = -one

do j = jsta , jend
do i = ista , iend

  BlankingFlag = .false.

  ! the nodes inside the blanking region are skipped
  if ( nblke /= 0 ) then

    do nb = 1 , nblke
      
      if ( i > le_blk_ia(1,nb) .and. i < le_blk_ib(1,nb) .and. &  
           j > le_blk_ja(1,nb) .and. j < le_blk_jb(1,nb) .and. &
           k > le_blk_ka(1,nb) .and. k < le_blk_kb(1,nb) ) then

        BlankingFlag = .true.

      end if

    end do
  
  end if

  if ( BlankingFlag ) cycle

  IDFreeSurfaceProc = -1
  
  ! Check if the free surface is in this proc and if it
  ! finds it, FreeSurfaceWithinProc returns as .true.
  call find_z_free_surface ( z  (i,j,ksta:kend)       , &
                             phi(i,j,ksta:kend)       , &
                             FreeSurfaceWithinProc    , &
                             z_free_surface             &
                            )
  
  ! If the free surface is in this proc, then I assign the ID
  if ( FreeSurfaceWithinProc ) then

    IDFreeSurfaceProc = comm1d_csi(1)%myrank

  end if

  ! I communicate the values of FreeSurfaceWithinProc among the
  ! vertical processors contained in the comm1d_csi(1)%mpi_comm
  ! communicator. If there is one that's true, then all the procs
  ! set FreeSurfaceWithinProc = .true. 
  
  !call MPI_Allreduce( mpi_in_place            , & 
  !                    FreeSurfaceWithinProc   , & 
  !                    1                       , & 
  !                    mpi_logical             , & 
  !                    mpi_lor                 , &
  !                    comm1d_csi(1)%mpi_comm  , & 
  !                    ierr                      &
  !                  )


  if ( FreeSurfaceWithinProc ) then
 
    ! The processor that found the value broadcasts it
    !call MPI_Bcast( z_free_surface          , &
    !                1                       , &
    !                MPI_REAL_TYPE           , &
    !                IDFreeSurfaceProc       , &
    !                comm1d_csi(1)%mpi_comm  , &
    !                ierr                      &
    !              )
 
  end if
  
  if ( IDFreeSurfaceProc /= -1 ) then
  
    !if (z_free_surface < 0.044 .or. z_free_surface > 0.056 ) then
    !  print *, 'myid = ', myid , 'i,j,k = ',i,j,k, ' , z_free_surface = ', z_free_surface
    !end if  

    ! Local Water Depth (negative within the air phase)
    
    do k = ksta , kend
      
      h(i,j,k) = z_free_surface - z(i,j,k)

    end do
  
  else
    
    print *, 'free surface not found' 
    
    ! h = 0
    h(i,j,:) = zero
    
  end if
  
end do
end do

  ! Blanking region h boundary condition
  if (nblk /= 0) then
    
    do nb = 1, nblk
  
      !=============================================================================== 
      ! ξ - direction
      !=============================================================================== 
      !
      ! Front side of an obstacle
      !
      !             i = li_blk_ib(1,nb)
      !             |##############
      !             |##          ##
      !   o----o----o## OBSTACLE ##
      !  i-2  i-1   |##          ##     
      !             |##############
      !   

      i = li_blk_ia(1,nb)
        
      if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 ) then
                
        do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) ) 

          h(i,j,:) = inter1 * h(i-1,j,:) + &
                     inter2 * h(i-2,j,:)

          h(i+1,j,:) = h(i,j,:) 
        
        end do

      end if    
      
      ! 
      ! Rear side of an obstacle
      !  
      !                i = li_blk_ib(1,nb)
      !  ##############|        
      !  ##          ##|        
      !  ## OBSTACLE ##o----o----o
      !  ##          ##|   i+1  i+2       
      !  ##############|
      !   
      
      i = li_blk_ib(1,nb)

      if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 ) then
                
        do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) )
        
          h(i,j,:) = inter1 * h(i+1,j,:) + &
                     inter2 * h(i+2,j,:)

          h(i-1,j,:) = h(i,j,:) 
        
        end do
       
      end if
       
      !=============================================================================== 
      ! η - direction
      !=============================================================================== 
              
      ! Front side of an obstacle
      j = li_blk_ja(1,nb)
       
      if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then
                    
        do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) )

          h(i,j,:) = inter1 * h(i,j-1,:) + &
                     inter2 * h(i,j-2,:) 

          h(i,j+1,:) = h(i,j,:) 

        end do
        
      end if    
       
      j = li_blk_jb(1,nb)
        
      ! Rear side of an obstacle
      if ( blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then
                
        do i =max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) )
            
          h(i,j,:) = inter1 * h(i,j+1,:) + &
                     inter2 * h(i,j+2,:) 

          h(i,j-1,:) = h(i,j,:) 

        end do     
        
      end if

      !=============================================================================== 
      ! Corners
      !=============================================================================== 

      ! Corner 1
      i = li_blk_ia(1,nb)
      j = li_blk_ja(1,nb)

      if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
           blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

        ! Zero phi gradient
        h(i,j,:)   = inter1 * h(i-1,j-1,:) + & 
                     inter2 * h(i-2,j-2,:)    

        h(i+1,j+1,:) = h(i,j,:) 

      end if

      ! Corner 2
      i = li_blk_ib(1,nb)
      j = li_blk_ja(1,nb)

      if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
           blktype(3,nb,myzone) == 0 .and. j > j_mysta+1 ) then

        h(i,j,:)   = inter1 * h(i+1,j-1,:) + & 
                     inter2 * h(i+2,j-2,:)    

        h(i-1,j+1,:) = h(i,j,:) 

      end if

      ! Corner 3
      i = li_blk_ib(1,nb)
      j = li_blk_jb(1,nb)

      if ( blktype(2,nb,myzone) == 0 .and. i < i_myend-1 .and. &
           blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

        h(i,j,:)   = inter1 * h(i+1,j+1,:) + & 
                     inter2 * h(i+2,j+2,:)    

        h(i-1,j-1,:) = h(i,j,:) 

      end if

      ! Corner 4
      i = li_blk_ia(1,nb)
      j = li_blk_jb(1,nb)

      if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta+1 .and. &
           blktype(4,nb,myzone) == 0 .and. j < j_myend-1 ) then

        h(i,j,:)   = inter1 * h(i-1,j+1,:) + & 
                     inter2 * h(i-2,j+2,:)    

        h(i+1,j-1,:) = h(i,j,:) 


      end if


        !=============================================================================== 
        ! ζ - direction
        !=============================================================================== 

        !k = li_blk_kb(1,nb)
        
        !if (k <= k_myend) then
                
          !do j = li_blk_ja(1,nb) , li_blk_jb(1,nb)
          !do i = li_blk_ia(1,nb) , li_blk_ib(1,nb)

          !  phi(i,j,k) = inter1 * phi(i,j,min(k+1,k_myend)) + &
          !               inter2 * phi(i,j,min(k+2,k_myend))
          !end do
          !end do     
        
        !end if  
       
        !k = li_blk_ka(1,nb)
        
        !if (k >= k_mysta) then
                
          !do j = li_blk_ja(1,nb) , li_blk_jb(1,nb)
          !do i = li_blk_ia(1,nb) , li_blk_ib(1,nb)

            !phi(i,j,k) = inter1 * phi(i,j,max(k-1,k_mysta)) + &
            !             inter2 * phi(i,j,max(k-2,k_mysta))
          !end do
          !end do

        !end if       
                   
    end do
  end if


contains

include 'rhs_exchng3_4d.F90'
include 'find_z_free_surface.F90'

end subroutine calc_h
