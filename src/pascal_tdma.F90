module pascal_tdma

  ! ========================================================================================
  !
	!  DESCRIPTION:
	!  -----------
  !  
  !  This module contains the PasCal Tridiagonal Matrix Algorithm (TDMA)
  !  solve tridiagonal systems in parallel using MPI. It was retrieved from
  !  the work Kim, K. H., Kang, J. H., Pan, X., & Choi, J. I. (2021). PaScaL_TDMA: A 
  !  library of parallel and scalable solvers for massive tridiagonal systems. Computer 
  !  Physics Communications, 260, 107722. 
  !	
  !  Unlike the origival version by Kim et al, this version is capable to solve 
  !  tridiagonal systems with blanked regions. It was adapted for that purpose, to be
  !  able to deal with nodes near the free surface.
  ! 
	!  AUTHORSHIP:
	!  ----------
	!  
	!  Written by Jorge Sandoval, University of Edinburgh / Universidad Catolica de Chile
	!  jsandoval@ed.ac.uk / jcsandov@uc.cl
	!  
	!  First update: 17/01/2024
	!  Last update : 17/01/2024
	!
  ! ========================================================================================

	use precision
	use global_mpi

	! ========================================================================================
	! 	DERIVED TYPES 
	! ========================================================================================

	type :: ptdma_plan_single

  	! ========================================================================================
	  !	
	  !  DESCRIPTION:
	  !  -----------
    !  Execution plan for a single tridiagonal system of equations.
    !  It uses the MPI_Igather function to build the reduced tridiagonal system of equations
    !  to a specified MPI process.
    !
  	! ========================================================================================

  	integer :: ptdma_world ! Single dimensional subcommunicator to assemble data for the 
  												 ! reduced TDMA
    integer :: n_row_rt    ! Number of rows of a reduced tridiagonal system after MPI_Gather

    integer :: gather_rank !Destination rank of MPI_Igather
    integer :: myrank      !Current rank ID in the communicator of ptdma_world
    integer :: nprocs      !Communicator size of ptdma_world

    ! Coefficient arrays after reduction, a: lower, b: diagonal, c: upper, d: rhs.
    ! The orginal dimension (n) is reduced to (2)
    real (kind = rdf) , allocatable, dimension(:) :: A_rd, B_rd, C_rd, D_rd
    !
    ! Coefficient arrays after transpose of a reduced system, a: lower, b: diagonal, 
    ! c: upper, d: rhs
    !
    ! The reduced dimension (2) changes to (2*np) after transpose.
    ! Coefficient arrays, a: lower, b: diagonal, c: upper, d: rhs
    real (kind = rdf) , allocatable, dimension(:) :: A_rt, B_rt, C_rt, D_rt   
    
  end type ptdma_plan_single

  contains

  subroutine PaScaL_TDMA_plan_single_create(plan, myrank, nprocs, mpi_world, gather_rank)

  	! ========================================================================================
	  !	
	  !  DESCRIPTION:
	  !  -----------
  	!  Create a plan for a single tridiagonal system of equations.
  	!
  	!  VARIABLES:
  	!  ----------  
  	!  plan        Plan for a single tridiagonal system of equations
  	!  myrank      Rank ID in mpi_world
  	!  nprocs      Number of MPI process in mpi_world
  	!  mpi_world   Communicator for MPI_Gather and MPI_Scatter of a reduced system
  	!  gather_rank Target rank where all coefficients are gathered into
  	!
  	! ========================================================================================

    implicit none

    type(ptdma_plan_single), intent(inout)  :: plan
    integer, intent(in)     :: myrank, nprocs, mpi_world, gather_rank

    integer :: nr_rd  ! Number of rows of a reduced tridiagonal system per process, 2
    integer :: nr_rt  ! Number of rows of a reduced tridiagonal system after MPI_Gather

    nr_rd = 2
    nr_rt = nr_rd*nprocs

    plan%myrank = myrank
    plan%nprocs = nprocs
    plan%gather_rank = gather_rank
    plan%ptdma_world = mpi_world
    plan%n_row_rt = nr_rt

    allocate( plan%A_rd(1:nr_rd), &
    					plan%B_rd(1:nr_rd), &
    					plan%C_rd(1:nr_rd), &
    					plan%D_rd(1:nr_rd) )
    
    allocate( plan%A_rt(1:nr_rt), &
    					plan%B_rt(1:nr_rt), &
    					plan%C_rt(1:nr_rt), &
    					plan%D_rt(1:nr_rt) )
        
  end subroutine PaScaL_TDMA_plan_single_create

	subroutine PaScaL_TDMA_plan_single_destroy( plan )

	  ! ========================================================================================
	  !	
	  !  DESCRIPTION:
	  !  -----------
		!  Deallocate the allocated arrays in the defined plan_single .
    !	 
  	!  VARIABLES:
  	!  ----------  
    !  plan        Plan for a single tridiagonal system of equations
    !
	  ! ========================================================================================

  	implicit none

    type(ptdma_plan_single), intent(inout)  :: plan

    deallocate(plan%A_rd, plan%B_rd, plan%C_rd, plan%D_rd)
    deallocate(plan%A_rt, plan%B_rt, plan%C_rt, plan%D_rt)

  end subroutine PaScaL_TDMA_plan_single_destroy

	subroutine PaScaL_TDMA_single_solve(plan, A, B, C, D, n_row)

	  ! ========================================================================================
	  !
	  !  DESCRIPTION:
	  !  -----------
    !	 Solve a single tridiagonal system of equations in paraller using MPI.
    !
  	!  VARIABLES:
  	!  ----------  
    !	 plan  : Plan for a single tridiagonal system of equation
    !	 A     : Coefficients in lower diagonal elements
    !	 B     : Coefficients in diagonal elements
    !	 C     : Coefficients in upper diagonal elements
    !	 D     : Coefficients in right-hand side terms
    !	 n_row : # of rows in each process, size of a tridiagonal matrix N divided by nprocs
    !
    !	 This version is also capable of solving traidiagonal matrices with null rows. I 
    !  implemented this feature by setting a tolerance on the matrices entries when the 
    !  coefficientes of the Thomas algorithm are calculated. The reduced tridiagonal system, 
    !  which results from the assembly of the first and last rows of each local matrix in 
    !  process, is solved using MultiGTSV instead of the original routine provided by Kim et
    !  al.
    !
		!  Adapted by Jorge Sandoval, University of Edinburgh / Universidad Catolica de Chile
		!  jsandoval@ed.ac.uk / jcsandov@uc.cl
		!  
		!  First update: 18/12/2023
		!  Last update : 02/01/2024
	  !
	  ! ========================================================================================
    
    use mpi

    implicit none

    type(ptdma_plan_single), intent(inout)   :: plan
    real(kind = rdf) , intent(inout)   :: A(1:n_row), B(1:n_row), C(1:n_row), D(1:n_row)
    integer, intent(in) :: n_row

    ! Temporary variables for computation and parameters for MPI functions
    real(kind = rdf) :: r
    integer :: i
    integer :: request(4), ierr

    real(kind=rdf), parameter :: tol = 1.0E-12
 		integer :: precision_number

    precision_number = kind(tol)

    if ( plan%nprocs == 1 ) then
      
      ! Solve the tridiagonal systems directly when nprocs = 1. 
			call tdma_single(A, B, C, D, n_row)
			return

		end if

		! Reduction step : elimination of lower diagonal elements
		
		if ( abs(B(1)) > tol ) then
			
			A(1) = A(1)/B(1)
			D(1) = D(1)/B(1)
			C(1) = C(1)/B(1)

		else

			A(1) = zero ! 0.0_rdf
			D(1) = zero ! 0.0_rdf
			C(1) = zero ! 0.0_rdf

		end if	

		if ( abs(B(2)) > tol ) then

			A(2) = A(2)/B(2)
			D(2) = D(2)/B(2)
			C(2) = C(2)/B(2)

		else

			A(2) = zero ! 0.0_rdf
			D(2) = zero ! 0.0_rdf
			C(2) = zero ! 0.0_rdf

		end if	

		do i=3,n_row

			if ( abs( B(i)-A(i)*C(i-1) ) > tol ) then

				r    =  one / (B(i)-A(i)*C(i-1))
				D(i) =  r*(D(i)-A(i)*D(i-1))
				C(i) =  r*C(i)
				A(i) = -r*A(i)*A(i-1)

			else

				D(i) = zero ! 0.0_rdf
				C(i) = zero ! 0.0_rdf
				A(i) = zero ! 0.0_rdf

			end if	

		end do
		
		! Reduction step : elimination of upper diagonal elements
		do i=n_row-2,2,-1

			D(i) = D(i)-C(i)*D(i+1)
			A(i) = A(i)-C(i)*A(i+1)
			C(i) =-C(i)*C(i+1)
		
		end do

		if ( abs( one - A(2)*C(1) ) > tol ) then

			r    = one / ( one - A(2)*C(1) )
			D(1) =  r*(D(1)-C(1)*D(2))
			A(1) =  r*A(1)
			C(1) = -r*C(1)*C(2)

		else

			D(1) = zero ! 0.0_rdf
			A(1) = zero ! 0.0_rdf
			C(1) = zero ! 0.0_rdf

		end if	

		! Construct a reduced tridiagonal system of equations per each rank. Each process has 
		! two reduced rows.
		plan%A_rd(1) = A(1) ; plan%A_rd(2) = A(n_row)
		
		if( abs(B(1)) > tol ) then
			plan%B_rd(1) = one !1.0_rdf 
		else
			plan%B_rd(1) = zero ! 0.0_rdf 
		end if

		if( abs( B(n_row) ) > tol ) then
			plan%B_rd(2) = one  ! 1.0_rdf
		else	
			plan%B_rd(2) = zero ! 0.0_rdf
		end if

		plan%C_rd(1) = C(1) ; plan%C_rd(2) = C(n_row)
		plan%D_rd(1) = D(1) ; plan%D_rd(2) = D(n_row)

		! Gather the coefficients of the reduced tridiagonal system to a defined rank, 
		! plan%gather_rank.


		call MPI_Igather(plan%A_rd, 2, MPI_REAL_TYPE, &
										 plan%A_rt, 2, MPI_REAL_TYPE, &
										 plan%gather_rank, plan%ptdma_world, request(1), ierr)
		
		call MPI_Igather(plan%B_rd, 2, MPI_REAL_TYPE, &
										 plan%B_rt, 2, MPI_REAL_TYPE, &
										 plan%gather_rank, plan%ptdma_world, request(2), ierr)
				
		call MPI_Igather(plan%C_rd, 2, MPI_REAL_TYPE, &
									   plan%C_rt, 2, MPI_REAL_TYPE, &
									   plan%gather_rank, plan%ptdma_world, request(3), ierr)
				
		call MPI_Igather(plan%D_rd, 2, MPI_REAL_TYPE, &
									   plan%D_rt, 2, MPI_REAL_TYPE, &
									   plan%gather_rank, plan%ptdma_world, request(4), ierr)


		call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

		! Solve the reduced tridiagonal system on plan%gather_rank.
		if( plan%myrank == plan%gather_rank ) then

			call tdma_single(plan%A_rt,plan%B_rt,plan%C_rt,plan%D_rt, plan%n_row_rt)
			
			!call MultiGTSV( plan%n_row_rt , & ! Size of the reduced system: 2 * nprocs
			!								1             , & ! # of rhs
			!								plan%A_rt     , & ! lower diagonal
			!								plan%B_rt     , & ! main  diagonal
 			!								plan%C_rt     , & ! upper diagonal
			!								plan%D_rt       & ! rhs
			!							)       

		end if

		! Scatter the solutions to each rank.

		select case (precision_number)

			case(4) ! Single precision case

				call MPI_Iscatter(plan%D_rt, 2, MPI_REAL_TYPE, &
													plan%D_rd, 2, MPI_REAL_TYPE, &
													plan%gather_rank, plan%ptdma_world, request(1), ierr)

			case(8) ! Double precision case

				call MPI_Iscatter(plan%D_rt, 2, MPI_DOUBLE_PRECISION, &
													plan%D_rd, 2, MPI_DOUBLE_PRECISION, &
													plan%gather_rank, plan%ptdma_world, request(1), ierr)

		end select

		call MPI_Waitall(1, request, MPI_STATUSES_IGNORE, ierr)

		! Update solutions of the modified tridiagonal system with the solutions of the reduced 
		! tridiagonal system.
		D(1 )    = plan%D_rd(1)
		D(n_row) = plan%D_rd(2)

		do i=2,n_row-1
			D(i) = D(i)-A(i)*D(1)-C(i)*D(n_row)
		end do

  end subroutine PaScaL_TDMA_single_solve


	subroutine tdma_single(a, b, c, d, n1)
		
	  ! ========================================================================================
	  !
	  !  DESCRIPTION:
	  !  -----------
		!  Solve a single tridiagonal system of equations using the Thomas algorithm.
		! 
  	!  VARIABLES:
  	!  ----------  
		!  a  : Coefficients in lower diagonal elements
		!  b  : Coefficients in diagonal elements
		!  c  : Coefficients in upper diagonal elements
		!  d  : Coefficients in the right-hand side terms
		!  n1 : # of rows in each process, dimension of tridiagonal matrix N divided by nprocs
		!
	  ! ========================================================================================

    implicit none

    integer, intent(in) :: n1
    real(kind = rdf) , intent(inout) :: a(n1), b(n1), c(n1), d(n1)
    
    integer :: i
    real(kind = rdf) :: r
    real(kind = rdf), parameter :: tol = 1.0E-12

    if ( abs ( b(1) ) > tol ) then

	    d(1)=d(1)/b(1)
  	  c(1)=c(1)/b(1)

  	else

  		d(1) = zero ! 0.0_rdf
  		c(1) = zero ! 0.0_rdf

  	end if

    do i=2,n1

    	if ( abs( b(i)-a(i)*c(i-1) ) > tol ) then

        r = one / ( b(i) - a(i)*c(i-1) )

        d(i)=r*(d(i)-a(i)*d(i-1))
        c(i)=r*c(i)
    
      else

      	d(i) = zero ! 0.0_rdf
      	c(i) = zero ! 0.0_rdf

     	end if

    end do

    do i=n1-1,1,-1
        d(i)=d(i)-c(i)*d(i+1)
    enddo

	end subroutine tdma_single


end module pascal_tdma
