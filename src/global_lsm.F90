module global_lsm
!variables globales para level set method

use precision
implicit none

real(kind = rdf), dimension(:), allocatable :: phi_zero, phi_static
!real(kind = rdf), dimension(:), allocatable :: sgndf,sgndf_n
real(kind = rdf), dimension(:), allocatable :: phi,phi_n
real(kind = rdf), dimension(:), allocatable :: h, hn
real(kind = rdf), dimension(:,:), allocatable :: phi_gradient
real(kind = rdf), dimension(:,:), allocatable :: h_gradient
! flag variable for identifying the nodes to be extrapolated using
! the zero-shear condition
real(kind = rdf), dimension(:), allocatable :: rsign

real(kind = rdf) :: epslsm !parametro malla epsilon
real(kind = rdf) :: epsReinitialisation !reini convergence threshold
real(kind = rdf) :: thetainc
real(kind = rdf) :: drho, dmu,diss_w,diss_a
real(kind = rdf) :: FrLSM, WeLSM
real(kind = rdf) :: IGitermax, IGtimestep
integer :: numiter_lsm
integer :: IG,IGstep2
integer :: RNitermax, RNfreq
integer :: phi_outputiter
integer :: k_surface
real(kind = rdf) :: RNtimestep
integer, dimension(:,:), allocatable :: btype_lsm
real (kind = rdf) :: delti_lsm
integer :: isnarrow, orderAD, nobstacles, bc_extrapolation_order
integer :: sussman_correction_method
real (kind = rdf) :: narrowcoef

! Order of the derivatives at the boundaries

integer :: OrderLSAdvectionBoundaries
integer :: OrderReinitialisationBoundaries

! switch to activate the call instruction to the levelsetmethod
! subroutine (it is read at init_LSM from ind3dmg_LSM.dat)

logical :: call_levelsetmethod 

! switch to call reinitilisation method (either Sussman or Hybrid)
logical :: call_reinitialisation 

! switch to activate hybrid reinitilisation
logical :: hybrid_reinitialisation

! to use eno2 approximation of the derivative at the boundaries 
! in the reinitialisation step.
logical :: ENOBCReinitialisation

! Logical that controls if I apply the Ausas Global Mass Correction
! Algorithm, or I just use the geometric redistancing in the 
! Geometric reinitialisation step
logical :: GlobalMassCorrection

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! global variables for Least-Square Method (lsqm)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

integer :: sweep_lsqm ! # of neighbour nodes to be considered for 
                      ! lsqm matrix sytems (usually three)

real(kind = rdf) :: radius_lsqm 

! size of the coefficient-matrix system

integer, parameter:: mms = 75 !111 ! 111 (75+36) is the experiment with the NDBC 
                               ! 75 ! 72 without coc,  75 for coc ! changed to 57 as Hokkaido code
integer, parameter:: nms = 36 ! changed to 27 to replicate Hokkaido code


! ConvergenceToleranceGeomReini Is the convergence tolerance of the neighbour 
! correction step (step 4). Is the global volume difference between the initial 
! volume enclosed by the neighbours and the volume once phi is corrected using 
! the algorithm. It's suggested to be third or fourth order accurate (dx^3)

real(kind = rdf) :: ConvergenceToleranceGeomReini

! Flag variable to activate/deactivate the total volume computation of the 
! water phase. 

logical :: TotalVolumeComputation

! limit ghost velocities to the maximum in the fluid domain
logical :: zero_pressure_fs

! limit ghost velocities to the maximum in the fluid domain
logical :: limit_ghost_velocities

! Big phi value to cutoff some calculations (it depends on the problem geometry)
real(kind = rdf) :: BigPhi

! Values to save specific iterations in the code
integer :: n_save_LSdebug_tsteps ! total of time steps to be saved
integer :: current_LSdebug_tstep ! it's been rewritten every time I save a debug solution
integer :: current_LSdebug_counter ! it's been rewritten every time I save a debug solution
integer, dimension(20) :: save_LSdebug_tsteps = -1 ! list with tsteps to be saved (max 20)

contains

real (kind = rdf) function heaviside(valor,epsi)

implicit none
real (kind = rdf),intent(in) :: valor
real (kind = rdf), intent(in) ::epsi

        if( valor < -epsi) then
               heaviside = zero
        else if(valor >epsi ) then
               heaviside = one
        else     
               heaviside = one/two * (one + valor/epsi + one/pi*sin(pi*valor/epsi))
        end if
end function heaviside

real (kind = rdf) function deltaf(valor,epsi)

implicit none
real (kind = rdf), intent(in) :: valor
real (kind = rdf), intent(in) :: epsi

        if (abs(valor) > epsi ) then
                deltaf = zero
        else
                deltaf = one/(two*epsi) * (one + cos(pi*valor/epsi))
        end if 
        

end function deltaf

!¿me conviene poner el ratio como parametro? Por ahora no lo hare y asumire que la variable es global (ahora es global)
!---------------------------------------------------------
real (kind = rdf) function rhoLSM(valor)

implicit none 
real (kind = rdf) , intent(in) :: valor

rhoLSM = drho + (one-drho)*heaviside(valor,epslsm)

end function rhoLSM
!----------------------------------------------------------
real (kind = rdf) function muLSM(valor)

implicit none 
real (kind = rdf) , intent(in) :: valor

muLSM = dmu + (one-dmu)*heaviside(valor,epslsm)

end function muLSM
!-----------------------------------------------------------
real (kind = rdf) function dissLSM(valor)

implicit none 
real (kind = rdf) , intent(in) :: valor

dissLSM = diss_a + (diss_w-diss_a)*heaviside(valor,epslsm)

end function dissLSM
!---------------------------------------------------
subroutine integral3d(gfunction,i_coord,j_coord,k_coord,deltax,deltay,deltaz,valint)

        integer,intent (in) :: i_coord, j_coord, k_coord
        real (kind = rdf),intent(in),dimension(i_coord-1:i_coord+1,&
        j_coord-1:j_coord+1,k_coord-1:k_coord+1) :: gfunction
        real (kind = rdf) , intent(in) :: deltax,deltay,deltaz
      
        real (kind = rdf) , intent(out) :: valint 

        !local variables
        real (kind = rdf) :: cof1, cof2,cof3,cof4
        
        cof1= 1.0_rdf/1728.0_rdf ; cof2 = 10.0_rdf/1728.0_rdf
        cof3 = 100.0_rdf/1728.0_rdf ; cof4 = 1000.0_rdf/1728.0_rdf

        valint = cof1 * gfunction(i_coord - 1, j_coord - 1, k_coord - 1) + cof2 * gfunction (i_coord - 1, j_coord, k_coord - 1) + &
                 cof1 * gfunction(i_coord - 1, j_coord + 1, k_coord - 1) + cof2 * gfunction(i_coord, j_coord - 1, k_coord - 1) + &
                 cof3 * gfunction(i_coord, j_coord, k_coord - 1)  + cof2 * gfunction(i_coord, j_coord + 1, k_coord - 1) + & 
                 cof1 * gfunction(i_coord + 1, j_coord - 1, k_coord - 1) + cof2 * gfunction(i_coord + 1, j_coord, k_coord - 1) + &
                 cof1 * gfunction(i_coord + 1, j_coord + 1, k_coord - 1) + cof2 * gfunction(i_coord - 1, j_coord - 1, k_coord) + &
                 cof3 * gfunction(i_coord - 1, j_coord, k_coord) + cof2 * gfunction(i_coord - 1, j_coord + 1, k_coord) + &
                 cof3 * gfunction(i_coord, j_coord - 1, k_coord) + cof4 * gfunction(i_coord, j_coord, k_coord) + &
                 cof3 * gfunction(i_coord, j_coord + 1, k_coord) + cof2 * gfunction(i_coord + 1, j_coord - 1, k_coord) + &
                 cof3 * gfunction(i_coord +1, j_coord, k_coord) + cof2 * gfunction(i_coord + 1, j_coord + 1, k_coord) + &
                 cof1 * gfunction(i_coord - 1, j_coord - 1, k_coord + 1) + cof2 * gfunction(i_coord - 1, j_coord, k_coord + 1) + &
                 cof1 * gfunction(i_coord - 1, j_coord + 1, k_coord + 1) + cof2 * gfunction(i_coord, j_coord - 1, k_coord + 1) + &
                 cof3 * gfunction(i_coord , j_coord , k_coord + 1) + cof2 * gfunction(i_coord, j_coord + 1, k_coord + 1) + &
                 cof1 * gfunction(i_coord + 1, j_coord - 1, k_coord + 1) + cof2 * gfunction(i_coord + 1, j_coord , k_coord + 1) + &
                 cof1 * gfunction(i_coord + 1, j_coord + 1, k_coord + 1)

        valint = deltax * deltay * deltaz * valint
                      

end subroutine integral3d


function GetLambdaMassCorrection ( PhiLocal         ,  &
                                   SignPhi0Local    ,  &
                                   NormGradPhiLocal      )
   implicit none

   real (kind = rdf) :: PhiLocal, SignPhi0Local, NormGradPhiLocal
   real (kind = rdf) :: GetLambdaMassCorrection

   real (kind = rdf) :: num, denom

   num   =   - deltaf(PhiLocal, epslsm) * SignPhi0Local * ( one - NormGradPhiLocal )
   denom = ( ( deltaf(PhiLocal, epslsm) ) ** 2 ) * NormGradPhiLocal 

   if ( abs(denom) < eps_sims ) then
      GetLambdaMassCorrection = zero
   else
      GetLambdaMassCorrection = num / denom
   end if

end function GetLambdaMassCorrection

end module global_lsm
