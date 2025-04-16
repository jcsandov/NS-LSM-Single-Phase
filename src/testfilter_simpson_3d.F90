function testfilter_simpson_3d ( var , aj ) result( var_f )

! ------------------------------------------------------------------------------
! DESCRIPTION:
!
! This routine calculates a filtered variable var over a ±1 domain around
! (i,j,k) in the three computational directions. 
!
! To calculate the filtered variable, this is integrated over the 
! neighbourhood domain. In curvilinear coordinates, the change of 
! variables theorem applies: 
! https://math.mit.edu/~larsh/teaching/F2007/handouts/changeofvariables.pdf
!
!                     
!  ∫∫∫ f(x,y,z) dxdydz = ∫∫∫ f( x(ξ,η,ζ) , y(ξ,η,ζ) , z(ξ,η,ζ)) * (1/J) dξdηdζ 
!   Ωx                    Ωξ                                                   
!
! Where J is the Jacobian of the transformation between from 
! (ξ,η,ζ) --> (x,y,z)
!                                                                _
! This transformation allows us to calculate a filtered variable f as
!
!
!              ∫∫∫ f( x(ξ,η,ζ) , y(ξ,η,ζ) , z(ξ,η,ζ)) * (1/J) dξdηdζ
!  _            Ωξ                                                  
!  f(x,y,z) = -------------------------------------------------------
!                                ∫∫∫  (1/J) dξdηdζ  
!                                 Ωξ                                                    
!
!
! The integrals of the numerator and denominator can be calculated using 
! the Simpson's formula for uniform grids. The coefficients for weighting
! nodal values are described in "A Multidimensional Analogue of the 
! Simpson's Formula of Integral" (Fujii, 2010), available here:
!
! https://arxiv.org/abs/0908.3533  
!
! ------------------------------------------------------------------------------

!input: var
!output: var_f

implicit none
 
! allocatable variablesI
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent (in)  :: var, aj
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: var_f

!local
real (kind = rdf) :: r_simp, weight

! Simpson's rule coefficients with the solver precision
real (kind = rdf), parameter :: c64 = 64.0_rdf
real (kind = rdf), parameter :: c16 = 16.0_rdf
real (kind = rdf), parameter :: c4  =  4.0_rdf

! Local weights p and m correspondigly mean +1 and -1 in that index
! so, for example, w_p0m means w(i+1,j,k-1) 
real (kind = rdf) :: w_000

real (kind = rdf) :: w_p00
real (kind = rdf) :: w_m00
real (kind = rdf) :: w_0p0
real (kind = rdf) :: w_0m0
real (kind = rdf) :: w_00p
real (kind = rdf) :: w_00m

real (kind = rdf) :: w_pp0
real (kind = rdf) :: w_p0p
real (kind = rdf) :: w_0pp
real (kind = rdf) :: w_mm0
real (kind = rdf) :: w_m0m
real (kind = rdf) :: w_0mm
real (kind = rdf) :: w_pm0
real (kind = rdf) :: w_mp0
real (kind = rdf) :: w_p0m
real (kind = rdf) :: w_m0p
real (kind = rdf) :: w_0mp
real (kind = rdf) :: w_0pm

real (kind = rdf) :: w_ppp
real (kind = rdf) :: w_mmm
real (kind = rdf) :: w_ppm
real (kind = rdf) :: w_pmm
real (kind = rdf) :: w_mpp
real (kind = rdf) :: w_mmp
real (kind = rdf) :: w_mpm
real (kind = rdf) :: w_pmp

!indexes
integer :: i_mysta , i_myend
integer :: j_mysta , j_myend
integer :: k_mysta , k_myend

integer :: i, j, k

integer :: ista , iend , jsta , jend , ksta , kend
! iplus, iminus
integer :: im   , ip   , jm   , jp   , km   , kp

! Physical boundaries within the processor
ista = il ; jsta = jl ; ksta = kl 
iend = iu ; jend = ju ; kend = ku 

if ( myback  == mpi_proc_null )  ista = il + igp  
if ( myleft  == mpi_proc_null )  jsta = jl + jgp  
if ( mydown  == mpi_proc_null )  ksta = kl + kgp  

if ( myfront == mpi_proc_null )  iend = iu - igp 
if ( myright == mpi_proc_null )  jend = ju - jgp 
if ( myup    == mpi_proc_null )  kend = ku - kgp 

!-------------------------------------------------------------------------------------------

i_mysta = il + igp
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp
j_myend = ju - jgp
k_myend = ku - kgp

! processes on the domain boundaries
! 

im = 1 ; jm = 1 ; km = 1
ip = 1 ; jp = 1 ; kp = 1

if ( myback  == mpi_proc_null )  i_mysta = il + igp + 1 ; im = 0
if ( myleft  == mpi_proc_null )  j_mysta = jl + jgp + 1 ; jm = 0
if ( mydown  == mpi_proc_null )  k_mysta = kl + kgp + 1 ; km = 0

if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1 ; ip = 0
if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1 ; jp = 0
if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1 ; kp = 0

! var_f initialisation
var_f = var

! Interior nodes ± 1 to have a node for the next filtering
do k = k_mysta - km , k_myend + kp
do j = j_mysta - jm , j_myend + jp
do i = i_mysta - im , i_myend + ip

    ! Local weights
    w_000  =  one  /  aj ( i   , j   , k   )
                                                                      
    w_p00  =  one  /  aj ( i+1 , j   , k   )
    w_m00  =  one  /  aj ( i-1 , j   , k   )
    w_0p0  =  one  /  aj ( i   , j+1 , k   )
    w_0m0  =  one  /  aj ( i   , j-1 , k   )
    w_00p  =  one  /  aj ( i   , j   , k+1 )
    w_00m  =  one  /  aj ( i   , j   , k-1 )
                     
    w_pp0  =  one  /  aj ( i+1 , j+1 , k   )
    w_p0p  =  one  /  aj ( i+1 , j   , k+1 )
    w_0pp  =  one  /  aj ( i   , j+1 , k+1 )
    w_mm0  =  one  /  aj ( i-1 , j-1 , k   )
    w_m0m  =  one  /  aj ( i-1 , j   , k-1 )
    w_0mm  =  one  /  aj ( i   , j-1 , k-1 )
    w_pm0  =  one  /  aj ( i+1 , j-1 , k   )
    w_mp0  =  one  /  aj ( i-1 , j+1 , k   )
    w_p0m  =  one  /  aj ( i+1 , j   , k-1 )
    w_m0p  =  one  /  aj ( i-1 , j   , k+1 )
    w_0mp  =  one  /  aj ( i   , j-1 , k+1 )
    w_0pm  =  one  /  aj ( i   , j+1 , k-1 )
                                            
    w_ppp  =  one  /  aj ( i+1 , j+1 , k+1 )
    w_mmm  =  one  /  aj ( i-1 , j-1 , k-1 )
    w_ppm  =  one  /  aj ( i+1 , j+1 , k-1 )
    w_pmm  =  one  /  aj ( i+1 , j-1 , k-1 )
    w_mpp  =  one  /  aj ( i-1 , j+1 , k+1 )
    w_mmp  =  one  /  aj ( i-1 , j-1 , k+1 )
    w_mpm  =  one  /  aj ( i-1 , j+1 , k-1 )
    w_pmp  =  one  /  aj ( i+1 , j-1 , k+1 )

    ! Numerator
    r_simp = c64 *      w_000 * var ( i   , j   , k   )           &
             +                                                    &  
             c16 *  (   w_p00 * var ( i+1 , j   , k   )      +    & 
                        w_m00 * var ( i-1 , j   , k   )      +    & 
                        w_0p0 * var ( i   , j+1 , k   )      +    & 
                        w_0m0 * var ( i   , j-1 , k   )      +    & 
                        w_00p * var ( i   , j   , k+1 )      +    & 
                        w_00m * var ( i   , j   , k-1 )   )       &
             +                                                    &     
             c4  *  (   w_pp0 * var ( i+1 , j+1 , k   )      +    &
                        w_p0p * var ( i+1 , j   , k+1 )      +    & 
                        w_0pp * var ( i   , j+1 , k+1 )      +    &
                        w_mm0 * var ( i-1 , j-1 , k   )      +    & 
                        w_m0m * var ( i-1 , j   , k-1 )      +    & 
                        w_0mm * var ( i   , j-1 , k-1 )      +    & 
                        w_pm0 * var ( i+1 , j-1 , k   )      +    & 
                        w_mp0 * var ( i-1 , j+1 , k   )      +    &
                        w_p0m * var ( i+1 , j   , k-1 )      +    & 
                        w_m0p * var ( i-1 , j   , k+1 )      +    & 
                        w_0mp * var ( i   , j-1 , k+1 )      +    & 
                        w_0pm * var ( i   , j+1 , k-1 )   )       &
             +                                                    &
                        w_ppp * var ( i+1 , j+1 , k+1 )      +    &
                        w_mmm * var ( i-1 , j-1 , k-1 )      +    & 
                        w_ppm * var ( i+1 , j+1 , k-1 )      +    & 
                        w_pmm * var ( i+1 , j-1 , k-1 )      +    &
                        w_mpp * var ( i-1 , j+1 , k+1 )      +    &
                        w_mmp * var ( i-1 , j-1 , k+1 )      +    & 
                        w_mpm * var ( i-1 , j+1 , k-1 )      +    & 
                        w_pmp * var ( i+1 , j-1 , k+1 )   

    ! Denominator
    weight = c64 *      w_000           &
             +                          &  
             c16 *  (   w_p00      +    & 
                        w_m00      +    & 
                        w_0p0      +    & 
                        w_0m0      +    & 
                        w_00p      +    & 
                        w_00m   )       &
             +                          &     
             c4  *  (   w_pp0      +    &
                        w_p0p      +    & 
                        w_0pp      +    &
                        w_mm0      +    & 
                        w_m0m      +    & 
                        w_0mm      +    & 
                        w_pm0      +    & 
                        w_mp0      +    &
                        w_p0m      +    & 
                        w_m0p      +    & 
                        w_0mp      +    & 
                        w_0pm   )       &
             +                          &
                        w_ppp      +    &
                        w_mmm      +    & 
                        w_ppm      +    & 
                        w_pmm      +    &
                        w_mpp      +    &
                        w_mmp      +    & 
                        w_mpm      +    & 
                        w_pmp   

    var_f(i,j,k) = r_simp / weight            
                          
end do
end do
end do


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!    
! ####### ######   #####  #######  #####  
! #       #     # #       #       #      
! #####   #     # #  #### #####    #####  
! #       #     # #     # #             # 
! ####### ######   #####  #######  #####  
!                                         
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!            i1,jm,km ---E11--- im,jm,km    
!                /|(8)           /|(7)                            
!        E12----/ |             /-|------E10                                   
!            i1,j1,km---E9---im,j1,km                   
!              |(5)           |(6)|                                   
!              | E8       E3  |   E7                                
!             E5  |       |   E6  |                                
!              | i1,jm,k1-----|-im,jm,k1           
!          E4--|-/(4)         |  /(3)                 
!              |/             | / -------E2
! ζ,k      i1,j1,k1---E1----im,j1,k1
! ^   η,j     (1)            (2)   
! |  7             
! | /
! |/
! +-----------> ξ,i

! Edge 1  : 1 - 2 --> i free ; j = 1  , k = 1
! Edge 1  : 2 - 3 --> j free ; i = im , k = 1
! Edge 3  : 3 - 4 --> i free ; j = jm , k = 1
! Edge 4  : 4 - 1 --> j free ; i = 1  , k = 1
!
! Edge 5  : 1 - 5 --> k free ; i = 1  , j = 1
! Edge 6  : 2 - 6 --> k free ; i = im , j = 1
! Edge 7  : 3 - 7 --> k free ; i = im , j = jm
! Edge 8  : 4 - 8 --> k free ; i = 1  , j = jm
!
! Edge 9  : 5 - 6 --> i free ; j = 1  , k = km
! Edge 10 : 6 - 7 --> j free ; i = im , k = km
! Edge 11 : 7 - 8 --> i free ; j = jm , k = km
! Edge 12 : 8 - 5 --> j free ; i = 1  , k = km
!
! i free : E1 , E3 , E9  , E11
! j free : E2 , E4 , E10 , E12
! k free : E5 , E6 , E7  , E8
!
! - - - - - - - - - - - - - - - - - - - - - - - - 

!Dirección i

! Edge 1
if (       mydown  == mpi_proc_null &
     .and. myleft  == mpi_proc_null   ) then
    
    k = ksta
    j = jsta

    do i = i_mysta, i_myend
    
        w_m00  = one / aj ( i-1 , j , k )
        w_000  = one / aj ( i   , j , k )
        w_p00  = one / aj ( i+1 , j , k )

        r_simp =       w_m00 * var ( i-1 , j , k ) + &
                  c4 * w_000 * var ( i   , j , k ) + &
                       w_p00 * var ( i+1 , j , k )
        
        weight =       w_m00 + &
                  c4 * w_000 + &
                       w_p00

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 9
if (       myup    == mpi_proc_null &
     .and. myleft  == mpi_proc_null   ) then
    
    k = kend
    j = jsta

    do i = i_mysta, i_myend
    
        w_m00  = one / aj ( i-1 , j , k )
        w_000  = one / aj ( i   , j , k )
        w_p00  = one / aj ( i+1 , j , k )

        r_simp =       w_m00 * var ( i-1 , j , k ) + &
                  c4 * w_000 * var ( i   , j , k ) + &
                       w_p00 * var ( i+1 , j , k )
        
        weight =       w_m00 + &
                  c4 * w_000 + &
                       w_p00

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 3
if (       mydown  == mpi_proc_null &
     .and. myright == mpi_proc_null   ) then
    
    k = ksta
    j = jend

    do i = i_mysta, i_myend
    
        w_m00  = one / aj ( i-1 , j , k )
        w_000  = one / aj ( i   , j , k )
        w_p00  = one / aj ( i+1 , j , k )

        r_simp =       w_m00 * var ( i-1 , j , k ) + &
                  c4 * w_000 * var ( i   , j , k ) + &
                       w_p00 * var ( i+1 , j , k )
        
        weight =       w_m00 + &
                  c4 * w_000 + &
                       w_p00

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 11
if (       myup    == mpi_proc_null &
     .and. myright == mpi_proc_null   ) then
    
    k = kend
    j = jend

    do i = i_mysta, i_myend
    
        w_m00  = one / aj ( i-1 , j , k )
        w_000  = one / aj ( i   , j , k )
        w_p00  = one / aj ( i+1 , j , k )

        r_simp =       w_m00 * var ( i-1 , j , k ) + &
                  c4 * w_000 * var ( i   , j , k ) + &
                       w_p00 * var ( i+1 , j , k )
        
        weight =       w_m00 + &
                  c4 * w_000 + &
                       w_p00

        var_f(i,j,k) = r_simp / weight

    end do

end if        


!Dirección j

! Edge 4
if (       myback  == mpi_proc_null &
     .and. mydown  == mpi_proc_null   ) then
    
    i = ista
    k = ksta

    do j = j_mysta , j_myend

        w_0m0  = one / aj ( i , j-1 , k )
        w_000  = one / aj ( i , j   , k )
        w_0p0  = one / aj ( i , j+1 , k )

        r_simp =       w_0m0 * var ( i , j-1 , k ) + &
                  c4 * w_000 * var ( i , j   , k ) + &
                       w_0p0 * var ( i , j+1 , k )


        weight =       w_0m0 + &
                  c4 * w_000 + &
                       w_0p0

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 2
if (       myfront == mpi_proc_null &
     .and. mydown  == mpi_proc_null   ) then
    
    i = iend
    k = ksta

    do j = j_mysta , j_myend

        w_0m0  = one / aj ( i , j-1 , k )
        w_000  = one / aj ( i , j   , k )
        w_0p0  = one / aj ( i , j+1 , k )

        r_simp =       w_0m0 * var ( i , j-1 , k ) + &
                  c4 * w_000 * var ( i , j   , k ) + &
                       w_0p0 * var ( i , j+1 , k )


        weight =       w_0m0 + &
                  c4 * w_000 + &
                       w_0p0

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 12
if (       myback  == mpi_proc_null &
     .and. myup    == mpi_proc_null   ) then
    
    i = ista
    k = kend

    do j = j_mysta , j_myend

        w_0m0  = one / aj ( i , j-1 , k )
        w_000  = one / aj ( i , j   , k )
        w_0p0  = one / aj ( i , j+1 , k )

        r_simp =       w_0m0 * var ( i , j-1 , k ) + &
                  c4 * w_000 * var ( i , j   , k ) + &
                       w_0p0 * var ( i , j+1 , k )


        weight =       w_0m0 + &
                  c4 * w_000 + &
                       w_0p0

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 10
if (       myfront == mpi_proc_null &
     .and. myup    == mpi_proc_null   ) then
    
    i = iend
    k = kend

    do j = j_mysta , j_myend

        w_0m0  = one / aj ( i , j-1 , k )
        w_000  = one / aj ( i , j   , k )
        w_0p0  = one / aj ( i , j+1 , k )

        r_simp =       w_0m0 * var ( i , j-1 , k ) + &
                  c4 * w_000 * var ( i , j   , k ) + &
                       w_0p0 * var ( i , j+1 , k )


        weight =       w_0m0 + &
                  c4 * w_000 + &
                       w_0p0

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Direccion k
! Edge 5
if (       myback == mpi_proc_null &
     .and. myleft == mpi_proc_null   ) then
    
    i = ista
    j = jsta

    do k = k_mysta , k_myend
    
        w_00m  = one / aj ( i , j , k-1 )
        w_000  = one / aj ( i , j , k   )
        w_00p  = one / aj ( i , j , k+1 )    

        r_simp =       w_00m * var ( i , j , k-1 )  +  &
                  c4 * w_000 * var ( i , j , k   )  +  & 
                       w_00p * var ( i , j , k+1 ) 

        weight =       w_00m +  &
                  c4 * w_000 +  & 
                       w_00p 

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 6
if (       myfront == mpi_proc_null &
     .and. myleft  == mpi_proc_null   ) then
    
    i = iend
    j = jsta

    do k = k_mysta , k_myend
    
        w_00m  = one / aj ( i , j , k-1 )
        w_000  = one / aj ( i , j , k   )
        w_00p  = one / aj ( i , j , k+1 )    

        r_simp =       w_00m * var ( i , j , k-1 )  +  &
                  c4 * w_000 * var ( i , j , k   )  +  & 
                       w_00p * var ( i , j , k+1 ) 

        weight =       w_00m +  &
                  c4 * w_000 +  & 
                       w_00p 

        var_f(i,j,k) = r_simp / weight

    end do

end if        

! Edge 7
if (       myfront == mpi_proc_null &
     .and. myright == mpi_proc_null    ) then
    
    i = iend
    j = jend

    do k = k_mysta , k_myend
    
        w_00m  = one / aj ( i , j , k-1 )
        w_000  = one / aj ( i , j , k   )
        w_00p  = one / aj ( i , j , k+1 )    

        r_simp =       w_00m * var ( i , j , k-1 )  +  &
                  c4 * w_000 * var ( i , j , k   )  +  & 
                       w_00p * var ( i , j , k+1 ) 

        weight =       w_00m +  &
                  c4 * w_000 +  & 
                       w_00p 

        var_f(i,j,k) = r_simp / weight

    end do

end if     

! Edge 8
if (       myback  == mpi_proc_null & 
     .and. myright == mpi_proc_null   ) then
    
    i = ista
    j = jend

    do k = k_mysta , k_myend
    
        w_00m  = one / aj ( i , j , k-1 )
        w_000  = one / aj ( i , j , k   )
        w_00p  = one / aj ( i , j , k+1 )    

        r_simp =       w_00m * var ( i , j , k-1 )  +  &
                  c4 * w_000 * var ( i , j , k   )  +  & 
                       w_00p * var ( i , j , k+1 ) 

        weight =       w_00m +  &
                  c4 * w_000 +  & 
                       w_00p 

        var_f(i,j,k) = r_simp / weight

    end do

end if        


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!  ######  ####### #     # #     # ######     #    ######  #     # 
!  #     # #     # #     # # #   # #     #  #   #  #     #   # #   
!  ######  #     # #     # #  #  # #     # #     # ######     #    
!  #     # #     # #     # #   # # #     # ####### #   #      #    
!  ######  #######  #####  #     # ######  #     # #     #    #    
!                                                                  
!  #######    #     #####  #######  #####  
!  #        #   #  #       #       #       
!  #####   #     # #       #####    #####  
!  #       ####### #       #             # 
!  #       #     #  #####  #######  #####  
!                                        
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


!Normal i

if ( myback  == mpi_proc_null ) then
    
    i = ista

    do k = k_mysta , k_myend
    do j = j_mysta , j_myend
    
        w_000  = one / aj ( i , j   , k   )
        w_00p  = one / aj ( i , j   , k+1 )
        w_0p0  = one / aj ( i , j+1 , k   )
        w_00m  = one / aj ( i , j   , k-1 )
        w_0m0  = one / aj ( i , j-1 , k   )
        w_0pp  = one / aj ( i , j+1 , k+1 )
        w_0mm  = one / aj ( i , j-1 , k-1 )
        w_0mp  = one / aj ( i , j-1 , k+1 )
        w_0pm  = one / aj ( i , j+1 , k-1 )

        r_simp = c16 * w_000 * var ( i , j   , k   ) + &
                  c4 * w_00p * var ( i , j   , k+1 ) + &
                  c4 * w_0p0 * var ( i , j+1 , k   ) + &
                  c4 * w_00m * var ( i , j   , k-1 ) + &
                  c4 * w_0m0 * var ( i , j-1 , k   ) + &
                       w_0pp * var ( i , j+1 , k+1 ) + &
                       w_0mm * var ( i , j-1 , k-1 ) + &
                       w_0mp * var ( i , j-1 , k+1 ) + &
                       w_0pm * var ( i , j+1 , k-1 )
        
        weight = c16 * w_000 + &
                  c4 * w_00p + &
                  c4 * w_0p0 + &
                  c4 * w_00m + &
                  c4 * w_0m0 + &
                       w_0pp + &
                       w_0mm + &
                       w_0mp + &
                       w_0pm

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if        

if ( myfront == mpi_proc_null )  then
    
    i = iend

    do k = k_mysta , k_myend
    do j = j_mysta , j_myend
    
        w_000  = one / aj ( i , j   , k   )
        w_00p  = one / aj ( i , j   , k+1 )
        w_0p0  = one / aj ( i , j+1 , k   )
        w_00m  = one / aj ( i , j   , k-1 )
        w_0m0  = one / aj ( i , j-1 , k   )
        w_0pp  = one / aj ( i , j+1 , k+1 )
        w_0mm  = one / aj ( i , j-1 , k-1 )
        w_0mp  = one / aj ( i , j-1 , k+1 )
        w_0pm  = one / aj ( i , j+1 , k-1 )

        r_simp = c16 * w_000 * var ( i , j   , k   ) + &
                  c4 * w_00p * var ( i , j   , k+1 ) + &
                  c4 * w_0p0 * var ( i , j+1 , k   ) + &
                  c4 * w_00m * var ( i , j   , k-1 ) + &
                  c4 * w_0m0 * var ( i , j-1 , k   ) + &
                       w_0pp * var ( i , j+1 , k+1 ) + &
                       w_0mm * var ( i , j-1 , k-1 ) + &
                       w_0mp * var ( i , j-1 , k+1 ) + &
                       w_0pm * var ( i , j+1 , k-1 )
        
        weight = c16 * w_000 + &
                  c4 * w_00p + &
                  c4 * w_0p0 + &
                  c4 * w_00m + &
                  c4 * w_0m0 + &
                       w_0pp + &
                       w_0mm + &
                       w_0mp + &
                       w_0pm

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if        

!Normal j

if ( myleft  == mpi_proc_null ) then
    
    j = jsta

    do k = k_mysta , k_myend
    do i = i_mysta , i_myend

        w_000  = one / aj ( i   , j , k   )
        w_00p  = one / aj ( i   , j , k+1 )
        w_p00  = one / aj ( i+1 , j , k   )
        w_00m  = one / aj ( i   , j , k-1 )
        w_m00  = one / aj ( i-1 , j , k   )
        w_p0p  = one / aj ( i+1 , j , k+1 )
        w_m0m  = one / aj ( i-1 , j , k-1 )
        w_m0p  = one / aj ( i-1 , j , k+1 )
        w_p0m  = one / aj ( i+1 , j , k-1 )
    
        r_simp = c16 * w_000 * var ( i   , j , k   ) + & 
                  c4 * w_00p * var ( i   , j , k+1 ) + & 
                  c4 * w_p00 * var ( i+1 , j , k   ) + & 
                  c4 * w_00m * var ( i   , j , k-1 ) + & 
                  c4 * w_m00 * var ( i-1 , j , k   ) + & 
                       w_p0p * var ( i+1 , j , k+1 ) + & 
                       w_m0m * var ( i-1 , j , k-1 ) + & 
                       w_m0p * var ( i-1 , j , k+1 ) + & 
                       w_p0m * var ( i+1 , j , k-1 ) 

        weight = c16 * w_000 + & 
                  c4 * w_00p + & 
                  c4 * w_p00 + & 
                  c4 * w_00m + & 
                  c4 * w_m00 + & 
                       w_p0p + & 
                       w_m0m + & 
                       w_m0p + & 
                       w_p0m 

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if        

if ( myright == mpi_proc_null ) then
    
    j = jend

    do k = k_mysta , k_myend
    do i = i_mysta , i_myend

        w_000  = one / aj ( i   , j , k   )
        w_00p  = one / aj ( i   , j , k+1 )
        w_p00  = one / aj ( i+1 , j , k   )
        w_00m  = one / aj ( i   , j , k-1 )
        w_m00  = one / aj ( i-1 , j , k   )
        w_p0p  = one / aj ( i+1 , j , k+1 )
        w_m0m  = one / aj ( i-1 , j , k-1 )
        w_m0p  = one / aj ( i-1 , j , k+1 )
        w_p0m  = one / aj ( i+1 , j , k-1 )
    
        r_simp = c16 * w_000 * var ( i   , j , k   ) + & 
                  c4 * w_00p * var ( i   , j , k+1 ) + & 
                  c4 * w_p00 * var ( i+1 , j , k   ) + & 
                  c4 * w_00m * var ( i   , j , k-1 ) + & 
                  c4 * w_m00 * var ( i-1 , j , k   ) + & 
                       w_p0p * var ( i+1 , j , k+1 ) + & 
                       w_m0m * var ( i-1 , j , k-1 ) + & 
                       w_m0p * var ( i-1 , j , k+1 ) + & 
                       w_p0m * var ( i+1 , j , k-1 ) 

        weight = c16 * w_000 + & 
                  c4 * w_00p + & 
                  c4 * w_p00 + & 
                  c4 * w_00m + & 
                  c4 * w_m00 + & 
                       w_p0p + & 
                       w_m0m + & 
                       w_m0p + & 
                       w_p0m 

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if

!Normal k

if ( mydown  == mpi_proc_null ) then
    
    k = ksta

    do j = j_mysta , j_myend
    do i = i_mysta , i_myend

        w_000  = one / aj ( i   , j   , k )
        w_p00  = one / aj ( i+1 , j   , k )
        w_0p0  = one / aj ( i   , j+1 , k )
        w_m00  = one / aj ( i-1 , j   , k )
        w_0m0  = one / aj ( i   , j-1 , k )
        w_pp0  = one / aj ( i+1 , j+1 , k )
        w_mm0  = one / aj ( i-1 , j-1 , k )
        w_pm0  = one / aj ( i+1 , j-1 , k )
        w_mp0  = one / aj ( i-1 , j+1 , k )

        r_simp = c16 * w_000 * var ( i   , j   , k ) + &
                  c4 * w_p00 * var ( i+1 , j   , k ) + &
                  c4 * w_0p0 * var ( i   , j+1 , k ) + &
                  c4 * w_m00 * var ( i-1 , j   , k ) + &
                  c4 * w_0m0 * var ( i   , j-1 , k ) + &
                       w_pp0 * var ( i+1 , j+1 , k ) + &
                       w_mm0 * var ( i-1 , j-1 , k ) + &
                       w_pm0 * var ( i+1 , j-1 , k ) + &
                       w_mp0 * var ( i-1 , j+1 , k ) 

        weight = c16 * w_000 + &
                  c4 * w_p00 + &
                  c4 * w_0p0 + &
                  c4 * w_m00 + &
                  c4 * w_0m0 + &
                       w_pp0 + &
                       w_mm0 + &
                       w_pm0 + &
                       w_mp0 

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if        

if ( myup    == mpi_proc_null ) then
    
    k = kend

    do j = j_mysta , j_myend
    do i = i_mysta , i_myend

        w_000  = one / aj ( i   , j   , k )
        w_p00  = one / aj ( i+1 , j   , k )
        w_0p0  = one / aj ( i   , j+1 , k )
        w_m00  = one / aj ( i-1 , j   , k )
        w_0m0  = one / aj ( i   , j-1 , k )
        w_pp0  = one / aj ( i+1 , j+1 , k )
        w_mm0  = one / aj ( i-1 , j-1 , k )
        w_pm0  = one / aj ( i+1 , j-1 , k )
        w_mp0  = one / aj ( i-1 , j+1 , k )

        r_simp = c16 * w_000 * var ( i   , j   , k ) + &
                  c4 * w_p00 * var ( i+1 , j   , k ) + &
                  c4 * w_0p0 * var ( i   , j+1 , k ) + &
                  c4 * w_m00 * var ( i-1 , j   , k ) + &
                  c4 * w_0m0 * var ( i   , j-1 , k ) + &
                       w_pp0 * var ( i+1 , j+1 , k ) + &
                       w_mm0 * var ( i-1 , j-1 , k ) + &
                       w_pm0 * var ( i+1 , j-1 , k ) + &
                       w_mp0 * var ( i-1 , j+1 , k ) 

        weight = c16 * w_000 + &
                  c4 * w_p00 + &
                  c4 * w_0p0 + &
                  c4 * w_m00 + &
                  c4 * w_0m0 + &
                       w_pp0 + &
                       w_mm0 + &
                       w_pm0 + &
                       w_mp0 

        var_f(i,j,k) = r_simp / weight
            
    end do
    end do

end if        

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!        ######  #          #    #     # #    # ### #     #  #####  
!        #     # #         # #   ##    # #   #   #  ##    # #     # 
!        #     # #        #   #  # #   # #  #    #  # #   # #       
!        ######  #       #     # #  #  # ###     #  #  #  # #  #### 
!        #     # #       ####### #   # # #  #    #  #   # # #     # 
!        #     # #       #     # #    ## #   #   #  #    ## #     # 
!        ######  ####### #     # #     # #    # ### #     #  #####  
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                                                            
if ( nblk /= 0 ) then

    do nb = 1,nblk

        ! Let's zero all the variables in the interior region of the blanking zone
        do k = max( k_mysta , li_blk_ka(1,nb)+1 ) , min( k_myend , li_blk_kb(1,nb)-1 )
        do j = max( j_mysta , li_blk_ja(1,nb)+1 ) , min( j_myend , li_blk_jb(1,nb)-1 )
        do i = max( i_mysta , li_blk_ia(1,nb)+1 ) , min( i_myend , li_blk_ib(1,nb)-1 )

            var_f(i,j,k) = var(i,j,k)

        end do
        end do
        end do

        i = li_blk_ia(1,nb)

        do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) )
        do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) )
    
            w_000  = one / aj ( i , j   , k   )
            w_00p  = one / aj ( i , j   , k+1 )
            w_0p0  = one / aj ( i , j+1 , k   )
            w_00m  = one / aj ( i , j   , k-1 )
            w_0m0  = one / aj ( i , j-1 , k   )
            w_0pp  = one / aj ( i , j+1 , k+1 )
            w_0mm  = one / aj ( i , j-1 , k-1 )
            w_0mp  = one / aj ( i , j-1 , k+1 )
            w_0pm  = one / aj ( i , j+1 , k-1 )
    
            r_simp = c16 * w_000 * var ( i , j   , k   ) + &
                      c4 * w_00p * var ( i , j   , k+1 ) + &
                      c4 * w_0p0 * var ( i , j+1 , k   ) + &
                      c4 * w_00m * var ( i , j   , k-1 ) + &
                      c4 * w_0m0 * var ( i , j-1 , k   ) + &
                           w_0pp * var ( i , j+1 , k+1 ) + &
                           w_0mm * var ( i , j-1 , k-1 ) + &
                           w_0mp * var ( i , j-1 , k+1 ) + &
                           w_0pm * var ( i , j+1 , k-1 )
            
            weight = c16 * w_000 + &
                      c4 * w_00p + &
                      c4 * w_0p0 + &
                      c4 * w_00m + &
                      c4 * w_0m0 + &
                           w_0pp + &
                           w_0mm + &
                           w_0mp + &
                           w_0pm
    
            var_f(i,j,k) = r_simp / weight
          
        end do
        end do


        i = li_blk_ib(1,nb)

        do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) )
        do j = max( j_mysta , li_blk_ja(1,nb) ) , min( j_myend , li_blk_jb(1,nb) )

            w_000  = one / aj ( i , j   , k   )
            w_00p  = one / aj ( i , j   , k+1 )
            w_0p0  = one / aj ( i , j+1 , k   )
            w_00m  = one / aj ( i , j   , k-1 )
            w_0m0  = one / aj ( i , j-1 , k   )
            w_0pp  = one / aj ( i , j+1 , k+1 )
            w_0mm  = one / aj ( i , j-1 , k-1 )
            w_0mp  = one / aj ( i , j-1 , k+1 )
            w_0pm  = one / aj ( i , j+1 , k-1 )
    
            r_simp = c16 * w_000 * var ( i , j   , k   ) + &
                      c4 * w_00p * var ( i , j   , k+1 ) + &
                      c4 * w_0p0 * var ( i , j+1 , k   ) + &
                      c4 * w_00m * var ( i , j   , k-1 ) + &
                      c4 * w_0m0 * var ( i , j-1 , k   ) + &
                           w_0pp * var ( i , j+1 , k+1 ) + &
                           w_0mm * var ( i , j-1 , k-1 ) + &
                           w_0mp * var ( i , j-1 , k+1 ) + &
                           w_0pm * var ( i , j+1 , k-1 )
            
            weight = c16 * w_000 + &
                      c4 * w_00p + &
                      c4 * w_0p0 + &
                      c4 * w_00m + &
                      c4 * w_0m0 + &
                           w_0pp + &
                           w_0mm + &
                           w_0mp + &
                           w_0pm
    
            var_f(i,j,k) = r_simp / weight

        end do
        end do



        j = li_blk_ja(1,nb)

        do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
        do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 
    
            w_000  = one / aj ( i   , j , k   )
            w_00p  = one / aj ( i   , j , k+1 )
            w_p00  = one / aj ( i+1 , j , k   )
            w_00m  = one / aj ( i   , j , k-1 )
            w_m00  = one / aj ( i-1 , j , k   )
            w_p0p  = one / aj ( i+1 , j , k+1 )
            w_m0m  = one / aj ( i-1 , j , k-1 )
            w_m0p  = one / aj ( i-1 , j , k+1 )
            w_p0m  = one / aj ( i+1 , j , k-1 )
        
            r_simp = c16 * w_000 * var ( i   , j , k   ) + & 
                      c4 * w_00p * var ( i   , j , k+1 ) + & 
                      c4 * w_p00 * var ( i+1 , j , k   ) + & 
                      c4 * w_00m * var ( i   , j , k-1 ) + & 
                      c4 * w_m00 * var ( i-1 , j , k   ) + & 
                           w_p0p * var ( i+1 , j , k+1 ) + & 
                           w_m0m * var ( i-1 , j , k-1 ) + & 
                           w_m0p * var ( i-1 , j , k+1 ) + & 
                           w_p0m * var ( i+1 , j , k-1 ) 
    
            weight = c16 * w_000 + & 
                      c4 * w_00p + & 
                      c4 * w_p00 + & 
                      c4 * w_00m + & 
                      c4 * w_m00 + & 
                           w_p0p + & 
                           w_m0m + & 
                           w_m0p + & 
                           w_p0m 
    
            var_f(i,j,k) = r_simp / weight

        end do
        end do

        j = li_blk_jb(1,nb)

        do k = max( k_mysta , li_blk_ka(1,nb) ) , min( k_myend , li_blk_kb(1,nb) ) 
        do i = max( i_mysta , li_blk_ia(1,nb) ) , min( i_myend , li_blk_ib(1,nb) ) 

            w_000  = one / aj ( i   , j , k   )
            w_00p  = one / aj ( i   , j , k+1 )
            w_p00  = one / aj ( i+1 , j , k   )
            w_00m  = one / aj ( i   , j , k-1 )
            w_m00  = one / aj ( i-1 , j , k   )
            w_p0p  = one / aj ( i+1 , j , k+1 )
            w_m0m  = one / aj ( i-1 , j , k-1 )
            w_m0p  = one / aj ( i-1 , j , k+1 )
            w_p0m  = one / aj ( i+1 , j , k-1 )
        
            r_simp = c16 * w_000 * var ( i   , j , k   ) + & 
                      c4 * w_00p * var ( i   , j , k+1 ) + & 
                      c4 * w_p00 * var ( i+1 , j , k   ) + & 
                      c4 * w_00m * var ( i   , j , k-1 ) + & 
                      c4 * w_m00 * var ( i-1 , j , k   ) + & 
                           w_p0p * var ( i+1 , j , k+1 ) + & 
                           w_m0m * var ( i-1 , j , k-1 ) + & 
                           w_m0p * var ( i-1 , j , k+1 ) + & 
                           w_p0m * var ( i+1 , j , k-1 ) 
    
            weight = c16 * w_000 + & 
                      c4 * w_00p + & 
                      c4 * w_p00 + & 
                      c4 * w_00m + & 
                      c4 * w_m00 + & 
                           w_p0p + & 
                           w_m0m + & 
                           w_m0p + & 
                           w_p0m 
    
            var_f(i,j,k) = r_simp / weight

        end do
        end do 

    end do

end if

 
end function testfilter_simpson_3d

