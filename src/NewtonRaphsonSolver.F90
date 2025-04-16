module NewtonRaphsonSolver

! TO DO: Have a look over all this routine, as I'm not sure it's doing what it's meant to.
! Probably there are parameters or steps that are not working properly.
! I'm gonna add an option to control.dat to perform the geometric redistancing by distancing
! only for the moment (JS, 18 / Jun / 2024)

use precision
use global_mpi 
use DataTypes
use TetrahedronMethods
use global_lsm, only: ConvergenceToleranceGeomReini

implicit none

contains

subroutine NewtonRaphson( ntetrahedra, KTetrahedraList, C_initial_guest, C_i, &
                          TotalWaterVolume )

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This routine computes three interpolation coefficients using a trilinear interpolation 
   ! described in 
   ! 
   ! Delandmeter, P., & Sebille, E. V. (2019). The Parcels v2. 0 Lagrangian framework: new 
   ! field interpolation schemes. Geoscientific Model Development, 12(8), 3571-3584.
   !
   ! It solves a three-equations non-linear system using a globally-convergent Newton-Rhapson 
   ! method, whose implementation is thoroughly described in
   
   ! Numerical receipes in Fortran, Vol I, by W. H. Press, 2nd ed., p379
   
   ! The routine returns three interpolation coefficients between 0 and 1 contained in the
   ! variable coeff
   !
   ! Jorge Sandoval, UoE/PUC. Edinburgh, July, 2022.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   implicit none
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! input-variables
   ! ------------------
   ! KTetrahedraList  : is the list with all the tetrahedrons obtained by the triangulation step
   !                   it's used for computing the volume enclosed by the modified phi 
   !                   distribution
   ! C_initial_guest : is the initial value of C for the solver
   ! 
   ! ------------------
   ! output-variables
   ! ------------------
   !
   ! C_i             : is the converged solution for C
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   integer, intent(in) :: ntetrahedra
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
   REAL(kind = rdf), intent(in)  :: C_initial_guest
   REAL(kind = rdf), intent(inout) :: C_i
   real (kind = rdf) , intent(in) :: TotalWaterVolume

   ! local variables for Newton-Rhapson subroutine

   REAL(kind = rdf) :: C_initial_guest_vec(1)
   REAL(kind = rdf) :: TotalTetrahedraVolume

   integer :: nmax,MAXITS,nzones
   parameter(nmax=2000000,MAXITS=20) !200)
   REAL(KIND = RDF):: rst(3)      
   logical::check      
   integer:: K, m, its,itss

   ! Initialisation of C_i as the initial guest given to the method
   C_initial_guest_vec(1) = C_initial_guest
   C_i                    = C_initial_guest
   
   ! its initialisation
   its = 1

   ! Total Volume enclosed by the linear reconstruction of ϕ, ϕh* (this is ϕ distribution
   ! after the advection step )

   TotalTetrahedraVolume = TotalWaterVolume


   call newt( ntetrahedra , KTetrahedraList , TotalTetrahedraVolume, &
              C_initial_guest_vec , 1 , check , MAXITS , its , itss    )

   ! If the solver converges, the new value of C_i is the updated value
   ! of C_initial_guest after the Newton - Rapshon method
   
   ! check == false --> the method CONVERGED and there's nothing to check

   if( its < MAXITS .and. .not.(check)) then
      C_i = C_initial_guest_vec(1)
      if ( myid == root ) print *, 'NewtonRhapson converged in ', its, ' its'
   else ! I drop the correction
      if ( myid == root ) print *, 'NewtonRhapson DIDNT converged. It ran ', its, ' its'
      C_i = zero
   end if

end subroutine NewtonRaphson

!-------------------------------------------------------------------------------
!
!       A globally convergent Newton-Raphson method (from Numerical
!       receipes in Fortran, Vol I, by  W. H. Press, 2nd ed., p379 ) 
!       Parameter: NP--dimension of the sysytem, MAXITS--maximum 
!       number of iterations within which convergence shall be 
!       obtained if it is supposed. 
!
!       * newt, fdjac, fmin, lnsrch, ludcmp, lubksb were obtained 
!         directly from the reference presented above.
!
!       * funcv and fdjac1 were defined specifically for the equations 
!         to be solved, but the guidelines and in-out parameters are 
!         described in the same reference.
!
!       * fmin was modified and was converted from a function to a subroutine.
!         it was because the use of external variables was generating 
!         compilation issues because of the references (Jorge Sandoval,
!         Edinburgh, November 8th, 2021)
!
!-------------------------------------------------------------------------------

SUBROUTINE newt(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, x,n,check,MAXITS,its,itss)
  
  integer, intent(in) :: ntetrahedra
  type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
  REAL(kind = rdf), intent(in) :: TotalTetrahedraVolume

  integer :: its, itss
  INTEGER n,nn,NP
  LOGICAL check
  REAL(KIND = RDF) x(n),fvec,TOLF,TOLMIN,TOLX,STPMX

  PARAMETER (NP=40,TOLF=1.e-12,TOLMIN=1.e-12,TOLX=1.e-12,   &
       &  STPMX=100.)
  COMMON /newtv/ fvec(NP),nn
  SAVE /newtv/
!C     USES fdjac,fmin,lnsrch,lubksb,ludcmp

  INTEGER i,j,indx(NP),MAXITS
  REAL(KIND = RDF) d,den,f,fold,stpmax,sum_aux,temp,test,fjac(NP,NP),g(NP),p(NP), &
       &  xold(NP),fmin_result

!  EXTERNAL fmin

  nn=n

  ! f and fmin_result initialisation
  f = zero
  fmin_result = zero

  ! fmin calls funcv internally, so this call update fvec value 
  call fmin( ntetrahedra , KTetrahedraList , TotalTetrahedraVolume , x , fmin_result )

  f = fmin_result

  test = zero

  do i = 1,n
     if( abs(fvec(i)) > test ) test = abs(fvec(i))
  end do

  if(test < 0.01_rdf*TOLF .and. abs(fvec(1)) < ConvergenceToleranceGeomReini ) then
     check=.false.
     return
  endif

  sum_aux = zero
  
  do i=1,n    
     sum_aux = sum_aux + x(i)**2
  end do

  stpmax = STPMX * max( sqrt(sum_aux) , float(n) )

  do its = 1 , MAXITS
     
     !write(*,*) its

     ! Jacobian computation
     call fdjac(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, n,x,fvec,NP,fjac)     

     ! If the Jacobian is about zero, we exit the Newton-Rhapson method
     
     sum_aux = zero
     
     do i = 1,n
      do j = 1,n
         sum_aux = sum_aux + fjac(j,i)
      end do
     end do

     if ( abs( sum_aux ) < ConvergenceToleranceGeomReini ) then 

      print *, ' '
      print *, 'sum_aux = ', sum_aux
      print *, 'fvec    = ', fvec(1)
      print *, 'zero Jacobian' 
      return

     end if

     call funcv(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, n,x,fvec)

     do i=1,n
        sum_aux = zero
        do j=1,n
           sum_aux = sum_aux + fjac(j,i)*fvec(j)
        end do
        g(i) = sum_aux
     end do
     
     do i=1,n
        xold(i) = x(i)
     end do

     fold = f
     
     do i=1,n
        p(i)=-fvec(i)
     end do

     ! Solving the linear system J * δx = F

     call ludcmp(fjac,n,NP,indx,d,itss)
     call lubksb(fjac,n,NP,indx,p)
     call lnsrch(ntetrahedra, KTetrahedraList,TotalTetrahedraVolume, n,xold,fold,g,p,x,f,stpmax,check)

     test=zero

     do i=1,n
        if( abs(fvec(i)) > test ) test = abs(fvec(i))
     end do

     if( test < TOLF .and. abs(fvec(1)) < ConvergenceToleranceGeomReini )then
        check=.false.
        return
     endif

     if(check) then ! spourious convergence
        print *, 'spourious convergence, fvec = ', fvec(1)
        test = zero
        den = max( f , one/two * n)
        
        do i=1,n
           temp = abs(g(i)) * max( abs(x(i)) , one )/den
           if(temp > test) test = temp
        end do

        if(test < TOLMIN)then
           check=.true.
        else
           check=.false.
        endif
        
        return
     
     endif

     test = zero
     
     do i=1,n
        temp = ( abs( x(i) - xold(i) ) )/max( abs(x(i)) , one )
        if(temp > test) test = temp
     end do

     if(test < TOLX .and. abs(fvec(1)) < ConvergenceToleranceGeomReini ) return
  
  end do

  if ( myid == root )  print *, 'MAXITS exceeded in newt, fvec = ', fvec(1)

END SUBROUTINE newt


SUBROUTINE fdjac(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, n,x,fvec,np,df)

  integer, intent(in) :: ntetrahedra
  type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
  REAL(kind = rdf), intent(in) :: TotalTetrahedraVolume

  INTEGER n,np,NMAX
  REAL(KIND = RDF) df(np,np),fvec(n),x(n),EPS
  !PARAMETER (NMAX=40,EPS=1.e-4)
  PARAMETER (NMAX=40,EPS=sqrt(eps_sims) )

!    USES funcv
  INTEGER i,j
  REAL(KIND = RDF) h,temp,f(NMAX)

  do j=1,n

     temp=x(j)

     h=EPS*abs(temp)

     !if(h.eq.0.) h=EPS
     if( h < eps_sims ) h=EPS

     x(j)=temp+h

     h=x(j)-temp
     
     call funcv(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume,n,x,f)
     
     x(j)=temp
     
     do i=1,n
        df(i,j)=(f(i)-fvec(i))/h
     end do
  
  end do
  return

END SUBROUTINE fdjac


SUBROUTINE fmin(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, x, fmin_result)

  integer, intent(in) :: ntetrahedra
  type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
  REAL(kind = rdf), intent(in) :: TotalTetrahedraVolume

  INTEGER n,NP
  REAL(KIND = RDF) , intent(out) :: fmin_result
  REAL(KIND = RDF) x(*),fvec
  
  PARAMETER (NP=40)
  COMMON /newtv/ fvec(NP),n
  SAVE /newtv/
!CU    USES funcv
  INTEGER i
  REAL(KIND = RDF) sum_aux

  fmin_result = zero

  call funcv(ntetrahedra, KTetrahedraList,TotalTetrahedraVolume,n,x,fvec)
  sum_aux=0.
  do i=1,n
     sum_aux=sum_aux+fvec(i)**2
  end do
  fmin_result=0.5*sum_aux
  return
END SUBROUTINE fmin

SUBROUTINE lnsrch(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, n,xold,fold,g,p,x,f,stpmax,check)

  integer, intent(in) :: ntetrahedra
  type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
  REAL(kind = rdf), intent(in) :: TotalTetrahedraVolume

  INTEGER n
  LOGICAL check
  REAL(KIND = RDF) f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF,TOLX, fmin_result
  PARAMETER (ALF=1.e-4,TOLX=1.e-7)
  !EXTERNAL func
!CU    USES func
  INTEGER i
  REAL(KIND = RDF) a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum_aux,temp, &
       &  test,tmplam
  check=.false.
  sum_aux=0.
  do i=1,n
     sum_aux=sum_aux+p(i)*p(i)
  end do
  sum_aux=sqrt(sum_aux)
  if(sum_aux > stpmax)then
     do i=1,n
        p(i)=p(i)*stpmax/sum_aux
     end do
  endif
  slope=0.
  do i=1,n
     slope=slope+g(i)*p(i)
  end do
  !if (slope.ge.0.) pause 'roundoff problem in lnsrch'
  test=0.
  do i=1,n
     temp=abs(p(i))/max(abs(xold(i)),1.)
     if(temp > test)test=temp
  end do
  alamin=TOLX/test
  alam=1.

  alam2=2.
  f2=2.
  fold2=fold

1     continue
  do i=1,n
     x(i)=xold(i)+alam*p(i)
  end do

  fmin_result = zero

  f = zero
  fmin_result = zero

  call fmin(ntetrahedra, KTetrahedraList,TotalTetrahedraVolume, x,fmin_result)
  f = fmin_result

  if(alam < alamin)then
     do i=1,n
        x(i)=xold(i)
     end do
     check=.true.
     return
  else if(f.le.fold+ALF*alam*slope)then
     return
  else
     if(alam.eq.1.)then
        tmplam=-slope/(2.*(f-fold-slope))
     else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if(a.eq.0.)then
           tmplam=-slope/(2.*b)
        else
           disc=b*b-3.*a*slope
           if(disc < 0.) then
              tmplam=0.5*alam
           else if (b.le.0.) then
              tmplam=(-b+sqrt(disc))/(3.*a)
           else
              tmplam=-slope/(b+sqrt(disc))
           endif
!           C              tmplam=(-b+sqrt(disc))/(3.*a)
        endif
        if(tmplam > .5*alam)tmplam=.5*alam
     endif
  endif
  alam2=alam
  f2=f
!  fold2=fold
  alam=max(tmplam,.1*alam)
  goto 1
END SUBROUTINE lnsrch

SUBROUTINE ludcmp(a,n,np,indx,d,itss)
  integer :: itss
  INTEGER n,np,indx(n),NMAX
  REAL(KIND = RDF) d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL(KIND = RDF) aamax,dum,sum_aux,vv(NMAX)
  itss=0
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)) > aamax) aamax=abs(a(i,j))
     end do
     !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
     if (aamax.eq.0.) itss=1
     if (aamax.eq.0.) return
     vv(i)=1./aamax
  end do
  do j=1,n
     do i=1,j-1
        sum_aux=a(i,j)
        do k=1,i-1
           sum_aux=sum_aux-a(i,k)*a(k,j)
        end do
        a(i,j)=sum_aux
     end do
     aamax=0.
     do i=j,n
        sum_aux=a(i,j)
        do k=1,j-1
           sum_aux=sum_aux-a(i,k)*a(k,j)
        end do
        a(i,j)=sum_aux
        dum=vv(i)*abs(sum_aux)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     end do
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        end do
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        end do
     endif
  end do
  return
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  REAL(KIND = RDF) a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL(KIND = RDF) sum_aux
  ii=0
  do i=1,n
     ll=indx(i)
     sum_aux=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum_aux=sum_aux-a(i,j)*b(j)
        end do
     else if (sum_aux.ne.0.) then
        ii=i
     endif
     b(i)=sum_aux
  end do
  do i=n,1,-1
     sum_aux=b(i)
     do j=i+1,n
        sum_aux=sum_aux-a(i,j)*b(j)
     end do
     b(i)=sum_aux/a(i,i)
  end do
  return
END SUBROUTINE lubksb


SUBROUTINE funcv(ntetrahedra, KTetrahedraList, TotalTetrahedraVolume, n, C_i, fvec)
  
   ! This is the left hand side of the equation 
   ! funcv(x_i) = 0 <--> ΔV( ϕh, ϕh* + C * ξh) = 0

   integer, intent(in) :: ntetrahedra
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: KTetrahedraList
   REAL(kind = rdf), intent(in) :: TotalTetrahedraVolume


   integer          :: n ! dummy argument in this case because the
                         ! the incognito variable is scalar
   REAL(kind = rdf) :: C_i(1)
   REAL(kind = rdf) :: fvec(1)

   ! dummy variables for Marching Tetrahedra

   REAL(kind = rdf), dimension(4,3)    :: VerticesCoordinates
   REAL(kind = rdf), dimension (4)     :: PhiAux
   REAL(kind = rdf), dimension (2,3,3) :: VerticesIsosurfaces
   integer                             :: ntriangles
   REAL(kind = rdf)                    :: IsosurfaceArea
   REAL(kind = rdf)                    :: WaterVolumeAux
   REAL(kind = rdf), dimension(4)      :: VerticesDistances
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! local variables
   REAL(kind = rdf) :: TotalVolumeAux , SDAux , CAux
   
   ! Global volume after communication

   real( kind = rdf ) :: TotalVolumeGlobal

   integer :: K

   ! loop over all the tetrahedra set to compute the total volume enclosed
   ! by the modified phi distribution

   ! Total Volume enclosed by the new phi distribution initialisation
   TotalVolumeAux = zero

   do K = 1, ntetrahedra

      ! I just consider the volume inside internal tetrahedrons to
      ! avoid considering some volumes more than once
      if ( .not. KTetrahedraList(K)%InternalTetrahedron ) cycle

      VerticesCoordinates(1,:) = (/ KTetrahedraList(K)%v1%x, &
                                    KTetrahedraList(K)%v1%y, &
                                    KTetrahedraList(K)%v1%z    /)  
   
      VerticesCoordinates(2,:) = (/ KTetrahedraList(K)%v2%x, &
                                    KTetrahedraList(K)%v2%y, &
                                    KTetrahedraList(K)%v2%z    /)  
      
      VerticesCoordinates(3,:) = (/ KTetrahedraList(K)%v3%x, &
                                    KTetrahedraList(K)%v3%y, &
                                    KTetrahedraList(K)%v3%z    /)  
   
      VerticesCoordinates(4,:) = (/ KTetrahedraList(K)%v4%x, &
                                    KTetrahedraList(K)%v4%y, &
                                    KTetrahedraList(K)%v4%z    /)  
   
      ! PhiAux =  ϕh* + C * ξh or 0 if the corrections produces a sign change 
   
      SDAux = KTetrahedraList(K)%v1%SDistanceFreeSurface 
      CAux  = C_i(1) * KTetrahedraList(K)%v1%xiI

      PhiAux(1) = SDAux + CAux / two * abs( sign( one , SDAux ) + sign( one , SDAux + CAux) )


      SDAux = KTetrahedraList(K)%v2%SDistanceFreeSurface 
      CAux  = C_i(1) * KTetrahedraList(K)%v2%xiI

      PhiAux(2) = SDAux + CAux / two * abs( sign( one , SDAux ) + sign( one , SDAux + CAux) )


      SDAux = KTetrahedraList(K)%v3%SDistanceFreeSurface 
      CAux  = C_i(1) * KTetrahedraList(K)%v3%xiI

      PhiAux(3) = SDAux + CAux / two * abs( sign( one , SDAux ) + sign( one , SDAux + CAux) )


      SDAux = KTetrahedraList(K)%v4%SDistanceFreeSurface 
      CAux  = C_i(1) * KTetrahedraList(K)%v4%xiI

      PhiAux(4) = SDAux + CAux / two * abs( sign( one , SDAux ) + sign( one , SDAux + CAux) )


      call MarchingTetrahedron(  VerticesCoordinates    , &
                                 PhiAux                 , &
                                 VerticesIsosurfaces    , &
                                 ntriangles             , &
                                 IsosurfaceArea         , &
                                 WaterVolumeAux         , &
                                 VerticesDistances           )


      ! Update Total Volume enclosed by the modified phi distribution
      TotalVolumeAux = TotalVolumeAux + WaterVolumeAux

   end do ! loop over the whole tetrahedra set
   
   ! Use mpi_allreduce to calculate the volume difference

   call mpi_allreduce( TotalVolumeAux, TotalVolumeGlobal , &
                       1 , MPI_REAL_TYPE , mpi_sum , mpi_comm_world, ierr )


   ! fvec(C_i) = ΔV(C_i) = ΔV( ϕh, ϕh* + C * ξh)
   fvec(1) = TotalTetrahedraVolume - TotalVolumeGlobal 

   !print *, ' '
   !print *, 'VolumeKTetrahedra AUX = ', TotalVolumeAux

   return

END SUBROUTINE funcv


end module NewtonRaphsonSolver