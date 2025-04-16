subroutine get_trilinear_interp_coefficients(cell, point, dis, initial_guest_newt, coeff)

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
   ! Jorge Sandoval, UoE/PUC. Edinburgh, October 18th, 2021.
   ! j.sandoval@ed.ac.uk / jcsandov@uc.cl
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   implicit none

   !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! input-variables
   ! ------------------
   ! cell: (x,y,z) coordinates of the eight nodes that compose the cell
   ! point: (x,y,z) coordinates of the point to be interpolated
   ! dis: distance between local i1, im, j1, km, k1, km face and point
   ! initial_guest_newt: three initial values for the interpolation coefficients before
   !                     Newton-Rhapson method
   !
   ! ------------------
   ! output-variables
   ! ------------------
   !
   ! coeff: three trilinear interpolation coefficients 
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! NOTES:
   ! Nodes nomenclature for interpolation subroutines
   ! v1 = v000, v2 = v100, v3 = v110, v4 = v010
   ! v5 = v001, v6 = v101, v7 = v111, v8 = v011
   !
   ! where v101 means vertex coordinates at (i+1,j,k+1)
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real (kind = rdf), intent(in)  :: cell(8,3), point(3), dis(6), initial_guest_newt(3) 
   real (kind = rdf), intent(out) :: coeff(3)

   ! local variables for Newton-Rhapson subroutine

   real (kind = rdf) :: p(3),v000(3),v100(3),v010(3),v001(3),v101(3), & 
                         v011(3), v110(3),v111(3)

   integer :: nmax,MAXITS,nzones
   parameter(nmax=2000000,MAXITS=200)
   real:: rst(3)      
   logical::check      
   integer:: m, its,itss
   
   !-----------------------------------------------------------------------------------------
   ! The COMMON block, a piece of shared memory in the computer, is another method for 
   ! passing information between program units. Data stored in a COMMON block is not passed 
   ! between program units via argument lists, but through the COMMON statement near the 
   ! beginning of each program unit. 
   ! For further information: https://www.obliquity.com/computer/fortran/common.html

      common/intp/v000,v100,v010,v001,v101,v011, &
        &   v110,v111,p

   !-----------------------------------------------------------------------------------------

   rst = initial_guest_newt

   call newt(rst,3,check,MAXITS,its,itss)
   
   if(its<MAXITS)then
      do m=1,3
         coeff(m)=rst(m)
      end do
   else
   
   ! When Newton-Rhapson doesn't converge, the coefficiente is asigned as the portion of 
   ! the total distance (denominator is the total distance between opposite cell faces)

      do m = 1,3
         coeff(m)=dis(2*m-1)/(dis(2*m-1)+dis(2*m))
      end do

   end if ! its

end subroutine get_trilinear_interp_coefficients

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
!
!-------------------------------------------------------------------------------

SUBROUTINE newt(x,n,check,MAXITS,its,itss)
  
  !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
  integer :: its, itss
  INTEGER n,nn,NP
  LOGICAL check
  REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX

  PARAMETER (NP=40,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,   &
       &  STPMX=100.)
  COMMON /newtv/ fvec(NP),nn
  SAVE /newtv/
!C     USES fdjac,fmin,lnsrch,lubksb,ludcmp

  INTEGER i,j,indx(NP),MAXITS
  REAL d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP), &
       &  xold(NP),fmin_result

!  EXTERNAL fmin

  nn=n

  ! f and fmin_result initialisation
  
  f = zero
  fmin_result = zero

  call fmin(x, fmin_result)

  f=fmin_result

  test=0.
  do i=1,n
     if(abs(fvec(i)).gt.test)test=abs(fvec(i))
  end do
  if(test.lt..01*TOLF) then
     check=.false.
     return
  endif

  sum=0.
  do i=1,n    
     sum=sum+x(i)**2
  end do

  stpmax=STPMX*max(sqrt(sum),float(n))

  do its=1,MAXITS
     write(*,*) its
!     call fdjac(n,x,fvec,NP,fjac)     
     call fdjac1(n,x,fjac(1:n,1:n))
     call funcv(n,x,fvec)

!     print *,fjac(1:n,1:n)
     do i=1,n
        sum=0.
        do j=1,n
           sum=sum+fjac(j,i)*fvec(j)
        end do
        g(i)=sum
     end do
     do i=1,n
        xold(i)=x(i)
     end do
     fold=f
     do i=1,n
        p(i)=-fvec(i)
     end do
!     write(*,*) '********'
     call ludcmp(fjac,n,NP,indx,d,itss)
!     write(*,*) '*******1'
     call lubksb(fjac,n,NP,indx,p)
!     write(*,*) '*******2'
     call lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
!     write(*,*) '*******3'
     test=0.
     do i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
     end do
     if(test.lt.TOLF)then
        check=.false.
        return
     endif
     if(check)then
        test=0.
        den=max(f,.5*n)
        do i=1,n
           temp=abs(g(i))*max(abs(x(i)),1.)/den
           if(temp.gt.test)test=temp
        end do
        if(test.lt.TOLMIN)then
           check=.true.
        else
           check=.false.
        endif
        return
     endif
     test=0.
     do i=1,n
        temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
        if(temp.gt.test)test=temp
     end do
     if(test.lt.TOLX)return
  end do
!c      pause 'MAXITS exceeded in newt'
END SUBROUTINE newt


SUBROUTINE fdjac(n,x,fvec,np,df)
  INTEGER n,np,NMAX
  REAL df(np,np),fvec(n),x(n),EPS
  PARAMETER (NMAX=40,EPS=1.e-4)
!    USES funcv
  INTEGER i,j
  REAL h,temp,f(NMAX)
  do j=1,n
     temp=x(j)
     h=EPS*abs(temp)
     if(h.eq.0.)h=EPS
     x(j)=temp+h
     h=x(j)-temp
     call funcv(n,x,f)
     x(j)=temp
     do i=1,n
        df(i,j)=(f(i)-fvec(i))/h
     end do
  end do
  return
END SUBROUTINE fdjac


SUBROUTINE fmin(x, fmin_result)
  INTEGER n,NP
  REAL , intent(out) :: fmin_result
  REAL x(*),fvec
  
  PARAMETER (NP=40)
  COMMON /newtv/ fvec(NP),n
  SAVE /newtv/
!CU    USES funcv
  INTEGER i
  REAL sum

  fmin_result = zero

  call funcv(n,x,fvec)
  sum=0.
  do i=1,n
     sum=sum+fvec(i)**2
  end do
  fmin_result=0.5*sum
  return
END SUBROUTINE fmin

SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
  INTEGER n
  LOGICAL check
  REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF,TOLX, fmin_result
  PARAMETER (ALF=1.e-4,TOLX=1.e-7)
  !EXTERNAL func
!CU    USES func
  INTEGER i
  REAL a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp, &
       &  test,tmplam
  check=.false.
  sum=0.
  do i=1,n
     sum=sum+p(i)*p(i)
  end do
  sum=sqrt(sum)
  if(sum.gt.stpmax)then
     do i=1,n
        p(i)=p(i)*stpmax/sum
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
     if(temp.gt.test)test=temp
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

  call fmin(x,fmin_result)
  f = fmin_result

  if(alam.lt.alamin)then
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
           if(disc.lt.0.) then
              tmplam=0.5*alam
           else if (b.le.0.) then
              tmplam=(-b+sqrt(disc))/(3.*a)
           else
              tmplam=-slope/(b+sqrt(disc))
           endif
!           C              tmplam=(-b+sqrt(disc))/(3.*a)
        endif
        if(tmplam.gt..5*alam)tmplam=.5*alam
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
  REAL d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  REAL aamax,dum,sum,vv(NMAX)
  itss=0
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     end do
     !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
     if (aamax.eq.0.) itss=1
     if (aamax.eq.0.) return
     vv(i)=1./aamax
  end do
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
     end do
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
        dum=vv(i)*abs(sum)
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
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  end do
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do
  return
END SUBROUTINE lubksb

SUBROUTINE funcv(n,rst,fvec)
  
  !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4

  integer :: i,n
  real (kind = rdf) ::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
  common/intp/v000,v100,v010,v001,v101,v011, &
       &   v110,v111,p
  real::rst(n),fvec(n)
  do i=1,n
     fvec(i)=v000(i)*(1-rst(1))*(1-rst(2))*(1-rst(3))+   &
          &  v100(i)*rst(1)*(1-rst(2))*(1-rst(3))+ &
          &  v010(i)*(1-rst(1))*rst(2)*(1-rst(3))+ &
          &  v001(i)*(1-rst(1))*(1-rst(2))*rst(3)+ &
          &  v101(i)*rst(1)*(1-rst(2))*rst(3)+     &
          &  v011(i)*(1-rst(1))*rst(2)*rst(3)+     &
          &  v110(i)*rst(1)*rst(2)*(1-rst(3))+     &
          &  v111(i)*rst(1)*rst(2)*rst(3)-p(i)
  end do
  return
END SUBROUTINE funcv

SUBROUTINE fdjac1(n,rst,fjac)

  !integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4
  integer :: n
  real (kind = rdf) ::p(3),v000(3),v100(3),v010(3),v001(3),v101(3),v011(3),v110(3),v111(3)
  common/intp/v000,v100,v010,v001,v101,v011, &
       &   v110,v111,p
  real::rst(n),fjac(n,n)

  fjac(1,1)=-v000(1)*(1-rst(2))*(1-rst(3))+  &
       &  v100(1)*(1-rst(2))*(1-rst(3))-  &
       &  v010(1)*rst(2)*(1-rst(3))-   &
       &  v001(1)*(1-rst(2))*rst(3)+   &
       &  v101(1)*(1-rst(2))*rst(3)-   &
       &  v011(1)*rst(2)*rst(3)+ &
       &  v110(1)*rst(2)*(1-rst(3))+   &
       &  v111(1)*rst(2)*rst(3)
  fjac(2,1)=-v000(2)*(1-rst(2))*(1-rst(3))+  &
       &  v100(2)*(1-rst(2))*(1-rst(3))-  &
       &  v010(2)*rst(2)*(1-rst(3))-   &
       &  v001(2)*(1-rst(2))*rst(3)+   &
       &  v101(2)*(1-rst(2))*rst(3)-   &
       &  v011(2)*rst(2)*rst(3)+ &
       &  v110(2)*rst(2)*(1-rst(3))+   &
       &  v111(2)*rst(2)*rst(3)
  fjac(3,1)=-v000(3)*(1-rst(2))*(1-rst(3))+  &
       &  v100(3)*(1-rst(2))*(1-rst(3))-  &
       &  v010(3)*rst(2)*(1-rst(3))-   &
       &  v001(3)*(1-rst(2))*rst(3)+   &
       &  v101(3)*(1-rst(2))*rst(3)-   &
       &  v011(3)*rst(2)*rst(3)+ &
       &  v110(3)*rst(2)*(1-rst(3))+   &
       &  v111(3)*rst(2)*rst(3)
  fjac(1,2)=-v000(1)*(1-rst(1))*(1-rst(3))-  &
       &  v100(1)*rst(1)*(1-rst(3))+   &
       &  v010(1)*(1-rst(1))*(1-rst(3))-  &
       &  v001(1)*(1-rst(1))*rst(3)-   &
       &  v101(1)*rst(1)*rst(3)+ &
       &  v011(1)*(1-rst(1))*rst(3)+   &
       &  v110(1)*rst(1)*(1-rst(3))+   &
       &  v111(1)*rst(1)*rst(3)
  fjac(2,2)=-v000(2)*(1-rst(1))*(1-rst(3))-  &
       &  v100(2)*rst(1)*(1-rst(3))+   &
       &  v010(2)*(1-rst(1))*(1-rst(3))-  &
       &  v001(2)*(1-rst(1))*rst(3)-   &
       &  v101(2)*rst(1)*rst(3)+ &
       &  v011(2)*(1-rst(1))*rst(3)+   &
       &  v110(2)*rst(1)*(1-rst(3))+   &
       &  v111(2)*rst(1)*rst(3)
  fjac(3,2)=-v000(3)*(1-rst(1))*(1-rst(3))-  &
       &  v100(3)*rst(1)*(1-rst(3))+   &
       &  v010(3)*(1-rst(1))*(1-rst(3))-  &
       &  v001(3)*(1-rst(1))*rst(3)-   &
       &  v101(3)*rst(1)*rst(3)+ &
       &  v011(3)*(1-rst(1))*rst(3)+   &
       &  v110(3)*rst(1)*(1-rst(3))+   &
       &  v111(3)*rst(1)*rst(3)
  fjac(1,3)=-v000(1)*(1-rst(1))*(1-rst(2))-  &
       &  v100(1)*rst(1)*(1-rst(2))-   &
       &  v010(1)*(1-rst(1))*rst(2)+   &
       &  v001(1)*(1-rst(1))*(1-rst(2))+  &
       &  v101(1)*rst(1)*(1-rst(2))+   &
       &  v011(1)*(1-rst(1))*rst(2)-   &
       &  v110(1)*rst(1)*rst(2)+ &
       &  v111(1)*rst(1)*rst(2)
  fjac(2,3)=-v000(2)*(1-rst(1))*(1-rst(2))-  &
       &  v100(2)*rst(1)*(1-rst(2))-   &
       &  v010(2)*(1-rst(1))*rst(2)+   &
       &  v001(2)*(1-rst(1))*(1-rst(2))+  &
       &  v101(2)*rst(1)*(1-rst(2))+   &
       &  v011(2)*(1-rst(1))*rst(2)-   &
       &  v110(2)*rst(1)*rst(2)+ &
       &  v111(2)*rst(1)*rst(2)
  fjac(3,3)=-v000(3)*(1-rst(1))*(1-rst(2))-  &
       &  v100(3)*rst(1)*(1-rst(2))-   &
       &  v010(3)*(1-rst(1))*rst(2)+   &
       &  v001(3)*(1-rst(1))*(1-rst(2))+  &
       &  v101(3)*rst(1)*(1-rst(2))+   &
       &  v011(3)*(1-rst(1))*rst(2)-   &
       &  v110(3)*rst(1)*rst(2)+ &
       &  v111(3)*rst(1)*rst(2)
END SUBROUTINE fdjac1