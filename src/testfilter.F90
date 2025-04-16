subroutine testfilter(var,varf)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Test filter of a var(i,j,k) 
  ! using the trapezoidal (1) or simpson (2) filters
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  

!input: var
!output: varf

implicit none
 
! allocatable variables
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent (in)  :: var
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent (out) :: varf

! local array
real (kind = rdf), dimension(27) :: weight, varfilter
real (kind = rdf) :: avgvar, avgvol
integer ::i,j,k,m, filtertype
integer :: i_mysta,i_myend,j_mysta,j_myend,k_mysta,k_myend

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!var(:,:,:) = var1(:,:,:)

i_mysta = il + igp 
j_mysta = jl + jgp
k_mysta = kl + kgp

i_myend = iu - igp        
j_myend = ju - jgp
k_myend = ku - kgp


! interior nodes only; processes on the domain boundaries
! 
if (myback == mpi_proc_null)  i_mysta = il + igp + 1
if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

if (myfront == mpi_proc_null) i_myend = iu - igp - 1
if (myright == mpi_proc_null) j_myend = ju - jgp - 1
if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

filtertype = 1 ! trapezoidal
! filtertype = 2 ! simpson

varf = var
  
! simpson test filter
do k=k_mysta,k_myend
do j=j_mysta,j_myend
do i=i_mysta,i_myend
        varfilter(1) = var(i,j,k)
        varfilter(2) = var(i+1,j,k)
        varfilter(3) = var(i-1,j,k)
        varfilter(4) = var(i,j+1,k)
        varfilter(5) = var(i,j-1,k)
        varfilter(6) = var(i,j,k+1)
        varfilter(7) = var(i,j,k-1)
        varfilter(8) = var(i+1,j+1,k)
        varfilter(9) = var(i+1,j-1,k)
        varfilter(10)= var(i-1,j+1,k)
        varfilter(11)= var(i-1,j-1,k)
        varfilter(12)= var(i+1,j,k+1)
        varfilter(13)= var(i+1,j,k-1)
        varfilter(14)= var(i-1,j,k+1)
        varfilter(15)= var(i-1,j,k-1)
        varfilter(16)= var(i,j+1,k+1)
        varfilter(17)= var(i,j+1,k-1)
        varfilter(18)= var(i,j-1,k+1)
        varfilter(19)= var(i,j-1,k-1)
        varfilter(20)= var(i+1,j+1,k+1)
        varfilter(21)= var(i+1,j+1,k-1)
        varfilter(22)= var(i+1,j-1,k+1)
        varfilter(23)= var(i+1,j-1,k-1)
        varfilter(24)= var(i-1,j+1,k+1)
        varfilter(25)= var(i-1,j+1,k-1)
        varfilter(26)= var(i-1,j-1,k+1)
        varfilter(27)= var(i-1,j-1,k-1)

        weight(1)  = sixtyfour / aj(i,j,k)
        weight(2)  = sixteen / aj(i+1,j,k)
        weight(3)  = sixteen / aj(i-1,j,k)
        weight(4)  = sixteen / aj(i,j+1,k)
        weight(5)  = sixteen / aj(i,j-1,k)
        weight(6)  = sixteen / aj(i,j,k+1)
        weight(7)  = sixteen / aj(i,j,k-1)
        weight(8)  = four / aj(i+1,j+1,k)
        weight(9)  = four / aj(i+1,j-1,k)
        weight(10) = four / aj(i-1,j+1,k)
        weight(11) = four / aj(i-1,j-1,k)
        weight(12) = four / aj(i+1,j,k+1)
        weight(13) = four / aj(i+1,j,k-1)
        weight(14) = four / aj(i-1,j,k+1)
        weight(15) = four / aj(i-1,j,k-1)
        weight(16) = four / aj(i,j+1,k+1)
        weight(17) = four / aj(i,j+1,k-1)
        weight(18) = four / aj(i,j-1,k+1)
        weight(19) = four / aj(i,j-1,k-1)
        weight(20) = one / aj(i+1,j+1,k+1)
        weight(21) = one / aj(i+1,j+1,k-1)
        weight(22) = one / aj(i+1,j-1,k+1)
        weight(23) = one / aj(i+1,j-1,k-1)
        weight(24) = one / aj(i-1,j+1,k+1)
        weight(25) = one / aj(i-1,j+1,k-1)
        weight(26) = one / aj(i-1,j-1,k+1)
        weight(27) = one / aj(i-1,j-1,k-1)
        
        avgvar = zero
        avgvol = zero
        do m=1,27
                avgvar = avgvar + weight(m)*varfilter(m)
                avgvol = avgvol + weight(m)
        end do
        varf(i,j,k) = avgvar / avgvol
end do  
end do  
end do  
 
end subroutine testfilter
