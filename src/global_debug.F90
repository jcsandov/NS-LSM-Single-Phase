module global_debug

!modulo para debugear las variables de las rutinas
!

use precision
implicit none

character (len = 256) :: OUTPUTdir
real (kind = rdf), dimension(10000000) :: debug_values  !valores para debugear de libre uso
integer, dimension(100000000) :: d_values_int
contains
!-----------------------------------------------

subroutine outputD1_real(varD,filename)

!recibe el vector en 1d como cuando estan definidas en init

use global, only :le_idx_a, le_idx_b, le_idx, &
                  li_ia, li_ib,       &
                  li_ja, li_jb,       &
                  li_ka, li_kb    
use global_mpi, only : myid
implicit none


!rutina que recibe una variable real (definida en cada nodo de la malla) 
!y cada procesador la escribe sin incluir los ghost point
        real (kind = rdf), dimension(le_idx_a(1):le_idx_b(1)),intent(in) :: varD
        real (kind = rdf), dimension(li_ia(1):li_ib(1),li_ja(1):li_jb(1),li_ka(1):li_kb(1)) :: varD_strip
        character (len = 256) :: filenameOUT
        character (len = *) :: filename
        integer:: ideb,jdeb,kdeb,ldeb
        integer:: myunit, offset
        
        do kdeb = li_ka(1), li_kb(1)
        do jdeb = li_ja(1), li_jb(1)
        do ideb = li_ia(1), li_ib(1)
                ldeb = le_idx(ideb,jdeb,kdeb,1)
                varD_strip(ideb,jdeb,kdeb) = varD(ldeb)
        end do
        end do
        end do
        offset = 60
        myunit = myid + offset

	
        write(filenameOUT,fmt='(2a,a,a,i2.2)' ) &
                 trim(OUTPUTdir), 'debugfiles/', &
                 trim(filename),'.',myunit -offset

        open  (unit = myunit, file = trim(filenameOUT), form = 'unformatted')
                write(unit=myunit) varD_strip                
     
        close(myunit)



end subroutine outputD1_real
!-----------------------------------------------

subroutine outputD3_real(varD,filename)

!recibe el vector en 3d como cuando estan dentro de las rutinas

use global, only : le_ia, le_ib, &
                   le_ja, le_jb, &
                   le_ka, le_kb, &
                   le_idx_a, le_idx_b
                   
use global_mpi, only : myid
implicit none


!rutina que recibe una variable real (definida en cada nodo de la malla) 
!y cada procesador la escribe sin incluir los ghost point
        real (kind = rdf), dimension(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)),intent(in) :: varD
        real (kind = rdf), dimension(le_idx_a(1):le_idx_b(1)) :: varD1 !transformacion a 1 dimension
        character (len = *) :: filename
        integer:: ideb,jdeb,kdeb,ldeb

        !paso el vector a 1d
        ldeb = 0
        do kdeb = le_ka(1),le_kb(1)
        do jdeb = le_ja(1),le_jb(1)
        do ideb = le_ia(1),le_ib(1)
                ldeb = ldeb + 1
                varD1(ldeb) = varD(ideb,jdeb,kdeb) 
        end do
        end do
        end do
        
        !aplico rutina 1d

        call outputD1_real(varD1,filename)



end subroutine outputD3_real
!-----------------------------------------------

subroutine gfunctionint(x_val,y_val,z_val,valfun)

implicit none
real (kind = rdf), intent(in) :: x_val,y_val,z_val
real (kind = rdf), intent(out):: valfun

valfun = x_val*y_val*z_val*sin(x_val*y_val*z_val)
!valfun = (x_val+1.0)**3*(y_val+1.0)**3*(z_val+1.0_rdf)**3


end subroutine gfunctionint


subroutine modifyq_daf(il,iu,jl,ju,kl,ku,igp,jgp,kgp,x,y,z,q)
!modifico q (pensando dentro de rutina solver daf)

Implicit none

integer, intent(in) :: igp,jgp,kgp
integer :: i_mysta_d, j_mysta_d, k_mysta_d, i_myend_d, j_myend_d, k_myend_d
integer :: i_d,j_d,k_d
integer, intent(in) :: il,jl,kl,iu,ju,ku
real (kind = rdf),dimension(1:4,il:iu,jl:ju,kl:ku) :: q
real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: x,y,z

i_mysta_d = il + igp
j_mysta_d = jl + jgp
k_mysta_d = kl + kgp

i_myend_d = iu - igp
j_myend_d = ju - jgp
k_myend_d = ku - kgp


do k_d = k_mysta_d, k_myend_d	
	do j_d = j_mysta_d, j_myend_d
		do i_d = i_mysta_d, i_myend_d
			q(1,i_d,j_d,k_d) = x(i_d,j_d,k_d)**2 + y(i_d,j_d,k_d)**2 + z(i_d,j_d,k_d)**2
		end do
	end do	
end do


end subroutine modifyq_daf

subroutine modifyvar(var)

!Por su forma puede modificar 3D y 1D. 
!Cuidado que solo modifica parte interior no gp hay que aplicarselo a la variable justo despues
use global, only : le_ia, le_ib, &
                   le_ja, le_jb, &
                   le_ka, le_kb, &
                   igp, jgp, kgp, &
                   x, y, z 
!modifico q pensado en forma 3D

Implicit none

integer :: i_mysta_d, j_mysta_d, k_mysta_d, i_myend_d, j_myend_d, k_myend_d
integer :: i_d,j_d,k_d
real (kind = rdf),dimension(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)), intent(inout) :: var
integer :: il,iu, jl,ju,kl,ku
real (kind = rdf), dimension(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)) :: xx,yy,zz

il = le_ia(1)
jl = le_ja(1)
kl = le_ka(1)
iu = le_ib(1)
ju = le_jb(1)
ku = le_kb(1)

call transform3D_to_1D(x,xx)
call transform3D_to_1D(y,yy)
call transform3D_to_1D(z,zz)

i_mysta_d = il + igp(1)
j_mysta_d = jl + jgp(1)
k_mysta_d = kl + kgp(1)

i_myend_d = iu - igp(1)
j_myend_d = ju - jgp(1)
k_myend_d = ku - kgp(1)


do k_d = k_mysta_d, k_myend_d	
	do j_d = j_mysta_d, j_myend_d
		do i_d = i_mysta_d, i_myend_d
			!var(i_d,j_d,k_d) = xx(i_d,j_d,k_d) + yy(i_d,j_d,k_d) + zz(i_d,j_d,k_d)
			var(i_d,j_d,k_d) = xx(i_d,j_d,k_d)**0.1 + yy(i_d,j_d,k_d)**0.1 + zz(i_d,j_d,k_d)**0.1      
		end do
	end do	
end do


end subroutine modifyvar


subroutine transform3D_to_1D(var1D,var3D)

use global, only : le_ia, le_ib, &
                   le_ja, le_jb, &
                   le_ka, le_kb
Implicit none
real(kind = rdf), dimension (le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)), intent(in) :: var1D
real(kind = rdf), dimension (le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)), intent(out) :: var3D

var3D = var1D

end subroutine transform3D_to_1D

subroutine dirRead_debug()

  open(unit = 42, file = "directory", form = "formatted")
	read(unit = 42, fmt = "(a)") OUTPUTdir
  close(unit = 42)

end subroutine dirRead_debug

subroutine index_outputD3(idx_i, idx_j, idx_k, varD, varD_value, whichproc, printscreen)

  !Rutina que plotea el valor de cualquier variable en tres dimensiones con el indice pedido
  !sin tener que preocuparse cual procesador es el que realmente tiene la variable
  
  use global, only :le_ia, le_ib, &
                    le_ja, le_jb, &
                    le_ka, le_kb, &
                    gi_ia, gi_ib, &
                    gi_ja, gi_jb, &
                    gi_ka, gi_kb
                    
  use global_mpi, only : myid
                    
  implicit none
  integer :: idx_i, idx_j, idx_k
  integer ::i_local, j_local, k_local
  real (kind = rdf), dimension(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1)),intent(in) :: varD
  
  real (kind = rdf),optional :: varD_value  
  integer, optional :: whichproc
  logical, optional :: printscreen
  !print *, gi_ia(1), gi_ib(1), myid
  !print *, gi_ja(1), gi_jb(1), myid
  !print *, gi_ka(1), gi_kb(1), myid
  
  !Encontrar nodo en dominio
  
  if( (gi_ia(1) .le. idx_i) .AND. (idx_i .le. gi_ib(1)) .AND. &
      (gi_ja(1) .le. idx_j) .AND. (idx_j .le. gi_jb(1)) .AND. &
      (gi_ka(1) .le. idx_k) .AND. (idx_k .le. gi_kb(1)) ) then
      
      i_local = idx_i - gi_ia(1) + 1
      j_local = idx_j - gi_ja(1) + 1
      k_local = idx_k - gi_ka(1) + 1
      
      if(present(varD_value)) varD_value = varD(i_local, j_local, k_local)
      if(present(printscreen)) then
        if(printscreen) print *, varD(i_local, j_local, k_local), myid
      end if
      if(present(whichproc)) whichproc = myid
  end if    
  
end subroutine index_outputD3

  
end module global_debug
