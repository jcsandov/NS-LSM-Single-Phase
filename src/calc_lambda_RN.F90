subroutine calc_lambda_RN(sgndf_zero,normagrad_zero,sgndf_RN,lambda_RN,timeStep)

!calcula lambda con una aproximación de primer orden.
!basado en sussmann fatemi

use global_LSM, only : heaviside, deltaf, sussman_correction_method

implicit none

real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: sgndf_zero,normagrad_zero
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (in):: sgndf_RN
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) :: num, denom, numerator, denominator
real (kind = rdf) , intent (in):: timeStep
real (kind = rdf) ,dimension (il:iu,jl:ju,kl:ku) , intent (out):: lambda_RN

!local

real (kind = rdf) :: volcell
real (kind = rdf) :: Loperator
real (kind = rdf) :: volnumerator,voldenominator
real (kind = rdf) :: aux, aux2


!index
integer :: i_sta,i_end
integer :: j_sta,j_end
integer :: k_sta,k_end

!-------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


!ahora calculo lambda con la aproximación del volumen como en el centro
i_sta = il + igp
j_sta = jl + jgp
k_sta = kl + kgp

i_end = iu - igp        
j_end = ju - jgp
k_end = ku - kgp

if (sussman_correction_method == 0) then
        do i = i_sta, i_end
        do j = j_sta, j_end
        do k = k_sta, k_end
        
                volcell = dc**(-1) * de**(-1) * dz**(-1)
                Loperator = deltaf(sgndf_zero(i,j,k),epslsm) *   (timeStep**(-1) * (sgndf_RN(i,j,k) - sgndf_zero(i,j,k)))
                volnumerator = volcell * Loperator
                voldenominator = volcell * deltaf(sgndf_zero(i,j,k),epslsm)**2 * normagrad_zero(i,j,k)
                
                if(voldenominator == zero) then
                        lambda_RN(i,j,k) = zero
                else
                
                        lambda_RN(i,j,k) = -volnumerator/voldenominator
                end if                        
                !if(myid == 0 .and. i==32 .and. j==31 .and. k==42) print *, '--------------------------'
                !if(myid == 0 .and. i==32 .and. j==31 .and. k==42) print *, volnumerator, voldenominator
                !if(myid == 0 .and. i==32 .and. j==31 .and. k==42) print *, deltaf(sgndf_zero(i,j,k),epslsm), normagrad_zero(i,j,k)
        end do
        end do
        end do

else if (sussman_correction_method == 1) then

! TO DO: calcular los integrandos para toda la malla del proc y luego usarla para el filtro es sumamente caro. Si esto funciona, hay 
! que camibarlos por unos for anidados, como lo del modelo de Kang, para aprovechar la localidad espacial de la lectura de memoria.
                
        num = zero
        denom = zero
        numerator = zero
        denominator = zero

        do i = il, iu
        do j = jl, ju
        do k = kl, ku
                num(i,j,k) = deltaf(sgndf_zero(i,j,k),epslsm) *   (timeStep**(-1.0_rdf) * (sgndf_RN(i,j,k) - sgndf_zero(i,j,k)))
                denom(i,j,k) = deltaf(sgndf_zero(i,j,k),epslsm)**two * normagrad_zero(i,j,k)
        end do
        end do
        end do

        call testfilter(num, numerator)
        call testfilter(denom, denominator)

        do i = i_sta, i_end
        do j = j_sta, j_end
        do k = k_sta, k_end


                if (abs(denominator(i,j,k)).le.1E-10) then
                        lambda_RN(i,j,k) = zero
                else
                        lambda_RN(i,j,k) = -numerator(i,j,k)/denominator(i,j,k)

                        !if (myid == root) then
                                !print *, '++++++++++++++++++++++++++++++++++++++++++++++++'
                                !print *, 'normagrad_zero = ', normagrad_zero(i,j,k)
                                !print *, 'deltaf(sgndf_zero(i,j,k),epslsm)**2.0 = ', deltaf(sgndf_zero(i,j,k),epslsm)**2.0
                                !print *, '-----------------------------------------'
                                !print *, 'num = ', num(i,j,k)
                                !print *, 'denom = ', denom(i,j,k)
                                !print *, '-----------------------------------------'
                                !print *, 'numerator = ', numerator(i,j,k)
                                !print *, 'denominator = ', denominator(i,j,k)                        
                                !print *, ' '
                                !print *, '-----------------------------------------'
                                !print *, 'lambda simpson = ', lambda_RN(i,j,k)
                                !print *, 'lambda point constant = ', aux 
                                !print *, '|lambda_simpson - lambda_point| = ', aux2
                                !print *, ' '
                        !end if
                end if
        

        end do
        end do
        end do
end if

end subroutine calc_lambda_RN
