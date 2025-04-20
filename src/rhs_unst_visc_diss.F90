
subroutine rhs_unst_visc_diss
  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Calculate second-order accurate unsteady terms and add
   ! current viscous and dissipation terms to rh

   ! Unlike previous codes we will add viscous, dissipation, and
   ! unsteady terms in the calling routine not in this one. 

   ! input
   !     aj(ijk)
   !     q(4,ijk)
   !     qn(4,ijk)
   !     qnm1(4,ijk)
   !     delti
   !     rh(4,ijk)

   ! output
   !     rh(4,ijk)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! local dummy reals
   !

   real (kind = rdf) :: dtaj
   real (kind = rdf) :: dq2dt, dq3dt, dq4dt

   !real (kind = rdf), dimension( il:iu , jl:ju , kl:ku ) :: dq2dt, dq3dt, dq4dt
   character (len = 256) :: debugname

   dq2dt = zero
   dq3dt = zero
   dq4dt = zero

   ! Add time derivative, viscous and dissipation terms to rh

   !    do k = 2, km - 1
   !    do j = 2, jm - 1
   !    do i = 2, im - 1


   !aca agregare los terminos de fuente S del Level set method (gravity forces.,surface forces y tension superficial)

!   do k = k_mysta , k_myend
!   do j = j_mysta , j_myend
!   do i = i_mysta , i_myend
!
!      gforce(4,i,j,k) =  rsign(i,j,k) * ( one / FrLSM**2 ) /  aj(i,j,k) 
!                              
!   end do
!   end do
!   end do
        
                !call outputD3_real(gforce(4,:,:,:),'gforce') 



!
!   print *, 'At rhs_unst_visc_diss :'
!   print *, 'max diss1 = ', maxval( diss(1,:,:,:) )  
!   print *, 'min diss1 = ', minval( diss(1,:,:,:) )  
!   print *, ' '
!
!   if ( abs(maxval( diss(1,:,:,:) )) > ten .or. &
!        abs(minval( diss(1,:,:,:) )) > ten ) then
!
!      print *, 'DISS ERROR AT rhs_unst_visc_diss'
!      stop
!
!   end if
!
!
!   print *, 'At rhs_unst_visc_diss :'
!   print *, 'max visc1 = ', maxval( visc(1,:,:,:) )  
!   print *, 'min visc1 = ', minval( visc(1,:,:,:) )  
!   print *, 'max visc2 = ', maxval( visc(2,:,:,:) )  
!   print *, 'min visc2 = ', minval( visc(2,:,:,:) )  
!   print *, 'max visc3 = ', maxval( visc(3,:,:,:) )  
!   print *, 'min visc3 = ', minval( visc(3,:,:,:) )  
!   print *, ' '

   !if ( abs(maxval( visc(1,:,:,:) )) > ten .or. &
   !     abs(minval( visc(1,:,:,:) )) > ten .or. &  
   !     abs(maxval( visc(2,:,:,:) )) > ten .or. &  
   !     abs(minval( visc(2,:,:,:) )) > ten .or. &  
   !     abs(maxval( visc(3,:,:,:) )) > ten .or. &  
   !     abs(minval( visc(3,:,:,:) )) > ten  ) then
   !
   !   print *, ' '
   !   print *, 'VISC ERROR AT rhs_unst_visc_diss'
   !   print *, ' '
   !   stop
   !
   !end if

   do k = k_mysta , k_myend
   do j = j_mysta , j_myend
   do i = i_mysta , i_myend

      dtaj = two * delti * aj(i,j,k)

      !dq2dt = e_source*(three*q(2,i,j,k)-four*qn(2,i,j,k)+qnm1(2,i,j,k)) / dtaj *rhoLSM(phi_n(i,j,k))
      !dq3dt = e_source*(three*q(3,i,j,k)-four*qn(3,i,j,k)+qnm1(3,i,j,k)) / dtaj *rhoLSM(phi_n(i,j,k))
      !dq4dt = e_source*(three*q(4,i,j,k)-four*qn(4,i,j,k)+qnm1(4,i,j,k)) / dtaj *rhoLSM(phi_n(i,j,k))

      dq2dt = e_source * ( three * q(2,i,j,k) - four * qn(2,i,j,k) + qnm1(2,i,j,k) ) / dtaj 
      dq3dt = e_source * ( three * q(3,i,j,k) - four * qn(3,i,j,k) + qnm1(3,i,j,k) ) / dtaj 
      dq4dt = e_source * ( three * q(4,i,j,k) - four * qn(4,i,j,k) + qnm1(4,i,j,k) ) / dtaj 

      rh(2,i,j,k) = rh(2,i,j,k) + dq2dt
      rh(3,i,j,k) = rh(3,i,j,k) + dq3dt
      rh(4,i,j,k) = rh(4,i,j,k) + dq4dt

   end do
   end do
   end do

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

   !write(debugname, fmt ='(a,i6.6)') 'dq2dt', iteraciontiempo
   !call outputD3_real( dq2dt, debugname ) 

   !write(debugname, fmt ='(a,i6.6)') 'dq3dt', iteraciontiempo
   !call outputD3_real( dq3dt, debugname ) 

   !write(debugname, fmt ='(a,i6.6)') 'dq4dt', iteraciontiempo
   !call outputD3_real( dq4dt, debugname ) 

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


end subroutine rhs_unst_visc_diss









