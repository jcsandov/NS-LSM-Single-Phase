subroutine FreeSurfacePressureExtrapolation( pWaterNode                 , &
                                             PhiAirNode                 , &
                                             PhiWaterNode               , &
                                             VelocityGradientAirNode    , &
                                             VelocityGradientWaterNode  , &
                                             nvecAirNode                , &
                                             nvecWaterNode              , &
                                             pExtp                        &
                                           )


   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! INPUT/OUTPUT variables

   real ( kind = rdf ), intent(in) :: pWaterNode, PhiAirNode, PhiWaterNode
   
   real ( kind = rdf ), dimension(3,3), intent(in) :: VelocityGradientAirNode , &
                                                      VelocityGradientWaterNode
   
   real ( kind = rdf ), dimension(3)  , intent(in) :: nvecAirNode , &
                                                      nvecWaterNode
   
   
   real ( kind = rdf ), intent(out) :: pExtp
   
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local variables

   real ( kind = rdf ), dimension(3,3) :: VelocityGradient_fs 
   real ( kind = rdf ), dimension(3)   :: nvec_fs   
   real ( kind = rdf )                 :: dui_dxj_ni_nj , pfs
   real ( kind = rdf ) :: tol

   integer :: i,j,k

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   tol = 1.0E-10
   pExtp = zero

   ! If the free surface is close to the AIR   node --> PhiWaterNode >> PhiAirNode 
   ! If the free surface is close to the WATER node --> PhiAirNode   >> PhiWaterNode 
   ! That's why the weight are "exchanged"
   
   !print *, ' '
   !print *, 'PhiAirNode = ', PhiAirNode
   !print *, 'PhiWaterNode = ', PhiWaterNode
   !print *, 'VelocityGradientAirNode = '
   !print *, VelocityGradientAirNode(1,:)
   !print *, VelocityGradientAirNode(2,:)
   !print *, VelocityGradientAirNode(3,:)
   !print *, 'VelocityGradientWaterNode = '
   !print *, VelocityGradientWaterNode(1,:)
   !print *, VelocityGradientWaterNode(2,:)
   !print *, VelocityGradientWaterNode(3,:)
   !print *, ' '


   VelocityGradient_fs = (     abs(PhiAirNode  ) * VelocityGradientWaterNode     &
                           +   abs(PhiWaterNode) * VelocityGradientAirNode   )   &
                           / ( abs(PhiAirNode) + abs(PhiWaterNode) )
   
   nvec_fs = (     abs(PhiAirNode  ) * nvecWaterNode     &
               +   abs(PhiWaterNode) * nvecAirNode   )   &
               / ( abs(PhiAirNode) + abs(PhiWaterNode) )
   
   
   ! Computation of ∂ui/∂xj * ni * nj at the free surface
   
   dui_dxj_ni_nj = zero
   
   do j = 1,3
   do i = 1,3
               
      dui_dxj_ni_nj = dui_dxj_ni_nj + VelocityGradient_fs(i,j) * nvec_fs(i) * nvec_fs(j) 
   
   end do
   end do
   
   ! Normal Dynamic Boundary Condition p_fs = -2/Re * ( ∂ui / ∂xj * ni * nj )_fs
   !pfs = -two / ren * dui_dxj_ni_nj
   pfs = zero !-two / ren * dui_dxj_ni_nj
   
   if ( PhiWaterNode < tol ) then
      pExtp = pfs
      return
   end if

   pExtp = ( (   abs( PhiAirNode ) + abs( PhiWaterNode ) ) * pfs                     &
               - abs( PhiAirNode ) * pWaterNode ) / ( abs( PhiWaterNode ) )   
   
end subroutine FreeSurfacePressureExtrapolation