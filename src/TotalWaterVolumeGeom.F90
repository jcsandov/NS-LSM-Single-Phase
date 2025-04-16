subroutine TotalWaterVolumeGeom(phizero, nnodes, ntetrahedra, NodesList, TetrahedraList, InterfaceNodesID)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! This subroutine computes the total volume of fluid given a spatial phi distribution
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   use DataTypes
   use precision
   use TetrahedronMethods


   implicit none

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Input/Output arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   real (kind = rdf), target, dimension(il:iu,jl:ju,kl:ku), intent(in) :: phizero
   integer, intent(in) :: nnodes, ntetrahedra
   type(pnode), dimension(1:nnodes), intent(in) :: NodesList
   type(tetrahedron), dimension(1:ntetrahedra), intent(in) :: TetrahedraList
   integer, dimension(il:iu,jl:ju,kl:ku), intent(in) :: InterfaceNodesID 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! Local Variables
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   logical :: FileExist
   real(kind = rdf) :: TotalWaterVolumeVar, VolumeAux   
   real(kind = rdf), dimension(3) :: vertA, vertB, vertC, vertD, vertE, vertF, vertG, vertH
   real(kind = rdf), dimension(3) :: AC, AF, AH, BE, CH, DB, ED, FC, GA, GB, GC, GD, GE, GF, GH, HF

   ! Auxiliary counters for distance computation
   integer :: ivol, jvol, kvol
   integer :: K

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend


   ! loop boundaries definition
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp
            
   ! processes on the domain boundaries
   
   if (myback  == mpi_proc_null)  i_mysta = il + igp + 1
   if (myleft  == mpi_proc_null)  j_mysta = jl + jgp + 1
   if (mydown  == mpi_proc_null)  k_mysta = kl + kgp + 1
   
   if (myfront == mpi_proc_null)  i_myend = iu - igp - 1
   if (myright == mpi_proc_null)  j_myend = ju - jgp - 1
   if (myup    == mpi_proc_null)  k_myend = ku - kgp - 1



! I look over the whole domain. If ϕ>0 and I'm on a pnode, then the volume is the one obtained by TetrahedronMethods
! if ϕ(i,j,k)>0 and all my neighbours ϕ(i±1, j±1, k±1)>0, then the volume is the cell volume

TotalWaterVolumeVar = zero

do ivol = i_mysta-1, i_myend
do jvol = j_mysta-1, j_myend
do kvol = k_mysta-1, k_myend

   ! we check if I'm in a one-phase cell

   if ( phizero( ivol   , jvol   , kvol   ) > 0 .and. &
        phizero( ivol+1 , jvol   , kvol   ) > 0 .and. &
        phizero( ivol+1 , jvol+1 , kvol   ) > 0 .and. &
        phizero( ivol   , jvol+1 , kvol   ) > 0 .and. &
        phizero( ivol   , jvol   , kvol+1 ) > 0 .and. &
        phizero( ivol+1 , jvol   , kvol+1 ) > 0 .and. &
        phizero( ivol+1 , jvol+1 , kvol+1 ) > 0 .and. &
        phizero( ivol   , jvol+1 , kvol+1 ) > 0          ) then

      !   NODE DISTRIBUTION OF THE COMPUTATIONAL CELL
      !
      !    i,j+1,k+1 ----i+1,j+1,k+1    
      !     /|(H)           /|(G)                            
      !    / |             / |                                  
      ! i,j,k+1-------i+1,j,k+1                   
      !   |(E)           |(F)|                                   
      !   |  |           |   |                                
      !   |  |           |   |                                
      !   |  i,j+1,k-----|-i+1,j+1,k           
      !   | /(D)         |  /(C)                 
      !   |/             | /
      ! i,j,k-----------i+1,j,k
      !  (A)               (B)

      ! 8 vertices coordinates
      vertA(:) = (/ x( ivol   , jvol   , kvol   ), y( ivol   , jvol   , kvol   ), z( ivol   , jvol   , kvol   ) /)
      vertB(:) = (/ x( ivol+1 , jvol   , kvol   ), y( ivol+1 , jvol   , kvol   ), z( ivol+1 , jvol   , kvol   ) /)
      vertC(:) = (/ x( ivol+1 , jvol+1 , kvol   ), y( ivol+1 , jvol+1 , kvol   ), z( ivol+1 , jvol+1 , kvol   ) /)
      vertD(:) = (/ x( ivol   , jvol+1 , kvol   ), y( ivol   , jvol+1 , kvol   ), z( ivol   , jvol+1 , kvol   ) /)
      vertE(:) = (/ x( ivol   , jvol   , kvol+1 ), y( ivol   , jvol   , kvol+1 ), z( ivol   , jvol   , kvol+1 ) /)
      vertF(:) = (/ x( ivol+1 , jvol   , kvol+1 ), y( ivol+1 , jvol   , kvol+1 ), z( ivol+1 , jvol   , kvol+1 ) /)
      vertG(:) = (/ x( ivol+1 , jvol+1 , kvol+1 ), y( ivol+1 , jvol+1 , kvol+1 ), z( ivol+1 , jvol+1 , kvol+1 ) /)
      vertH(:) = (/ x( ivol   , jvol+1 , kvol+1 ), y( ivol   , jvol+1 , kvol+1 ), z( ivol   , jvol+1 , kvol+1 ) /)

      ! edges vectors
      AC = vertC - vertA
      AF = vertF - vertA
      AH = vertH - vertA

      BE = vertE - vertB

      CH = vertH - vertC

      DB = vertB - vertD

      ED = vertD - vertE

      FC = vertC - vertF

      GA = vertA - vertG
      GB = vertB - vertG
      GC = vertC - vertG
      GD = vertD - vertG
      GE = vertE - vertG
      GF = vertF - vertG
      GH = vertH - vertG

      HF = vertF - vertH

      ! volume formula for a general hexahedronn obtained from
      ! Davies, D. E., & Salmond, D. J. (1985). Calculation of the volume of a general  
      ! hexahedron for flow predictions. AIAA journal, 23(6), 954-956.

      VolumeAux =   -one / twelve * &
                     (    Dot_Product( (GA+GB) , ( CrossProduct(DB,AC) ) ) &
                       +  Dot_Product( (GA+GE) , ( CrossProduct(BE,AF) ) ) &
                       +  Dot_Product( (GA+GD) , ( CrossProduct(ED,AH) ) ) &
                       +  Dot_Product( GF      , ( CrossProduct(HF,GE) ) ) &
                       +  Dot_Product( GH      , ( CrossProduct(CH,GD) ) ) &
                       +  Dot_Product( GC      , ( CrossProduct(FC,GB) ) ) &
                     )

      TotalWaterVolumeVar = TotalWaterVolumeVar + VolumeAux 

      if (VolumeAux < eps_sims) then
        print *, 'VolumeAux = ', VolumeAux
      end if

   end if


end do
end do
end do

! All the rest of the fluid volume is contributed by the one contained within tetrahedra

do K = 1, ntetrahedra
  
  if ( TetrahedraList(K)%WaterVolume  < - eps_sims ) then
    print *, 'TetrahedraList(K)%WaterVolume = ', TetrahedraList(K)%WaterVolume
  end if

   TotalWaterVolumeVar = TotalWaterVolumeVar + TetrahedraList(K)%WaterVolume
end do

! Global mass file writing

inquire(file = "GlobalMassGeom.txt", exist = FileExist)

if (FileExist) then
  open( 12, file = "GlobalMassGeom.txt", status = "old", &
            position = "append", action = "write"       )
else
  open( 12, file = "GlobalMassGeom.txt", status = "new", & 
            action = "write")
end if

write(12, *) TotalWaterVolumeVar
close(12)


end subroutine TotalWaterVolumeGeom
