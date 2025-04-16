subroutine WriteTotalWaterVolume( GlobalVolume )

   real ( kind = rdf ), intent(in) :: GlobalVolume
   logical :: FileExist 


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ! GLOBAL VOLUME WRITING FILE
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   ! Global mass file writing
            
   inquire(file = "GlobalVolume.txt", exist = FileExist)
            
   if (FileExist) then
      open( 12 , file     = "GlobalVolume.txt" , status = "old",   &
                 position = "append"           , action = "write"    )
   else
      
      open( 12 , file   = "GlobalVolume.txt", status = "new"      , & 
                 action = "write"                                    )
   end if
            
   write(12, *) GlobalVolume
   
   close(12)
      

end subroutine WriteTotalWaterVolume