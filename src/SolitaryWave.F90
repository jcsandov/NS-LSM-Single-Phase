module SolitaryWave
   
   use precision
   
   implicit none

   ! Solitary wave parameters

   real ( kind = rdf ), parameter :: hdepth    = one                                   , & 
                                     Hwave     = two/ten                               , &
                                     Steepness = Hwave / hdepth                        , &
                                     Froude    = 0.128_rdf                             , &
                                     Celerity  = one / Froude * sqrt( hdepth + Hwave ) 

   contains

   function EtaWave(x,t)

      implicit none
      real ( kind = rdf ) :: x,t
      real ( kind = rdf ) :: EtaWave
      real ( kind = rdf ) :: aux 

      !                   _                          _
      !                  |      ( 3 H )    (x - Ct)   |
      ! η(x,t) = H sech2 | sqrt (-----) *  --------   |
      !                  |_     ( 4 h )       h      _|
      !


      aux = ( x - Celerity * t ) / hdepth * sqrt( ( three * Hwave ) / ( four * hdepth ) ) 

      EtaWave = Hwave * ( one / cosh( aux ) )**2

   end function EtaWave

   function uwave(x,z,t)

      real ( kind = rdf ) :: x,z,t
      real ( kind = rdf ) :: uwave
      real ( kind = rdf ) :: e        ! ϵ : Wave Steepness    
      real ( kind = rdf ) :: eta2H    ! η/H 
      real ( kind = rdf ) :: z2h      ! z/h

      e     = Steepness
      eta2H = EtaWave(x,t) / Hwave
      z2h   = z / hdepth      

      uwave =               eta2H * ( e + 3 * e**2 * ( one/six - one_half * z2h**2 ) ) &
               - ( eta2H * e )**2 * ( seven/four - nine/four * z2h**2                )

      uwave = uwave * sqrt( hdepth / Froude**2 )

   end function uwave


   function wwave(x,z,t)

      real ( kind = rdf ) :: x,z,t
      real ( kind = rdf ) :: wwave
      real ( kind = rdf ) :: e        ! ϵ = H/h : Wave Steepness    
      real ( kind = rdf ) :: eta2H    ! η/H 
      real ( kind = rdf ) :: eta2hd   ! η/h 
      real ( kind = rdf ) :: X2h      ! X/h 
      real ( kind = rdf ) :: z2h      ! z/h

      e      = Steepness
      eta2H  = EtaWave(x,t) / Hwave
      eta2hd = EtaWave(x,t) / hdepth      
      X2h    = ( x - Celerity * t ) / hdepth
      z2h    = z / hdepth      

      wwave = sqrt( three * e ) * z2h * eta2hd * tanh( sqrt( three * e / four ) * X2h ) * &
              ( one + one_half * e * ( one - seven * eta2H - (one - three * eta2H) * z2h**2 ) )   

      wwave = wwave * sqrt( hdepth / Froude**2 )

   end function wwave

   subroutine SolitaryWaveFlowField(x,z,t,p,u,w) 

      real ( kind = rdf ) :: x , z , t
      real ( kind = rdf ) :: ufs, wfs, umod , ufsmod
      real ( kind = rdf ) , intent(out) :: p , u , w 

      u = uwave( x , z , t )
      w = wwave( x , z , t )

      umod = sqrt(u**2 + w**2)

      ufs = uwave( x , hdepth + EtaWave(x,t) , t )
      wfs = wwave( x , hdepth + EtaWave(x,t) , t )

      ufsmod = sqrt(ufs**2 + wfs**2)

      ! Pressure is updated using non-dimensional Bernoulli's equation

      ! z + Fr^2 * p + Fr^2 * |u|^2 / 2 = constant

      !
      !              h + η - z      | u( x , h+η , t ) |^2 - | u( x , h+η , t ) |^2     
      ! p(x,z,t) =  ------------ +  -------------------------------------------------  
      !                Fr^2                                2
      !

      p = ( hdepth + EtaWave(x,t) - z ) / Froude**2 + &
          one_half * ( ufsmod**2 - umod**2 )

   end subroutine SolitaryWaveFlowField

end module SolitaryWave
