
c-----------------------------------------------------------------------
   Menter SST(original 1994) Turbulence model constants
c
c-----------------------------------------------------------------------
      module turbKW
      real*8 mu, kappa, a1, CDES1, CDES2,  Cd1, Cd2
      real*8 alp1, alp2, beta1, beta2, sigk1, sigk2, sigw1, sigw2, tenpowerminustwenty
      parameter (
     &  mu                     = 0.09
     &  kappa                   = 0.41
     &  a1                      = 0.31
     &  CDES1                   = 0.78
     &  CDES2                   = 0.61
     &  Cd1                     = 20
     &  Cd2                     = 3

     &  alp1                    = 5./9.
     &  beta1                   = 0.075
     & 	 sigk1                   = 0.85
     & 	 sigw1                   = 0.5
     & 	 alp2                    = 0.44
     & 	 beta2                   = 0.0828
     & 	 sigk2                   = 1
     &	  sigw2                   = 0.856
     & 	 tenpowerminustwenty     = 1.0d-20
     &     )

      end module
