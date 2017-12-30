
c-----------------------------------------------------------------------
!   Menter SST(revised 2003) Turbulence model constants
c
c-----------------------------------------------------------------------
      module turbKW
      real*8 mu, kappa, a1, CDES1, CDES2,  Cd1, Cd2, tenpowerminusten
      real*8 alp1, alp2, beta1, beta2, sigk1, sigk2, sigw1, sigw2, tenpowerminustwenty
      parameter (
     &  CmuKW                     = 0.09 !!This variable is \beta^* in the NASA document
     &  kappa                   = 0.41
     &  a1                      = 0.31
     &  CDES1                   = 0.78
     &  CDES2                   = 0.61
     &  Cd1                     = 20
     &  Cd2                     = 3

     &  alp1                    = 5.0/9.0
     &  beta1                   = 0.075
     & 	 sigk1                   = 0.85
     & 	 sigw1                   = 0.5
     & 	 alp2                    = 0.44
     & 	 beta2                   = 0.0828
     & 	 sigk2                   = 1
     &	 sigw2                   = 0.856
     & 	 tenpowerminustwenty     = 1.0d-20
     &   gam1                    = 0.55316666d0
     &   gam2                    = 0.44035466d0
     &   tenpowerminusten     = 1.0d-10
     &     )

      end module
