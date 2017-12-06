!!! modifications start at line 168 now
!!! all constants defined. all arguments for switching functions defined
!!! except CD_w. Everything else still remains to be done
!!! note that the original k-w formulation involved density being conserved along 
!!! with k and omega in consitituent equations. Here, I need to make sure that is the case.
!!! As such, this function is not multi-phase friendly yet
      subroutine elm3omega   (kay,	omega, 
      &     dwall,	gradV,
      &     srcRat1,src1,	srcJac)
      !
      !.... Data declaration
      !
      include "common.h"
      real*8	kay(npro),		omega(npro),
      &          dwall(npro)
      real*8	gradV(npro,3,3)
      
      real*8	rmu(npro),
      &          rho(npro)
      real*8	srcRat(npro,2),
      &		src(npro,2),		srcJac(npro,4),
      &          srcRat1(npro),          src1(npro)
      !
      integer advdiff
      integer	e
      real*8  fct1
      real*8  tmp1,   tmp2  ,  tmp0
      !
      real*8  k,    kInv,   kq,     omg,   y,    omgInv,   omgP2,
      &          ss,   mut, mut_k, mut_omg,  rat,      jac
      
      
      real*8  Ce1,   Ce2,  omgP2Inv,   CmuKE,  sigmaKE,       nu,
      &          kk,  kqInv,     nuInv,     Rey,     Ret,    Rey_k,
      &          Ret_k,        Ret_omg,     ff1,      f2,    kkInv,
      &          RetInv,       fmukeInv,  fmukeP2Inv
      
      real*8  fmuke,  fmuke_k,  fmuke_omg
      
      !...    Addings due to Durbins correction
      
      real*8  T1st, T2nd, T2ndInv, two3rdCmuInv, tri8thq,
      &          ssq, ssqInv, two3rdq, pk1, pk2,
      &          pk, pk_k, pk_omg, Tscale,Tinv,
      &          Tscale_k, Tscale_omg, ff1_k, ff1_omg, f2_k,f2_omg,
      &          ff1Inv, f2Inv
      
      
      
      !
      !	include "elmf.h"
      !
      !.... Compute src and its jacobians
      !
      !.... set Menter SST  k-omega Model constants
      rho(:)=datmat(1,1,1)
      rmu(:)=datmat(1,2,1)
      vinf(:)datmat(1,4,2)!freestream velocity
      chln(:)datmat(1,4,3)!characteristic velocity
      
      Cmu			= 0.09
      kappa			= 0.41
      a1			= 0.31
      CDES1			= 0.78
      CDES2			= 0.61
      Cd1			= 20
      Cd2			= 3
      
      alp1			= 5./9.
      beta1			= 0.075
      sigk1			= 0.85
      sigw1			= 0.5
      alp2			= 0.44
      beta2			= 0.0828
      sigk2			= 1
      sigw2			= 0.856
      tenpowerminustwenty	= 1.0d-20
      
      
      
      
      
      !Ce1             = 1.44
      !Ce2             = 1.92
      !CmuKE           = 0.09
      !sigmaKE         = 1.3
      
      !two3rdCmuInv    = 2./3./CmuKE
      !tri8thq         = SQRT(3./8.)
      !two3rdq         = SQRT(2./3.)
      !
      advdiff = 0
      if(advdiff.eq.0)then ! not advection-diffusion
      do e = 1, npro
      
         nuInv           = rho(e)/rmu(e)
      !   mutden          = 
      !   muInv           = rho(el)
         k		    = abs(kay(e))
         if (k.lt.1.e-32) k=0
         omg		    = abs(omega(e))
         y		    = dwall(e)
      !--------------patch
      !            if(k.gt.0.45) k=0.45
      !            if(omg.gt.2158) omg=2158
      !--------------------------------
      
      !
         kInv	    = 0
         kq	            = sqrt(k)
         kqInv           = 0
         kkInv           = 0
      
      !.... limiting k.  Instead of saying k.ne.0
      
         if ( k .gt.1.e-32 ) then
            kInv         = 1. / k
            kqInv        = 1./sqrt(k)
            kkInv        = kInv*kInv
         endif
      
         kk              =   k * k
      
         omgP2           = omg * omg
         omgInv	    = 0
         omgP2Inv        = 0
      
      
      !.....  limiting omega.  Instead of saying omg.ne.0
      
         if ( omg .gt.1.e-32) then
            omgInv        = 1. / omg
            omgP2Inv      = omgInv*omgInv
         endif
      
      !
         ss		= gradV(e,1,1) ** 2
      &			+ gradV(e,2,2) ** 2
      &			+ gradV(e,3,3) ** 2
      &			+ 0.5
      &			* ( (gradV(e,2,3) + gradV(e,3,2)) ** 2
      &			  + (gradV(e,3,1) + gradV(e,1,3)) ** 2
      &			  + (gradV(e,1,2) + gradV(e,2,1)) ** 2 )
      !
         ssq                 = sqrt(ss)
         ssqInv              = 0
      
         if(ss.ne.0) ssqInv  = 1./sqrt(ss)
      
         Rey                 = kq *      y * nuInv
         Ret                 = kk * omgInv * nuInv
         RetInv              = 0
      
      !...     limitng Ret so that it does not get 'nan' error
      
         if(Ret.lt.1.d100.AND.Ret.gt.zero) RetInv=1./Ret
      
         Rey_k    =  0.5 * y * kqInv  * nuInv
         Ret_k    =  2.  * k * omgInv * nuInv
         Ret_omg  = -kk  * omgP2Inv   * nuInv
      
      
      
      
      
      
      
      
      
      
      
      
      
      
         arg11den = Cmu*omg*y
         arg11    = kq / arg11den
         arg12num = 500.*nu
         arg12den = y * y * omg
         arg12    = arg12num / arg12den
         arg13num = 4.*rho(e)*sigw2*k
         argCD_kw = 2*(`rho(e)*sigw2/omg)*!!!!!!!!!!!!!!!!!!!!!!!gradk*gradw**********************************
         CD_kw = max(argCD_kw,tenpowerminustwenty)	
         arg13den = CD_kw*y*y
         arg13  = arg13num/arg13den 
         mutden	    = max(a1*omg,F2*ssq)
         mut		    = rho(e) * a1 *k / mutden 
       
         arg1 = min(max(arg11,arg12),arg13)
      
         F1   = tanh(arg1 ** 4)
      
         arg2 = max(2*arg11,arg12)
       
         F2   = tanh(arg2 **2)
      
         if (F1.gt.1) then  !limiting case for F1
            F1 = 1.0
         else if (F1.lt.1e-32) then
            F1 = 0.0
         endif 
      
         if (F2.gt.1) the   !limiting case for F2
            F2 = 1.0
         else if (F2.lt.1e-32) then
            F2 = 0.0
         endif 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of current modification from k-e to k-w
       
      !   tmp1     = exp(-0.0165*Rey)
      
      !   fmuKE    = (1. -tmp1) ** 2 * (1.+20.5*RetInv) ! fmu of Lam-Bremhorst
      
      !   fmuKEInv = 0.0
      !   fmuKEP2Inv = 0.0
      !....   limiting fmuKE.  fmuke max ~ 1, and it could get very small near the wall.
      
         if(fmuKe.gt.1e-32) then
            fmuKEInv = 1./fmuke
            fmuKEP2Inv = fmuKEInv*fmuKEInv
         endif
      
         fmuKE_k  = 0.033*(1.-tmp1)*(1.+20.5*RetInv)*Rey_k*tmp1
      &                -20.5 *(1.-tmp1)**2 * Ret_k*RetInv**2
      
         fmuKE_omg= -20.5*(1.-tmp1)**2* Ret_omg*RetInv**2
      
         ff1      =  1. + ( 0.05*fmuKEInv) ** 3 ! f1 as in Lam-Bremhorst
         f2       =  1. - exp(- Ret ** 2 ) ! f2 as in Lam-Bremhorst
      
         ff1Inv=zero
         f2Inv=zero
         if(ff1.gt.1.0e-32)    ff1Inv=1./ff1
         if(f2.gt.1.0e-32)     f2Inv =1./f2
      
      
         ff1_k    = -0.000375*fmuKE_k  *fmuKEInv**4
         ff1_omg  = -0.000375*fmuKE_omg*fmuKEInv**4
      
         f2_k     = 2.* Ret * Ret_k   * exp(-Ret**2)
         f2_omg   = 2.* Ret * Ret_omg * exp(-Ret**2)
      
         T1st = k * omgInv      ! 1st time scale
         T2nd = fmuKEInv*two3rdCmuInv*tri8thq*ssqInv ! 2nd time scale
      
      !...     Depending on the choice of T (to limit k growth near stagnation)
      
         if (T1st.lt.T2nd) then
      
            Tscale      = T1st
            Tinv        =  0
      
      !                 if(Tscale.ne.0)  Tinv=1./Tscale
      
      !...       Limiting time scale s.t. Tinv**4 not go over 1.e160
      
            if(Tscale.gt.1.e-32)  Tinv=1./Tscale
      !ccccTHIS IS WHAT WAS HERE WHEN I GOT IT FROM JE HOON
      !cccc                 Tscale_k    = fmuKE_k*k*omgInv+fmuKE*omgInv
      !cccc                 Tscale_omg  = fmuKE_omg*k*omgInv -fmuKE*k*omgP2Inv
            Tscale_k= omgInv ! Acusolve's choice
            Tscale_omg= -k*omgP2Inv ! AcuSolve's choice
      
         else
      
            Tscale      = T2nd
            Tinv        =  0
      
            if(Tscale.gt.1.e-32) Tinv=1./Tscale
      !ccccTHIS IS WHAT WAS HERE WHEN I GOT IT FROM JE HOON
      !cccc            Tscale_k   = 0
      !cccc            Tscale_omg =0
            Tscale_k= two3rdCmuInv*tri8thq*ssqInv*
      &                   (-fmuKEP2Inv*fmuKE_k)  ! acusolve's choice
            Tscale_omg= two3rdCmuInv*tri8thq*ssqInv*
      &                   (-fmuKEP2Inv*fmuKE_omg)  ! acusolve's choice
      
         endif
      
      !...      After Limiting all values, i feel its unnecessary to limit jacobians
      !...      since they are already limited.
      !...      Two other routines which defines these quantities are the same
      
      
      
      
      !...   recall that tmp1= exp(-0.0165*Rey)
      
      
      
         mut  = rho(e)*CmuKE*fmuKE* k *Tscale
      
         mut_k = rho(e)*CmuKE*(
      &                 fmuKE_k*k*Tscale
      &               + fmuKE    *Tscale
      &               + fmuKE  *k*Tscale_k
      &                            )
      
         mut_omg  = rho(e) * CmuKE *k*(
      &        fmuKE_omg*Tscale+
      &        fmuKE    *Tscale_omg
      &                            )
      
         tmp0        = 2 * ss
      
      
         pk1     = mut*tmp0
         pk2     = rho(e)*two3rdq*k*sqrt(ss)
      
      
         if (pk1.lt.pk2) then
            pk     = pk1
            pk_k   = mut_k  * tmp0
            pk_omg = mut_omg * tmp0
         else
            pk   = pk2
            pk_k = kInv*pk2
            pk_omg =0
         endif
      
      
         src(e,1)    = pk    - rho(e) * omg
         src(e,2)    = ff1 * Ce1 * pk * Tinv
      &                   -rho(e)* Ce2 * f2 * omg * Tinv
      
         srcRat(e,1) = -kInv*(  pk    - rho(e) * omg   )
         srcRat(e,2) = -omgInv*(
      &                    ff1 * Ce1 * pk * Tinv
      &                   -rho(e)* Ce2 * f2 * omg * Tinv
      &                         )
      
         srcJac(e,1) = -(pk_k)  ! jacobian for k PDE alone
         srcJac(e,3) = -(pk_omg-rho(e)) ! d(Fsrck)/d(omega)
         srcJac(e,2) = -(
      &            ff1_k * Ce1 * pk    * Tinv
      &           +ff1   * Ce1 * pk_k  * Tinv
      &           -ff1   * Ce1 * pk    * Tscale_k * Tinv**2
      &           -rho(e)  * Ce2 * f2_k  * omg * Tinv
      &           +rho(e)  * Ce2 * f2    * omg * Tscale_k * Tinv**2 ! do not touch
      &                      ) ! d(Fsrce)/d(k)
         srcJac(e,4) = -(
      &           ff1_omg * Ce1 * pk   * Tinv
      &           +ff1     * Ce1 * pk_omg * Tinv
      &        -ff1     * Ce1 * pk   * Tscale_omg * Tinv**2
      &           -rho(e)    * Ce2 * f2_omg     * omg * Tinv
      &           -rho(e)    * Ce2 * f2       * Tinv
      &        +rho(e)    * Ce2 * f2* omg * Tscale_omg * Tinv**2
      &              ) ! jacobian for omega PDE alone
      
      
      
      
      enddo
      !
      !.... Ensure positivity of srcJac
      !
      do e = 1, npro
      
      
         srcJac(e,1)		= max( srcJac(e,1), srcRat(e,1), 0.d0 )
         srcJac(e,4)		= max( srcJac(e,4), srcRat(e,2), 0.d0 )
      
      !            write(777,*) e,   srcJac(e,1),   srcJac(e,4)
      
      !
         tmp1		= min( srcJac(e,1) * srcJac(e,4),
      &				       (srcJac(e,1)-srcRat(e,1)) *
      &				       (srcJac(e,4)-srcRat(e,2)) )
         tmp2		= srcJac(e,2) * srcJac(e,3)
      
      
         if ( tmp2 .gt. tmp1 ) then
            tmp2		= sqrt(tmp1/tmp2)
            srcJac(e,2)	= tmp2 * srcJac(e,2)
            srcJac(e,3)	= tmp2 * srcJac(e,3)
         endif
      !	    srcJac(e,2)	= 0
      !	    srcJac(e,3)	= 0
      enddo
      if(isclr.eq.1)then        ! kay
         srcrat1 = srcrat(:,1)
         src1 = src(:,1)
      else if (isclr.eq.2) then ! omega
         srcrat1 = srcrat(:,2)
         src1 = src(:,2)
      endif
      
      else ! Advection-diffusion
      ! Advection-diffusion case
      srcrat1 = zero
      src1 = zero
      srcjac = zero
      endif
      !
      !.... Compute viscosity
      !
      !	do e = 1, npro
      !	    viscTot(e,1) = rmu(e) + xmuT(e)
      !
      !	    viscTot(e,2) = rmu(e) + xmuT(e)
      !     &                   / (SigmaKE)
      !	enddo
      !
      !.... Compute PDE residual
      !
      !	do e = 1, nElems
      !	    pdeRes(e,1)	= dens
      !     1			* ( masFct      * tkeTd(e)
      !     2		      	  + velK(e,1)   * gradK(e,1)
      !     3		      	  + velK(e,2)   * gradK(e,2)
      !     4		      	  + velK(e,3)   * gradK(e,3) )
      !     5		      	- src(e,1)
      !     6		      	- diffK(e)
      !	    pdeRes(e,2)	= dens
      !     1			* ( masFct      * tomgTd(e)
      !     2		      	  + velomg(e,1) * gradomg(e,1)
      !     3		      	  + velomg(e,2) * gradomg(e,2)
      !     4		      	  + velomg(e,3) * gradomg(e,3) )
      !     5		      	- src(e,2)
      !     6		      	- diffomg(e)
      !	enddo
      !
      !.... end of fElm3KomgCoef()
      !
      return
      end
