!!! all constants defined in TurbKW. This function is not multi-phase friendly yet


      subroutine elm3komega ( kay ,	omega ,
     &     dwall,	gradV,
     &     srcRat1 , src1 ,	srcJac)
      !
      !.... Data declaration
      !
      use turbKW
      include "common.h"
!!INPUTS
      real*8	kay(npro),		omega(npro),
     &          dwall(npro),  gradV(npro,nsd,nsd)
!!OUTPUTS
      real*8	srcRat1(npro) ,  src1(npro) ,	srcJac(npro,4)
!!below set of definitions are for e3qvarkwSclr subroutine
      real*8   yl(npro,nshl,ndof),    shp(npro,nshl),
     &         shgl(npro,nsd,nshl),   xl(npro,nenl,nsd),
     &         dxidx(npro,nsd,nsd),   WdetJ(npro),
     &         gradK(npro,nsd),    gradW(npro,nsd) !!check v
!!LOCALS
      real*8	rmu(npro),
     &          rho(npro)  !!check v

      real*8	srcRat(npro,2),
     &		src(npro,2)
      !
      integer advdiff
      integer	e , n ! n is possibly redundant. clean up required
!      real*8  fct1
!      real*8  tmp1,   tmp2  ,  tmp0
      !
      real*8  k,  kq,     omg,   y, mut, mut_k,
     &      mut_w, nu, nuInv, F1, F2, P1


      real*8  alpha, theta, chi, psi, ss, ssq

      real*8  alp, beta, sigk, sigw, gamma, P1, Pk, Pw

      real*8  Pk_k, Pk_w, Pw_k, Pw_w

      real*8  kInv,omgInv

      !
      !	include "elmf.h" <--- Don't know what this does
      !
      !.... Compute src and its jacobians
      !
      !.... set Menter SST (2003 paper)  k-omega Model constants
      rho(:)=datmat(1,1,1)
      rmu(:)=datmat(1,2,1)
!****      vinf(:)=datmat(1,4,2)!freestream velocity
!****      chln(:)=datmat(1,4,3)!characteristic velocity




      advdiff = 0
      if(advdiff . eq . 0) then ! not advection-diffusion
        call e3qvarkwSclr  (yl,       shgl,         xl,
     &                        gradK,   gradW,  dxidx,        WdetJ )

      do e = 1, npro

         nuInv           = rho(e)/rmu(e)  !possibly redundant. Need to Check.
         nu              = rmu(e)/rho(e)
      !   mutden          =
      !   muInv           = rho(el)
         k  		    = abs(kay(e))
         if (k . lt . 1.e-32) kInv=0
         omg		    = abs(omega(e))
         if (omg . lt . 1.e-32) omgInv=0
         y	  	    = dwall(e)
      !--------------patch
      !            if(k.gt.0.45) k=0.45
      !            if(omg.gt.2158) omg=2158
      !--------------------------------

      ! alpha = dui/dxj*dui/dxj
          alpha = gradV(e,1,1) ** 2 + gradV(e,1,2) ** 2
     &  + gradV(e,1,3) ** 2 + gradV(e,2,1) ** 2
     &  + gradV(e,2,2) ** 2 + gradV(e,2,3) ** 2
     &  + gradV(e,3,1) ** 2 + gradV(e,3,2) ** 2
     &  + gradV(e,3,3) ** 2

      ! theta = dui/dxj*duj/dxi
         theta = gradV(e,1,1) ** 2 + gradV(e,2,2) ** 2
     &  + gradV(e,3,3) ** 2 + 2*gradV(e,2,1)*gradV(e,1,2)
     &  + 2*gradV(e,3,1)*gradV(e,1,3)
     &  + 2*gradV(e,3,2)*gradV(e,2,3)

      ! chi = duk/dxk*dum/dxm = div(u)^2
         chi = (gradV(e,1,1) + gradV(e,2,2)
     &     + gradV(e,3,3)) ** 2

      ! psi = duk/dxk = div(u) = sqrt(chi)
         psi = (gradV(e,1,1) + gradV(e,2,2)
     &     + gradV(e,3,3))

         ss		= gradV(e,1,1) ** 2
     &			+ gradV(e,2,2) ** 2
     &			+ gradV(e,3,3) ** 2
     &			+ 0.5
     &			* ( (gradV(e,2,3) + gradV(e,3,2)) ** 2
     &			  + (gradV(e,3,1) + gradV(e,1,3)) ** 2
     &			  + (gradV(e,1,2) + gradV(e,2,1)) ** 2 )
  !
         ssq                 = sqrt(ss)


       delKdelW = gradK(e,1)*gradW(e,1)  +  gradK(e,2)*gradW(e,2)
     &              gradK(e,3)*gradW(e,3)
!! getblendfunc defined at the end of this subrountine Scroll down
       call getblendfunc (delKdelW, k, omg, y, rho(e), nu, F1, F2 )

       alp = F1 * alp1 + (1 - F1) * alp2
       beta = F1 * beta1 + (1 - F1) * beta2
       sigk = F1 * sigk1 + (1 - F1) * sigk2
       sigw = F1 * sigw1 + (1 - F1) * sigw2
       gamma = F1 * gam1 + (1 - F1) * gam2


       if (a1 * omg . gt . F2 * ssq) then
         mut = k / omg
         mut_k = 1.0 / omg  ! dmut/dk
         mut_w = - k / (omg ** 2) ! dmut/dw
       else
         mut = a1 * k / F2 * ssq
         mut_k = a1 / F2 * ssq
         mut_w = 0.0
       endif

       P1 = mut * ( alpha + theta + 2.0/3.0 * chi)
     &   - 2.0/3.0 * rho(e) * k * psi !! complete production term for k


       if (P1 . lt . 10 * CmuKW * rho(e) * omg * k) then
         Pk = P1 - CmuKW * rho(e) * omg * k !! complete source term for k
!! dPk/dk
         Pk_k = ( alpha + theta - 2.0/3.0 * chi ) * mut_k
     &       - 2.0/3.0 * rho(e) * psi - CmuKW * rho(e) * omg
!! dPk/dw
         Pk_w = ( alpha + theta - 2.0/3.0 * chi ) * mut_w
     &       - CmuKW * rho(e) * k


!! complete source term for w given below:
         Pw = P1 * gamma * rho(e) / mut - beta * rho(e) * ( omg ** 2 )
     &       + 2 * (1.0-F1) * rho(e) * sigw2 * delKdelW / omg
!! dPw/dk
         Pw_k = 2 * gamma * ( rho(e) ** 2 ) * k * psi * mut_k /
     &     ( mut ** 2 ) - 2 * gamma * ( rho(e) ** 2 ) * psi / mut
!! dPw/dw
         Pw_w = 2 * gamma * ( rho(e) ** 2 ) * k * psi * mut_w /
     &     ( mut ** 2 ) - 2 * beta * rho(e) * omg - 2 * ( 1.0 - F1 ) *
     &     rho(e) * sigw2 * delKdelW / ( omg ** 2 )
       else
!! production term (limiting case)
         Pk = 9 * CmuKW * rho(e) * omg * k
!! dPk/dk
         Pk_k = 9 * CmuKW * rho(e) * omg
!! dPk/dw
         Pk_w = 9 * CmuKW * rho(e) * k


!! complete source term for w given below:
         Pw = 10 * CmuKW * omg * k * gamma * ( rho(e) ** 2 ) / mut
     &     - beta * rho(e) * ( omg ** 2 )
     &       + 2 * ( 1.0 - F1 ) * rho(e) * sigw2 * delKdelW / omg
!! dPw/dk
         Pw_k = 10 * CmuKW * omg * gamma * ( rho(e) ** 2 ) / mut
     &     - 10 * CmuKW * omg * k * gamma * ( rho(e) ** 2 ) /
     &     ( mut ** 2 ) * mut_k
!! dPw/dw
         Pw_w = 10 * CmuKW * k * gamma * ( rho(e) ** 2 ) / mut
     &     - 10 * CmuKW * omg * k * gamma * ( rho(e) ** 2 ) /
     &     ( mut ** 2 ) * mut_w  -  2 * beta * rho(e) * omg
     &      - 2 * (1.0-F1) * rho(e) * sigw2 * delKdelW / ( omg ** 2 )
        endif



         src(e,1)    = Pk
         src(e,2)    = Pw

         srcRat(e,1) = -kInv*(  Pk   )
         srcRat(e,2) = -omgInv*(  Pw  )

         srcJac(e,1) = Pk_k  ! jacobian for k PDE alone
         srcJac(e,3) = Pk_w ! d(Fsrck)/d(omega)
         srcJac(e,2) = Pw_k ! d(Fsrce)/d(k)
         srcJac(e,4) = Pw_w ! jacobian for omega PDE alone




      enddo
      !
      !.... Ensure positivity of srcJac
      !
      do e = 1, npro


!         srcJac(e,1)		= max( srcJac(e,1), srcRat(e,1), 0.d0 )
!         srcJac(e,4)		= max( srcJac(e,4), srcRat(e,2), 0.d0 )

         srcJac(e,1)		= max( srcJac(e,1), 0.d0 ) !not sure if srcrat
         ! condition is to be used in the max argument
         srcJac(e,4)		= max( srcJac(e,4), 0.d0 )

      !            write(777,*) e,   srcJac(e,1),   srcJac(e,4)
!****Not sure if artifact of the elm3keps routine given below should be used
      !
!         tmp1		= min( srcJac(e,1) * srcJac(e,4),
!     &				       (srcJac(e,1)-srcRat(e,1)) *
!     &				       (srcJac(e,4)-srcRat(e,2)) )
!         tmp2		= srcJac(e,2) * srcJac(e,3)


!         if ( tmp2 .gt. tmp1 ) then
!            tmp2		= sqrt(tmp1/tmp2)
!            srcJac(e,2)	= tmp2 * srcJac(e,2)
!            srcJac(e,3)	= tmp2 * srcJac(e,3)
!         endif
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
      return
      end



!! Calculates blending functions for an element. Has to be called inside of
!! an   (e =1:npro) loop
      subroutine getblendfunc (delKdelW, kay, omg, dwl, rho, nu, F1, F2 )
        use turbKW ! access to KW model constants
        include "common.h"
!! note: all imputs are values of properties at the element level (example: rho
!!  = rho(e), kay = kay(e))  in this subroutine. Keep that in mind
!INPUTS
        double precision kay, omg, nu, dwl, rho, delKdelW
!OUTPUTS
        double precision F1, F2
!LOCALS
        double precision kq
        double precision arg11, arg11den, arg12, arg12num, arg12den
        double precision arg13, arg13num, arg13den, argCD_kw, CD_kw
        double precision arg1, arg21, arg22, arg2


        kq = sqrt(kay)


        arg11den = CmuKW*omg*dwl
        arg11    = kq / arg11den


        arg12num = 500.*nu
        arg12den = dwl * dwl * omg
        arg12    = arg12num / arg12den


        arg13num = 4.*rho(e)*sigw2*k
        argCD_kw = 2*(rho(e)* sigw2 / omg)  *   delKdelW
        CD_kw = max(argCD_kw,tenpowerminustwenty)
        arg13den = CD_kw*dwl*dwl
        arg13  = arg13num/arg13den


        arg1 = min(max(arg11 , arg12) , arg13 )

        F1 = tanh ( arg ** 4)
        if (F1 . lt . 0.0) then
          F1 = 0.0
        endif
        if (F1 . gt . 1.0) then
          F1 = 1.0
        endif

        arg21 = 2*kq / arg11den

!          arg22num = 500.*nu
!          arg22den = dwl * dwl * omg
        arg22    = arg12num / arg12den

        arg2 = max(arg21 , arg22)

        F2 = tanh( arg2 * arg2)

        if (F2 . lt . 0.0) then
          F2 = 0.0
        endif
        if (F2 . gt . 1.0) then
          F2 = 1.0
        endif

      return
      end
