!!  place in /commmon/
!!  generates gradV for easy use/retrieval. (Assembled variable)
         subroutine gradVgen(yl, gradV)

           include "common.h"

            real*8  gradV(npro,nsd,nsd),  xl(npro,nenl,nsd),
      &          shgl(npro,nsd,nshl),     dxidx(npro,nsd,nsd),
      &           shg(npro,nshl,nsd),     WdetJ(npro),
      &           yl(npro,nshl,ndof)
            integer i


            call e3metric( xl,         shgl,        dxidx,
      &               shg,        WdetJ)

            do n = 1, nshl
c
c          du_i/dx_j
c
c           i j   indices match array where V is the velocity (u in the notes)
                 gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * yl(:,n,2)
                 gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * yl(:,n,3)
                 gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * yl(:,n,4)
c
                 gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * yl(:,n,2)
                 gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * yl(:,n,3)
                 gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * yl(:,n,4)
c
                 gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * yl(:,n,2)
                 gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * yl(:,n,3)
                 gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * yl(:,n,4)

            enddo
