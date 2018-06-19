      subroutine bspcoe
     x  (nospl,kmax1,kmax2,kmax3,csp,bscx,bscy,bscz,ww1,ww2,ww3)

c**********************************************************************
c     
c     dl_poly subroutine to calculate B-spline coefficients for 
c     Euler exponential splines.
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith july 1998
c     
c     wl
c     2001/05/30 12:40:00
c     1.2
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      complex*16 ccc

      dimension csp(mxspl)
      complex*16 bscx(kmaxd),bscy(kmaxe),bscz(kmaxf)
      complex*16 ww1(kmaxd),ww2(kmaxe),ww3(kmaxf)

#ifdef VAMPIR
      call VTBEGIN(129, ierr)
#endif
c     calculate B-splines at knots

        csp(1)=0.d0
        csp(2)=1.d0
        
        do k=3,nospl
          
          csp(k)=0.d0
          
          do j=k,2,-1
            
            csp(j)=(dble(j-1)*csp(j)+dble(k-j+1)*csp(j-1))/dble(k-1)
            
          enddo
          
        enddo
        
c     calculate B-spline coefficients

      do i=0,kmax1-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww1(mod(i*k,kmax1)+1)

        enddo

        bscx(i+1)=ww1(mod(i*(nospl-1),kmax1)+1)/ccc

      enddo

      do i=0,kmax2-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww2(mod(i*k,kmax2)+1)

        enddo

        bscy(i+1)=ww2(mod(i*(nospl-1),kmax2)+1)/ccc

      enddo

      do i=0,kmax3-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww3(mod(i*k,kmax3)+1)

        enddo

        bscz(i+1)=ww3(mod(i*(nospl-1),kmax3)+1)/ccc

      enddo

#ifdef VAMPIR
      call VTEND(129, ierr)
#endif
      return
      end
