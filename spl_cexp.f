      subroutine spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

c***********************************************************************
c     
c     dl-poly routine to create complex exponential arrays for
c     b-splines
c     
c     copyright daresbury laboratory 1998
c     
c     author w smith oct 1998
c     
c     wl
c     2002/05/31 14:44:50
c     1.3
c     Exp
c     
c***********************************************************************
      
      implicit real*8(a-h,o-z)
      
      complex*16 ww1(ndiv1),ww2(ndiv2),ww3(ndiv3)

      data tpi/6.283185307179586d0/

#ifdef VAMPIR
      call VTBEGIN(158, ierr)
#endif
c     initialise complex exponential factors
      
      ww1(1)=(1.d0,0.d0)

      do i=1,ndiv1/2

        arg=(tpi/dble(ndiv1))*dble(i)
        ww1(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww1(ndiv1+1-i)=conjg(ww1(i+1))

      enddo
      
      ww2(1)=(1.d0,0.d0)

      do i=1,ndiv2/2

        arg=(tpi/dble(ndiv2))*dble(i)
        ww2(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww2(ndiv2+1-i)=conjg(ww2(i+1))

      enddo
      
      ww3(1)=(1.d0,0.d0)

      do i=1,ndiv3/2

        arg=(tpi/dble(ndiv3))*dble(i)
        ww3(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww3(ndiv3+1-i)=conjg(ww3(i+1))

      enddo

#ifdef VAMPIR
      call VTEND(158, ierr)
#endif
      return
      end
