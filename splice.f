      subroutine splice
     x      (idnode,mxnode,natms,listme,listot,xxx,yyy,zzz,buffer)
c
c*********************************************************************
c     
c     dl_poly subroutine for splicing together coordinate arrays
c     across a number of processors during shake algorithm
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1993
c     author    - w. smith       march 1993
c
c     second version of splice
c
c     wl
c     2001/05/30 12:40:24
c     1.4
c     Exp
c
c*********************************************************************
c

#include "dl_params.inc"
      
      dimension listme(mxatms),listot(mxatms)
      dimension xxx(natms),yyy(natms),zzz(natms)
      dimension buffer(mxbuff)
#ifdef VAMPIR
      call VTBEGIN(127, ierr)
#endif
c
c     check buffer size

      if(mxbuff.lt.6*natms) call error(idnode,190)
c
c     load initial transfer buffers


      j=3*natms
      n3=3*natms

      do i=1,natms

         if(listot(i).gt.0)then

            if(listme(i).gt.0)then

               buffer(j+1)=xxx(i)
               buffer(j+2)=yyy(i)
               buffer(j+3)=zzz(i)

            else

               buffer(j+1)=0.d0
               buffer(j+2)=0.d0
               buffer(j+3)=0.d0

            endif

            j=j+3

         endif

      enddo

      lastot=j-n3

c
c     splice constraint coordinates

      if(lastot.gt.0) call gdsum(buffer(n3+1),lastot,buffer(1))

c
c     reconstitute coordinate arrays

      j=n3

      do i=1,natms

         if(listot(i).gt.0)then

            xxx(i)=buffer(j+1)/dble(listot(i))
            yyy(i)=buffer(j+2)/dble(listot(i))
            zzz(i)=buffer(j+3)/dble(listot(i))

            j=j+3

         endif

      enddo

#ifdef VAMPIR
      call VTEND(127, ierr)
#endif
      return
      end
