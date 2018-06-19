      subroutine bspgen
     x  (nospl,natms,xxx,yyy,zzz,bspx,bspy,bspz,bsdx,bsdy,bsdz)

c***********************************************************************
c
c     dl_poly subroutine to calculate B-splines for SPME method
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

      dimension xxx(mxatms),bspx(mxspme,mxspl),bsdx(mxspme,mxspl)
      dimension yyy(mxatms),bspy(mxspme,mxspl),bsdy(mxspme,mxspl)
      dimension zzz(mxatms),bspz(mxspme,mxspl),bsdz(mxspme,mxspl)

#ifdef VAMPIR
      call VTBEGIN(130, ierr)
#endif
c     construct B-splines

      do i=1,natms

        bsdx(i,1)=1.d0
        bsdy(i,1)=1.d0
        bsdz(i,1)=1.d0
        bsdx(i,2)=-1.d0
        bsdy(i,2)=-1.d0
        bsdz(i,2)=-1.d0
        bspx(i,1)=xxx(i)-int(xxx(i))
        bspy(i,1)=yyy(i)-int(yyy(i))
        bspz(i,1)=zzz(i)-int(zzz(i))
        bspx(i,2)=1.d0-xxx(i)+int(xxx(i))
        bspy(i,2)=1.d0-yyy(i)+int(yyy(i))
        bspz(i,2)=1.d0-zzz(i)+int(zzz(i))

      enddo

      do k=3,nospl

        do i=1,natms

          bspx(i,k)=0.d0
          bspy(i,k)=0.d0
          bspz(i,k)=0.d0

        enddo

        do j=k,2,-1

          if(k.eq.nospl)then
            
            do i=1,natms
              
              bsdx(i,j)=bspx(i,j)-bspx(i,j-1)
              bsdy(i,j)=bspy(i,j)-bspy(i,j-1)
              bsdz(i,j)=bspz(i,j)-bspz(i,j-1)
              
            enddo
            
          endif
          
          do i=1,natms
            
            aaa=xxx(i)+dble(j-1)-int(xxx(i))
            bbb=yyy(i)+dble(j-1)-int(yyy(i))
            ccc=zzz(i)+dble(j-1)-int(zzz(i))
            bspx(i,j)=(aaa*bspx(i,j)+(dble(k)-aaa)*bspx(i,j-1))/
     x        dble(k-1)
            bspy(i,j)=(bbb*bspy(i,j)+(dble(k)-bbb)*bspy(i,j-1))/
     x        dble(k-1)
            bspz(i,j)=(ccc*bspz(i,j)+(dble(k)-ccc)*bspz(i,j-1))/
     x        dble(k-1)

          enddo

        enddo

        if(k.eq.nospl)then
          
          do i=1,natms
            
            bsdx(i,1)=bspx(i,1)
            bsdy(i,1)=bspy(i,1)
            bsdz(i,1)=bspz(i,1)
            
          enddo
          
        endif

        do i=1,natms

          bspx(i,1)=(xxx(i)-int(xxx(i)))*bspx(i,1)/dble(k-1)
          bspy(i,1)=(yyy(i)-int(yyy(i)))*bspy(i,1)/dble(k-1)
          bspz(i,1)=(zzz(i)-int(zzz(i)))*bspz(i,1)/dble(k-1)

        enddo

      enddo

#ifdef VAMPIR
      call VTEND(130, ierr)
#endif
      return
      end

