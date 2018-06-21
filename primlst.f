      subroutine primlst
     x  (idnode,mxnode,natms,imcon,rprim,lentry,list,
     x  cell,xxx,yyy,zzz,xdf,ydf,zdf)

c     
c*************************************************************************************
c     
c     dlpoly routine to split interaction list into primary and secondary
c     neighbours for use with multiple timestep method
c     
c     copyright daresbury laborartory
c     
c     author - t. forester february 1993
c     
c     wl
c     2000/01/18 14:05:52
c     1.3
c     Exp
c     
c************************************************************************************

#include "dl_params.inc"
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension cell(9)

#ifdef VAMPIR
      call VTBEGIN(82, ierr)
#endif

      rprim2 = rprim*rprim
      ii = 0

      do i = 1+idnode,natms,mxnode

        ii = ii + 1

        do j = 1,lentry(ii)

          k = iabs(list(ii,j))
          xdf(j) = xxx(i) - xxx(k)
          ydf(j) = yyy(i) - yyy(k)
          zdf(j) = zzz(i) - zzz(k)

        enddo           
c     
c     apply minimum image convention
        
        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)

c     assign atoms as primary or secondary


        do j = 1,lentry(ii)

c     calculate interatomic distance

          rsq = xdf(j)**2+ydf(j)**2+zdf(j)**2

          if(rsq.lt.rprim2)then
            
            list(ii,j) = -iabs(list(ii,j))
            
c     compile primary neighbour list array  : -ve indices

          else

            list(ii,j) = iabs(list(ii,j))
c     compile secondary neighbour list array : +ve indices
            
          endif

        enddo

      enddo

#ifdef VAMPIR
      call VTEND(82, ierr)
#endif
      return

      end
