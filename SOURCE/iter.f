! 20 OCT 05 - IUCHI - COMMENT OUT ERROR, DEBFAC 
!
      subroutine iter
     x   (iflag,idnode,mxnode,natms,imcon,
     x    dipx,dipy,dipz,efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,emux,emuy,emuz,
     x    efddmurecx,efddmurecy,efddmurecz,polr2,toler,sor_omega,
     x    buffer)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating the induced-dipole
c
c     parallel replicated data version
c
c     copyright - voth group 2003
c     author    - c. j. burnham sept 2003
c
!     Last updated: 7 July 2006 by S. Iuchi
!
c***********************************************************************
c
      
#include "dl_params.inc"

      dimension polr2(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      real(8), intent(in) :: sor_omega
      dimension buffer(mxbuff)


c
c     no external field

      extx=0.d0
      exty=0.d0
      extz=0.d0
c
c
!!!      debfac=4.8d0
      deltadip=0.d0

      oldfac=sor_omega
      fac=1.d0-oldfac

      ndip=0

!     set up atoms numbers for nodes

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

      do i=iatm0,iatm1

         if(polr2(i).gt.toler) then

            ndip=ndip+1

            oldmux=dipx(i)
            oldmuy=dipy(i)
            oldmuz=dipz(i)

            dipx(i)=oldfac*oldmux
     x           +fac*polr2(i)*(efieldkx(i)+emux(i)+efddmurecx(i)
     x           +efdcrecx(i)+extx)/r4pie0
            dipy(i)=oldfac*oldmuy
     x           +fac*polr2(i)*(efieldky(i)+emuy(i)+efddmurecy(i)
     x           +efdcrecy(i)+exty)/r4pie0
            dipz(i)=oldfac*oldmuz
     x           +fac*polr2(i)*(efieldkz(i)+emuz(i)+efddmurecz(i)
     x           +efdcrecz(i)+extz)/r4pie0

            deltadip=deltadip+(
     x         (dipx(i)-oldmux)**2
     x        +(dipy(i)-oldmuy)**2
     x        +(dipz(i)-oldmuz)**2)

         endif
      enddo

      if (mxnode.gt.1) then
         buffer(1) = deltadip
         buffer(2) = ndip
         call gdsum(buffer(1), 2, buffer(3))
         deltadip = buffer(1)
         ndip = buffer(2)

         call merge(idnode,mxnode,natms,mxbuff,dipx,
     x              dipy,dipz,buffer)
      end if

!!!      deltadip=dsqrt(deltadip/dble(ndip))/debfac
      if (ndip.eq.0) ndip=1
      deltadip = sqrt( deltadip / dble(ndip) ) 

      if(deltadip.le.toler) iflag=1

      return
      end
