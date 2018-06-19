      subroutine dcp_0
     x  (dipmas,dpms,polr2,dipx,dipy,dipz,vdxx,vdyy,vdzz,
     x  engdke,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,
     x  idnode,mxnode,natms,imcon,tstep,buffer)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     modified for dipole degrees of freedom - t. yan dec 2003
c     
c     wl
c     2000/01/18 14:05:47
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"

      dimension polr2(mxatms),dpms(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension vdxx(mxatms),vdyy(mxatms),vdzz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      
      dimension buffer(mxbuff)

      logical lfirst
      data lfirst/.true./
      save lfirst

#ifdef VAMPIR
      call VTBEGIN(32, ierr)
#endif
c     
c     set kinetic energy accumulator
      
      engdke=0.d0
c     
c     block indices
      
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c
c     assign dipole mass

      if (lfirst) then

        do i=iatm0,iatm1

          if ( polr2(i) > 1.d-6) then
             dpms(i)=dipmas/polr2(i)*r4pie0
          else
             dpms(i) = 0.d0
          endif

        enddo

        lfirst = .false.

      endif
c     
c     move atoms by leapfrog algorithm
      
      do i=iatm0,iatm1
c     
c     store old dipole velocities
        
        udxx=vdxx(i)
        udyy=vdyy(i)
        udzz=vdzz(i)
c     
c     calculate new dipole velocities

        fdxx=-dipx(i)/polr2(i)*r4pie0+efieldkx(i)+efdcrecx(i)+
     x       emux(i)+efddmurecx(i)
        fdyy=-dipy(i)/polr2(i)*r4pie0+efieldky(i)+efdcrecy(i)+
     x       emuy(i)+efddmurecy(i)
        fdzz=-dipz(i)/polr2(i)*r4pie0+efieldkz(i)+efdcrecz(i)+
     x       emuz(i)+efddmurecz(i)

c        fdxx=-dipx(i)/polr2(i)+(efieldkx(i)+efdcrecx(i)+
c     x       emux(i)+efddmurecx(i))/r4pie0
c        fdyy=-dipy(i)/polr2(i)+(efieldky(i)+efdcrecy(i)+
c     x       emuy(i)+efddmurecy(i))/r4pie0
c        fdzz=-dipz(i)/polr2(i)+(efieldkz(i)+efdcrecz(i)+
c     x       emuz(i)+efddmurecz(i))/r4pie0

        if ( dpms(i) > 0.d0 ) then
           vdxx(i)=vdxx(i)+(tstep/dpms(i))*fdxx
           vdyy(i)=vdyy(i)+(tstep/dpms(i))*fdyy
           vdzz(i)=vdzz(i)+(tstep/dpms(i))*fdzz
        else
           vdxx(i) = 0.d0
           vdyy(i) = 0.d0
           vdzz(i) = 0.d0
        endif
c     
c     update dipoles
        
        dipx(i)=dipx(i)+tstep*vdxx(i)
        dipy(i)=dipy(i)+tstep*vdyy(i)
        dipz(i)=dipz(i)+tstep*vdzz(i)
c     
c     calculate kinetic energy
        
        engdke=engdke+(dpms(i)/8.d0)*
     x    ((vdxx(i)+udxx)**2+(vdyy(i)+udyy)**2+(vdzz(i)+udzz)**2)

      enddo
c     
c     global exchange of dipoles and dipole velocities
      
      if(mxnode.gt.1)then
        
        nbuff=mxbuff
        call gdsum(engdke,1,buffer)
        call merge(idnode,mxnode,natms,nbuff,dipx,dipy,dipz,buffer)
!VB        call merge(idnode,mxnode,natms,nbuff,vdxx,vdyy,vdzz,buffer)
        
      endif
#ifdef VAMPIR
      call VTEND(32, ierr)
#endif
      return
      end
