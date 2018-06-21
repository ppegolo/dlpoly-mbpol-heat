      subroutine spme_for
     x  (idnode,mxnode,nospl,natms,kmax1,kmax2,kmax3,rvolm,
     x  epsq,chge,txx,tyy,tzz,fxx,fyy,fzz,bspx,bspy,bspz,
     x  bsdx,bsdy,bsdz,rcell,qqq,buffer)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using smoothed particle mesh ewald method
c     
c     parallel replicated data version (part 1)
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith oct 1998
c     
c     part 1 - reciprocal space terms (fourier part)
c
c     wl
c     2002/05/31 14:44:50
c     1.3
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension chge(mxatms),rcell(9),fff(3)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms),buffer(mxbuff)
      dimension bspx(mxspme,mxspl),bspy(mxspme,mxspl),bspz(mxspme,mxspl)
      dimension bsdx(mxspme,mxspl),bsdy(mxspme,mxspl),bsdz(mxspme,mxspl)
      complex*16 qqq(kmaxd,kmaxe,kmaxf)

#ifdef VAMPIR
      call VTBEGIN(159, ierr)
#endif
      fac=-2.d0*rvolm*r4pie0/epsq

c     set up atom numbers for nodes

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

c     calculate forces

      do i=iatm0,iatm1

        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0

        do l=1,nospl
          
          ll=int(tzz(i))-l+2
          if(ll.gt.kmax3)ll=1
          if(ll.lt.1)ll=ll+kmax3
          
          do k=1,nospl
            
            kk=int(tyy(i))-k+2
            if(kk.gt.kmax2)kk=1
            if(kk.lt.1)kk=kk+kmax2
            
            do j=1,nospl
              
              jj=int(txx(i))-j+2
              if(jj.gt.kmax1)jj=1
              if(jj.lt.1)jj=jj+kmax1
              
              qsum=real(qqq(jj,kk,ll))
              bdx=qsum*bsdx(i,j)*bspy(i,k)*bspz(i,l)*dble(kmax1)
              bdy=qsum*bspx(i,j)*bsdy(i,k)*bspz(i,l)*dble(kmax2)
              bdz=qsum*bspx(i,j)*bspy(i,k)*bsdz(i,l)*dble(kmax3)
              
              fxx(i)=fxx(i)+fac*chge(i)*(bdx*rcell(1)+bdy*rcell(2)+
     x          bdz*rcell(3))
              fyy(i)=fyy(i)+fac*chge(i)*(bdx*rcell(4)+bdy*rcell(5)+
     x          bdz*rcell(6))
              fzz(i)=fzz(i)+fac*chge(i)*(bdx*rcell(7)+bdy*rcell(8)+
     x          bdz*rcell(9))
              
            enddo
            
          enddo
          
        enddo

      enddo

c     remove COM drift arising from SPME approximations

      fff(1)=0.d0
      fff(2)=0.d0
      fff(3)=0.d0

      do i=iatm0,iatm1

        fff(1)=fff(1)+fxx(i)
        fff(2)=fff(2)+fyy(i)
        fff(3)=fff(3)+fzz(i)

      enddo

      if(mxnode.gt.1)call gdsum(fff,3,buffer)

      fff(1)=fff(1)/dble(natms)
      fff(2)=fff(2)/dble(natms)
      fff(3)=fff(3)/dble(natms)

      do i=iatm0,iatm1

        fxx(i)=fxx(i)-fff(1)
        fyy(i)=fyy(i)-fff(2)
        fzz(i)=fzz(i)-fff(3)

      enddo
#ifdef VAMPIR
      call VTEND(159, ierr)
#endif
      return
      end
