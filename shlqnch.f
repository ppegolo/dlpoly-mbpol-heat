      subroutine shlqnch
     x     (idnode,mxnode,ntshl,temp,listshl,weight,vxx,vyy,vzz,buffer)
      
c     
c*********************************************************************
c     
c     dl_poly subroutine for quenching the internal bond energies
c     in ions defined by shell model
c     
c     copyright - daresbury laboratory 1994
c     author w.smith july  1994
c     
c     wl
c     2001/05/30 12:40:23
c     1.4
c     Exp
c
c*********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension listshl(mxshl,3)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension weight(mxatms),buffer(mxbuff)

#ifdef VAMPIR
      call VTBEGIN(62, ierr)
#endif
c     
c     permitted core-shell internal kinetic energy 
      
c      pke=3.d0*boltz
      pke=boltz*temp*1.d-4
c
c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     
c     calculate core and shell velocities from total momentum
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1
        
        i=listshl(m,2)
        j=listshl(m,3)

        rmu=(weight(i)*weight(j))/(weight(i)+weight(j))
        
        dvx=vxx(j)-vxx(i)
        dvy=vyy(j)-vyy(i)
        dvz=vzz(j)-vzz(i)

        scl=sqrt(pke/(rmu*(dvx*dvx+dvy*dvy+dvz*dvz)))

        tmx=weight(i)*vxx(i)+weight(j)*vxx(j)
        tmy=weight(i)*vyy(i)+weight(j)*vyy(j)
        tmz=weight(i)*vzz(i)+weight(j)*vzz(j)
        
        vxx(i)=tmx/(weight(i)+weight(j))-scl*rmu*dvx/weight(i)
        vxx(j)=tmx/(weight(i)+weight(j))+scl*rmu*dvx/weight(j)
        vyy(i)=tmy/(weight(i)+weight(j))-scl*rmu*dvy/weight(i)
        vyy(j)=tmy/(weight(i)+weight(j))+scl*rmu*dvy/weight(j)
        vzz(i)=tmz/(weight(i)+weight(j))-scl*rmu*dvz/weight(i)
        vzz(j)=tmz/(weight(i)+weight(j))+scl*rmu*dvz/weight(j)

      enddo

      if(mxnode.gt.1) call shlmerge
     x     (idnode,mxnode,ntshl,listshl,vxx,vyy,vzz,buffer)

#ifdef VAMPIR
      call VTEND(62, ierr)
#endif
      return

      end

