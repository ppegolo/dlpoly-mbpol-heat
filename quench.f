      subroutine quench
     x  (imcon,idnode,mxnode,natms,nscons,tolnce,
     x  lashap,lishap,listcon,listme,listot,
     x  buffer,cell,dxt,dyt,dzt,uxx,uyy,uzz,vxx,vyy,vzz,
     x  weight,xxt,xxx,yyt,yyy,zzt,zzz)
c     
c*********************************************************************
c     
c     dl_poly subroutine for quenching the bond energies in the 
c     initial structure of a molecule defined by constraints
c     
c     copyright - daresbury laboratory 1992
c     author w.smith november 1992
c     
c     wl
c     2000/01/18 14:05:53
c     1.3
c     Exp
c     
c*********************************************************************
c     
      
#include "dl_params.inc"
      
      logical test
      
      dimension listme(mxatms),listot(mxatms)
      dimension lashap(mxproc)
      dimension lishap(mxlshp),listcon(mxcons,3)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension buffer(mxbuff),cell(9),weight(mxatms)

#ifdef VAMPIR
      call VTBEGIN(60, ierr)
#endif
c     
c     calculate bond vectors
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
        dxt(k)=xxx(i)-xxx(j)
        dyt(k)=yyy(i)-yyy(j)
        dzt(k)=zzz(i)-zzz(j)
        
      enddo
      
      call images(imcon,0,1,nscons,cell,dxt,dyt,dzt)
      
c     
c     normalise bond vectors
      
      do k=1,nscons
        
        ddd=sqrt(dxt(k)**2+dyt(k)**2+dzt(k)**2)
        
        dxt(k)=dxt(k)/ddd
        dyt(k)=dyt(k)/ddd
        dzt(k)=dzt(k)/ddd
        
      enddo
      
c     
c     start of quenching cycle
      
      do icyc=1,mxshak
        
c     
c     initialise velocity correction arrays
        
        do i=1,natms
          
          uxx(i)=0.d0
          uyy(i)=0.d0
          uzz(i)=0.d0
          
        enddo
        
        
c     
c     calculate velocity corrections and error
        
        esig=0.d0
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          vvv=dxt(k)*(vxx(i)-vxx(j))+dyt(k)*(vyy(i)-vyy(j))+
     x      dzt(k)*(vzz(i)-vzz(j))
          
          esig=max(esig,abs(vvv))
          
          ww1=weight(j)*vvv/(weight(i)+weight(j))
          ww2=weight(i)*vvv/(weight(i)+weight(j))
          
          uxx(i)=uxx(i)-ww1*dxt(k)
          uyy(i)=uyy(i)-ww1*dyt(k)
          uzz(i)=uzz(i)-ww1*dzt(k)
          uxx(j)=uxx(j)+ww2*dxt(k)
          uyy(j)=uyy(j)+ww2*dyt(k)
          uzz(j)=uzz(j)+ww2*dzt(k)
          
        enddo
        
        test=(esig.lt.tolnce)
        
        if(mxnode.gt.1)call gstate(test)
        
        if(test)go to 100
        
c     
c     transport velocity adjustments to other nodes
        
        if(mxnode.gt.1)then
          
          call shmove
     x      (idnode,mxnode,natms,lashap,lishap,uxx,uyy,uzz,
     x      xxt,yyt,zzt,buffer)
          
        endif
        
c     
c     update velocities
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          vxx(i)=vxx(i)+uxx(i)/dble(listme(i))
          vyy(i)=vyy(i)+uyy(i)/dble(listme(i))
          vzz(i)=vzz(i)+uzz(i)/dble(listme(i))
          vxx(j)=vxx(j)+uxx(j)/dble(listme(j))
          vyy(j)=vyy(j)+uyy(j)/dble(listme(j))
          vzz(j)=vzz(j)+uzz(j)/dble(listme(j))
          
        enddo
        
      enddo
      
c     
c     error exit if quenching fails
      
      call error(idnode,70)
      
  100 continue
      
c     
c     splice velocity arrays across nodes
      
      if(mxnode.gt.1) call splice
     x  (idnode,mxnode,natms,listme,listot,vxx,vyy,vzz,buffer)

#ifdef VAMPIR
      call VTEND(60, ierr)
#endif
      return
      end
