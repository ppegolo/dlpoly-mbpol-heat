      subroutine shlfrc
     x   (idnode,imcon,mxnode,ntshl,engshl,virshl,listshl,cell,
     x   fxx,fyy,fzz,prmshl,xxx,yyy,zzz,xdf,ydf,zdf,stress,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating shell model spring energy and 
c     force terms in molecular dynamics.
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith        july 1994
c     
c     wl
c     2000/01/18 14:05:55
c     1.5
c     Exp
c
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension listshl(mxshl,3)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension prmshl(mxtshl),buffer(mxbuff),stress(9)
      
c     
c     define spring potential function and derivative
c     using the parameters in array prmshl
      
      vharm(r,k)=0.5d0*prmshl(k)*r**2

#ifdef VAMPIR
      call VTBEGIN(26, ierr)
#endif
c
c     check adequate workspace is available

      if(mxxdf.lt.mxshl)call error(idnode,423)

c
c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     
c     calculate core-shell separation vectors
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1
c     
c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        
c     
c     components of bond vector
        
        xdf(m)=xxx(i)-xxx(j)
        ydf(m)=yyy(i)-yyy(j)
        zdf(m)=zzz(i)-zzz(j)
        
      enddo
      
c     
c     periodic boundary condition
      
      call images(imcon,0,1,m,cell,xdf,ydf,zdf)

c     
c     zero shell energy and virial accumulators
      
      engshl=0.d0
      virshl=0.d0
      
c     
c     loop over all specified core-shell units
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1

c
c     index of potential parameters

        kk=listshl(m,1)

c     
c     core-shell separation
        
        rij=sqrt(xdf(m)**2+ydf(m)**2+zdf(m)**2)
        
c     
c     calculate scalar constant terms
        
        omega=vharm(rij,kk)
        gamma=prmshl(kk)

c     
c     calculate spring energy and virial
        
        engshl=engshl+omega
        virshl=virshl+gamma*rij*rij
        
c     
c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        
c
c     calculate spring forces

        ffx = -gamma*xdf(m)
        ffy = -gamma*ydf(m) 
        ffz = -gamma*zdf(m)

        fxx(i) = fxx(i) + ffx
        fyy(i) = fyy(i) + ffy
        fzz(i) = fzz(i) + ffz
        
        fxx(j) = fxx(j) - ffx
        fyy(j) = fyy(j) - ffy
        fzz(j) = fzz(j) - ffz

#ifdef STRESS
c     
c     calculate stress tensor
        
        stress(1) = stress(1) + xdf(m)*ffx
        stress(2) = stress(2) + xdf(m)*ffy
        stress(3) = stress(3) + xdf(m)*ffz
        
        stress(5) = stress(5) + ydf(m)*ffy
        stress(6) = stress(6) + ydf(m)*ffz
        
        stress(9) = stress(9) + zdf(m)*ffz
#endif        
      enddo
      
#ifdef STRESS
c     
c     complete stress tensor

      stress(4) = stress(2)
      stress(7) = stress(3)
      stress(8) = stress(6)
#endif      
c     
c     sum contributions to potential and virial
      
      if(mxnode.gt.1) then
        
        buffer(3)=engshl
        buffer(4)=virshl
      
        call gdsum(buffer(3),2,buffer)

        engshl=buffer(3)
        virshl=buffer(4)
        
      endif
      
#ifdef VAMPIR
      call VTEND(26, ierr)
#endif
      return
      end
