      subroutine hkgen
     x  (idnode,nhko,nlatt,alpha,drewd,rcut,ahk,hon,fon,dhn)

c***********************************************************************
c     
c     dl_poly subroutine for generating convergence function
c     arrays for hautman klein ewald method (up to order 3 only)
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith february 2000
c     
c     wl
c     2001/06/12 13:10:08
c     1.1
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"
      
      dimension ahk(0:mxhko)
      dimension hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko),fon(mxegrd,0:7)

#ifdef VAMPIR
      call VTBEGIN(170, ierr)
#endif
      if(nhko.gt.mxhko)call error(idnode,332)

c     define effective cutoff

      ecut=rcut*dble(2*nlatt+1)

c     define grid resolution for potential arrays
      
      drewd=ecut/dble(mxegrd-4)

c     calculate HKE coefficients

      ahk(0)=1.d0

      do i=1,nhko

        ahk(i)=-0.25d0*ahk(i-1)*dble(2*i*(2*i-1))/dble(i*i)

      enddo

c     generate convergence function arrays

      do i=1,mxegrd

        hon(i,0)=0.d0
        hon(i,1)=dble(i-1)*drewd
        hon(i,2)=(2.d0*alpha/sqrpi)*exp(-(alpha*hon(i,1))**2)

      enddo

c     generate error function and derivatives by recursion

      do k=300,1,-1

        den=1.d0/dble(2*k-1)
        fac=(2.d0*alpha**2)**(k-1)

        do i=1,mxegrd

          hon(i,0)=den*(hon(i,0)*hon(i,1)**2+fac*hon(i,2))

        enddo

        if(k.le.2*nhko+2)then

          do i=1,mxegrd

            fon(i,k-1)=hon(i,0)

          enddo

        endif

      enddo
        
c     zeroth order function
c     note: hon(1,0)=2.d0*alpha/sqrpi

      do i=1,mxegrd

        hon(i,0)= fon(i,0)
        dhn(i,0)=-fon(i,1)

      enddo

      if(nhko.eq.0)then

        sss=dble(mxegrd-1)*drewd
        aaa=abs(1.d0-hon(mxegrd,nhko)*sss)
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)
c          call error(idnode,486)

        endif
#ifdef VAMPIR
        call VTEND(170, ierr)
#endif
        return

      endif

c     first order function
c     note: hon(1,1)=8.d0*alpha**3/(3.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,1)=-(2.d0*fon(i,1)-fon(i,2)*ss2)
        dhn(i,1)= (4.d0*fon(i,2)-fon(i,3)*ss2)
          
      enddo
        
      if(nhko.eq.1)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)
c          call error(idnode,486)

        endif
#ifdef VAMPIR
        call VTEND(170, ierr)
#endif
        return

      endif

c     second order function
c     note: hon(1,2)=64.d0*alpha**5/(45.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,2)=(8.d0*fon(i,2)+ss2*(-8.d0*fon(i,3)+ss2*fon(i,4)))/9.d0
        dhn(i,2)=(-24.d0*fon(i,3)+ss2*(12.d0*fon(i,4)-ss2*fon(i,5)))
     x          /9.d0

      enddo

      if(nhko.eq.2)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)
c          call error(idnode,486)

        endif
#ifdef VAMPIR
        call VTEND(170, ierr)
#endif
        return
        
      endif

c     third order function (enough for anyone!)
c     note: hon(1,3)=768.d0*alpha**7/(14175.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,3)=-(48.d0*fon(i,3)+ss2*(-72.d0*fon(i,4)+ss2*(
     x    18.d0*fon(i,5)-ss2*fon(i,6))))/225.d0
        dhn(i,3)= (192.d0*fon(i,4)+ss2*(-144.d0*fon(i,5)+ss2*(
     x    24.d0*fon(i,6)-ss2*fon(i,7))))/225.d0
        
      enddo

      if(nhko.eq.3)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)
c          call error(idnode,486)

        endif
c$$$        if(idnode.eq.0)then
c$$$          do k=0,nhko
c$$$            
c$$$            do i=2,mxegrd-1
c$$$              
c$$$              sss=dble(i-1)*drewd
c$$$              ssa=sss*alpha
c$$$              ddh=0.5d0*(hon(i+1,k)-hon(i-1,k))/(drewd*sss)
c$$$              ssh=hon(i,k)*sss**(2*k+1)
c$$$              ssd=dhn(i,k)
c$$$              write(*,'(1p,5e14.6)')ssa,ssh,ssd,ddh,ddh-ssd
c$$$              
c$$$            enddo
c$$$            
c$$$            write(*,'(a)')'&'
c$$$            
c$$$          enddo
c$$$        endif
#ifdef VAMPIR
        call VTEND(170, ierr)
#endif
        return
        
      endif
      return
      end
