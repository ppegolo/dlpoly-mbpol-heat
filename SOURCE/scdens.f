      subroutine scdens 
     x  (idnode,imcon,mxnode,ntpvdw,natms,rvdw,dlrpot,engmet,virden,
     x  ilist,lentry,list,lstvdw,ltpvdw,ltype,buffer,cell,rho,vvv,
     x  rsqdf,xdf,xxx,ydf,yyy,zdf,zzz,elrcm,vlrcm)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating local density in metals
c     using the verlet neighbour list and sutton-chen potentials
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith june 1995
c     
c     wl
c     2002/08/01 14:46:27
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension ilist(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension lstvdw(mxvdw),ltype(mxatms),ltpvdw(mxvdw)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),elrcm(0:mxsvdw),vlrcm(0:mxsvdw)
      dimension cell(9),rho(mxatms),buffer(mxbuff)
      
#ifdef VAMPIR
      call VTBEGIN(110, ierr)
#endif
c     
c     initialise energy accumulator
      
      engmet=0.d0
      virden=0.d0
c     
c     initialise density array
      
      do i=1,natms

        rho(i)=0.d0

      enddo
      
c     
c     calculate local atomic density
      
      ii=0
      
c     
c     outer loop over atoms
      
      do i=idnode+1,natms,mxnode
        
        ii=ii+1
        
c     
c     calculate interatomic distances
        
        do k=1,lentry(ii)
          
          j=list(ii,k)
          ilist(k) = j
          
          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)
          
        enddo
        
c     
c     periodic boundary conditions
        
        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
        
c     
c     square of distances
        
        do k=1,lentry(ii)
          
          rsqdf(k) = xdf(k)**2+ydf(k)**2+zdf(k)**2
          
        enddo
        
c     
c     calculate contributions to local density

        call denloc
     x    (i,lentry(ii),ntpvdw,rvdw,dlrpot,ilist,ltype,lstvdw,ltpvdw,
     x    rsqdf,rho,vvv)
        
      enddo
      
c     
c     global sum of local atomic densities
      
      call gdsum(rho,natms,buffer)
      
c     
c     calculate square root of density and contribution to energy
      
      do i=1,natms
        
        if(rho(i).gt.0.d0)then
          
          rhosqr=sqrt(rho(i)+elrcm(ltype(i)))
          engmet=engmet+rhosqr
          rho(i)=1.d0/rhosqr
          virden=virden+vlrcm(ltype(i))/rhosqr

        endif

      enddo

      engmet=-engmet/dble(mxnode)
      
#ifdef VAMPIR
      call VTEND(110, ierr)
#endif
      return
      end


