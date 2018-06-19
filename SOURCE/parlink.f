      subroutine parlink
     x  (newlst,natms,idnode,mxnode,imcon,rcut,delr,
     x  lct,link,lexatm,nexatm,lentry,list,lstfrz,
     x  cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on link-cell method.
c     frozen atoms taken into account
c     
c     to be used with the link version of exclude :exclude_link
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester september 1993.
c     
c     wl
c     2001/06/12 12:54:50
c     1.8
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lchk,newlst,linc,newjob,lfrzi,ldo
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension lexatm(msatms,mxexcl),nexatm(msatms),lstfrz(mxatms)
      dimension cell(9),rcell(9),celprp(10)
      dimension nix(508),niy(508),niz(508)
      dimension lct(mxcell),link(mxatms)
      dimension buffer(mxbuff)
      
      save newjob
      data newjob/.true./
      data nix/0,1,0,0,-1,1,0,-1,1,0,-1,1,-1,1,2,0,0,-2,2,-1,1,0,-2,2,0,
     x0,-1,1,0,-1,1,-2,2,-2,2,-1,1,-1,1,-1,1,-2,2,0,-2,2,0,-2,2,-2,2,
     x-1,1,-2,2,-2,2,-1,1,-2,2,-2,2,3,0,0,-3,3,-1,1,0,-3,3,0,0,-1,1,0,
     x-1,1,-3,3,-3,3,-1,1,-1,1,-1,1,-3,3,-2,2,0,-3,3,0,0,-2,2,0,-2,2,
     x-3,3,-3,3,-2,2,-1,1,-3,3,-3,3,-1,1,-1,1,-2,2,-2,2,-1,1,-2,2,-3,3,
     x-3,3,-2,2,-2,2,-2,2,-3,3,0,-3,3,0,-3,3,-3,3,-1,1,-3,3,-3,3,-1,1,
     x-3,3,-3,3,-2,2,-3,3,-3,3,-2,2,-3,3,-3,3,4,0,0,-4,4,-1,1,0,-4,4,0,
     x0,-1,1,0,-1,1,-4,4,-4,4,-1,1,-1,1,-1,1,-4,4,-2,2,0,-4,4,0,0,-2,2,
     x0,-2,2,-4,4,-4,4,-2,2,-1,1,-4,4,-4,4,-1,1,-1,1,-2,2,-2,2,-1,1,-2,
     x2,-4,4,-4,4,-2,2,-2,2,-2,2,-4,4,-3,3,0,-4,4,0,0,-3,3,0,-3,3,-4,4,
     x-4,4,-3,3,-1,1,-4,4,-4,4,-1,1,-1,1,-3,3,-3,3,-1,1,-3,3,-4,4,-4,4,
     x-3,3,-2,2,-4,4,-4,4,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,-4,4,-4,4,-3,3,
     x-3,3,-3,3,-4,4,0,-4,4,0,-4,4,-4,4,-1,1,-4,4,-4,4,-1,1,-4,4,-4,4,
     x-2,2,-4,4,-4,4,-2,2,-4,4,-4,4,-3,3,-4,4,-4,4,-3,3,5,0,0,-5,5,-1,
     x1,0,-5,5,0,0,-1,1,0,-1,1,-5,5,-5,5,-1,1,-1,1,-1,1,-5,5,-2,2,0,-5,
     x5,0,0,-2,2,0,-2,2,-5,5,-5,5,-2,2,-1,1,-5,5,-5,5,-1,1,-1,1,-2,2,
     x-2,2,-1,1,-2,2,-5,5,-5,5,-2,2,-2,2,-2,2,-5,5,-3,3,0,-5,5,0,0,-3,
     x3,0,-3,3,-5,5,-5,5,-3,3,-1,1,-5,5,-5,5,-1,1,-1,1,-3,3,-3,3,-1,1,
     x-3,3,-5,5,-5,5,-3,3,-2,2,-5,5,-5,5,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,
     x-5,5,-5,5,-3,3,-3,3,-3,3/
      data niy/  0,0,1,0,1,1,-1,0,0,1,-1,-1,1,1,0,2,0,1,1,2,2,-2,0,0,2,
     x-1,0,0,1,-2,-2,-1,-1,1,1,2,2,-1,-1,1,1,2,2,-2,0,0,2,-2,-2,2,2,-2,
     x-2,-1,-1,1,1,2,2,-2,-2,2,2,0,3,0,1,1,3,3,-3,0,0,3,-1,0,0,1,-3,-3,
     x-1,-1,1,1,3,3,-1,-1,1,1,2,2,3,3,-3,0,0,3,-2,0,0,2,-3,-3,-2,-2,2,
     x2,3,3,-3,-3,-1,-1,1,1,3,3,-2,-2,-1,-1,1,1,2,2,-3,-3,-2,-2,2,2,3,
     x3,-2,-2,2,2,3,3,-3,0,0,3,-3,-3,3,3,-3,-3,-1,-1,1,1,3,3,-3,-3,3,3,
     x-3,-3,-2,-2,2,2,3,3,-3,-3,3,3,0,4,0,1,1,4,4,-4,0,0,4,-1,0,0,1,-4,
     x-4,-1,-1,1,1,4,4,-1,-1,1,1,2,2,4,4,-4,0,0,4,-2,0,0,2,-4,-4,-2,-2,
     x2,2,4,4,-4,-4,-1,-1,1,1,4,4,-2,-2,-1,-1,1,1,2,2,-4,-4,-2,-2,2,2,
     x4,4,-2,-2,2,2,3,3,4,4,-4,0,0,4,-3,0,0,3,-4,-4,-3,-3,3,3,4,4,-4,
     x-4,-1,-1,1,1,4,4,-3,-3,-1,-1,1,1,3,3,-4,-4,-3,-3,3,3,4,4,-4,-4,
     x-2,-2,2,2,4,4,-3,-3,-2,-2,2,2,3,3,-4,-4,-3,-3,3,3,4,4,-3,-3,3,3,
     x4,4,-4,0,0,4,-4,-4,4,4,-4,-4,-1,-1,1,1,4,4,-4,-4,4,4,-4,-4,-2,-2,
     x2,2,4,4,-4,-4,4,4,-4,-4,-3,-3,3,3,4,4,0,5,0,1,1,5,5,-5,0,0,5,-1,
     x0,0,1,-5,-5,-1,-1,1,1,5,5,-1,-1,1,1,2,2,5,5,-5,0,0,5,-2,0,0,2,-5,
     x-5,-2,-2,2,2,5,5,-5,-5,-1,-1,1,1,5,5,-2,-2,-1,-1,1,1,2,2,-5,-5,
     x-2,-2,2,2,5,5,-2,-2,2,2,3,3,5,5,-5,0,0,5,-3,0,0,3,-5,-5,-3,-3,3,
     x3,5,5,-5,-5,-1,-1,1,1,5,5,-3,-3,-1,-1,1,1,3,3,-5,-5,-3,-3,3,3,5,
     x5,-5,-5,-2,-2,2,2,5,5,-3,-3,-2,-2,2,2,3,3,-5,-5,-3,-3,3,3,5,5,-3,
     x-3,3,3/
      data niz/0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,2,0,0,0,0,1,1,1,1,2,2,2,
     x2,1,1,1,1,1,1,1,1,2,2,2,2,0,0,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,
     x2,2,2,0,0,3,0,0,0,0,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,3,3,3,3,0,0,
     x0,0,2,2,2,2,3,3,3,3,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
     x3,3,2,2,2,2,2,2,2,2,3,3,3,3,0,0,3,3,3,3,1,1,1,1,3,3,3,3,3,3,3,3,
     x2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,0,0,4,0,0,0,0,1,1,1,1,4,4,4,4,1,
     x1,1,1,1,1,1,1,4,4,4,4,0,0,0,0,2,2,2,2,4,4,4,4,1,1,1,1,1,1,1,1,2,
     x2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,0,0,0,0,3,
     x3,3,3,4,4,4,4,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,2,
     x2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,4,
     x4,4,4,0,0,4,4,4,4,1,1,1,1,4,4,4,4,4,4,4,4,2,2,2,2,4,4,4,4,4,4,4,
     x4,3,3,3,3,4,4,4,4,4,4,4,4,0,0,5,0,0,0,0,1,1,1,1,5,5,5,5,1,1,1,1,
     x1,1,1,1,5,5,5,5,0,0,0,0,2,2,2,2,5,5,5,5,1,1,1,1,1,1,1,1,2,2,2,2,
     x2,2,2,2,5,5,5,5,5,5,5,5,2,2,2,2,2,2,2,2,5,5,5,5,0,0,0,0,3,3,3,3,
     x5,5,5,5,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,2,2,2,2,
     x2,2,2,2,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,5,5,5,5/
#ifdef VAMPIR
      call VTBEGIN(10, ierr)
#endif
      if(newlst)then
        
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)
     x    call error(idnode,300)
        lchk = .true.
        ibig = 0
c     
c     zero link arrays

        do i = 1,natms
          link(i)=0
        enddo
c     
c     construct pair force neighbour list
        
        do i=1,msatms
          
          lentry(i)=0
          
        enddo
c     
c     real space cut off 
        
        rcsq = (rcut+delr)**2
c     
c     create mock cell vector for non-periodic system
        
        if(imcon.eq.0.or.imcon.eq.6) then
c     
c     find maximum x,y,z postions
          
          xm=0.d0
          ym=0.d0
          zm=0.d0
          
          do i =1,natms
            
            xm = max(xm,abs(xxx(i)))
            ym = max(ym,abs(yyy(i)))
            zm = max(zm,abs(zzz(i)))
            
          enddo
          
          if(imcon.eq.0) then

            cell(1) = max(2.d0*xm+rcut+delr,3.d0*(rcut+delr))
            cell(5) = max(2.d0*ym+rcut+delr,3.d0*(rcut+delr))
            cell(2) =0.d0
            cell(3) =0.d0
            cell(4) =0.d0
            cell(6) =0.d0
            cell(7) =0.d0
            cell(8) =0.d0

          endif

          cell(9) = max(2.d0*zm+rcut+delr,3.d0*(rcut+delr),cell(9))
          
        endif
        
        call dcell(cell,celprp)
        call invert(cell,rcell,det)
c     
c     ratio of link cell length to cut off diameter - max value is 5
        
c     irat = nint((rcut+delr)/rlink)
c     irat = min(max(irat,1),5)
        
        irat = 1
        
c     
c     number of subcells
        
  100   if (irat.eq.1) then 
          
          nsbcll = 14
          
        elseif(irat.eq.2) then
          
          nsbcll = 63
          
        elseif(irat.eq.3) then
          
          nsbcll = 156
          
        elseif(irat.eq.4) then
          
          nsbcll = 307
          
        elseif(irat.eq.5) then
          
          nsbcll = 508
          
        endif
        
        ilx = int(max(1.d0,celprp(7)*dble(irat)/(rcut+delr)))
        ily = int(max(1.d0,celprp(8)*dble(irat)/(rcut+delr)))
        ilz = int(max(1.d0,celprp(9)*dble(irat)/(rcut+delr)))
c     
c     check there are enough link cells
        
        ilx = max(ilx,2*irat+1)
        ily = max(ily,2*irat+1)
        ilz = max(ilz,2*irat+1)
        linc = .false.
        if(irat.eq.6) call error(idnode,305)
        if(linc) goto 100
        
        ncells = ilx*ily*ilz
        if(ncells.gt.mxcell) call error(idnode,392)
        
c     
c     calculate link cell indices
        
        do i = 1,ncells
          
          lct(i)=0
          
        enddo
c     
c     link-cell cutoff for reduced space
        
        xdc = dble(ilx)
        ydc = dble(ily)
        zdc = dble(ilz)
        
c     
c     reduced space coordinates
        if(newjob) then
          
          newjob=.false.
          call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
          nbuff = mxbuff
          if(mxnode.gt.1)  call merge
     x      (idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
          
        endif
        
        do i = 1,natms
          
          tx=xxx(i)
          ty=yyy(i)
          tz=zzz(i)

          uxx(i)=(rcell(1)*tx+rcell(4)*ty+rcell(7)*tz)+0.5d0
          uyy(i)=(rcell(2)*tx+rcell(5)*ty+rcell(8)*tz)+0.5d0
          uzz(i)=(rcell(3)*tx+rcell(6)*ty+rcell(9)*tz)+0.5d0
          
        enddo
        
c     
c     link neighbours 
        
        do i = 1,natms
          
          ix = min(int(xdc*uxx(i)),ilx-1)
          iy = min(int(ydc*uyy(i)),ily-1)
          iz = min(int(zdc*uzz(i)),ilz-1)
          
          icell = 1+ix+ilx*(iy+ily*iz)
          
          j = lct(icell)
          lct(icell)=i
          link(i)=j
          
        enddo
        
c     
c     set control variables for loop over subcells
        
        ix=1
        iy=1
        iz=1
c     
c     primary loop over subcells
        
        do ic = 1,ncells
          
          ii=lct(ic)
          if(ii.gt.0) then
c     
c     secondary loop over subcells
            
            ik=0
            
            do kc = 1,nsbcll
              
              i=ii
  
              cx = 0.d0
              cy = 0.d0
              cz = 0.d0
              jx=ix+nix(kc)
              jy=iy+niy(kc)
              jz=iz+niz(kc)
c     
c     minimum image convention
              
              if(jx.gt.ilx) then
                
                jx = jx-ilx
                cx = 1.d0
                
              elseif(jx.lt.1) then
                
                jx = jx+ilx
                cx =-1.d0
                
              endif
              
              if(jy.gt.ily) then
                
                jy = jy-ily
                cy = 1.d0
                
              elseif(jy.lt.1) then
                
                jy = jy+ily
                cy =-1.d0
                
              endif
              
              if(jz.gt.ilz) then
                
                jz = jz-ilz
                cz = 1.d0
                
              elseif(jz.lt.1) then
                
                jz = jz+ilz
                cz =-1.d0
                
              endif
c     
c     index of neighbouring cell
              
              jc =jx+ilx*((jy-1)+ily*(jz-1))
              j=lct(jc)
c     
c     ignore if empty
              
              if(j.gt.0) then
                
  200           continue
c     
c     test if site is of interest to this node
                
                if(mod(i-1,mxnode).eq.idnode) then

c     
c     i's index for this processor
                  ik = ((i-1)/mxnode)+1
c
c     test if i is a frozen atom

                  lfrzi=(lstfrz(i).ne.0)
                  
                  if(ic.eq.jc) j=link(i)
                  if(j.gt.0) then
                    
  300               continue

c     
c     test of frozen atom pairs
                    
                    ldo=.true.
                    if(lfrzi) ldo=(lstfrz(j).eq.0)
                    
                    if(ldo) then
c     
c     distance in real space : minimum image applied
                      
                      sxd = uxx(j)-uxx(i)+cx
                      syd = uyy(j)-uyy(i)+cy
                      szd = uzz(j)-uzz(i)+cz
                      
                      xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                      yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                      zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                      
                      if(imcon.eq.6)then

                        rsq = xd**2+yd**2
                        
                      else
                        
                        rsq = xd**2+yd**2+zd**2

                      endif
c     
c     test of distance
                      if(rcsq.gt.rsq) then
c     
c     test for excluded atom 
c     
                        linc = .true.
                        do ixl =1,nexatm(ik)
                          
                          if(lexatm(ik,ixl).eq.j) linc=.false.
                          
                        enddo
                        
                        if(linc) then
                          
                          lentry(ik)=lentry(ik)+1
                          
                          if(lentry(ik).gt.mxlist) then
                            
                            ibig = max(ibig,lentry(ik))
                            lchk = .false.
                            
                          else
                            
                            list(ik,lentry(ik))=j
                            
                          endif
                          
                        endif
                        
                      endif
                      
                    endif
                    
                    j=link(j)
                    if(j.ne.0) goto 300
                    
                  endif
                  
                endif
                
                j=lct(jc)
                i=link(i)
                
                if(i.ne.0) goto 200
                
              endif
              
            enddo
            
          endif
          
          ix=ix+1
          if(ix.gt.ilx) then
            
            ix=1
            iy=iy+1
            
            if(iy.gt.ily) then
              
              iy=1
              iz=iz+1
              
            endif
            
          endif
          
        enddo
c     
c     terminate job if neighbour list array exceeded
        
        if(mxnode.gt.1) call gstate(lchk)
        
        if(.not.lchk) then
          call gimax(ibig,1,idum)
          if(idnode.eq.0) then
            write(nrite,*) ' mxlist must be >=  ',ibig
            write(nrite,*) ' mxlist is currenty ',mxlist
          endif
          call error(idnode,106)
        endif
        
      endif
#ifdef VAMPIR
      call VTEND(10, ierr)
#endif
      return
      end
