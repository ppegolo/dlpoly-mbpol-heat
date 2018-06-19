      subroutine parlinkneu
     x  (newlst,lneut,natms,nneut,idnode,mxnode,imcon,rcut,delr,
     x  lentry,lct,link,lstneu,list,lstfrz,neulst,
     x  cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on link-cell method with neutral groups
c     frozen atoms taken into account
c     
c     to be used with the link version of exclude :excludeneu_link
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1996
c     author    - t. forester january 1996.
c     
c     wl
c     2000/01/18 14:05:51
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical lchk,newlst,linc,newjob,lfrzi,ldo,swop,lneut
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension lstfrz(mxatms)
      dimension cell(9),rcell(9),celprp(10)
c$$$  dimension nix(27),niy(27),niz(27)
      dimension nix(14),niy(14),niz(14)
      dimension lct(mxcell),link(mxatms),lstneu(mxatms)
      dimension buffer(mxbuff)
      
      save newjob
      data newjob/.true./
c$$$  data nix/0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,
c$$$  x  0,1,1,0,-1,-1,-1,0,1/
c$$$  data niy/0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,
c$$$  x  0,0,1,1,1,0,-1,-1,-1/
c$$$  data niz/0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,
c$$$  x  -1,-1,-1,-1,-1,-1,-1,-1,-1/

      data nix/0,1,0,0,-1,1,0,-1,1,0,-1,1,-1,1/
      data niy/ 0,0,1,0,1,1,-1,0,0,1,-1,-1,1,1/
      data niz/0,0,0,1,0,0,1,1,1,1,1,1,1,1/

#ifdef VAMPIR
      call VTBEGIN(13, ierr)
#endif
      lchk = .true.
      ibig = 0
      if(newlst)then
        
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)
     x    call error(idnode,300)
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
c     ratio of link cell length to cut off diameter 
        
        irat = 1
c     
c     number of subcells
        
        nsbcll = 14
        
        ilx = max(int(celprp(7)*dble(irat)/(rcut+delr)),3)
        ily = max(int(celprp(8)*dble(irat)/(rcut+delr)),3)
        ilz = max(int(celprp(9)*dble(irat)/(rcut+delr)),3)
c     
c     check there are enough link cells
        
        linc =.false.
        if(ilx.lt.2*irat+1) linc = .true.
        if(ily.lt.2*irat+1) linc = .true.
        if(ilz.lt.2*irat+1) linc = .true.
        if(linc) call error(idnode,305)
        
        ncells = ilx*ily*ilz
        if(ncells.gt.mxcell) then
          
          if(idnode.eq.0) write(nrite,*) 'mxcell must be >= ',ncells
          call  error(idnode,392)
          
        endif
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
                
                ineu = lstneu(i)
                ik = 0
c     
c     i's  group index for this processor
                
                if(mod(ineu-1,mxnode).eq.idnode)
     x            ik = ((ineu-1)/mxnode)+1
c     
c     test if i is a frozen atom
                
                lfrzi=(lstfrz(i).ne.0)
                
                if(ic.eq.jc) j=link(i)
                if(j.gt.0) then
                  
  300             continue
                  
                  jneu = lstneu(j)
c     
c     swop tests for switching of group indices,
c     ldo for 'doing' interaction
                  
                  swop = .false.
                  ldo = (ik.gt.0)
                  jneua = jneu
                  ineua = ineu
                  ika =ik
c     
c     keep only Brode-Ahlrichs pairs
                  
                  if(jneua.ge.ineua) then
                    
                    if(jneua-ineua.gt.nneut/2) then 
                      
                      swop=(mod(jneu-1,mxnode).eq.idnode)
                      if(swop) then 
                        ldo = ((nneut+ineua-jneua).le.(nneut-1)/2)
                      else
                        ldo=.false.
                      endif
                      
                    endif
                    
                  elseif(nneut+jneua-ineua.gt.(nneut-1)/2) then
                    
                    swop=(mod(jneu-1,mxnode).eq.idnode)
                    if(swop) then
                      ldo = ((ineua-jneua).le.nneut/2)
                    else
                      ldo=.false.
                    endif
                    
                  endif

                  if(swop.and.ldo) then
                    jneua = ineu
                    ineua = jneu
                    ika = ((jneu-1)/mxnode)+1
                  endif
c     
c     test of frozen atom pairs
                  
                  if(lfrzi.and.ldo) ldo=(lstfrz(j).eq.0)
c     
c     check we haven't already included this group in the list ...
                  
                  if(ldo) then
                    do jneua1 = 1,min(lentry(ika),mxlist)
                      if(list(ika,jneua1).eq.jneua) ldo=.false.
                      if(.not.ldo) goto 23
                    enddo
   23               continue
                  endif

                  if(ldo) then
                    
c     
c     distance in real space : minimum image applied
                    
                    sxd = uxx(j)-uxx(i)+cx
                    syd = uyy(j)-uyy(i)+cy
                    szd = uzz(j)-uzz(i)+cz
                    
                    xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                    yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                    zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                    
                    rsq = xd*xd+yd*yd+zd*zd
                    
c     
c     test of distance
                    if(rsq.lt.rcsq) then
                      
                      lentry(ika)=lentry(ika)+1
                      if(lentry(ika).gt.mxlist) then
                        
                        ibig = max(ibig,lentry(ika))
                        lchk = .false.
                        
                      else
                        
                        list(ika,lentry(ika))=jneua
                        
                      endif
                      
                    endif
                    
                  endif
                  
                  j=link(j)
                  if(j.ne.0) goto 300
                  
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
            write(nrite,*)'mxlist must be at least ',ibig
            write(nrite,*)'mxlist is currently ',mxlist
          endif
          call error(idnode,107)
        endif
c     
c     sort list into order ..
c     use link as a work array
        ik=0
        do i = 1+idnode,nneut,mxnode
          ik=ik+1
          do j = 1,lentry(ik)
            link(j) = list(ik,j)
          enddo
          call shellsort(lentry(ik),link)
c     
c     ensure Brode-Ahlrichs ordering
          
          i1 = lentry(ik)+1
          j1 = 0
          do j = 1,lentry(ik)
            if(link(j).ge.i) then
              j1 = j1+1
              list(ik,j1)=link(j)
              i1 = min(i1,j)
            endif
          enddo
          do j = 1,i1-1
            j1 = j1+1
            list(ik,j1)=link(j)
          enddo

        enddo
        
      endif

#ifdef VAMPIR
      call VTEND(13, ierr)
#endif
      return
      end



