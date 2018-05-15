      subroutine sysbook
     x  (loglnk,lneut,lshmov,lcnb,idnode,imcon,mxnode,natms,
     x  nneut,ngrp,nscons,ntangl,ntbond,ntcons,ntdihd,ntinv,
     x  ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,degfre,degrot,
     x  keybnd,keyang,lashap,lexsit,lexatm,lexatm2,lishap,listang,
     x  listbnd,listcon,listdih,listinv,listin,listme,listot,
     x  listyp,lstang,lstbnd,lstcon,lstdih,lstinv,lstfre,
     x  lstfrz,lstgtp,lstgst,lstrgd,lstcsit,nexatm,nexatm2,nexsit,
     x  numang,lstme,lstbod,lstout,lstshl,numbonds,numcon,
     x  numdih,numinv,numgsit,indpmf,lstpmf,numpmf,npmf,ind,
     x  lstpmt,listpm,numgrp,nummols,numsit,numteth,numshl,
     x  itest,index,kscons,msite,mconst,lsttet,listtet,
     x  listshl,neulst,lstneu,buffer,cell,gcmx,gcmy,gcmz,
     x  gmass,gxx,gyy,gzz,prmdih,q0,q1,q2,q3,rotinx,rotiny,
     x  rotinz,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x  accum,gaxs,rotmin)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine  defining global bookkeeping
c     arrays
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     
c     wl
c     2000/01/18 14:05:58
c     1.3
c     Exp
c     
c***********************************************************************
c     

      use ttm_forces,  only: lpolintra
      
#include "dl_params.inc"
      
      logical loglnk,safe,lneut,lshmov,lcnb
      
      dimension lstcsit(2*mxcons),lstinv(mxtinv,4)
      dimension lstneu(mxatms),neulst(mxneut)
      dimension prmdih(mxtdih,mxpdih),weight(mxatms)
      dimension numcon(mxtmls),numang(mxtmls),numbonds(mxtmls)
      dimension lstang(mxtang,3),lexatm(msatms,mxexcl),nexatm(msatms)
      dimension lexatm2(msatms,mxexcl),nexatm2(msatms)
      dimension lstbnd(mxtbnd,2),lexsit(mxsite,mxexcl),nexsit(mxsite)
      dimension keybnd(mxtbnd),lstme(mxatms),listshl(mxshl,3)
      dimension lstcon(mxtcon,2),lstdih(mxtdih,4),cell(9)
      dimension listbnd(mxbond,3),listdih(mxdihd,5),numshl(mxtmls)
      dimension listcon(mxcons,3),listang(mxangl,4),lstfrz(mxatms)
      dimension listme(mxatms),lashap(mxproc),listin(mxatms)
      dimension lishap(mxlshp),listot(mxatms),keyang(mxtang)
      dimension nummols(mxtmls),numsit(mxtmls),numdih(mxtmls)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp)
      dimension lstgtp(mxgrp),listyp(mxungp),numgrp(mxtmls)
      dimension lstgst(mxungp,mxngp),lstshl(mxtshl,2)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension lstbod(mxatms),lstout(mxatms),listinv(mxinv,5)
      dimension lsttet(mxteth),listtet(msteth,2),numteth(mxtmls)
      dimension npmf(2),numpmf(mxtmls),indpmf(mxspmf),numinv(mxtmls)
      dimension lstpmf(mxspmf,mspmf),listpm(mxpmf),lstpmt(mxpmf)
      dimension buffer(mxbuff),gmass(mxungp)
      dimension itest(mxtmls),index(mxtmls),kscons(0:mxproc-1)
      dimension msite(mxtmls),mconst(mxtmls)
      dimension ind(mxgrp,3),gaxs(mxungp,9)
      dimension rotmin(mxungp),accum(mxungp)
#ifdef VAMPIR
      call VTBEGIN(5, ierr)
#endif
c     
c     neutral group bookkeeping: sites must be listed consecutively
      
      if(lneut) then
        
        if(lstneu(1).ne.1) call error(idnode,230)
        
        neulst(1) = 1
        nneut = 1
        
        do i = 2,natms
          
          safe=.false.
          if(lstneu(i).eq.lstneu(i-1)) safe =.true.
          if(lstneu(i).eq.lstneu(i-1)+1) then
            
            safe=.true.
            nneut = nneut + 1
            if(nneut.gt.mxneut-1) call error(idnode,220)
            neulst(nneut) = i
            
          endif
          
          if(.not.safe) call error(idnode,230)
          
        enddo
        
        neulst(nneut+1)=natms+1
        
      endif
c     
c     rigid body bookkeeping 
      
      call quatbook
     x  (idnode,imcon,mxnode,natms,ngrp,ntpmls,ntfree,
     x  degfre,degrot,numgsit,numgrp,nummols,numsit,
     x  lstme,lstfrz,listyp,lstfre,listin,lstgtp,ind,
     x  lstgst,lstrgd,lstbod,xxx,yyy,zzz,weight,q0,q1,
     x  q2,q3,cell,buffer,gxx,gyy,gzz,gcmx,gcmy,gcmz,
     x  xxt,yyt,zzt,txx,tyy,tzz,gmass,rotinx,rotiny,
     x  rotinz,accum,gaxs,rotmin)

c     
c     construct list of excluded pair interactions
      
      if(lneut) then
 
        call exclude
     x    (idnode,mxnode,natms,ntpmls,keybnd,keyang,lexatm,lexsit,
     x    lstang,lstbnd,lstcon,lstdih,lstinv,nexatm,nexsit,numang,
     x    numbonds,numcon,numdih,numinv,numsit,numgrp,listyp,lstgst,
     x    numgsit,numshl,lstshl,prmdih)

        call excludeneu
     x    (idnode,mxnode,nneut,lexatm,lexsit,
     x    neulst,nexatm,nexsit,nummols,numsit)
        
      elseif(.not.lneut) then
        
        call exclude
     x    (idnode,mxnode,natms,ntpmls,keybnd,keyang,lexatm,lexsit,
     x    lstang,lstbnd,lstcon,lstdih,lstinv,nexatm,nexsit,numang,
     x    numbonds,numcon,numdih,numinv,numsit,numgrp,listyp,lstgst,
     x    numgsit,numshl,lstshl,prmdih)

        if(loglnk) then

          call exclude_link
     x      (idnode,mxnode,ntpmls,lexatm,lexsit,
     x       nexatm,nexsit,nummols,numsit,lexatm2,nexatm2)

        else

          ! commenting excluded_link out is results in considering charge
          ! contributions to polarizability from all non-excluded atoms
          ! (regardless of inter/intramolecular designations).
          if(lpolintra.eqv..false.)then
            call exclude_link
     x        (idnode,mxnode,ntpmls,lexatm,lexsit,
     x         nexatm,nexsit,nummols,numsit,lexatm2,nexatm2)
          endif

          call exclude_atom
     x      (idnode,mxnode,natms,ntpmls,lexatm,
     x      lexsit,nexatm,nexsit,nummols,numsit)

        endif

      endif
c     
c     construct interaction lists for bonded forces

      call intlist
     x  (lshmov,idnode,mxnode,natms,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntteth,ntshl,ntpmf,nspmf,lashap,lishap,
     x  listang,listbnd,listcon,listdih,listinv,listshl,listin,
     x  listme,listot,lstang,lstbnd,lstcon,lstdih,lstinv,lstfrz,
     x  lsttet,listtet,lstshl,numang,numbonds,numcon,numdih,numinv,
     x  nummols,numsit,numteth,numshl,numpmf,npmf,indpmf,lstpmf,
     x  listpm,lstpmt,itest,index,kscons,msite,mconst)

      lcnb=.false.
      if(ntcons.gt.0.and.ngrp.gt.0) then
        
        call passquat
     x    (lcnb,idnode,mxnode,natms,ngrp,nscons,ntpmls,listin,
     x    listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp,numgsit)
        
      endif
#ifdef VAMPIR
      call VTEND(5, ierr)
#endif
      return
      end
