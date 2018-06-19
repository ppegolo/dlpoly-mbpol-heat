      subroutine engforce
     x (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm,
     x  lmetal,keyens,cell,celprp,newlst,delr,tstep,vxx,vyy,
     x  vzz,xold,yold,zold,engcp3,epsq,rcut,vircp3,lexatm,
     x  nexatm,noxatm,lentry,list,lstfrz,xxx,yyy,zzz,xdf,ydf,
     x  zdf,flx,fly,flz,chge,stresl,lct,link,uxx,uyy,uzz,
     x  buffer,lneut,lms,nneut,neulst,lstneu,lnsq,lgofr,lzeql,
     x  keyfce,multt,nstep,nstbgr,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,rprim,ilist,fpx,fpy,fpz,ggg,
     x  rdf,rsqdf,vvv,loglnk,kmax1,kmax2,kmax3,nhko,nlatt,
     x  ntpvdw,nospl,fplan,bplan,alpha,drewd,jlist,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,polr,polr2,rho,ewlbuf,toler,
     x  ckc,cks,clm,elc,els,emc,ckr,skr,ercp,ntptbp,rcut3b,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,slm,ntbond,ntpfbp,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,lstout,xxt,yyt,zzt,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,stress,
     x  rcuttb,engtbp,virtbp,listin,lst,lsttbp,ltptbp,prmtbp,
     x  rcutfb,engfbp,virfbp,lstfbp,ltpfbp,prmfbp,rcut4b,
     x  engbnd,virbnd,keybnd,listbnd,prmbnd,xdab,ydab,zdab,
     x  ntangl,engang,virang,keyang,listang,prmang,xdbc,
     x  ydbc,zdbc,ntdihd,engdih,virdih,keydih,listdih,
     x  prmdih,xdcd,ydcd,zdcd,ntshl,engshl,virshl,listshl,
     x  newjob,ntinv,enginv,virinv,keyinv,listinv,prminv,prmshl,
     x  ntteth,keytet,listtet,engtet,virtet,prmtet,xxs,yys,zzs,
     x  keyfld,engfld,virfld,prmfld,weight,elrc,virlrc,
     x  ntpmls,nummols,natms_r,natms_f,lttm,nttm2,listttm2,
     x  nthole,lacs,lads,ascc,ascd,keyres,
     x  uuu,gammattm)
c     
c***********************************************************************
c     
c     calculate potential energy and force
c
c***********************************************************************
c     
#include "dl_params.inc"
      
      logical lzeql,lgofr,lpolar,lthole,lcp,lttm,lms(mxneut)
      logical newjob,newlst,lneut,loglnk,lnsq,lmetal,lacs,lads

      complex*16 qqq,ww1,ww2,ww3,bscx,bscy,bscz
      
      dimension cell(9),celprp(10),elrcm(2),stresl(9)
      dimension stress(9),vlrcm(2),nauxfft(4)

      dimension numtyp(mxsvdw),ltype(mxatms),ilist(mxxdf),jlist(mxxdf)
      dimension prmvdw(mxvdw,mxpvdw),ltpvdw(mxvdw),lstvdw(mxvdw)
      dimension dens(mxsvdw)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension xold(msatms),yold(msatms),zold(msatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist),lstfrz(mxatms)
      dimension lexatm(msatms,mxexcl),nexatm(msatms)
      dimension nexatm2(msatms),lexatm2(msatms,mxexcl)
      dimension noxatm(msatms),lct(mxcell)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension buffer(mxbuff)
      dimension neulst(mxneut),lstneu(mxatms),lstout(mxatms)
      dimension flx(mxatms),fly(mxatms),flz(mxatms)
      dimension fpx(mxatms),fpy(mxatms),fpz(mxatms)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension rdf(mxrdf,mxvdw),rsqdf(mxxdf)
      dimension key1(kmaxd),key2(kmaxe),key3(kmaxf)
      dimension polr(mxatms),polr2(mxatms),potcc(mxatms)
      dimension ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld)
      dimension ckr(mxewld),skr(mxewld)
      dimension elc(mxewld,0:1),els(mxewld,0:1),ewlbuf(mxebuf)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc)
      dimension erc(mxegrd),fer(mxegrd),rho(mxatms)
      dimension ercp(mxegrd,0:3)
      dimension csp(mxspl),ffttable(mxftab)
      dimension qqc(kmaxd,kmaxe,kmaxf)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension bspx(mxspme,mxspl),bspy(mxspme,mxspl)
      dimension bspz(mxspme,mxspl),bsdx(mxspme,mxspl)
      dimension bsdy(mxspme,mxspl),bsdz(mxspme,mxspl)
      dimension qqq(kmaxd,kmaxe,kmaxf)
      dimension bscx(kmaxd),bscy(kmaxe),bscz(kmaxf)
      dimension ww1(kmaxd),ww2(kmaxe),ww3(kmaxf)
      dimension ahk(0:mxhko),zzn(mxxdf),zzd(mxxdf),sss(mxxdf)
      dimension hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko)
      dimension pp(2*mxhko),znp(mxhke,0:2*mxhko)
      dimension zgc(0:2*mxhko),zgs(0:2*mxhko),crn(0:mxhko,0:mxhko)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms),rcut4b(mxfbp)
      dimension lst(mxcell),link(mxatms),rcut3b(mxtbp)
      dimension lsttbp(mxtbp),ltptbp(mxtbp),prmtbp(mxtbp,mxptbp)
      dimension listin(mxatms)
      dimension lstfbp(mxfbp),ltpfbp(mxfbp),prmfbp(mxfbp,mxpfbp)
      dimension keybnd(mxtbnd),listbnd(mxbond,3)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension prmbnd(mxtbnd,mxpbnd)
      dimension keyang(mxtang),listang(mxangl,4),prmang(mxtang,mxpang)
      dimension xdbc(msbad),ydbc(msbad),zdbc(msbad)
      dimension keydih(mxtdih),listdih(mxdihd,5),prmdih(mxtdih,mxpdih)
      dimension xdcd(msbad),ydcd(msbad),zdcd(msbad)
      dimension keyinv(mxtinv),listinv(mxinv,5),prminv(mxtinv,mxpinv)
      dimension xxs(mxatms),yys(mxatms),zzs(mxatms)
      dimension prmtet(mxteth,mxpbnd)
      dimension listtet(msteth,2),keytet(mxteth)
      dimension listshl(mxshl,3),prmshl(mxtshl)
      dimension prmfld(mxfld),weight(mxatms)
      dimension listttm2(mxatms),nummols(mxatms)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif
c
c     for ttm2 water model, use fictitious degree of freedom

      if (lttm) then

        natms=natms_f
c
c     assign coordinates to the M sites jcp 116(2002)5115 Eq. (A5)

        mttm2=listttm2(nttm2)
        do i=1,nummols(ntpmls)
          mttm2=mttm2+1
          ioxy=listttm2(3*i-2)
          ih1=listttm2(3*i-1)
          ih2=listttm2(3*i)
          xh1o=xxx(ih1)-xxx(ioxy)
          yh1o=yyy(ih1)-yyy(ioxy)
          zh1o=zzz(ih1)-zzz(ioxy)
          xh2o=xxx(ih2)-xxx(ioxy)
          yh2o=yyy(ih2)-yyy(ioxy)
          zh2o=zzz(ih2)-zzz(ioxy)
          xxx(mttm2)=xxx(ioxy)+(xh1o+xh2o)*gammattm/2.d0
          yyy(mttm2)=yyy(ioxy)+(yh1o+yh2o)*gammattm/2.d0
          zzz(mttm2)=zzz(ioxy)+(zh1o+zh2o)*gammattm/2.d0
        enddo

      endif
c
c     initialize energy and virial accumulators for bonds,angles
c     dihedrals, inversions, short range potential and electrostatics
c     constraint virial,  com virial and field terms
      
      engbnd = 0.d0
      virbnd = 0.d0
      engang = 0.d0
      virang = 0.d0
      engdih = 0.d0
      virdih = 0.d0
      enginv = 0.d0
      virinv = 0.d0
      engtbp = 0.d0
      virtbp = 0.d0
      engfbp = 0.d0
      virfbp = 0.d0
      engsrp = 0.d0
      virsrp = 0.d0
      engcpe = 0.d0
      vircpe = 0.d0
      vircon = 0.d0
      vircom = 0.d0
      engfld = 0.d0
      virfld = 0.d0
      engshl = 0.d0
      virshl = 0.d0
      shlke  = 0.d0
      engtet = 0.d0
      virtet = 0.d0
      virpmf = 0.d0
      
      do i = 1,9
        stress(i) = 0.d0
      enddo
c
c     reset sutton chen long range corrections (constant pressure only)

      if(lmetal) then
        if(keyens.ge.4.and.keyens.le.7) then
          call lrcmetal
     x      (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x      lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm)

        endif
      endif
c     
c     initialise the force arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0

      enddo
c
c     initialize electric field and dipole force arrays
c     if polarizability is used with car-parrinello

      if (lpolar .and. lcp) then

        do i=1,natms

          efieldkx(i)=0.d0
          efieldky(i)=0.d0
          efieldkz(i)=0.d0

          efdcrecx(i)=0.d0
          efdcrecy(i)=0.d0
          efdcrecz(i)=0.d0

          efddmurecx(i)=0.d0
          efddmurecy(i)=0.d0
          efddmurecz(i)=0.d0

          emux(i)=0.d0
          emuy(i)=0.d0
          emuz(i)=0.d0

        enddo

      endif
c     
c     calculate volume of simulation cell
      
      if(imcon.ne.0.and.imcon.ne.6)then
        
        call dcell(cell,celprp)
        volm=celprp(10)
        if(imcon.eq.4)then

          volm=0.5d0*celprp(10)

        elseif(imcon.eq.5)then
        
          volm=0.5d0*celprp(10)

        elseif(imcon.eq.7)then
        
          volm=0.5d0*celprp(10)

        endif

      else
        
        volm=0.d0
        
      endif
c     
c     test for updating of Verlet list
      
c      call vertest
c     x  (newlst,idnode,mxnode,natms,delr,tstep,vxx,vyy,vzz,
c     x  xold,yold,zold)

      newlst = .true.
c     
c     set up nonbonded interaction (verlet) list
      
      newlst=((newjob).or.(newlst))

      if (newlst) then
c     
c     coulombic accumulators for interactions outside rcut
        
        engcp3 = 0.d0
        vircp3 = 0.d0

        if(.not.lneut) then
          
          if(lnsq) then 
c     
c     calculate distant interactions explicitly
            
            call parlst_nsq
     x        (newlst,natms,idnode,mxnode,imcon,engcp3,
     x        epsq,rcut,vircp3,lexatm,nexatm,noxatm,lentry,
     x        list,lstfrz,cell,xxx,yyy,zzz,xdf,ydf,zdf,flx,
     x        fly,flz,chge,stresl)
            
          elseif(loglnk) then
c     
c     ignore real space distant interactions
            
            call parlink
     x        (newlst,natms,idnode,mxnode,imcon,rcut,delr,
     x        lct,link,lexatm,nexatm,lentry,list,lstfrz,
     x        cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)

            
          else
            
            call parlst
     x        (newlst,natms,idnode,mxnode,imcon,rcut,delr,lexatm,
     x        nexatm,noxatm,lentry,list,lstfrz,cell,xxx,yyy,zzz,
     x        xdf,ydf,zdf)
            
          endif
          
        else
          
          if(.not.loglnk) then 

            call parneulst
     x        (newlst,lneut,lms,nneut,idnode,mxnode,imcon,rcut,delr,
     x        lentry,list,lstfrz,neulst,cell,xxx,yyy,zzz,xdf,ydf,zdf)
          
          else

            call parlinkneu
     x        (newlst,lneut,natms,nneut,idnode,mxnode,imcon,rcut,delr,
     x        lentry,lct,link,lstneu,list,lstfrz,neulst,
     x        cell,xxx,yyy,zzz,uxx,uyy,uzz,buffer)

          endif

        endif

        if (lttm) then
c
c     remove intramolecular m-site interactions for ttm2 water

          ii=0
          do iatm=idnode+1,natms,mxnode

            ii=ii+1

            mttm1=listttm2(1)
            mttm2=listttm2(nttm2)
            ki=1
 1005       continue
            do k=ki,lentry(ii)
              jatm=list(ii,k)
              if (iatm.gt.mttm2 .and. jatm.le.mttm2
     x            .and. jatm.ge.mttm1) then
                idum=iatm-mttm2
                jdum1=listttm2(3*idum-2)
                jdum2=listttm2(3*idum-1)
                jdum3=listttm2(3*idum)
                if (jatm.eq.jdum1 .or. jatm.eq.jdum2
     x              .or. jatm.eq.jdum3) then
                  do l=k,lentry(ii)-1
                    list(ii,l)=list(ii,l+1)
                  enddo
                  lentry(ii)=lentry(ii)-1
                  ki=k
                  goto 1005
                endif
              endif
              if (jatm.gt.mttm2 .and. iatm.le.mttm2
     x            .and. iatm.ge.mttm1) then
                idum=jatm-mttm2
                jdum1=listttm2(3*idum-2)
                jdum2=listttm2(3*idum-1)
                jdum3=listttm2(3*idum)
                if (iatm.eq.jdum1 .or. iatm.eq.jdum2
     x              .or. iatm.eq.jdum3) then
                  do l=k,lentry(ii)-1
                    list(ii,l)=list(ii,l+1)
                  enddo
                  lentry(ii)=lentry(ii)-1
                  ki=k
                  goto 1005
                endif
              endif

            enddo

          enddo

c       ii=0
c        do i=idnode+1,natms,mxnode
c        ii=ii+1
c        write(nrite,*)'i,lentry(ii)',i,lentry(ii)
c        do m=1,lentry(ii)
c        write(nrite,*)i,list(ii,m)
c        enddo
c       enddo
c	stop

        endif

      endif
c     
c     calculate pair forces, including coulombic forces

      do i = 1,9
        stress(i) = stresl(i)
      enddo

      if(lnsq) then
c     
c     multiple timestep - all-pairs
        
        call multiple_nsq
     x    (lnsq,lgofr,lzeql,newlst,idnode,imcon,keyfce,
     x    multt,mxnode,natms,nstep,nstbgr,nsteql,numrdf,
     x    delr,dlrpot,engcpe,engsrp,engcp3,epsq,rcut,
     x    rprim,rvdw,vircpe,virsrp,vircp3,ilist,lentry,
     x    list,lstvdw,ltpvdw,ltype,buffer,cell,chge,flx,fly,
     x    flz,fxx,fyy,fzz,fpx,fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,
     x    xxx,ydf,yyy,zdf,zzz,stress)
        
      elseif(.not.lneut) then         
c     
c     single timestep
        
        if(multt.eq.1) then

          call forces
     x  (lmetal,loglnk,lgofr,lzeql,idnode,imcon,keyfce,kmax1,
     x  kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,natms,nstbgr,nstep,
     x  nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,engcpe,
     x  engsrp,epsq,rcut,rvdw,vircpe,virsrp,volm,ilist,jlist,
     x  lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,buffer,cell,chge,polr,polr2,
     x  toler,ckc,cks,clm,elc,els,emc,ckr,skr,ercp,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,ggg,rdf,rsqdf,slm,vvv,
     x  xdf,xxx,ydf,yyy,zdf,zzz,stress,rho,elrcm,vlrcm,ewlbuf,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,lttm,nttm2,listttm2,lads,lacs,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  nthole,ascc,ascd,efieldkx,efieldky,efieldkz,
     x  efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,keyres)

        else

          call multiple
     x      (loglnk,lgofr,lzeql,newlst,idnode,imcon,keyfce,nlatt,
     x      kmax1,kmax2,kmax3,nhko,multt,mxnode,natms,nstep,nstbgr,
     x      nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,
     x      engcpe,engsrp,epsq,rcut,rprim,rvdw,vircpe,virsrp,volm,
     x      ilist,jlist,lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,
     x      nauxfft,key1,key2,key3,buffer,cell,chge,ckc,cks,clm,elc,
     x      els,emc,ems,enc,ens,erc,fer,flx,fly,flz,fpx,fpy,fpz,
     x      fxx,fyy,fzz,ggg,rdf,rsqdf,slm,stress,vvv,xdf,xxx,ydf,yyy,
     x      zdf,zzz,ewlbuf,csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,
     x      bsdy,bsdz,qqq,bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,
     x      zzn,zzd,sss,hon,dhn,pp,znp,zgc,zgs,crn)

        endif

      elseif(lneut) then
c     
c     neutral groups
        
        if(multt.eq.1) then
          
          call forcesneu
     x      (lgofr,lzeql,idnode,imcon,keyfce,mxnode,natms,
     x      nneut,nstbgr,nstep,nsteql,numrdf,dlrpot,engcpe,
     x      engsrp,epsq,rcut,delr,rvdw,vircpe,virsrp,ilist,jlist,
     x      lentry,lexatm,list,lstfrz,link,lstout,lstvdw,
     x      ltype,nexatm,neulst,buffer,cell,chge,fxx,fyy,
     x      fzz,ggg,rdf,rsqdf,stress,vvv,xxx,yyy,zzz,xdf,
     x      ydf,zdf,txx,tyy,tzz,xxt,yyt,zzt,uxx,uyy,uzz)
          
        else
          
          call multipleneu
     x      (lgofr,lzeql,newlst,idnode,imcon,keyfce,multt,
     x      mxnode,natms,nneut,nstbgr,nstep,nsteql,numrdf,
     x      delr,dlrpot,engcpe,engsrp,epsq,rprim,rcut,rvdw,
     x      vircpe,virsrp,ilist,jlist,lentry,lexatm,list,
     x      lstout,lstvdw,ltype,nexatm,neulst,link,lstfrz,
     x      buffer,cell,chge,fxx,fyy,fzz,flx,fly,flz,fpx,
     x      fpy,fpz,ggg,rdf,rsqdf,vvv,xdf,xxx,ydf,yyy,zdf,
     x      zzz,stress,txx,tyy,tzz,xxt,yyt,zzt,uxx,
     x      uyy,uzz)
          
        endif
        
      endif
c     
c     add in long range corrections to energy and pressure
      
      engsrp = engsrp + elrc + elrcm(1)
      virsrp = virsrp + virlrc + vlrcm(1)
c     
c     calculate three body forces

      if (ntptbp.gt.0) call thbfrc
     x  (idnode,mxnode,natms,imcon,rcuttb,engtbp,virtbp,
     x  listin,ltype,lst,lct,link,lsttbp,ltptbp,prmtbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut3b,stress,buffer)
c     
c     calculate four body forces

      if (ntpfbp.gt.0) call fbpfrc
     x  (idnode,mxnode,natms,imcon,rcutfb,engfbp,virfbp,
     x  listin,ltype,lst,lct,link,lstfbp,ltpfbp,prmfbp,cell,
     x  xxx,yyy,zzz,fxx,fyy,fzz,rcut4b,stress,buffer)
c     
c     calculate bond forces
      
      if (ntbond.gt.0) call bndfrc
     x  (idnode,imcon,mxnode,ntbond,engbnd,virbnd,keybnd,listbnd,
     x  cell,fxx,fyy,fzz,prmbnd,xxx,yyy,zzz,xdab,ydab,zdab,stress,
     x  buffer)
c
c     calculate valence angle forces
      
      if (ntangl.gt.0) call angfrc
     x  (idnode,imcon,mxnode,ntangl,engang,virang,keyang,listang,
     x  cell,fxx,fyy,fzz,prmang,xxx,yyy,zzz,xdab,ydab,zdab,xdbc,
     x  ydbc,zdbc,stress,buffer)
c     
c     calculate dihedral forces
      
      if (ntdihd.gt.0) call dihfrc
     x  (idnode,imcon,mxnode,ntdihd,keyfce,dlrpot,epsq,engcpe,
     x  engdih,engsrp,rcut,vircpe,virdih,virsrp,keydih,listdih,
     x  ltype,lstvdw,buffer,cell,chge,fxx,fyy,fzz,prmdih,xxx,yyy,
     x  zzz,xdab,ydab,zdab,xdbc,ydbc,zdbc,xdcd,ydcd,zdcd,vvv,ggg,
     x  stress)
c     
c     calculate inversion forces
      
      if (ntinv.gt.0) call invfrc
     x  (idnode,imcon,mxnode,ntinv,enginv,virinv,keyinv,listinv,
     x  buffer,cell,fxx,fyy,fzz,prminv,xxx,yyy,zzz,xdab,ydab,zdab,
     x  xdbc,ydbc,zdbc,xdcd,ydcd,zdcd,stress)

c     
c     calculate tethered atom forces
      
      if(ntteth.gt.0) call tethfrc
     x  (idnode,mxnode,imcon,natms,nstep,ntteth,keytet,listtet,
     x  engtet,virtet,buffer,cell,fxx,fyy,fzz,prmtet,xxx,yyy,
     x  zzz,xxs,yys,zzs,xdab,ydab,zdab,stress)
c
c     calculate shell model forces
      
      if (ntshl.gt.0) call shlfrc
     x   (idnode,imcon,mxnode,ntshl,engshl,virshl,listshl,cell,
     x   fxx,fyy,fzz,prmshl,xxx,yyy,zzz,xdab,ydab,zdab,stress,
     x   buffer)
c     
c     external field
      
      if(keyfld.gt.0) call extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld,cell,
     x  chge,prmfld,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,weight)
c     
c     global summation of force arrays (basic replicated data strategy)

      if(mxnode.gt.1) then
#ifdef VAMPIR
      call VTBEGIN(28, ierr)
#endif
        j=0
        do i=1,natms

          buffer(j+1)=fxx(i)
          buffer(j+2)=fyy(i)
          buffer(j+3)=fzz(i)
          j=j+3

        enddo

        call gdsum(buffer(1),3*natms,buffer(3*natms+1))

        j=0
        do i=1,natms
        
          fxx(i)=buffer(j+1)
          fyy(i)=buffer(j+2)
          fzz(i)=buffer(j+3)
          j=j+3

        enddo
#ifdef VAMPIR
      call VTEND(28, ierr)
#endif
      endif


#ifdef STRESS      
      if(mxnode.gt.1) call gdsum(stress,9,buffer)
c
c     add long range correction to diagonal terms of stress tensor

      stress(1) = stress(1) - (virlrc+vlrcm(1))/3.d0
      stress(5) = stress(5) - (virlrc+vlrcm(1))/3.d0
      stress(9) = stress(9) - (virlrc+vlrcm(1))/3.d0
#endif

      if (lttm) then
c
c     for ttm2 water use real degree of freedom

        natms=natms_r
c
c     distribute force on the m-site to oxygen and hydrogen

        mttm2=listttm2(nttm2)
        do i=1,nummols(ntpmls)
          mttm2=mttm2+1
          ioxy=listttm2(3*i-2)
          ih1=listttm2(3*i-1)
          ih2=listttm2(3*i)
          fxx(ioxy)=fxx(ioxy)+(1.d0-gammattm)*fxx(mttm2)
          fyy(ioxy)=fyy(ioxy)+(1.d0-gammattm)*fyy(mttm2)
          fzz(ioxy)=fzz(ioxy)+(1.d0-gammattm)*fzz(mttm2)
          fxx(ih1)=fxx(ih1)+fxx(mttm2)*gammattm/2.d0
          fyy(ih1)=fyy(ih1)+fyy(mttm2)*gammattm/2.d0
          fzz(ih1)=fzz(ih1)+fzz(mttm2)*gammattm/2.d0
          fxx(ih2)=fxx(ih2)+fxx(mttm2)*gammattm/2.d0
          fyy(ih2)=fyy(ih2)+fyy(mttm2)*gammattm/2.d0
          fzz(ih2)=fzz(ih2)+fzz(mttm2)*gammattm/2.d0
        enddo

      endif
c
c     frozen atoms option

      call freeze(idnode,mxnode,natms,lstfrz,vxx,vyy,vzz,fxx,fyy,fzz)
c
c     total virial (excluding constraint virial and c.o.m virial)
c      for npt routines     note: virsrp already includes virlrc

      virtot = vircpe+virsrp+virbnd+virtbp+virfld+virang+virshl+virtet
c
c     total potential energy

      uuu=engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+engfbp
     x   +engshl+enginv

      return

      write(*,*)'in engforce:'
      write(*,*)'engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+
     x   +engfbp+engshl+enginv'
      write(*,'(10f8.3)')engsrp/engunit,engcpe/engunit,engbnd/engunit,
     x   engang/engunit,engdih/engunit,engfld/engunit,engtbp/engunit,
     x   engfbp/engunit,engshl/engunit,enginv/engunit
      write(*,*)'ucc =',engcpe/engunit
      write(*,*)'uintra =',(engbnd+engang)/engunit
      write(*,*)'in engforce',uuu/engunit

        stop 'in engforce'

      return
      end
