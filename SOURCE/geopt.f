      subroutine geopt
     x (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm,
     x  lmetal,keyens,cell,celprp,newlst,delr,tstep,vxx,vyy,
     x  vzz,xold,yold,zold,engcp3,epsq,rcut,vircp3,lexatm,
     x  nexatm,noxatm,lentry,list,lstfrz,xxx,yyy,zzz,xdf,ydf,
     x  zdf,flx,fly,flz,chge,stresl,lct,link,uxx,uyy,uzz,
     x  buffer,lneut,lms,nneut,neulst,lstneu,lnsq,lgofr,lzeql,
     x  keyfce,multt,nstep,nstbgr,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,rprim,vircpe,virsrp,ilist,fpx,fpy,fpz,ggg,
     x  rdf,rsqdf,vvv,loglnk,kmax1,kmax2,kmax3,nhko,nlatt,
     x  ntpvdw,nospl,fplan,bplan,alpha,drewd,jlist,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,polr,rho,ewlbuf,toler,
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
     x  ntangl,engang,virang,keyang,listang,prmang,xdbc,temp,
     x  ydbc,zdbc,ntdihd,engdih,virdih,keydih,listdih,
     x  prmdih,xdcd,ydcd,zdcd,ntshl,engshl,virshl,listshl,
     x  newjob,ntinv,enginv,virinv,keyinv,listinv,prminv,prmshl,
     x  ntteth,keytet,listtet,engtet,virtet,prmtet,xxs,yys,zzs,
     x  keyfld,engfld,virfld,prmfld,weight,elrc,virlrc,nstrun,
     x  ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,xxx1,yyy1,zzz1,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,
     x  ooo,p,xi,g,h,pcom,xicom,xt,ftol)

#include "dl_params.inc"
c
c     random number seed, dynamical factor and 
c     accept-reject ratio

      parameter(iseed=17,scaldyn=1.d-3,smax=1.5d0,accrej=0.5d0)

      logical lzeql,lgofr,lpolar,lthole,lcp,lms(mxneut)
      logical newjob,newlst,lneut,loglnk,lnsq,lmetal

      complex*16 qqq,ww1,ww2,ww3,bscx,bscy,bscz

      dimension cell(9),celprp(10),elrcm(2),eta(9),stresl(9)
      dimension stress(9),vlrcm(2),npmf(2),nauxfft(4)

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
      dimension polr(mxatms),potcc(mxatms)
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
      dimension p(6*mxatms),xi(6*mxatms),g(6*mxatms),h(6*mxatms)
      dimension pcom(6*mxatms),xicom(6*mxatms),xt(6*mxatms)
      dimension ooo(6*mxatms)
      dimension nummols(mxtmls),numsit(mxtmls),numbonds(mxtmls)
      dimension wgtsit(mxsite)
      dimension listin2(4*mxatms),listbnd2(mxatms),listbnd3(mxatms)
      dimension xxx1(mxatms),yyy1(mxatms),zzz1(mxatms)
      dimension xcm(mxatms),ycm(mxatms),zcm(mxatms)
      dimension vxcm(mxatms),vycm(mxatms),vzcm(mxatms)
      dimension xdab2(mxatms),ydab2(mxatms),zdab2(mxatms)

      dimension atype1(19),atype2(4),atype3(19),atype4(4)

      data atype1/'N ','C ','N ','C ','C ','CT','CT','CT','H5','H4',
     &            'H4','H1','H1','H1','H1','H1','HC','HC','HC'/
      data atype2/'N','OS','OS','OS'/
      data atype3/'N ','C ','N ','C ','C ','C','C','C','H','H',
     &            'H','H','H','H','H','H','H','H','H'/
      data atype4/'N','O','O','O'/

      character*2 atype1,atype2,atype3,atype4

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif
c
c     initialize arrays

      if (.not.lpolar) then
        do i=1,natms
          dipx(i)=0.d0
          dipy(i)=0.d0
          dipz(i)=0.d0
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
          potcc(i)=0.d0
          polr(i)=0.d0
        enddo
      endif
c
c     convergence criterion

      ftol = ftol*engunit/418.4d0
c
      nstep = 1
c
c     geometry optimization on basin hopping
c     beta value

      beta = 1.d0/temp/boltz

c     optimization cycle

      nopt = 0
c
c     scaling factor in angstrom

      scalmc = 1.d0
c
c     initialize the variable for the lowest energy value 

      potlowest = 1.d10
c
c     counter for acceptence

      nacc = 0
c
c     for polarizable forcefield, the very first step use iteration 
c     to get self-converged dipoles

      if (lpolar) lcp = .false.

  110   continue
c
c     system center-of-mass

      wtt = 0.d0
      xscm = 0.d0
      yscm = 0.d0
      zscm = 0.d0
      do i=1,natms
        wtt=wtt+weight(i)
        xscm=xscm+weight(i)*xxx(i)
        yscm=yscm+weight(i)*yyy(i)
        zscm=zscm+weight(i)*zzz(i)
      enddo
      xscm = xscm/wtt
      yscm = yscm/wtt
      zscm = zscm/wtt
      do i=1,natms
        xxx(i) = xxx(i) - xscm
        yyy(i) = yyy(i) - yscm
        zzz(i) = zzz(i) - zscm
      enddo

c
c     save initial coordinates and dipoles

      call pack(lpolar,natms,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
c
c     optmization

      nopt = nopt + 1

      call frprmn
     x (idnode,imcon,mxnode,natms,ntpatm,engunit,rvdw,volm,
     x  lstvdw,ltpvdw,ltype,numtyp,prmvdw,dens,elrcm,vlrcm,
     x  lmetal,keyens,cell,celprp,newlst,delr,tstep,vxx,vyy,
     x  vzz,xold,yold,zold,engcp3,epsq,rcut,vircp3,lexatm,
     x  nexatm,noxatm,lentry,list,lstfrz,xxx,yyy,zzz,xdf,ydf,
     x  zdf,flx,fly,flz,chge,stresl,lct,link,uxx,uyy,uzz,
     x  buffer,lneut,lms,nneut,neulst,lstneu,lnsq,lgofr,lzeql,
     x  keyfce,multt,nstep,nstbgr,nsteql,numrdf,dlrpot,engcpe,
     x  engsrp,rprim,vircpe,virsrp,ilist,fpx,fpy,fpz,ggg,
     x  rdf,rsqdf,vvv,loglnk,kmax1,kmax2,kmax3,nhko,nlatt,
     x  ntpvdw,nospl,fplan,bplan,alpha,drewd,jlist,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,polr,rho,ewlbuf,toler,
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
     x  p,xi,g,h,pcom,xicom,xt,ftol,fret,niter,nopt)
c
c     check geometry because pbc is used

c      do i=1,natms
c        if (dabs(xxx(i)).gt.cell(1)/2.d0 .or.
c     x      dabs(yyy(i)).gt.cell(1)/2.d0 .or.
c     x      dabs(zzz(i)).gt.cell(1)/2.d0) then
c          fret = 1.d20
c          goto 130
c        endif
c      enddo

      write(90,*)nopt,niter,fret/engunit
      do i=1,natms
        write(90,*),xxx(i),yyy(i),zzz(i)
      enddo
      write(91,'(i6,2x,f20.12,2x,f12.6,2x,I6)')nopt, 
     x          fret/engunit,scalmc,niter
  130 continue
c
c     store the lowest energy configuration

      if (potlowest .gt. fret) then
        potlowest = fret
        open(56,file='min_config.com')
        write(56,*)'#am1'
        write(56,*)
        write(56,'(2I8,f12.4)')nopt,niter,fret/engunit
        write(56,*)
        write(56,*)'0  1'

        open(55,file='new_config')
        write(55,'(2I8,f12.4)')nopt,niter,fret/engunit
        write(55,*)'        0         1    600000    0.1000000000E-02'
        write(55,'(3f20.10)')(cell(i),i=1,3)
        write(55,'(3f20.10)')(cell(i),i=4,6)
        write(55,'(3f20.10)')(cell(i),i=7,9)
        nat = 0
        do j=1,ntpmls
           do k=1,nummols(j)
             do l=1,numsit(j)
               nat = nat + 1
               if (j.eq.1) then
                 write(55,'(a,I18)')atype1(l),nat
                 write(56,'(a,3x,3f12.6)')atype3(l),
     x              xxx(nat),yyy(nat),zzz(nat)
               else
                 write(55,'(a,I18)')atype2(l),nat
                 write(56,'(a,3x,3f12.6)')atype4(l),
     x              xxx(nat),yyy(nat),zzz(nat)
               endif
               write(55,'(3f20.9)')xxx(nat),yyy(nat),zzz(nat)
             enddo
           enddo
        enddo

        close(55)
        close(56)

        write(88,*)nopt,niter,fret/engunit
        do i=1,natms
          write(88,*),xxx(i),yyy(i),zzz(i)
        enddo

      endif

      if(nopt.gt.nstrun) return

      if (nopt.eq.1) then
        fret0 = fret
        nacc = nacc+1
        goto 120
      endif
c
c     acceptence critetion

      duu = fret-fret0
      acc = dexp(-beta*duu)
      dum = ran2(iseed)

      if(fret/engunit.lt.300.d0)
     x write(93,'(i9,2x,5f18.6,2x,i7)')nopt,
     x   fret0/engunit,fret/engunit,acc,dum,scalmc,niter

      if (acc .lt. dum) then
c
c     restore previous coordinates and dipoles

        call unpack(lpolar,natms,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

      else

        fret0 = fret
        nacc = nacc+1

      endif

  120 continue
c
c    monte carlo trial move
c     randomly translate molecules

      iran=int(ran2(iseed)*ntpmls)+1
      jran=int(ran2(iseed)*nummols(i))+1
      tmov1 = scalmc*(0.5d0-ran2(iseed))
      tmov2 = scalmc*(0.5d0-ran2(iseed))
      tmov3 = scalmc*(0.5d0-ran2(iseed))
      iflag=0
      ia = 0
      do i=1,ntpmls
        do j=1,nummols(i)
          do k=1,numsit(i)
            ia = ia+1
            if (i.eq.iran .and. j.eq.jran) then
              xxx(ia) = xxx(ia) + tmov1
              yyy(ia) = yyy(ia) + tmov2
              zzz(ia) = zzz(ia) + tmov3
              iflag=1
            endif
          enddo
          if(iflag.eq.1) goto 161
        enddo
      enddo
  161 continue
c
c     randamly rotate molecules 

      call cenmas
     x  (idnode,imcon,mxnode,ntbond,listbnd,cell,
     x  natms,ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,
     x  xxx,yyy,zzz,vxx,vyy,vzz,xxx1,yyy1,zzz1,xdab,ydab,zdab,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,buffer)

      iran=int(ran2(iseed)*ntpmls)+1
      jran=int(ran2(iseed)*nummols(i))+1

      rand=ran2(iseed)*scalmc
      phi=2.d0*pi*rand
      rand=ran2(iseed)*scalmc
      csthta=1.0d0-2.0d0*rand
      rand=ran2(iseed)*scalmc
      chi=2.d0*pi*rand
      thta=dacos(csthta)
      snthta=dsin(thta)
      snphi=dsin(phi)
      csphi=dcos(phi)
      snchi=dsin(chi)
      cschi=dcos(chi)
      rxx=csthta*csphi*cschi-snphi*snchi
      rxy=-csthta*csphi*snchi-snphi*cschi
      rxz=snthta*csphi
      ryx=csthta*snphi*cschi+csphi*snchi
      ryy=-csthta*snphi*snchi+csphi*cschi
      ryz=snthta*snphi
      rzx=-snthta*cschi
      rzy=snthta*snchi
      rzz=csthta

      iflag=0
      ia = 0
      ic = 0
      do i=1,ntpmls
        do j=1,nummols(i)
          ic=ic+1
          do k=1,numsit(i)
            ia=ia+1
            if (i.eq.iran .and. j.eq.jran) then
              qxx=xxx(ia)-xcm(ic)
              qyy=yyy(ia)-ycm(ic)
              qzz=zzz(ia)-zcm(ic)
              xxx(ia)=xcm(ic)+qxx*rxx+qyy*rxy+qzz*rxz
              yyy(ia)=ycm(ic)+qxx*ryx+qyy*ryy+qzz*ryz
              zzz(ia)=zcm(ic)+qxx*rzx+qyy*rzy+qzz*rzz
              iflag=1
            endif
          enddo
          if(iflag.eq.1) goto 162
        enddo
      enddo
  162 continue
c
c    displace an atom randomly

cc      do i=1,natms
c      idatm = int(ran2(iseed)*natms)+1
c      tmov = scalmc*(0.5d0-ran2(iseed))
c      xxx(idatm) = xxx(idatm) + tmov
c      tmov = scalmc*(0.5d0-ran2(iseed))
c      yyy(idatm) = yyy(idatm) + tmov
c      tmov = scalmc*(0.5d0-ran2(iseed))
c      zzz(idatm) = zzz(idatm) + tmov
cc      enddo
c
c     swap two no3 molecules

      goto 150
  140 continue
      ia1 = int(ran2(iseed)*nummols(2))+1+nummols(1)
      ia2 = int(ran2(iseed)*nummols(2))+1+nummols(1)
      if (ia1 .eq. ia2) goto 140
      cumx1 = xcm(ia1)
      dumy1 = ycm(ia1)
      dumz1 = zcm(ia1)
      cumx2 = xcm(ia2)
      dumy2 = ycm(ia2)
      dumz2 = zcm(ia2)
      if (ia1.le.nummols(1)) then
        nt1=numsit(1)
      else
        nt1=numsit(2)
      endif
      if (ia2.le.nummols(1)) then
        nt2=numsit(1)
      else
        nt2=numsit(2)
      endif
      do i=(ia1-1)*nt1+1,ia1*nt1
        xxx(i) = xxx(i) - xcm(ia1) + xcm(ia2)
        yyy(i) = yyy(i) - ycm(ia1) + ycm(ia2)
        zzz(i) = zzz(i) - zcm(ia1) + zcm(ia2)
      enddo
      do i=(ia2-1)*nt1+1,ia2*nt1
        xxx(i) = xxx(i) - xcm(ia2) + xcm(ia1)
        yyy(i) = yyy(i) - ycm(ia2) + ycm(ia1)
        zzz(i) = zzz(i) - zcm(ia2) + zcm(ia1)
      enddo
  150 continue
c
c     dynamically scale the box size

      if (mod(nopt,10).eq.0) then
        dum = dble(nacc)/10.d0
        write(92,*)nopt,dum,scalmc
        if (dum .le. accrej) then
          scalmc = scalmc - scaldyn
        else 
          scalmc = scalmc + scaldyn
        endif
        nacc=0
      endif

      goto 110

      return
      end


      subroutine pack(lpolar,natms,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

#include "dl_params.inc"
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension ooo(6*mxatms)

      do i=1,natms
        ooo(3*i-2) = xxx(i)
        ooo(3*i-1) = yyy(i)
        ooo(3*i  ) = zzz(i)
      enddo
      if (lpolar) then
        n3a = 3*natms
        do i=1,natms
          ooo(n3a+3*i-2) = dipx(i)
          ooo(n3a+3*i-1) = dipy(i)
          ooo(n3a+3*i  ) = dipz(i)
        enddo
      endif

      return
      end


      subroutine unpack(lpolar,natms,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

#include "dl_params.inc"
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension ooo(6*mxatms)

      do i=1,natms
        xxx(i) = ooo(3*i-2)
        yyy(i) = ooo(3*i-1)
        zzz(i) = ooo(3*i)
      enddo
      if (lpolar) then
        n3a = 3*natms
        do i=1,natms
          dipx(i) = ooo(n3a+3*i-2)
          dipy(i) = ooo(n3a+3*i-1)
          dipz(i) = ooo(n3a+3*i)
        enddo
      endif

      return
      end

