      subroutine geopt_pt
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
     x  ntangl,engang,virang,keyang,listang,prmang,xdbc,temp,
     x  ydbc,zdbc,ntdihd,engdih,virdih,keydih,listdih,
     x  prmdih,xdcd,ydcd,zdcd,ntshl,engshl,virshl,listshl,
     x  newjob,ntinv,enginv,virinv,keyinv,listinv,prminv,prmshl,
     x  ntteth,keytet,listtet,engtet,virtet,prmtet,xxs,yys,zzs,
     x  keyfld,engfld,virfld,prmfld,weight,elrc,virlrc,nstrun,
     x  ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,xxx1,yyy1,zzz1,keyres,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,ascc,ascd,
     x  natms_r,natms_f,lttm,nttm2,listttm2,nthole,lacs,lads,
     x  ooo,p,xi,g,h,pcom,xicom,xt,ftol,npartem,temgap)
c
c     doing global minimum search using parallel tempering with
c     basin hopping Monte Carlo algorithm

#include "dl_params.inc"
c
c     random number seed, dynamical factor and 
c     accept-reject ratio

      parameter(iseed=17,scaldyn=.1d0,smax=1.5d0,accrej=0.5d0)

      logical lzeql,lgofr,lpolar,lthole,lcp,lttm,lms(mxneut)
      logical newjob,newlst,lneut,loglnk,lnsq,lmetal,lwrt,lacs,lads

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
      dimension p(6*mxatms),xi(6*mxatms),g(6*mxatms),h(6*mxatms)
      dimension pcom(6*mxatms),xicom(6*mxatms),xt(6*mxatms)
      dimension ooo(22,6*mxatms)
      dimension nummols(mxtmls),numsit(mxtmls),numbonds(mxtmls)
      dimension wgtsit(mxsite)
      dimension listin2(4*mxatms),listbnd2(mxatms),listbnd3(mxatms)
      dimension xxx1(mxatms),yyy1(mxatms),zzz1(mxatms)
      dimension xcm(mxatms),ycm(mxatms),zcm(mxatms)
      dimension vxcm(mxatms),vycm(mxatms),vzcm(mxatms)
      dimension xdab2(mxatms),ydab2(mxatms),zdab2(mxatms)
      dimension scalmc(10),nacc(10),uu0(10),beta(10),tempdum(10)
      dimension listttm2(mxatms)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif
c
c     disable pbc

c      imcon=0
c
c     set timestep

      nstep=1
c
c     write intermediate results

      lwrt = .true.
c
c     initialize storage array

      do i=1,21
        do j=1,6*mxatms
          ooo(i,j)=0.d0
        enddo
      enddo
c
c     initialize dynamical scaling factor and acceptance counter

      do i=1,npartem
        scalmc(i)=0.50d0
        nacc(i)=0
      enddo
c
c     initialize coordinates

      do ipt=1,npartem

        call pack(lpolar,natms,ipt,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
      enddo
c
c     convergence criterion

      ftol = ftol*engunit/418.4d0
c
c     monte carlo beta value 1/kT

      do i=1,npartem
        if (i.le.2) then
          tempdum(i) = temp + dble(i-1)*temgap
        else
          tempdum(i) = tempdum(i-1)**2/tempdum(i-2)
        endif
        beta(i) = 1.d0/boltz/tempdum(i)
      enddo
c
c     initialize the optimization cycle counter

      nopt = 0
c
c     initialize the variable for the lowest energy value 

      potlowest = 1.d10
c
c     start geometry optimization with parallel tempering and
c     basin hopping

  110   continue

      nopt = nopt + 1

      if(nopt.gt.nstrun) return

      do ipt=1,npartem
c
c     get coordinates

        call unpack(lpolar,natms,ipt,xxx,yyy,zzz,
     x            dipx,dipy,dipz,ooo)

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

        call pack(lpolar,natms,ipt+10,xxx,yyy,zzz,
     x            dipx,dipy,dipz,ooo)
c
c     geometry optimization

2002    continue
        call frprmn
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
     x  keyres,nthole,lacs,lads,ascc,ascd,ooo,
     x  p,xi,g,h,pcom,xicom,xt,ftol,fret,niter,nopt,ipt,lwrt)

!        write(101,'(i6,2x,f24.10,2x,2f12.4,2x,2I5)')nopt,
!     x       fret/engunit,scalmc(ipt),tempdum(ipt),niter,ipt
c
c     penelty

ccc        if(fret.lt.-100.*engunit) fret=1.d8
c
c     store the lowest energy configuration

        if (potlowest .gt. fret) then
          potlowest = fret
          call minconfig(nopt,ipt,niter,natms,ntpmls,nummols,numsit,
     x                   lttm,fret,engunit,cell,xxx,yyy,zzz)
        endif
        write(*,*)'optimized dipole:'
        do i=1,natms
        write(*,*)dipx(i)*4.8,dipy(i)*4.8,dipz(i)*4.8
        enddo
ccc	stop 'in geopt_pt'

        if(nopt.eq.1) uu0(ipt) = fret
c
c     acceptence critetion

        duu = fret-uu0(ipt)
        acc = dexp(-beta(ipt)*duu)
        dum = ran2(iseed)

        write(90+ipt,'(i9,2x,5f18.6,2x,i7)')nopt,
     x   uu0(ipt)/engunit,fret/engunit,acc,dum,scalmc(ipt),niter

        if (acc .gt. dum) then
c
c     accept the new configuration

          uu0(ipt) = fret
          nacc(ipt) = nacc(ipt)+1

        else

          call unpack(lpolar,natms,ipt+10,xxx,yyy,zzz,
     x            dipx,dipy,dipz,ooo)

        endif

        call pack(lpolar,natms,ipt,xxx,yyy,zzz,
     x            dipx,dipy,dipz,ooo)

      enddo  !end one cycle of parallel tempering

      if(nopt.ge.nstrun) return
c
c     swap two configurations

      if (npartem.gt.1) then

        nswap = int((ran2(iseed))*(npartem-1))+1

        bbb = beta(nswap)-beta(nswap+1)
        duu = uu0(nswap)-uu0(nswap+1)
        acc = dexp(bbb*duu)
        dum = ran2(iseed)

        write(102,'(3i6,2x,2f20.12,2x,2f12.6)')nopt,nswap,
     x   nswap+1,uu0(nswap)/engunit,uu0(nswap+1)/engunit,acc,dum

        ndum=0
        if (acc.ge.dum) then
c
c      swap!

          ndim=21
          nswap1=nswap+1

          call unpack(lpolar,natms,nswap,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
          call pack(lpolar,natms,ndim,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
          call unpack(lpolar,natms,nswap1,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
          call pack(lpolar,natms,nswap,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
          call unpack(lpolar,natms,ndim,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)
          call pack(lpolar,natms,nswap1,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

          udum=uu0(nswap)
          uu0(nswap)=uu0(nswap1)
          uu0(nswap1)=udum
          ndum=1

        endif

        write(102,'(3i6,2x,2f20.12,2x,2f12.6,2x,i6)')nopt,nswap,
     x   nswap+1,uu0(nswap)/engunit,uu0(nswap+1)/engunit,acc,dum,ndum

      endif
c
c     monte carlo trial move

      do ipt=1,npartem
c
c     unpack coordinates

        call unpack(lpolar,natms,ipt,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

        imove = int(ran2(iseed)+0.5d0)

        if (imove.eq.10) then
c
c     1. randomly translate molecules

          iran=int(ran2(iseed)*ntpmls)+1
          jran=int(ran2(iseed)*nummols(iran))+1
          tmov1 = scalmc(ipt)*(0.5d0-ran2(iseed))
          tmov2 = scalmc(ipt)*(0.5d0-ran2(iseed))
          tmov3 = scalmc(ipt)*(0.5d0-ran2(iseed))
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
  161     continue

        else
c
c     2. randamly rotate molecules 

          call cenmas
     x  (idnode,imcon,mxnode,ntbond,listbnd,cell,
     x  natms,ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,
     x  xxx,yyy,zzz,vxx,vyy,vzz,xxx1,yyy1,zzz1,xdab,ydab,zdab,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,buffer)

          iran=int(ran2(iseed)*ntpmls)+1
          jran=int(ran2(iseed)*nummols(iran))+1

          rand=ran2(iseed)*scalmc(ipt)/10.d0
          phi=2.d0*pi*rand
          rand=ran2(iseed)*scalmc(ipt)/10.d0
          csthta=1.0d0-2.0d0*rand
          rand=ran2(iseed)*scalmc(ipt)/10.d0
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
              if(iflag.eq.1) goto 163
            enddo
          enddo
  163     continue

        endif !end trial move of translation or rotation
c
c     swap two no3 molecules

        goto 150
  140   continue
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
  150   continue
c
c     dynamically scale the box size

        if (mod(nopt,10).eq.0) then
          dum = dble(nacc(ipt))/10.d0
          write(103,'(2i6,2x,3f10.4)')nopt,ipt,tempdum(ipt),
     x       dum,scalmc(ipt)
          if (dum .le. accrej) then
            scalmc(ipt) = scalmc(ipt)*(1.d0-scaldyn)
          else 
            scalmc(ipt) = scalmc(ipt)*(1.d0+scaldyn)
          endif
          nacc(ipt)=0
        endif

        call pack(lpolar,natms,ipt,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

      enddo  !end monte carlo trial move of parallel tempering
      goto 110

      return
      end


      subroutine pack(lpolar,natms,ndim,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

#include "dl_params.inc"
      logical lpolar
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension ooo(22,6*mxatms)

      do i=1,natms
        ooo(ndim,3*i-2) = xxx(i)
        ooo(ndim,3*i-1) = yyy(i)
        ooo(ndim,3*i  ) = zzz(i)
      enddo
      if (lpolar) then
        n3a = 3*natms
        do i=1,natms
          ooo(ndim,n3a+3*i-2) = dipx(i)
          ooo(ndim,n3a+3*i-1) = dipy(i)
          ooo(ndim,n3a+3*i  ) = dipz(i)
        enddo
      endif

      return
      end


      subroutine unpack(lpolar,natms,ndim,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

#include "dl_params.inc"
      logical lpolar
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension ooo(22,6*mxatms)

      do i=1,natms
        xxx(i) = ooo(ndim,3*i-2)
        yyy(i) = ooo(ndim,3*i-1)
        zzz(i) = ooo(ndim,3*i)
      enddo
      if (lpolar) then
        n3a = 3*natms
        do i=1,natms
          dipx(i) = ooo(ndim,n3a+3*i-2)
          dipy(i) = ooo(ndim,n3a+3*i-1)
          dipz(i) = ooo(ndim,n3a+3*i)
        enddo
      endif

      return
      end


      subroutine minconfig(nopt,ipt,niter,natms,ntpmls,nummols,numsit,
     x                   lttm,fret,engunit,cell,xxx,yyy,zzz)

#include "dl_params.inc"
      logical lttm
      dimension cell(9)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension nummols(mxtmls),numsit(mxtmls)
      dimension atype1(19),atype2(4),atype3(19),atype4(4)
      character*2 atype1,atype2,atype3,atype4
      data atype1/'OW','HW','HW','C ','C ','CT','CT','CT','H5','H4',
     &            'H4','H1','H1','H1','H1','H1','HC','HC','HC'/
      data atype2/'N','OS','OS','OS'/
      data atype3/'O ','H ','H ','C ','C ','C','C','C','H','H',
     &            'H','H','H','H','H','H','H','H','H'/
      data atype4/'N','O','O','O'/

      open(56,file='min_config.com')
      write(56,*)'#am1'
      write(56,*)
      write(56,'(3I6,f12.4)')nopt,ipt,niter,fret/engunit
      write(56,*)
      write(56,*)'0  1'

      open(55,file='new_config')
      write(55,'(3I6,f12.4)')nopt,ipt,niter,fret/engunit
      write(55,*)'        0         1    600000    0.1000000000E-02'
      write(55,'(3f20.10)')(cell(i),i=1,3)
      write(55,'(3f20.10)')(cell(i),i=4,6)
      write(55,'(3f20.10)')(cell(i),i=7,9)
      nat = 0
      ntpmls2=ntpmls
      if(lttm)ntpmls2=ntpmls-1
      do j=1,ntpmls2
        do k=1,nummols(j)
          do l=1,numsit(j)
            nat = nat + 1
            if (j.eq.1) then
              write(55,'(a,I18)')atype1(l),nat
              write(56,'(a,3x,3f12.6)')atype3(l),
     x           xxx(nat),yyy(nat),zzz(nat)
            else
              write(55,'(a,I18)')atype2(l),nat
              write(56,'(a,3x,3f12.6)')atype4(l),
     x           xxx(nat),yyy(nat),zzz(nat)
            endif
            write(55,'(3f20.9)')xxx(nat),yyy(nat),zzz(nat)
          enddo
        enddo
      enddo

      close(55)
      close(56)

cc      write(88,'(3I6,f20.12)')nopt,ipt,niter,fret/engunit
cc      do i=1,natms
cc        write(88,*)xxx(i),yyy(i),zzz(i)
cc      enddo

      return
      end
