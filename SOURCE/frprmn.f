      subroutine frprmn
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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

#include "dl_params.inc"
      parameter (itmax=500000,nrt=20,eps=1.d-20)
#include "comopt2.inc"
      dimension p(6*mxatms),xi(6*mxatms),g(6*mxatms),h(6*mxatms)
      dimension ooo(22,6*mxatms)
      logical lwrt

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif
c
c     for polarizable forcefield, the very first step use iteration
c     to get self-converged dipoles

      nrep=0
  202 continue

      if (lpolar .and. nopt.eq.1) lcp=.false.

      call engforce
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
     x  fp)

      if (lwrt) write(*,*)nopt,fp/engunit
      if (lwrt) write(88,*)nopt,fp/engunit
      do i=1,natms
        p(3*i-2) = xxx(i)
        p(3*i-1) = yyy(i)
        p(3*i  ) = zzz(i)
        xi(3*i-2) = -fxx(i)
        xi(3*i-1) = -fyy(i)
        xi(3*i  ) = -fzz(i)
cc        if (lwrt) write(88,*),xxx(i),yyy(i),zzz(i)
      enddo

      if (lpolar) then

        lcp = .true.

        n3a = 3*natms
        do i=1,natms

        fdxx=-dipx(i)/polr2(i)*r4pie0+efieldkx(i)+efdcrecx(i)+
     x       emux(i)+efddmurecx(i)
        fdyy=-dipy(i)/polr2(i)*r4pie0+efieldky(i)+efdcrecy(i)+
     x       emuy(i)+efddmurecy(i)
        fdzz=-dipz(i)/polr2(i)*r4pie0+efieldkz(i)+efdcrecz(i)+
     x       emuz(i)+efddmurecz(i)

          p(n3a+3*i-2) = dipx(i)
          p(n3a+3*i-1) = dipy(i)
          p(n3a+3*i  ) = dipz(i)
          xi(n3a+3*i-2) = -fdxx
          xi(n3a+3*i-1) = -fdyy
          xi(n3a+3*i  ) = -fdzz

        enddo
      endif

      if (lpolar) then
        n = 6*natms
      else
        n = 3*natms
      endif

      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue

      do 14 its=1,itmax

        niter=its
c
c     save coordinates and dipoles

        call pack(lpolar,natms,22,xxx,yyy,zzz,
     x          dipx,dipy,dipz,ooo)

        call linmin
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  p,xi,pcom,xicom,xt,fret)
c
c     save coordinates and dipoles

c        call pack(lpolar,natms,ipt,xxx,yyy,zzz,
c     x          dipx,dipy,dipz,ooo)

c
c     black hole?

        if(dabs(fret/fp).gt.5.) then
          call unpack(lpolar,natms,22,xxx,yyy,zzz,
     x              dipx,dipy,dipz,ooo)
          nrep=nrep+1
          if(nrep.gt.2) then
            fret=fp
            return
          endif
          goto 202
        endif

        do i=1,natms
          xxx(i) = p(3*i-2)
          yyy(i) = p(3*i-1)
          zzz(i) = p(3*i)
        enddo
        if (lpolar) then
          n3a = 3*natms
          do i=1,natms
            dipx(i) = p(n3a+3*i-2)
            dipy(i) = p(n3a+3*i-1)
            dipz(i) = p(n3a+3*i)
          enddo
        endif
 
        if (lwrt) then
          write(*,'(3i6,2x,3f24.14)')nopt,ipt,its,
     x      fp/engunit,fret/engunit,(fret-fp)/engunit

cc          if (mod(its,nrt).eq.0 .or. its.le.10) then
cc            write(88,*)nopt,its,fp/engunit
cc            do i=1,natms
cc              write(88,*)xxx(i),yyy(i),zzz(i)
cc            enddo
cc          endif
        endif

        if(2.*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+eps))
     x    return

        call engforce
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
     x  fp)

        do i=1,natms
          xi(3*i-2) = -fxx(i)
          xi(3*i-1) = -fyy(i)
          xi(3*i  ) = -fzz(i)
        enddo

        if (lpolar) then
          n3a = 3*natms
ccc	write(*,*)'in frprmn, step by step dipoles'
          do i=1,natms

            fdxx=(-dipx(i)/polr2(i)*r4pie0+efieldkx(i)+efdcrecx(i)+
     x         emux(i)+efddmurecx(i))
            fdyy=(-dipy(i)/polr2(i)*r4pie0+efieldky(i)+efdcrecy(i)+
     x         emuy(i)+efddmurecy(i))
            fdzz=(-dipz(i)/polr2(i)*r4pie0+efieldkz(i)+efdcrecz(i)+
     x         emuz(i)+efddmurecz(i))

            xi(n3a+3*i-2) = -fdxx
            xi(n3a+3*i-1) = -fdyy
            xi(n3a+3*i  ) = -fdzz

ccc        write(*,*)dipx(i)*4.8,dipy(i)*4.8,dipz(i)*4.8

          enddo
        endif

        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
ccc          dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.)return
        gam=dgg/gg
        if(nrep.gt.1)gam=gam/100.
        nrep=0
c-----------------------------------
c     adjust search step?
c-----------------------------------
ccc        if(its.gt.10 .and. dabs(gam).gt.1.)gam=gam/10.
ccc	write(*,*)gg,dgg,gam
c	if(gam.gt.10000.)then
c	do j=1,n,3
c	write(*,'i6,2x,3e20.12')1+(j-1)/3,xi(j),xi(j+1),xi(j+2)
c	enddo
c	stop
c	endif
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
        gam1=gam
14    continue
ccc      pause 'frprmn maximum iterations exceeded'
      return
      end

      subroutine linmin
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  p,xi,pcom,xicom,xt,fret)

#include "dl_params.inc"
      parameter (tol=1.d-4)
cu    uses dbrent,f1dim,mnbrak
      external f1dim
#include "comopt2.inc"
      dimension p(6*mxatms),xi(6*mxatms)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      if (lpolar) then
        n = 6*natms
      else
        n = 3*natms
      endif
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.

      call mnbrak
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax,xx,bx,fa,fx,fb,f1dim)

      fret=dbrent
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax,xx,bx,f1dim,tol,xmin)

      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      end

      function f1dim
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  x)

#include "dl_params.inc"
#include "comopt2.inc"

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      if (lpolar) then
        n = 6*natms
        n3a = 3*natms
      else
        n = 3*natms
      endif

      do 11 j=1,n
        xt(j)=pcom(j)+x*xicom(j)
11    continue

      do i=1,natms
        xxx(i) = xt(3*i-2)
        yyy(i) = xt(3*i-1)
        zzz(i) = xt(3*i)
      enddo
      if (lpolar) then
        do i=1,natms
          dipx(i) = xt(n3a+3*i-2)
          dipy(i) = xt(n3a+3*i-1)
          dipz(i) = xt(n3a+3*i)
        enddo
      endif

      call engforce
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
     x  uuu)

      f1dim=uuu

      return
      end


      subroutine mnbrak
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax,bx,cx,fa,fb,fc,func)

#include "dl_params.inc"
      parameter (gold=1.618034, glimit=2., tiny=1.d-20)
      external func
#include "comopt2.inc"

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      fa=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax)

      fb=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  bx)

      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+gold*(bx-ax)
      fc=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/
     x    (2.*dsign(max(dabs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+gold*(cx-bx)
          fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            fb=fc
            fc=fu
            fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
        else
          u=cx+gold*(cx-bx)
          fu=func
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      end


      function brent
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax,bx,cx,f,tol,xmin)

#include "dl_params.inc"
      parameter (itmax=100,cgold=.3819660,zeps=1.0d-10)
      external f
#include "comopt2.inc"

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.

      fx=f
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  x)

      fv=fx
      fw=fx
      do 11 iter=1,itmax
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+zeps
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          ppp=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) ppp=-ppp
          q=dabs(q)
          etemp=e
          e=d
          if(dabs(ppp).ge.dabs(.5*q*etemp).or.ppp.le.q*(a-x)
     x      .or.ppp.ge.q*(b-x)) goto 1
          d=ppp/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=dsign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(dabs(d).ge.tol1) then
          u=x+d
        else
          u=x+dsign(tol1,d)
        endif

        fu=f
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)

        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      end


      function dbrent
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  ax,bx,cx,f,tol,xmin)

#include "dl_params.inc"
      parameter (itmax=100,zeps=1.0d-10)
      logical ok1,ok2
      external f
#include "comopt2.inc"

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  x)
      fv=fx
      fw=fx
      dx=0.d0
      if(lpolar) n3a=3*natms
      do j=1,natms
        dx=dx+(fxx(j)*xicom(3*j-2) +
     x        fyy(j)*xicom(3*j-1) +
     x        fzz(j)*xicom(3*j))
        if (lpolar) then

          fdxx=-dipx(j)/polr2(j)*r4pie0+efieldkx(j)+efdcrecx(j)+
     x       emux(j)+efddmurecx(j)
          fdyy=-dipy(j)/polr2(j)*r4pie0+efieldky(j)+efdcrecy(j)+
     x       emuy(j)+efddmurecy(j)
          fdzz=-dipz(j)/polr2(j)*r4pie0+efieldkz(j)+efdcrecz(j)+
     x       emuz(j)+efddmurecz(j)

          dx=dx+fdxx*xicom(n3a+3*j-2) +
     x          fdyy*xicom(n3a+3*j-1) +
     x          fdzz*xicom(n3a+3*j)
        endif
      enddo
      dx=-dx
      dv=dx
      dw=dx
      do 11 iter=1,itmax
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+zeps
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
2       if(dabs(d).ge.tol1) then
          u=x+d
          fu=f
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
        else
          u=x+sign(tol1,d)
          fu=f
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
     x  lexatm2,nexatm2,key1,key2,key3,polr2,rho,ewlbuf,toler,
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
     x  pcom,xicom,xt,
     x  u)
          if(fu.gt.fx)goto 3
        endif
        du=0.d0
        do j=1,natms
          du=du+(fxx(j)*xicom(3*j-2) +
     x          fyy(j)*xicom(3*j-1) +
     x          fzz(j)*xicom(3*j))
          if (lpolar) then

          fdxx=-dipx(j)/polr2(j)*r4pie0+efieldkx(j)+efdcrecx(j)+
     x       emux(j)+efddmurecx(j)
          fdyy=-dipy(j)/polr2(j)*r4pie0+efieldky(j)+efdcrecy(j)+
     x       emuy(j)+efddmurecy(j)
          fdzz=-dipz(j)/polr2(j)*r4pie0+efieldkz(j)+efdcrecz(j)+
     x       emuz(j)+efddmurecz(j)

            du=du+fdxx*xicom(n3a+3*j-2) +
     x            fdyy*xicom(n3a+3*j-1) +
     x            fdzz*xicom(n3a+3*j)
          endif
        enddo
        du=-du
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      end

