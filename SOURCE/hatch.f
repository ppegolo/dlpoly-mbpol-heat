      subroutine hatch
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
     x  pcom,xicom,xt,nprt)
ccc     x  p,xi,g,h,pcom,xicom,xt,ftol,fret,niter,nopt,nprt)

#include "dl_params.inc"
#include "comopt2.inc"
ccc      dimension p(6*mxatms),xi(6*mxatms),g(6*mxatms),h(6*mxatms)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      open(nprt)
      write(nprt,*)fihatch(numtyp,mxsvdw),fihatch(ltype,mxatms)
      write(nprt,*)fihatch(ilist,mxxdf),fihatch(jlist,mxxdf)
      write(nprt,*)fhatch2(prmvdw,mxvdw,mxpvdw)
      write(nprt,*)fihatch(ltpvdw,mxvdw),fihatch(lstvdw,mxvdw)
      write(nprt,*)fhatch(dens,mxsvdw)
c 5
      write(nprt,*)fhatch(vxx,mxatms),fhatch(vyy,mxatms)
      write(nprt,*)fhatch(vzz,mxatms)
      write(nprt,*)fhatch(xold,msatms),fhatch(yold,msatms)
      write(nprt,*)fhatch(zold,msatms)
      write(nprt,*)fhatch(xxx,mxatms),fhatch(yyy,mxatms)
c 10
      write(nprt,*)fhatch(zzz,mxatms)
      write(nprt,*)fhatch(xdf,mxatms),fhatch(ydf,mxatms)
      write(nprt,*)fhatch(zdf,mxatms)
      write(nprt,*)fihatch(lentry,msatms)
      write(nprt,*)fihatch2(list,msatms,mxlist),fihatch(lstfrz,mxatms)
c 15
      write(nprt,*)fihatch2(lexatm,msatms,mxexcl)
      write(nprt,*)fihatch(nexatm,msatms)
      write(nprt,*)fihatch(nexatm2,msatms)
      write(nprt,*)fihatch2(lexatm2,msatms,mxexcl)
      write(nprt,*)fihatch(noxatm,msatms),fihatch(lct,mxcell)
c 20
      write(nprt,*)fhatch(fxx,mxatms),fhatch(fyy,mxatms)
      write(nprt,*)fhatch(fzz,mxatms)
      write(nprt,*)fhatch(chge,mxatms)
      write(nprt,*)fhatch(uxx,mxatms),fhatch(uyy,mxatms)
      write(nprt,*)fhatch(uzz,mxatms)
c 25
      write(nprt,*)fhatch(buffer,mxbuff)
      write(nprt,*)fihatch(neulst,mxneut),fihatch(lstneu,mxatms)
      write(nprt,*)fihatch(lstout,mxatms)
      write(nprt,*)fhatch(flx,mxatms),fhatch(fly,mxatms)
      write(nprt,*)fhatch(flz,mxatms)
c 30
      write(nprt,*)fhatch(fpx,mxatms),fhatch(fpy,mxatms)
      write(nprt,*)fhatch(fpz,mxatms)
      write(nprt,*)fhatch2(vvv,mxgrid,mxvdw)
      write(nprt,*)fhatch2(ggg,mxgrid,mxvdw)
      write(nprt,*)fhatch2(rdf,mxrdf,mxvdw)
c 35
      write(nprt,*)fhatch(rsqdf,mxxdf)
      write(nprt,*)fihatch(key1,kmaxd),fihatch(key2,kmaxe)
      write(nprt,*)fihatch(key3,kmaxf)
      write(nprt,*)fhatch(polr,mxatms),fhatch(potcc,mxatms)
      write(nprt,*)fhatch(ckc,mxewld),fhatch(cks,mxewld)
c 40
      write(nprt,*)fhatch(clm,mxewld),fhatch(slm,mxewld)
      write(nprt,*)fhatch(ckr,mxewld),fhatch(skr,mxewld)
      write(nprt,*)fhatch20(elc,mxewld,1),fhatch20(els,mxewld,1)
      write(nprt,*)fhatch(ewlbuf,mxebuf)
      write(nprt,*)fhatch20(emc,mxewld,kmaxb)
c 45
      write(nprt,*)fhatch20(ems,mxewld,kmaxb)
      write(nprt,*)fhatch20(enc,mxewld,kmaxc)
      write(nprt,*)fhatch20(ens,mxewld,kmaxc)
      write(nprt,*)fhatch(erc,mxegrd),fhatch(fer,mxegrd)
      write(nprt,*)fhatch(rho,mxatms)
c 50
      write(nprt,*)fhatch20(ercp,mxegrd,3)
      write(nprt,*)fhatch(csp,mxspl),fhatch(ffttable,mxftab)
      write(nprt,*)fhatch3(qqc,kmaxd,kmaxe,kmaxf)
      write(nprt,*)fhatch(txx,mxatms),fhatch(tyy,mxatms)
      write(nprt,*)fhatch(tzz,mxatms)
c 55
      write(nprt,*)fhatch2(bspx,mxspme,mxspl)
      write(nprt,*)fhatch2(bspy,mxspme,mxspl)
      write(nprt,*)fhatch2(bspz,mxspme,mxspl)
      write(nprt,*)fhatch2(bsdx,mxspme,mxspl)
      write(nprt,*)fhatch2(bsdy,mxspme,mxspl)
c 60
      write(nprt,*)fhatch2(bsdz,mxspme,mxspl)
      write(nprt,*)fhatch0(ahk,mxhko),fhatch(zzn,mxxdf)
      write(nprt,*)fhatch(zzd,mxxdf),fhatch(sss,mxxdf)
      write(nprt,*)fhatch20(hon,mxegrd,mxhko)
      write(nprt,*)fhatch20(dhn,mxegrd,mxhko)
c 65
      write(nprt,*)fhatch(pp,2*mxhko),fhatch20(znp,mxhke,2*mxhko)
      write(nprt,*)fhatch0(zgc,2*mxhko),fhatch0(zgs,2*mxhko)
      write(nprt,*)fhatch(dipx,mxatms),fhatch(dipy,mxatms)
      write(nprt,*)fhatch(dipz,mxatms)
      write(nprt,*)fhatch(emux,mxatms),fhatch(emuy,mxatms)
c 70
      write(nprt,*)fhatch(emuz,mxatms)
      write(nprt,*)fhatch(efieldkx,mxatms),fhatch(efieldky,mxatms)
      write(nprt,*)fhatch(efieldkz,mxatms)
      write(nprt,*)fhatch(efddmurecx,mxatms),fhatch(efddmurecy,mxatms)
      write(nprt,*)fhatch(efddmurecz,mxatms)
c 75
      write(nprt,*)fhatch(xxt,mxatms),fhatch(yyt,mxatms)
      write(nprt,*)fhatch(zzt,mxatms),fhatch(rcut4b,mxfbp)
      write(nprt,*)fihatch(lst,mxcell),fihatch(link,mxatms)
      write(nprt,*)fhatch(rcut3b,mxtbp)
      write(nprt,*)fihatch(lsttbp,mxtbp),fihatch(ltptbp,mxtbp)
c 80
      write(nprt,*)fhatch2(prmtbp,mxtbp,mxptbp)
      write(nprt,*)fihatch(listin,mxatms)
      write(nprt,*)fihatch(lstfbp,mxfbp),fihatch(ltpfbp,mxfbp)
      write(nprt,*)fhatch2(prmfbp,mxfbp,mxpfbp)
      write(nprt,*)fihatch(keybnd,mxtbnd),fihatch2(listbnd,mxbond,3)
c 85
      write(nprt,*)fhatch(xdab,msbad),fhatch(ydab,msbad)
      write(nprt,*)fhatch(zdab,msbad)
      write(nprt,*)fhatch2(prmbnd,mxtbnd,mxpbnd)
      write(nprt,*)fihatch(keyang,mxtang),fihatch2(listang,mxangl,4)
      write(nprt,*)fhatch2(prmang,mxtang,mxpang)
c 90
      write(nprt,*)fhatch(xdbc,msbad),fhatch(ydbc,msbad)
      write(nprt,*)fhatch(zdbc,msbad)
      write(nprt,*)fihatch(keydih,mxtdih),fihatch2(listdih,mxdihd,5)
      write(nprt,*)fhatch2(prmdih,mxtdih,mxpdih)
      write(nprt,*)fhatch(xdcd,msbad),fhatch(ydcd,msbad)
c 95
      write(nprt,*)fhatch(zdcd,msbad)
      write(nprt,*)fihatch(keyinv,mxtinv),fihatch2(listinv,mxinv,5)
      write(nprt,*)fhatch2(prminv,mxtinv,mxpinv)
      write(nprt,*)fhatch(xxs,mxatms),fhatch(yys,mxatms)
      write(nprt,*)fhatch(zzs,mxatms)
c 100
      write(nprt,*)fhatch2(prmtet,mxteth,mxpbnd)
      write(nprt,*)fihatch2(listtet,msteth,2),fihatch(keytet,mxteth)
      write(nprt,*)fihatch2(listshl,mxshl,3),fhatch(prmshl,mxtshl)
      write(nprt,*)fhatch(prmfld,mxfld),fhatch(weight,mxatms)
      write(nprt,*)fhatch(pcom,6*mxatms),fhatch(xicom,6*mxatms)
c 105
      write(nprt,*)fhatch(xt,6*mxatms)
ccc      write(nprt,*)fhatch(p,6*mxatms),fhatch(xi,6*mxatms)
ccc     x      fhatch(g,6*mxatms),fhatch(h,6*mxatms)

      close(nprt)
      return
      end

      function fhatch(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      fhatch=0.d0
      do i=1,n
        fhatch=fhatch+(dabs(a(i)))**1.2d0-a(i)
      enddo
      return
      end

      function fhatch0(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      fhatch=0.d0
      do i=0,n
        fhatch=fhatch+(dabs(a(i)))**1.2d0-a(i)
      enddo
      return
      end

      function fihatch(m,n)
      implicit real*8 (a-h,o-z)
      dimension m(n)
      fihatch=0
      do i=1,n
        fm=dble(m(i))
        fihatch=fihatch+(dabs(fm))**1.2d0-fm
      enddo
      return
      end

      function fhatch2(a,n1,n2)
      implicit real*8 (a-h,o-z)
      dimension a(n1,n2)
      fhatch2=0.d0
      do i=1,n1
        do j=1,n2
          fhatch=fhatch+(dabs(a(i,j)))**1.2d0-a(i,j)
        enddo
      enddo
      return
      end

      function fihatch2(m,n1,n2)
      implicit real*8 (a-h,o-z)
      dimension m(n1,n2)
      fihatch2=0.d0
      do i=1,n1
        do j=1,n2
          fm=dble(m(i,j))
          fihatch=fihatch+(dabs(fm))**1.2d0-fm
        enddo
      enddo
      return
      end

      function fhatch20(a,n1,n2)
      implicit real*8 (a-h,o-z)
      dimension a(n1,n2)
      fhatch2=0.d0
      do i=1,n1
        do j=0,n2
          fhatch=fhatch+(dabs(a(i,j)))**1.2d0-a(i,j)
        enddo
      enddo
      return
      end

      function fhatch3(a,n1,n2,n3)
      implicit real*8 (a-h,o-z)
      dimension a(n1,n2,n3)
      fhatch2=0.d0
      do i=1,n1
        do j=1,n2
          do k=1,n3
            fhatch=fhatch+(dabs(a(i,j,k)))**1.2d0-a(i,j,k)
          enddo
        enddo
      enddo
      return
      end

