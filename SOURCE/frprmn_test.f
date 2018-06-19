      subroutine frprmn_test
     x  (loglnk,lgofr,lzeql,idnode,imcon,keyfce,kmax1,
     x  kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,natms,nstbgr,nstep,
     x  nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,engcpe,
     x  engsrp,epsq,rcut,rvdw,vircpe,virsrp,volm,ilist,jlist,
     x  lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,buffer,cell,chge,polr,
     x  toler,ckc,cks,clm,elc,els,emc,ckr,skr,ercp,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,ggg,rdf,rsqdf,slm,vvv,
     x  xdf,xxx,ydf,yyy,zdf,zzz,stress,rho,elrcm,vlrcm,ewlbuf,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,weight,engunit,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,
     x  p,xi,g,h,pcom,xicom,xt,ftol,fret,niter)

#include "dl_params.inc"
      parameter (itmax=5000,nrt=20,eps=1.d-10)
#include "comopt.inc"

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      call forces_opt
     x  (loglnk,lgofr,lzeql,idnode,imcon,keyfce,kmax1,
     x  kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,natms,nstbgr,nstep,
     x  nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,engcpe,
     x  engsrp,epsq,rcut,rvdw,vircpe,virsrp,volm,ilist,jlist,
     x  lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,nauxfft,
     x  lexatm2,nexatm2,key1,key2,key3,buffer,cell,chge,polr,
     x  toler,ckc,cks,clm,elc,els,emc,ckr,skr,ercp,
     x  ems,enc,ens,erc,fer,fxx,fyy,fzz,ggg,rdf,rsqdf,slm,vvv,
     x  xdf,xxx,ydf,yyy,zdf,zzz,stress,rho,elrcm,vlrcm,ewlbuf,
     x  csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqq,
     x  bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,zzn,zzd,sss,hon,
     x  dhn,pp,znp,zgc,zgs,crn,lcp,weight,engunit,
     x  lpolar,lthole,athole,dipx,dipy,dipz,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,potcc,
     x  fp)

      write(*,*)fp/engunit

      write(88,*)fp/engunit
      do i=1,natms
        p(3*i-2) = xxx(i)
        p(3*i-1) = yyy(i)
        p(3*i  ) = zzz(i)
        xi(3*i-2) = fxx(i)
        xi(3*i-1) = fyy(i)
        xi(3*i  ) = fzz(i)
        write(88,'(3f9.3)'),xxx(i),yyy(i),zzz(i)
      enddo

      if (lpolar) then

        lcp = .true.

        n3a = 3*natms
        do i=1,natms

        fdxx=-dipx(i)/polr(i)*r4pie0+efieldkx(i)+efdcrecx(i)+
     x       emux(i)+efddmurecx(i)
        fdyy=-dipy(i)/polr(i)*r4pie0+efieldky(i)+efdcrecy(i)+
     x       emuy(i)+efddmurecy(i)
        fdzz=-dipz(i)/polr(i)*r4pie0+efieldkz(i)+efdcrecz(i)+
     x       emuz(i)+efddmurecz(i)

          p(n3a+3*i-2) = dipx(i)
          p(n3a+3*i-1) = dipy(i)
          p(n3a+3*i  ) = dipz(i)
          xi(n3a+3*i-2) = fdxx
          xi(n3a+3*i-1) = fdyy
          xi(n3a+3*i  ) = fdzz

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

	return
      end
