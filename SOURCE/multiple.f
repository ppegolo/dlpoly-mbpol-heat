      subroutine multiple
     x  (loglnk,lgofr,lzeql,newlst,idnode,imcon,keyfce,nlatt,
     x  kmax1,kmax2,kmax3,nhko,multt,mxnode,natms,nstep,nstbgr,
     x  nsteql,numrdf,nospl,fplan,bplan,alpha,dlrpot,drewd,
     x  engcpe,engsrp,epsq,rcut,rprim,rvdw,vircpe,virsrp,volm,
     x  ilist,jlist,lentry,lexatm,list,lstvdw,ltpvdw,ltype,nexatm,
     x  nauxfft,key1,key2,key3,buffer,cell,chge,ckc,cks,clm,elc,
     x  els,emc,ems,enc,ens,erc,fer,flx,fly,flz,fpx,fpy,fpz,
     x  fxx,fyy,fzz,ggg,rdf,rsqdf,slm,stress,vvv,xdf,xxx,ydf,yyy,
     x  zdf,zzz,ewlbuf,csp,qqc,txx,tyy,tzz,bspx,bspy,bspz,bsdx,
     x  bsdy,bsdz,qqq,bscx,bscy,bscz,ffttable,ww1,ww2,ww3,ahk,
     x  zzn,zzd,sss,hon,dhn,pp,znp,zgc,zgs,crn)

c***************************************************************************
c     
c     DL_POLY subroutine for multiple time step algorithm
c     reciprocal space calculated on long time steps.
c     
c     copyright daresbury laboratory
c     
c     author  t. forester,  may 1993
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ Ewald sum                         : ewald1,2,3,4
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ truncated and shifted coulombic   : coul4
c     = 10,11 ----- reaction field                    : coul3
c     = 12,13 ----- Smoothed Particle Mesh Ewald      : ewald[_spme,2,3,4]
c     = 14,15 ----- Hautman-Klein-Ewald               : hkewald1,2,3,4
c     
c     wl
c     2002/05/31 14:01:12
c     1.10
c     Exp
c     
c****************************************************************************
      
#include "dl_params.inc"
      
      logical newplst,newlst,lgofr,lzeql,lgr,loglnk,lewald,lspme
      logical lhke
      complex*16 qqq,ww1,ww2,ww3,bscx,bscy,bscz
      
      dimension ahk(0:mxhko),zzn(mxxdf),zzd(mxxdf),sss(mxxdf)
      dimension hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko)
      dimension pp(2*mxhko),znp(mxhke,0:2*mxhko)
      dimension zgc(0:2*mxhko),zgs(0:2*mxhko),crn(0:mxhko,0:mxhko)
      dimension key1(kmaxd),key2(kmaxe),key3(kmaxf)
      dimension lentry(msatms),list(msatms,mxlist),ilist(mxxdf)
      dimension ltype(mxatms),ltpvdw(mxvdw),jlist(mxxdf)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension lstvdw(mxvdw),nauxfft(4)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension chge(mxatms),ewlbuf(mxebuf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension flx(mxatms),fly(mxatms),flz(mxatms)
      dimension fpx(mxatms),fpy(mxatms),fpz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension buffer(mxbuff)
      dimension ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld)
      dimension elc(mxewld,0:1),els(mxewld,0:1)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc)
      dimension erc(mxegrd),fer(mxegrd)
      dimension rdf(mxrdf,mxvdw)
      dimension stress(9),stresl(9),stresp(9)
      dimension txx(mxatms), tyy(mxatms), tzz(mxatms)
      dimension csp(mxspl),ffttable(mxftab)
      dimension ww1(kmaxd),ww2(kmaxe),ww3(kmaxf)
      dimension bspx(mxspme,mxspl),bspy(mxspme,mxspl)
      dimension bspz(mxspme,mxspl),bsdx(mxspme,mxspl)
      dimension bsdy(mxspme,mxspl),bsdz(mxspme,mxspl)
      dimension qqc(kmaxd,kmaxe,kmaxf)
      dimension qqq(kmaxd,kmaxe,kmaxf)
      dimension bscx(kmaxd),bscy(kmaxe),bscz(kmaxf)
#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif

      save engcpl,engsrl,vircpl,virsrl,nstep0,numlsts,engcp1,vircp1,
     x  engsr1,virsr1,stresl,stresp      
      data numlsts/-1/
#ifdef VAMPIR
      call VTBEGIN(16, ierr)
#endif
      lhke=(keyfce/2.eq.7)
      lspme=(keyfce/2.eq.6)
      lewald=(keyfce/2.eq.1)
c     
c     set up primary and secondary neighbour lists if needed
      
      if(newlst) nstep0 = nstep
      
      newplst =(newlst).or.(mod(nstep-nstep0,multt).eq.0)
      
      if (newplst) then        
        
        numlsts = numlsts + 1

        call  primlst
     x    (idnode,mxnode,natms,imcon,rprim,lentry,list,
     x    cell,xxx,yyy,zzz,xdf,ydf,zdf)
        
      endif
      
      if(newplst.or.(mod(nstep-nstep0,multt).le.1)) then
c     
c     zero accumulators for secondary neighbour energies and virial
        
        engcpl = 0.d0
        vircpl = 0.d0
        engsrl = 0.d0
        virsrl = 0.d0
        
c     
c     zero secondary forces
        
        do i = 1,natms
          
          flx(i) = 0.d0
          fly(i) = 0.d0
          flz(i) = 0.d0
          
        enddo

        do i = 1,9
          stresl(i) = 0.d0
        enddo
        
      endif
c     
c     zero accumulators for total energies and virials
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
c     
c     zero force arrays
      
      do i = 1,natms
        
        fxx(i) = 0.d0
        fyy(i) = 0.d0
        fzz(i) = 0.d0
        
      enddo
c     
c     flag for accumulating rdfs
      
      lgr =.false.
      if(nstbgr.gt.0)lgr=(mod(numlsts,nstbgr).eq.0)
      lgr = (lgr.and.(newplst.and.lgofr))
      lgr = (lgr.and.((.not.lzeql).or.(nstep-nsteql.gt.0)))
      
      
      if(newplst.or.(mod(nstep-nstep0,multt).le.1)) then
c     
c     use  SECONDARY NEIGHBOURS*******************************************
        
c     
c     fourier contribution to coulombic forces
        
        if(lewald.or.lspme.or.lhke)then

          if(lewald)then
            
            call ewald1
     x        (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x        engac1,viracc,alpha,volm,epsq,cell,chge,xxx,
     x        yyy,zzz,flx,fly,flz,elc,emc,enc,els,ems,ens,
     x        ckc,cks,clm,slm,stresl,buffer,ewlbuf)
c     
c     hautman-klein-ewald method
            
          elseif(lhke)then
            
            call hkewald1
     x        (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x        engac1,viracc,alpha,epsq,cell,ahk,chge,xxx,
     x        yyy,zzz,flx,fly,flz,elc,emc,els,ems,ckc,cks,
     x        stresl,crn,pp,znp,zgc,zgs,buffer)
c
c     real space terms of hk-ewald

            call hkewald2
     x        (idnode,mxnode,nhko,nlatt,imcon,natms,engac2,
     x        virac2,drewd,rcut,epsq,cell,chge,ahk,zzn,zzd,
     x        hon,dhn,xxx,yyy,zzz,flx,fly,flz,stresl,xdf,ydf,
     x        zdf,sss,rsqdf)

            engac1=engac1+engac2
            viracc=viracc+virac2
            
          elseif(lspme) then
c
c     smoothed particle mesh ewald
            
            call ewald_spme
     x        (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,fplan,
     x         bplan,engac1,viracc,alpha,volm,epsq,nauxfft,
     x         key1,key2,key3,cell,chge,
     x         xxx,yyy,zzz,txx,tyy,tzz,flx,fly,flz,stresl,buffer,
     x         csp,ww1,ww2,ww3,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqc,
     x         qqq,bscx,bscy,bscz,ffttable)

          endif

          engcpl=engcpl+engac1
          vircpl=vircpl+viracc
c     
c     calculate corrections for intramolecular coulomb terms in 
c     Ewald sum
c     
c     note: if using link cells - have double counted excluded 
c     interactions use temporary adjustment to relative dielectric
c     constant

          eps = epsq
          if(loglnk) eps = epsq*2.0d0
          
c     
c     outer loop over atoms
          
          ii=0
          
          do i=idnode+1,natms,mxnode
            
            ii=ii+1
            
c     
c     calculate interatomic distances
            
            do k=1,nexatm(ii)
              
              j=lexatm(ii,k)
              jlist(k)=j
              
              xdf(k)=xxx(i)-xxx(j)
              ydf(k)=yyy(i)-yyy(j)
              zdf(k)=zzz(i)-zzz(j)
              
            enddo
            
c     
c     periodic boundary condition
            
            call images(imcon,0,1,nexatm(ii),cell,xdf,ydf,zdf)
            
c     
c     calculate correction terms
            
            if(lhke)then

              call hkewald3
     x          (i,ii,engacc,viracc,eps,nexatm,lexatm,chge,
     x          xdf,ydf,zdf,flx,fly,flz,stresl)
              
            else

              call ewald3
     x          (i,ii,engacc,viracc,alpha,eps,nexatm,
     x          lexatm,chge,xdf,ydf,zdf,flx,fly,flz,stresl)
              
            endif

            engac1=engac1+engacc
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          enddo
          
        endif

c     
c     outer loop over atoms  **** OUTER LOOP STARTS HERE *****
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     
c     calculate interatomic distances
          
          ik = 0
          
          do k=1,lentry(ii)
            
            j=list(ii,k)
            
            if(j.gt.0) then
              
              ik = ik + 1
              ilist(ik) = j
              xdf(ik)=xxx(i)-xxx(j)
              ydf(ik)=yyy(i)-yyy(j)
              zdf(ik)=zzz(i)-zzz(j)
              
            endif
            
          enddo
c     
c     periodic boundary conditions
          
          call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
c     
c     square of distance

          do k = 1,ik

            rsqdf(k) = xdf(k)**2 + ydf(k)**2 + zdf(k)**2

          enddo
c     
c     accumulate radial distribution functions
          
          if(lgr) call rdf0(i,ik,rcut,ilist,ltype,lstvdw,rdf,rsqdf)
          
c     
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1) then
            
            call srfrce
     x        (i,ik,engacc,viracc,rvdw,dlrpot,ilist,ltype,lstvdw,
     x        ltpvdw,rsqdf,xdf,ydf,zdf,flx,fly,flz,vvv,ggg,stresl)
            
            engsrl=engsrl+engacc
            virsrl=virsrl+viracc
            
          endif
c     
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)

          if (lewald.or.lspme)then

            call ewald2
     x        (i,ik,engacc,viracc,drewd,rcut,epsq,ilist,
     x        chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,erc,fer,stresl)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif (keyfce/2.eq.2) then
            
            call coul2(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif (keyfce/2.eq.3) then
            
            call coul0(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.4) then
            
            call coul4(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc

          elseif(keyfce/2.eq.5) then
            
            call coul3 (i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,flx,fly,flz,stresl)  
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          endif
          
        enddo
        
      endif

c     
c     use PRIMARY NEIGHBOURS 
      
c     
c     outer loop over atoms  **** OUTER LOOP STARTS HERE *****
      
      ii=0
      
      do i=idnode+1,natms,mxnode
        
        ii=ii+1
        
c     
c     calculate interatomic distances
        
        ik = 0
        
        do k=1,lentry(ii)
          
          j= -list(ii,k)
          
          if(j.gt.0) then
            
            ik = ik + 1
            ilist(ik) = j
            xdf(ik)=xxx(i)-xxx(j)
            ydf(ik)=yyy(i)-yyy(j)
            zdf(ik)=zzz(i)-zzz(j)
            
          endif
          
        enddo
c     
c     periodic boundary conditions
        
        call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
c     
c     square of distance

        do k = 1,ik

          rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2

        enddo
c     
c     accumulate radial distribution functions
        
        if (lgr) call rdf0(i,ik,rcut,ilist,ltype,lstvdw,rdf,rsqdf)
c     
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1) then
          
          call srfrce
     x      (i,ik,engacc,viracc,rvdw,dlrpot,ilist,ltype,lstvdw,
     x      ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
c     
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)
        
        if (lewald.or.lspme.or.lhke) then
          
          if(newplst.or.
     x      (mod(nstep-nstep0,multt).le.1)) then
            
            if(lhke)then

              call hkewald4
     x          (i,ik,engacc,viracc,engacl,viracl,
     x          rcut,epsq,ilist,chge,rsqdf,xdf,ydf,zdf,
     x          fxx,fyy,fzz,flx,fly,flz,stress,stresl)

            else

              call ewald4
     x          (i,ik,engacc,viracc,engacl,viracl,drewd,
     x          rcut,epsq,ilist,chge,rsqdf,xdf,ydf,zdf,
     x          fxx,fyy,fzz,flx,fly,flz,erc,fer,stress,stresl)

            endif
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            engcpl=engcpl+engacl
            vircpl=vircpl+viracl
            
          else
            
            call coul0(i,ik,engacc,viracc,rcut,epsq,
     x        ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
              
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          endif
          
        elseif (keyfce/2.eq.2) then
          
          call coul2(i,ik,engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif (keyfce/2.eq.3) then
          
          call coul0(i,ik,engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4) then
          
          call coul4(i,ik,engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc

        elseif(keyfce/2.eq.5) then
          
          call coul3 (i,ik,engacc,viracc,rcut,epsq,
     x      ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
      enddo
c     
c     END of PRIMARY neighbour loop
      
c     
c     counter for rdf statistics outside loop structure
      
      if(lgr) numrdf = numrdf+1
      
c     
c     add in secondary neighbour contributions to force, energy etc
      
      if(newplst) then
        
        do i = 1,natms
          
          fpx(i) = flx(i)
          fpy(i) = fly(i)
          fpz(i) = flz(i)
          
        enddo
        
        do i = 1,9
          stresp(i) = stresl(i)
        enddo
        
        engsr1 = engsrl
        virsr1 = virsrl
        engcp1 = engcpl
        vircp1 = vircpl

      endif
      
      if(mod(nstep-nstep0,multt).eq.1) then
        
        do i = 1,natms
          
          flx(i) = flx(i) - fpx(i)
          fly(i) = fly(i) - fpy(i)
          flz(i) = flz(i) - fpz(i)
          
        enddo
        
        do i = 1,9
          stresl(i) = stresl(i) - stresp(i)
        enddo

        virsrl = virsrl - virsr1
        engsrl = engsrl - engsr1
        vircpl = vircpl - vircp1
        engcpl = engcpl - engcp1
        
      endif
      
      ann = dble(mod(nstep-nstep0,multt))

      do i = 1,natms
        
        fxx(i) =  fpx(i) + flx(i)*ann + fxx(i) 
        fyy(i) =  fpy(i) + fly(i)*ann + fyy(i) 
        fzz(i) =  fpz(i) + flz(i)*ann + fzz(i) 
        
      enddo
      
      do i = 1,9
        stress(i) = stress(i) + stresp(i) + stresl(i)*ann
      enddo

      engsrp = engsr1 + engsrl*ann + engsrp
      virsrp = virsr1 + virsrl*ann + virsrp
      engcpe = engcp1 + engcpl*ann + engcpe
      vircpe = vircp1 + vircpl*ann + vircpe
c     
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1) then
        
        buffer(5)=engsrp
        buffer(6)=virsrp
        buffer(7)=engcpe
        buffer(8)=vircpe
        
        call gdsum(buffer(5),4,buffer(1))
        
        engsrp=buffer(5)
        virsrp=buffer(6)
        engcpe=buffer(7)
        vircpe=buffer(8)
        
      endif
      
#ifdef VAMPIR
      call VTEND(16, ierr)
#endif

      return
      end
