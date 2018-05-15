      subroutine result
     x  (lgofr,lpgr,lzden,cfgname,atmnam,unqatm,idnode,imcon,
     x  keyens,mxnode,natms,nzden,nstep,ntpatm,ntpvdw,numacc,
     x  numrdf,chip,chit,conint,rcut,tstep,volm,lstvdw,buffer,
     x  cell,dens,eta,fxx,fyy,fzz,ravval,rdf,ssqval,stkval,
     x  stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,zzz,xx0,yy0,
     x  zz0,zdens,xxs,yys,zzs,dipx,dipy,dipz,lpolar,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,chitd,conintd)
c     
c***********************************************************************
c     
c     dl_poly subroutine for writing simulation summary and
c     saving the restart data
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c     wl
c     2000/01/18 14:05:54
c     1.5
c     Exp
c
c***********************************************************************
c     
      
#include "dl_params.inc"

      logical lgofr,lpgr,lzden,check,lpolar,lcp,ldpts
      
      character*80 cfgname
      character*8  atmnam(mxatms),unqatm(mxsite)

      dimension lstvdw(mxvdw)
      dimension rdf(mxrdf,mxvdw),eta(9)
      dimension zdens(mxrdf,mxsvdw)
      dimension cell(9),dens(mxsvdw)
      dimension stpval(mxnstk),sumval(mxnstk),ssqval(mxnstk)
      dimension zumval(mxnstk),ravval(mxnstk),stkval(mxstak,mxnstk)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)
      dimension buffer(mxbuff)
      dimension xxs(mxatms),yys(mxatms),zzs(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension vdxx(mxatms),vdyy(mxatms),vdzz(mxatms)
      dimension xmsd(mxatms),ymsd(mxatms),zmsd(mxatms)
      
#ifdef VAMPIR
      call VTBEGIN(69, ierr)
#endif
c     
c     save restart data

      call revive
     x  (lgofr,lzden,lpolar,cfgname,atmnam,idnode,imcon,mxnode,
     x  natms,nstep,nzden,numacc,numrdf,chip,chit,chitd,conint,
     x  tstep,buffer,cell,fxx,fyy,fzz,ravval,rdf,ssqval,
     x  stkval,stpval,sumval,vxx,vyy,vzz,xxx,yyy,zumval,
     x  zzz,xx0,yy0,zz0,zdens,xxs,yys,zzs,eta,dipx,dipy,dipz,
     x  vdxx,vdyy,vdzz,xmsd,ymsd,zmsd,lcp,ldpts,conintd)
c
c     write final instantaneous induced dipoles

c        open(68)
c        do i=1,mxatms
c        write(68,*)dipx(i),dipy(i),dipz(i)
c        enddo
c        close(68)

c     
c     calculate final fluctuations

      do i=1,mxnstk

        ssqval(i)=sqrt(ssqval(i))

      enddo

c     
c     final averages and fluctuations
      
      call timchk(0,timelp)
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,/,1x,'run terminated after',i8,' steps.',
     x    ' final averages calculated over',i8,' steps.',/,/)") 
     x    nstep,numacc
        write(nrite,"(1x,130('-'),
     x    /,/,1x,'    step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',
     x    5x,'eng_vdw',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',5x,
     x    'eng_dih',5x,'eng_tet',/,1x,'time(ps)',5x,' eng_pv',4x,
     x    'temp_rot',5x,'vir_cfg',5x,'vir_vdw',5x,'vir_cou',5x,
     x    'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,
     x    1x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',
     x    5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',5x,'vir_pmf',
     x    7x,'press',/,/,
     x    1x,130('-'))")          

        write(nrite,'(1x,i8,1p,9e12.4,/,1x,0p,f8.3,1p,9e12.4,
     x    /,1x,0p,f8.2,1p,9e12.4)') 
     x    nstep,(sumval(i),i=1,9),
     x    dble(nstep)*tstep,(sumval(i),i=10,18),
     x    timelp,(sumval(i),i=19,27)
        write(nrite,"(/,1x,' r.m.s. ',1p,9e12.4,/,1x,'fluctn. ',
     x    1p,9e12.4,/,9x,9e12.4)") (ssqval(i),i=1,27)
        write(nrite,"(1x,130('-'))")

        if(numacc.gt.0) then
          iadd = 27
c     
c     write out estimated diffusion coefficients

          if(numacc.gt.0) then
            write(nrite,"(/,/,12x,'Approximate 3D Diffusion',
     x       '  coefficients (10^-9 m^2 / s)',/,/,12x,'atom',7x,' D ')")

            do i = 1,ntpatm
              iadd = iadd+1
              dc = (ravval(iadd)-sumval(iadd))/
     x          (3.d0*dble(numacc-min(mxnstk,numacc-1))*tstep)*10.d0
c              if(dc.lt.1d-10) dc = 0.d0
              write(nrite,'(12x,a8,1p,e13.4)') unqatm(i),dc
            enddo
          endif
#ifdef STRESS
c     
c     print out average pressure tensor

          write(nrite,"(/,/,16x,'Average pressure tensor',
     x      39x,'  m.s. fluctuations ',/)")

          do i = iadd,iadd+6,3
            write(nrite,'(9x,1p,3e12.4,24x,3e12.4)')
     x        (sumval(i+j),j = 1,3),(ssqval(i+j)**2,j = 1,3)
          enddo
          iadd = iadd+9

          write(nrite,'(/,12x,a,1p,e12.4)') 'trace/3. ',
     x     (sumval(iadd)+sumval(iadd-4)+sumval(iadd-8))/3.d0

          write(nrite,'(/,9x,a,24x,a)')
     x         '    average surface tension (mN / m)',
     x         'fluctuation of surface tension (mN / m)'
          surtendum= cell(9)*(2.d0*sumval(iadd)-sumval(iadd-4)
     x        -sumval(iadd-8))/4.d0*10.d0
          std1=ssqval(iadd)/dsqrt(dble(nstep))
          std2=ssqval(iadd-4)/dsqrt(dble(nstep))
          std3=ssqval(iadd-8)/dsqrt(dble(nstep))
          surtenflc= cell(9)*(2.d0*std1+std2+std3)/4.d0*10.d0
          write(nrite,'(/,21x,e12.4,48x,e12.4,/,/)')
     x        surtendum,surtenflc
#endif
c     
c     write out mean cell vectors for npt 

          if(keyens.gt.3.and.(keyens.le.7)) then

            write(nrite,"(/,/,17x,'Average cell vectors',
     x        41x,'r.m.s. fluctuations ',/)")

            do i = iadd,iadd+6,3
              write(nrite,'(3f20.10,9x,1p,3e12.4)')
     x          (sumval(i+j),j = 1,3),(ssqval(i+j),j = 1,3)
            enddo
            iadd = iadd+9
            
          endif
c     
c     write out remaining registers

          check = .false.
          do i = iadd+1,mxnstk

            if((sumval(i).ne.0.d0).or.(ssqval(i).ne.0.d0)) check= .true.          

          enddo

          if(check) then

            write(nrite,"(/,/,12x,
     x        'Remaining non-zero statistics registers ',
     x        /,/,12x,'Register',7x,'Average value',8x,'r.m.s. fluc.')")

            do i = iadd+1,mxnstk

              if((sumval(i).ne.0.d0).or.(ssqval(i).ne.0.d0))
     x          write(nrite,'(10x,i10,2f20.10)') i,sumval(i),ssqval(i)

            enddo

          endif

        endif
c     
c     print out sample of final configuration 
        
        write(nrite,"(/,/,1x,'sample of final configuration',/)")
        write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x    7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',7x,'fx(i)',7x,
     x    'fy(i)',7x,'fz(i)',/,/)")
        io=(natms+19)/20
        
        do i=1,natms,io
          
          write(nrite,"(1x,i6,1p,3e12.4,3e12.4,3e12.4)") 
     x      i,xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),
     x      fxx(i),fyy(i),fzz(i)
          
        enddo

      endif
c     
c     average volume

      avvol = sumval(19)
      if(imcon.eq.0.or.imcon.eq.6) then
        avvol = 4.d0*pi/3.d0*rcut**3
        volm = avvol
      endif
c     
c     calculate and print radial distribution functions

      if(lgofr.and.lpgr.and.(numrdf.gt.0)) then 
c     
c     scale densities for average volume

        do i = 1,ntpatm
          dens(i) = dens(i)*(volm/avvol)
        enddo

        call  rdf1
     x    (lpgr,cfgname,unqatm,idnode,mxnode,ntpatm,ntpvdw,numrdf,
     x    rcut,avvol,lstvdw,dens,rdf,buffer)

      endif

      if(lzden.and.lpgr.and.(nzden.gt.0)) then 

        zlen=abs(cell(3))+abs(cell(6))+abs(cell(9))

        call  zden1
     x    (lpgr,cfgname,unqatm,idnode,mxnode,ntpatm,nzden,
     x    avvol,zlen,zdens,buffer)

      endif

      if (imcon.eq.0) volm = 0.d0
c     
c     print final time check
      
      call timchk(1,timelp)
      
#ifdef VAMPIR
      call VTEND(69, ierr)
#endif
      return
      end
