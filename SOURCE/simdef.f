! 07 JAN 11 - VB - FP's CHANGES FROM pimdpol
! 06 MAR 06 - IUCHI - FLAG LUNFORMAT TO CONTROL FORMAT OF HISTORY FILE
! 16 NOV 05 - IUCHI - CALL INIT_PS_TYPE_DMS IF LDMS
! 16 NOV 05 - IUCHI - INTRODUCING FLAGS LDMS AND LTTMFOLD
! 14 NOV 05 - IUCHI - INTRODUCING FLAG LDECOMP
! 09 NOV 05 - IUCHI - INTRODUCING FLAG LTIP4P_GEO
! 07 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
! 07 NOV 05 - IUCHI - FLAG LPOLOFF FOR POLARIZATION TURN OFF
! 27 OCT 05 - IUCHI - COMMENT OUT MINOR BUG, LSHAKE
! 20 OCT 05 - IUCHI - REPLACE 4.8 BY EATD FROM MODULE
!
      subroutine simdef
     x  (lfcap,lgofr,lnsq,loptim,lpgr,ltraj,ltscal,lzeql,lzden,lpolar,
     x  lthole,sysname,idnode,mxnode,intsta,istraj,keyens,keyfce,lopt2,
     x  lfd,keyres,keytrj,kmax1,kmax2,kmax3,multt,nstack,nstbgr,nstbpo,
     x  nhko,nlatt,nstbts,nsteql,nstraj,nstrun,nospl,fplan,bplan,ftol,
     x  alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup,taut,
     x  athole,athole12,athole13,temp,timcls,timjob,tolnce,tstep,
     x  ffttable,
     x  lcp,ldpts,lads,lacs,lttm,lttm3,nthole,dipmas,diptmp,tautd,toler,
     x  npartem,temgap,
     x  delta,ascc,ascd,gammattm,laspc,sor_omega,aspc_k_param,
     x  aspc_iter_nstep)
      use multibead, only: bead_suffix
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading in the simulation control 
c     parameters
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith july 1992.
c     
c     modified
c     author   - t. forester      may  1993
c     amended  - t.forester       sept 1994 - dl_poly_1.1
c     amended  - t.forester       nov  1994 - macro version
c     amended  - w.smith          dec  1994 - t3d adaptation
c     
c     wl
c     2001/08/31 11:13:52
c     1.13
c     Exp
c     
!     Last updated: 6 March 2006 by S. Iuchi
!
c***********************************************************************
c     
! from module
      use format_history,  only: lunformat
      use ps_type_dms,     only: ldms, lttmfold, init_ps_type_dms,
     x                           fac_r, fac_theta
      use ttm_forces,      only: lpoloff, ltip4p_geo, ldecomp, lpolintra
      use unit_parameters, only: eatd

#include "dl_params.inc"
#include "fftw3.f"

      parameter (mega = 1000000)
      
      character*80 sysname
      character*256 record,record1
      character*20 blank

#ifdef FFTW
      FFTW_PLAN_TYPE fplan,bplan
#else
      integer fplan, bplan
#endif
      
      logical ltscal,lzeql,loptim,ltraj,lfcap,lgofr,lpolar,lthole,
     x  lpgr,lpres,safe,lnsq,lzden,lewald,lspme,lhke,lcp,ldpts,lopt2,
     x  lpartem,lttm,lttm3,lfd,lacs,lads,laspc!,lshake
      logical lstep,ltemp,lcut,ldelr,lprim,lforc,lens,lvdw,lrvdw,kill
      
      dimension cell(9),celprp(10),ffttable(mxftab)
      real(8), intent(out) :: sor_omega
      integer, intent(out) :: aspc_k_param
      integer, intent(out) :: aspc_iter_nstep
      
      data blank/'                    '/
#ifdef VAMPIR
      call VTBEGIN(3, ierr)
#endif
c     
c     intitialize system variables: temperature,pressure,ensemble key
c     force key, cutoff, primary cutoff, verlet shell width, relative
c     dielectric constant,timestep,temperature scaling flag, 
c     temp scaling interval
      nhko = 0
      nlatt = 0
      nsteql = 0
      nstrun = 0
      keyres = 0
      keyens = 0
      npartem = 1
      temgap = 0.d0
      taut = 0.d0
      tautd = 0.1d0
      ltscal = .false.
      lewald=.false.
      lspme=.false.
      lhke=.false.
      nstbts = 0
      lgofr = .false.
      nstbgr = 1
      lpgr = .false.
      lzeql =.true.
      loptim = .false.
      lopt2 = .false.
      lpartem = .false.
      nstbpo = 100
      nstack = mxstak
      intsta = 100
      ltraj = .false.
      nstraj = 0
      istraj = 1
      keytrj = 0
      alpha = 0.d0
      kmax1 = 0
      kmax2 = 0
      kmax3 = 0
      nospl = min(8,mxspl)
      lfcap = .false.
      fmax = 1000.d0
      keyfce = 0
      multt = 1
      tstep = 0.d0
      temp = 0.d0
      press = 0.d0
      rcut = 0.d0
      rprim = 0.d0
      rvdw = 0.d0
      delr = 0.d0
      epsq = 1.d0
      tolnce = 1.d-8
      quattol= 1.d-8
      timjob = 0.d0
      timcls = 0.d0
      athole = 1.d5
      athole12 = athole
      athole13 = athole
      ascc = 1.d5
      ascd = 1.d5
      gammattm = 0.d0
      nthole = 0
      dipmas = 0.d0
      toler = 1.d-8/eatd
      ftol = 1.d-10
      
      ltemp = .false.
      lstep = .false.
      lcut  = .false.
      ldelr = .false.
      lprim = .false.
      lforc = .false.
      lens = .false.
      lvdw = .false.
      lrvdw =.false.
      lpres =.false.
      kill = .false.
      lnsq = .false.
      lzden = .false.

      lpolar=.false.
      lthole=.false.
      lacs=.false.
      lads=.false.
      lcp=.false.
      ldpts=.false.
      lttm=.false.
      lttm3=.false.
      lfd=.false.
      lpolintra=.false.
!      lshake=.false.
      lpoloff=.false.    ! polarization turns off
      ltip4p_geo=.false. ! TIP4P geometry
      ldecomp=.false.    ! force docomposition
      ldms=.false.       ! DMS 
      lttmfold=.true.    ! method of charge assignment
      lunformat=.false.  ! HISTORY file format
      diptmp=1.d-5
      delta=1.d-6
      laspc=.false.
      sor_omega = 0.3d0
      aspc_k_param = 4
      aspc_iter_nstep = 6

c     
c     open the simulation input file
      
      if(idnode.eq.0) open(nread,file='CONTROL',status='old')
c     
c     read job title
      
c$$$  read(nread,'(a80)',end=1000) sysname
      call getrec(safe,idnode,mxnode,nread,record)
      if(.not.safe)go to 1000
      sysname=record(1:80)
      if(idnode.eq.0) then 
        
        write(nrite,"(3(1x,130('*'),/),1x,20('*'),5x,a80,5x,20('*'),/,
     x    3(1x,130('*'),/),/,/,1x,'SIMULATION CONTROL PARAMETERS',/)") 
     x    sysname
        
      endif
c     
c     read and process directives from CONTROL file
      
      do nrecs = 1,mega
        
c$$$  read(nread,'(a80)',end=1000) record
        call getrec(safe,idnode,mxnode,nread,record)
        if(.not.safe)go to 1000
        
        record(81:100)=blank(1:20)
c     
c     convert to lowercase and strip out leading blanks
        
        call lowcase(record,80)
        call strip(record,80)
        
        if(record(1:1).eq.'#'.or.record(1:1).eq.' ') then
c     
c     record is commented out
          
        elseif(record(1:5).eq.'steps') then
c     
c     number of timesteps
          
          nstrun = intstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'selected number of timesteps',3x,i10)") nstrun
          
        elseif(record(1:5).eq.'equil') then
c     
c     number of equilibration timesteps
          
          nsteql = intstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'equilibration period        ',3x,i10)") nsteql
          
        elseif(record(1:7).eq.'restart') then
c     
c     restart control
          
          record(1:73) = record(8:80)
          call strip(record,73)
          
          if(record(1:5).eq.'scale') then
            
            keyres = 2
            if(idnode.eq.0) write(nrite,
     x        "(/,1x,'scaled restart requested')")
            
          else
            
            keyres = 1
            if(idnode.eq.0) write(nrite,"(/,1x,'restart requested')")
            
          endif
          
        elseif(record(1:8).eq.'ensemble') then
c     
c     ensemble
      
          record1 = record
          record(1:72) = record(9:80)
          call strip(record,72)
          
          if(record(1:3).eq.'nve') then
            
            if(idnode.eq.0) write(nrite,
     x        "(/,1x,'microcanonical ensemble')")
            if(lens) then
              call error(idnode,-414)
              kill = .true.
            endif
            lens=.true.
            
          elseif(record(1:3).eq.'nvt') then
            
            record(1:69) = record(4:72)
            call strip(record,69)
            
            if(record(1:5).eq.'evans') then
              
              keyens = 1
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Evans Gaussian temperature constraints',
     x          ' in use')") 
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            elseif(record(1:3).eq.'ber') then
              
              keyens = 2
              taut = dblstr(record,69,idum)
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Berendsen thermostat',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4)")
     x          taut
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            elseif(record(1:6).eq.'hoover') then
              
              keyens = 3
              taut = dblstr(record,69,idum)
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Nose-Hoover ',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4)")
     x          taut
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            else
              
              kill = .true.
              if(idnode.eq.0) write(nrite,"(/,/,a80)") record1
              call error(idnode,-3)
              
            endif
            
          elseif(record(1:3).eq.'npt') then
            
            record(1:69) = record(4:72)
            call strip(record,69)
            
            if(record(1:3).eq.'ber') then
              
              keyens = 4
              taut = dblstr(record,69,idum)
              ilen = 69-idum 
              record(1:ilen) = record(idum:69)
              taup = dblstr(record,ilen,idum)
              
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Berendsen isotropic N-P-T',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4,
     x          /,1x,'barostat relaxation time       ',1p,e12.4)")
     x          taut,taup
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            elseif(record(1:6).eq.'hoover') then
              
              keyens = 5
              taut = dblstr(record,69,idum)
              ilen = 69-idum 
              record(1:ilen) = record(idum:69)
              taup = dblstr(record,ilen,idum)
              
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Nose-Hoover  (Melchionna) isotropic N-P-T ',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4,
     x          /,1x,'barostat relaxation time       ',1p,e12.4)")
     x          taut,taup
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            else
              
              kill = .true.
              if(idnode.eq.0) write(nrite,"(/,/,a80)") record1
              call error(idnode,-3)
              
            endif
            
          elseif(record(1:3).eq.'nst') then
            
            record(1:69) = record(4:72)
            call strip(record,69)
            
            if(record(1:3).eq.'ber') then
              
              keyens = 6
              taut = dblstr(record,69,idum)
              ilen = 69-idum 
              record(1:ilen) = record(idum:69)
              taup = dblstr(record,ilen,idum)
              
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Berendsen anisotropic N-P-T',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4,
     x          /,1x,'barostat relaxation time       ',1p,e12.4)")
     x          taut,taup
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
            elseif(record(1:6).eq.'hoover') then
              
              keyens = 7
              taut = dblstr(record,69,idum)
              ilen = 69-idum 
              record(1:ilen) = record(idum:69)
              taup = dblstr(record,ilen,idum)
              
              if(idnode.eq.0) write(nrite,
     x          "(/,1x,'Nose-Hoover (Melchionna) anisotropic N-P-T ',
     x          /,1x,'thermostat relaxation time     ',1p,e12.4,
     x          /,1x,'barostat relaxation time       ',1p,e12.4)")
     x          taut,taup
              if(lens) then
                call error(idnode,-414)
                kill = .true.
              endif
              lens=.true.
              
              
            else
              
              kill = .true.
              if(idnode.eq.0) write(nrite,"(/,/,a80)") record1
              call error(idnode,-3)
              
            endif

          elseif(record(1:3).eq.'pmf') then
            
            keyens = 8
            if(idnode.eq.0) write(nrite,
     x        "(/,1x,'potential of mean force calulation (NVE)')")
            if(lens) then
              call error(idnode,-414)
              kill = .true.
            endif
            lens=.true.
  
          else
          
            call error(idnode,-436)
            kill=.true.

          endif
          
        elseif(record(1:5).eq.'scale') then
          
          nstbts = intstr(record,80,idum)
          if(nstbts.gt.0) then
            ltscal =.true.
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'temperature scaling on' 
     x        /,1x,'temperature scaling interval',3x,i10)")
     x        nstbts
            
          endif
          
        elseif(record(1:3).eq.'rdf') then
          
          nstbgr = intstr(record,80,idum)
          
          if(nstbgr.gt.0) then
            
            lgofr =.true.
            if(idnode.eq.0)write(nrite,
     x        "(/,/,1x,'radial distribution functions on ',
     x        /,1x,'g(r) collection interval    ',3x,i10)")
     x        nstbgr
            
          endif
          
        elseif(record(1:4).eq.'zden') then
          
          lzden = .true.
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'Z density profile requested')")
          
        elseif(record(1:9).eq.'print rdf') then
          
          lpgr =.true.
          lpgr =(lgofr.and.lpgr)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'g(r) printing option         ',3x,l10)")lpgr
          
        elseif(record(1:7).eq.'collect') then
          
          lzeql =.false.
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'equilibration included in overall averages')")
          
        elseif(record(1:4).eq.'zero') then
          
          loptim =.true.
          ltemp = .true.
          ltscal = .false.
          if(temp.eq.0.d0) temp=1.d0
          
          if(idnode.eq.0) then
            write(nrite,
     x      "(/,1x,'zero K optimisation requested:')")
            if(temp.gt.10.d0) write(nrite,
     x      "(' temperature reset',1p,e12.4)") 10.d0
          endif

          temp = max(temp,10.d0)

        elseif(record(1:3).eq.'opt') then

          lopt2 = .true.
          ftol = dblstr(record,80,idum)
          write(nrite,*)'geometry optimization on basin hopping'
          write(nrite,*)'convergence criterion (energy in kcal/mol) = ', 
     x                  ftol

        elseif(record(1:6).eq.'partem') then

            lpartem = .true.
            npartem = intstr(record,80,idum)
            ilen = 80-idum
            record(1:ilen) = record(idum:80)
            temgap=dblstr(record,ilen,idum)

        elseif(record(1:5).eq.'print') then
          
          nstbpo = intstr(record,80,idum)
          nstbpo = max(nstbpo,1)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'data printing interval      ',3x,i10)") nstbpo
          
        elseif(record(1:5).eq.'stack') then
          
          nstack = intstr(record,80,idum)
c     
c     reset stack limit if too large
          
          nstack=min(nstack,mxstak)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'data stacking interval      ',3x,i10)") nstack
          
        elseif(record(1:5).eq.'stats') then
          
          intsta = intstr(record,80,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'statistics file interval    ',3x,i10)") intsta
          
        elseif(record(1:4).eq.'traj') then
          
          ltraj = .true.
          nstraj = intstr(record,80,idum)
          ilen = 80 - idum
          record(1:ilen) = record(idum:80)
          istraj = intstr(record,ilen,idum)
          istraj = max(istraj,1)
          jlen = ilen - idum
          record(1:jlen) = record(idum:ilen)
          keytrj = intstr(record,jlen,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'trajectory file option on  ',
     x      /,1x,'trajectory file start       ',3x,i10,
     x      /,1x,'trajectory file interval    ',3x,i10
     x      /,1x,'trajectory file info key    ',3x,i10)")
     x      nstraj,istraj,keytrj

        elseif(record(1:5).eq.'polar') then

          lpolar=.true.
          lpolintra=.true. ! default behavior is to include all
                           ! non-excluded charges in polarizability
                           ! calculation

        elseif(record(1:4).eq.'iter') then

          toler = dblstr(record,80,idum)

          if (idnode.eq.0) then
            write(nrite,"(/,1X,'polarizability based on iteration')")
            write(nrite,*)'convergence criterion (debye) = ',toler
          endif
c
c     convert from Debye to e.A

!          toler = toler/4.8d0
          toler = toler / eatd
        
        elseif(record(1:10).eq.'aspc_omega') then

          sor_omega = dblstr(record,80,idum)
          if(sor_omega.gt.1.d0 .or. sor_omega.lt.0.d0) then
              call error(idnode,-816)
          endif


          if (idnode.eq.0) then
            write(nrite,*) 'ASPC omega (SOR param) = ', sor_omega
          endif

        elseif(record(1:6).eq.'aspc_k') then

          aspc_k_param = intstr(record,80,idum)
          if(aspc_k_param.lt.0) then
              call error(idnode,-819)
          endif

          ! set the recommended sor_omega parameter
          !     See Kolafa, dx.doi.org/10.1002/jcc.10385
          if(aspc_k_param .eq. 2) then
              sor_omega = 4.0d0/7.0d0
              aspc_iter_nstep = aspc_k_param + 2 ! 4
          elseif(aspc_k_param .eq. 4) then
              sor_omega = 6.0d0/11.0d0
              aspc_iter_nstep = aspc_k_param + 2 ! 6
          endif


          if (idnode.eq.0) then
            write(nrite,*) 'ASPC prediction k val = ', aspc_k_param
          endif

       elseif(record(1:4).eq.'aspc') then

          laspc = .true.

          if (idnode.eq.0) then
            write(nrite,*) 'polarizability using ASPC  '
            write(nrite,*) ' (Always Stable Predictor Corrector) method'
          endif

          ! default parameter assumes prediction k == 4
          sor_omega = 6.0d0/11.0d0
          aspc_k_param = 4
          aspc_iter_nstep = 6 ! does iteration for first 6 step
                              ! to accumulate dipoles used for prediction

        elseif(record(1:3).eq.'car') then

          dipmas = dblstr(record,80,idum)
          lcp=.true.
          if (idnode.eq.0)
     x    write(nrite,"(/,1X,'polarizability based on Car-Parrinello
     x    '/,1X,'fictious mass of dipole   ',3x,1p,e12.4)")dipmas

        elseif(record(1:5).eq.'dtemp') then

          diptmp = dblstr(record,80,idum)
          ldpts = .true.
          if (idnode.eq.0)
     x    write(nrite,"(/,1X,'dipole couple to hoover thermostat 
     x    '/,1X,'dipole temperature   ',3x,1p,e12.4)")diptmp

        elseif(record(1:7).eq.'dhoover') then

          tautd = dblstr(record,80,idum)
          ldpts = .true.
          if (idnode.eq.0)
     x    write(nrite,"(/,1X,'thermostat relaxation time',
     x    3x,1p,e12.4)")tautd

        elseif(record(1:5).eq.'thole') then

          lthole=.true.
          if(idnode.eq.0) write(nrite,*) 'Thole scheme is used'

        elseif(record(1:5).eq.'exp-1') then

          nthole = 1
          if (idnode.eq.0)
     x    write(nrite,*)'use Thole exponential decay function'

        elseif(record(1:5).eq.'exp-3') then

          nthole = 3
          if (idnode.eq.0)
     x    write(nrite,*)'use Thole exponential-3 decay function'

        elseif(record(1:18).eq.'all charge smeared') then

          lacs=.true.

          if (idnode.eq.0)
     x    write(nrite,*)'all charges are smeared!'

        elseif(record(1:18).eq.'all dipole smeared') then

          lads=.true.
          if (idnode.eq.0)
     x    write(nrite,*)'all dipoles are smeared!'

        elseif(record(1:3).eq.'ttm') then

          lttm=.true.
          if (idnode.eq.0) then
!            if(.not.lshake)write(nrite,*)'TTM2-R water model is used!'
!            if(lshake)write(nrite,*)'TTM2-F water model is used!'
             write(nrite,*) 'TTM water model is used:'
          endif

       elseif(record(1:5).eq.'gamma') then

          gammattm = dblstr(record,80,idum)

       elseif(record(1:12).eq.'atholeonetwo') then

          athole12 = dblstr(record,80,idum)

       elseif(record(1:14).eq.'atholeonethree') then

          athole13 = dblstr(record,80,idum)

       elseif(record(1:6).eq.'athole') then

          athole = dblstr(record,80,idum)

       elseif(record(1:4).eq.'ascc') then

          ascc = dblstr(record,80,idum)

       elseif(record(1:4).eq.'ascd') then

          ascd = dblstr(record,80,idum)

       elseif(record(1:5).eq.'fac_r') then

          fac_r = dblstr(record,80,idum)

       elseif(record(1:6).eq.'fac_th') then

          fac_theta = dblstr(record,80,idum)

       elseif(record(1:5).eq.'vttm2') then

          lthole = .true.     !thole
          nthole = 3          !exp-3
          lacs = .true.       !all charge smeared
          lads = .true.       !all dipole smeared
          ldms = .true.       !dms
          lttm = .true.       !ttm
          lttm3 = .false.     !ttm3-f charge assignment
          lttmfold = .false.  !false -> newdms

          call init_ps_type_dms

          gammattm = 0.4267d0
          athole = 0.3d0
          athole_bonded = 0.3d0
          ascc = 0.2d0
          ascd = 0.2d0

          if (idnode.eq.0) then
             write(nrite,*) 'TTM water model is used:'
             write(nrite,*) 'TTM2 model '
             write(nrite,*) 'Thole scheme is used'
             write(nrite,*) 'using Thole exponential-3 decay function'
             write(nrite,
     x  "(/,1x,'Thole attenuating factor    ',3x,1p,e12.4)")athole
             write(nrite,*) 'all charges are smeared!'
             write(nrite,
     x  "(1x,'TTM2 charge-charge smearing factor ',3x,1p,e12.4)")ascc
             write(nrite,*) 'all dipoles are smeared!'
             write(nrite,
     x  "(1x,'TTM2 charge-dipole smearing factor ',3x,1p,e12.4)")ascd
             write(nrite,
     x            "(/,1x,'gamma factor ',3x,1p,e12.4)")gammattm
             write(nrite,*) 'Partridge and Schwenke DMS is used'
             write(nrite,*) 'New charge assignment is used'
          endif

       elseif(record(1:5).eq.'vttm3') then

          lthole = .true.       !thole
          nthole = 3            !exp-3 
          lacs = .true.         !all charge smeared
          lads = .true.         !all dipole smeared
          ldms = .true.         !dms
          lttm = .true.         !ttm
          lttmfold = .false.    !false -> newdms
          
          lttm3 = .true.        !ttm3f ad-hoc modification of dms
          fac_r = 0.5
          fac_theta = 0.012
          
          call init_ps_type_dms
          
          gammattm = 0.4600d0
          athole = 0.175d0
          athole_bonded = 0.175d0
          ascc = 0.175d0
          ascd = 0.175d0
          
          if (idnode.eq.0) then
             write(nrite,*) 'TTM water model is used:'
             write(nrite,*) 'TTM3 model '
             write(nrite,*) 'Thole scheme is used'
             write(nrite,*) 'using Thole exponential-3 decay function'
             write(nrite,
     x  "(/,1x,'Thole attenuating factor    ',3x,1p,e12.4)")athole
             write(nrite,*) 'all charges are smeared!'
             write(nrite,
     x  "(1x,'TTM3 charge-charge smearing factor ',3x,1p,e12.4)")ascc
             write(nrite,*) 'all dipoles are smeared!'
             write(nrite,
     x  "(1x,'TTM3 charge-dipole smearing factor ',3x,1p,e12.4)")ascd
             write(nrite,
     x            "(/,1x,'gamma factor ',3x,1p,e12.4)")gammattm
             write(nrite,*) 'Partridge and Schwenke DMS is used'
             write(nrite,*) 'New charge assignment is used'
             write(nrite,*) 'Using the TTM3-F modification of the dms'
             write(nrite,
     x            "(/,1x,'bond stretching adjustment    ',3x,1p,e12.4)")
     x            fac_r
             write(nrite,
     x            "(/,1x,'bending adjustment    ',3x,1p,e12.4)")
     x            fac_theta
          endif

       elseif(record(1:5).eq.'vttm4') then

          lacs = .true.
          lads = .true.
          ldms = .true.
          lttm = .true.
          lttm3 = .false.
          lttmfold = .false.

          lthole = .true.
          nthole = 4

          gammattm = 0.426706882d0

          athole = 0.055d0 ! non-bonded dipole-dipole
          athole12 = 0.626d0 ! bonded dipole-dipole
          athole13 = 0.626d0 ! bonded dipole-dipole

          ascc = 0.4d0
          ascd = 0.4d0

          call init_ps_type_dms

          if (idnode.eq.0) then
             write(nrite,*) 'TTM4 model'
             write(nrite,
     x  "(/,1x,'Thole attenuating factor (non-bonded) ',3x,1p,e12.4)")
     x       athole
             write(nrite,
     x  "(/,1x,'Thole attenuating factor (12) ',3x,1p,e12.4)")
     x       athole12
             write(nrite,
     x  "(/,1x,'Thole attenuating factor (13) ',3x,1p,e12.4)")
     x       athole13
             write(nrite,
     x  "(1x,'TTM4 charge-charge smearing factor ',3x,1p,e12.4)")ascc
             write(nrite,
     x  "(1x,'TTM4 charge-dipole smearing factor ',3x,1p,e12.4)")ascd
          endif

        elseif(record(1:8).eq.'polintra') then

          lpolintra=.true.

        elseif(record(1:6).eq.'poloff') then  ! polarization turn off

           lpoloff=.true.
           if(idnode.eq.0) then 
              write(nrite,'(/)')
              write(nrite,*) 'Turn off polarization part!'
           endif

        elseif(record(1:8).eq.'tip4pgeo') then  ! TIP4P geometry

           ltip4p_geo=.true.
           if(idnode.eq.0) then 
              write(nrite,'(/)')
              write(nrite,*) 'M site position: TIP4P'
           endif

        elseif(record(1:6).eq.'decomp') then ! output force decompositions

           ldecomp=.true.
           if(idnode.eq.0) then  
              write(nrite,'(/)')
              write(nrite,*) 'Decomposition forces are output'
           endif

        elseif(record(1:3).eq.'dms') then ! PS DMS

           ldms=.true.
           if(idnode.eq.0) then  
              write(nrite,'(/)')
              write(nrite,*) 'Partridge and Schwenke DMS is used'
           endif
           call init_ps_type_dms  

        elseif(record(1:7).eq.'newdms') then ! original charge assignment

           lttmfold=.false.
           if(idnode.eq.0) then  
              write(nrite,*) 'New charge assignment is used'
           endif

        elseif(record(1:6).eq.'finite') then

          lfd=.true.
          delta = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x  "(/,1x,'Finite difference test with displ.',3x,1p,e12.4)") delta

        elseif(record(1:8).eq.'unformat') then

          lunformat=.true.
          if(idnode.eq.0) then
             write(nrite,'(/)')
             write(nrite,*) 'HISTORY is written by unformatted style'
          endif

        elseif(record(1:5).eq.'ewald'.or. 
     x         record(1:4).eq.'spme'.or.
     x         record(1:3).eq.'hke') then
c     
c     read Ewald or HK-Ewald or SPM-Ewald sum parameters
          
          lhke=(record(1:3).eq.'hke')
          lspme=(record(1:4).eq.'spme')
          lewald=(record(1:5).eq.'ewald')
          if(lewald)keyfce=2
          if(lspme)keyfce=12
          if(lhke)keyfce=14
          if(idnode.eq.0) open(nconf,file='CONFIG'//bead_suffix)
          call getrec(safe,idnode,mxnode,nconf,record1)
          call getrec(safe,idnode,mxnode,nconf,record1)
          imcon = intstr(record1(11:11),10,idum)
          if(.not.lhke.and.(imcon.eq.0.or.imcon.eq.6)) then

            call error(idnode,-180)
            kill=.true.

          endif
          
          if(record(7:15).eq.'precision' .or.
     x       record(6:14).eq.'precision' .or.
     x       record(5:13).eq.'precision') then
            
            record(1:87)=record(14:100)
            eps = dblstr(record(1:1),87,idum)
            if(idnode.eq.0) write(nrite,
     x        "(/,1x,'Ewald sum  precision    ',7x,1p,e12.4)") eps
            
            if(lhke) then

              ilen=87-idum 
              record(1:ilen)=record(idum:87)
              nhko=intstr(record(1:1),ilen,idum)
              nhko=min(nhko,3)
              if(nhko.eq.0)nhko=1
              jlen=ilen-idum
              record(1:jlen)=record(idum:ilen)
              nlatt=intstr(record(1:1),jlen,idum)
              nlatt=min(nlatt,2)
              if(nlatt.eq.0)nlatt=1
            
            endif

            if(.not.lcut) then
              call error(idnode,-433)
              kill=.true.
            else
c     
c     retreive cell vectors

              call getrec(safe,idnode,mxnode,nconf,record1)
              cell(1) = dblstr(record1(1:1),20,idum)
              cell(2) = dblstr(record1(21:21),20,idum)
              cell(3) = dblstr(record1(41:41),20,idum)
              call getrec(safe,idnode,mxnode,nconf,record1)
              cell(4) = dblstr(record1(1:1),20,idum)
              cell(5) = dblstr(record1(21:21),20,idum)
              cell(6) = dblstr(record1(41:41),20,idum)
              call getrec(safe,idnode,mxnode,nconf,record1)
              cell(7) = dblstr(record1(1:1),20,idum)
              cell(8) = dblstr(record1(21:21),20,idum)
              cell(9) = dblstr(record1(41:41),20,idum)
c
c     compute alpha and the kmax

              if(lewald.or.lspme)then

                call dcell(cell,celprp)
                eps = min(abs(eps),0.5d0)
                tol = sqrt(abs(log(eps*rcut)))
                alpha = sqrt(abs(log(eps*rcut*tol)))/rcut
                tol1 = sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
                fac = 1.d0
                if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) 
     x            fac = 2.d0**(1.d0/3.d0)
                kmax1 = nint(0.25d0 + fac*celprp(1)*alpha*tol1/pi)
                kmax2 = nint(0.25d0 + fac*celprp(2)*alpha*tol1/pi)
                kmax3 = nint(0.25d0 + fac*celprp(3)*alpha*tol1/pi)

              elseif(lhke)then

                if(nhko.eq.0)then
                  if(eps.le.1.d-6)then
                    alpha=3.46d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=3.14d0/rcut
                  else
                    alpha=2.76d0/rcut
                  endif
                elseif(nhko.eq.1)then
                  if(eps.le.1.d-6)then
                    alpha=4.37d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=4.08d0/rcut
                  else
                    alpha=3.75d0/rcut
                  endif                
                elseif(nhko.eq.2)then
                  if(eps.le.1.d-6)then
                    alpha=5.01d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=4.74d0/rcut
                  else
                    alpha=4.44d0/rcut
                  endif
                elseif(nhko.eq.3)then
                  if(eps.le.1.d-6)then
                    alpha=5.55d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=5.28d0/rcut
                  else
                    alpha=5.00d0/rcut
                  endif
                endif
                alpha=alpha/dble(2*nlatt+1)
                if(abs(cell(9)).lt.1.d-8)cell(9)=1.d0
                call dcell(cell,celprp)
                tol=2.d0*alpha*sqrt(abs(log(eps*alpha)))
                tol1=2.d0*alpha*sqrt(abs(log(eps*alpha*tol)))
                kmax1 = nint(0.25d0 + 0.5d0*celprp(1)*tol1/pi)
                kmax2 = nint(0.25d0 + 0.5d0*celprp(2)*tol1/pi)
                kmax3=1

              endif

            endif

          else

            alpha = dblstr(record,80,idum)
            ilen = 80-idum 
            record(1:ilen) = record(idum:80)
            
            kmax1=intstr(record,ilen,idum)
            jlen = ilen-idum 
            record(1:jlen) = record(idum:ilen)
            
            kmax2=intstr(record,jlen,idum)
            ilen = jlen-idum 
            record(1:ilen) = record(idum:jlen)
            
            if(lhke)then

              kmax3=1
              nhko=intstr(record,ilen,idum)
              nhko=min(nhko,3)
              jlen=ilen-idum
              record(1:jlen)=record(idum:ilen)
              nlatt=intstr(record(1:1),jlen,idum)
              nlatt=min(nlatt,2)

            else

              kmax3=intstr(record,ilen,idum)

            endif
            
          endif

c     if spme double kmax and set to next power of 2, with current upper
c     limit of 512.

          if(lspme)then

             kmaxpow2 = 1
             do while (kmax1.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2 = kmaxpow2 * 2
             end do
             kmax1 = 2 * kmaxpow2

             kmaxpow2 = 1
             do while (kmax2.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2 = kmaxpow2 * 2
             end do
             kmax2 = 2 * kmaxpow2

             kmaxpow2 = 1
             do while (kmax3.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2 = kmaxpow2 * 2
             end do
             kmax3 = 2 * kmaxpow2

          endif

          if(idnode.eq.0) then

            close(nconf)

            if(lspme)then

              write(nrite,
     x          "(/,1x,'Electrostatics : SPME  ')")

              write(nrite,
     x          "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x          /,1x,'Ewald kmax1 kmax2 kmax3     ',3i5)") 
     x          alpha,kmax1/2,kmax2/2,kmax3/2
          
            elseif(lhke)then

              write(nrite,
     x          "(/,1x,'Electrostatics : Hautman-Klein-Ewald sum  ')")

              write(nrite,
     x          "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x          /,1x,'Ewald kmax1 kmax2              ',2i5)") 
     x          alpha,kmax1,kmax2

              write(nrite,
     x          "(1x,'HKE expansion order     ',7x,i10,
     x          /,1x,'HKE lattice control     ',7x,i10)") nhko,nlatt
          
            else

              write(nrite,
     x          "(/,1x,'Electrostatics : Ewald sum  ')")

              write(nrite,
     x          "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x          /,1x,'Ewald kmax1 kmax2 kmax3     ',3i5)") 
     x          alpha,kmax1,kmax2,kmax3
          
            endif

          endif

          if (lspme) then

c     Initialize fft tables
#if defined FFTW
      call dfftw_plan_dft_3d(fplan,kmaxd,kmaxe,kmaxf,
     x     qqq,qqq,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
      call dfftw_plan_dft_3d(bplan,kmaxd,kmaxe,kmaxf,
     x     qqq,qqq,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
#elif defined SGICRAY
             call zzfft3d( 0,kmaxd,kmaxe,kmaxf,1.d0,dummy,1,1,
     x                     dummy,1,1,ffttable,dummy,dummy )

#elif defined CRAY
             call ccfft3d( 0,kmaxd,kmaxe,kmaxf,1.d0,dummy,1,1,
     x                     dummy,1,1,ffttable,dummy,dummy )
#endif
          endif


          if(lspme)then

             if(kmax1.gt.kmaxd.or.kmax2.gt.kmaxe.or.kmax3.gt.kmaxf)then
            
                kill =.true.
                call error(idnode,-185)
            
             endif

          elseif(lhke)then

             if(kmax2.gt.kmaxb) then
            
                kill =.true.
                call error(idnode,-185)
            
             endif

           else

             if(kmax2.gt.kmaxb.or.kmax3.gt.kmaxc) then
            
                kill =.true.
                call error(idnode,-185)
            
             endif

          endif

          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          lforc=.true.
          
        elseif(record(1:6).eq.'distan') then
          
          keyfce = 4
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'Electrostatics : Distance dependent dielectric')")
          
          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          
          lforc=.true.
          
        elseif(record(1:4).eq.'coul') then
          
          keyfce = 6
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'Electrostatics : Coulombic potential')")
          
          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          
          lforc=.true.
          
        elseif(record(1:5).eq.'shift') then
          
          keyfce = 8
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'Electrostatics : Shifted Coulombic potential')")
          
          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          
          lforc=.true.
          
        elseif(record(1:8).eq.'reaction') then
          
          keyfce = 10
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'Electrostatics : Reaction field')")
          
          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          
          lforc=.true.
          
        elseif(record(1:3).eq.'cap') then
          
          lfcap = .true.
          fm = dblstr(record,80,idum)
          if(fm.gt.0.d0) fmax = fm
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'force capping :',16x,1p,e12.4,' kT/A')") fmax
          
        elseif(record(1:6).eq.'no vdw') then
          
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'short-range potential terms off')")
          lvdw=.true.
          
        elseif(record(1:7).eq.'no elec') then
          
          keyfce = 0
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'electrostatic potential terms off')")
          
          if(lforc) then
            call  error(idnode,-416)
            kill =.true.
          endif
          
          lforc=.true.
          
        elseif(record(1:4).eq.'mult') then
          
          multt = intstr(record,80,idum)
          multt = max(1,multt)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'multiple timestep interval  ',3x,i10)") multt
          
        elseif(record(1:8).eq.'timestep') then
          
          lstep = .true.
          tstep = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'simulation timestep         ',3x,1p,e12.4)") tstep
          
        elseif(record(1:4).eq.'temp') then
          
          ltemp = .true.
          temp = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'simulation temperature      ',3x,1p,e12.4)") temp
          
        elseif(record(1:4).eq.'pres') then
          
          press = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'simulation pressure (katm)  ',3x,1p,e12.4)") press
c     
c     convert from katm to internal units of pressure
          
          press = press/prsunt
          lpres=.true.
          
        elseif(record(1:4).eq.'prim') then
c     
c     primary cutoff
          
          lprim =.true.
          rprim = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'primary neighbour cut off   ',3x,1p,e12.4)") rprim
          
        elseif(record(1:3).eq.'cut') then
c     
c     cutoff
          
          lcut = .true.
          rcut = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'real space cut off          ',3x,1p,e12.4)") rcut
          
        elseif(record(1:4).eq.'rvdw') then
c     
c     cutoff for short range potentials
          
          rvdw = dblstr(record,80,idum)
          lrvdw = .true.
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'real space cut off (vdw)    ',3x,1p,e12.4)") rvdw
          
        elseif(record(1:4).eq.'delr')then
c     
c     Verlet shell width
          
          ldelr = .true.
          delr = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'border width of Verlet shell',3x,1p,e12.4)") delr
          
        elseif(record(1:3).eq.'eps') then
c     
c     relative dielectric constant
          
          epsq = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,/,1x,'relative dielectric constant',2x,1p,e12.4)") epsq
          
        elseif(record(1:5).eq.'shake') then
c     
c     tolerance for shake

!          lshake=.true.
          tolnce = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'tolerance for SHAKE         ',3x,1p,e12.4)") tolnce
          
        elseif(record(1:10).eq.'quaternion') then
c     
c     tolerance for quaternion integration
          
          quattol = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'tolerance for Quaternions   ',3x,1p,e12.4)") quattol
          
        elseif(record(1:8).eq.'job time') then
c     
c     time for simulation (in seconds)
          
          timjob = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'user allocated job time (s) ',3x,1p,e12.4)") timjob
          
        elseif(record(1:10).eq.'close time') then
c     
c     time for winding up a job (in seconds)
          
          timcls = dblstr(record,80,idum)
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'job closure time        (s) ',3x,1p,e12.4)") timcls
          
        elseif(record(1:9).eq.'all pairs') then
c     
c     full minimum image - N^2 interactions each timestep
          
          lnsq =.true.
          if(idnode.eq.0) write(nrite,
     x      "(/,1x,'All-pairs requested for electrostatics')")
c
c     set ewald_spme interpolation order

        elseif(record(1:5).eq.'nospl') then

          nospl = intstr(record,80,idum)
c     
c     close control file
          
        elseif(record(1:6).eq.'finish')then
          
          goto 2000
          
        else
          
          kill = .true.
          if(idnode.eq.0) write(nrite,"(/,/,a80)") record
          call error(idnode,-3)
          
        endif
        
      enddo
c     
c     more than mega records in CONTROL file
      
      if(idnode.eq.0)close (nread)
      call error(idnode,17)
c     
c     unexpected end of file
      
 1000 continue
      if(idnode.eq.0)close (nread)
      call error(idnode,53)
      
c     
c     safe termination of reading CONTROL
      
 2000 continue
      if(idnode.eq.0)close(nread)
      
      if(.not.lforc) then
        
        if(.not.lvdw) then
          
          kill = .true.
          call error(idnode,-383)
          
        endif
        
      else
c     
c     turn on short range forces
        
        if(.not.lvdw) keyfce = keyfce+1
        
      endif
c     
c     error checking 

      if (.not.lopt2 .and. lpartem) then

        kill = .true.
        write(nrite,*)'partem has to be specified with opt!'

      endif

      if (lpartem .and. npartem.gt.10) then

        kill = .true.
        write(nrite,*)'number of parallel tempering can not be
     x greater than 10!'

      endif
      
      if(.not.ltemp) then
        
        kill = .true.
        call error(idnode,-380)
        
      endif
      
      if(.not.lstep) then
        
        kill = .true.
        call error(idnode,-381)
        
      endif
      
      if(.not.lcut) then
        
        kill = .true.
        call error(idnode,-382)
        
      endif
c     
c     check if van der Waals cutoff set
      
      if(.not.lrvdw.and.mod(keyfce,2).eq.1) then
        
        if(lcut)then
          
          rvdw=rcut
          
        else
          
          kill = .true.
          call error(idnode,-402)

        endif      
        
      endif
      
      if(.not.ldelr) then
        
        kill = .true.
        call error(idnode,-384)
        
      endif
      
      if(multt.gt.1) then
        
        if(.not.lprim) then
          
          kill = .true.
          call error(idnode,-385)
          
        elseif(rprim.gt.rcut) then
          
          kill = .true.
          call error(idnode,-386)
          
        endif
        
      endif
c     
c     check settings in nvt ensemble
      
      if(keyens.ge.2.and.keyens.le.3) then

        if(taut.le.0.d0) then

          kill=.true.
          call error(idnode,-464)

        endif

      endif
c     
c     check settings in npt ensemble
      
      if(keyens.ge.4.and.keyens.le.7) then
        
        if(.not.lpres) then
          
          kill = .true.
          call error(idnode,-387)
          
        endif
c
c     check barostat and thermostat rates non zero

        if(taut.le.0.d0) then

          kill=.true.
          call error(idnode,-464)

        endif
        if(taup.le.0.d0) then

          kill=.true.
          call error(idnode,-466)

        endif

      endif
c     
c     check multiple timestep cutoffs are sensible
      
      if(multt.gt.1) then
        if(rcut-rprim.lt.delr) then
        
          kill = .true.
          call error(idnode,-398)
        
        endif
      endif
c
c     check rcut > rvdw (for verlet list constructor)

      if(rcut.lt.rvdw) then 
        
        kill = .true.
        call error(idnode,-400)
        
      endif
c
c     check spme is not being used with incorrect pbc

      if(lspme)then

        if(imcon.eq.0.or.imcon.eq.6)then

          kill=.true.
          call error(idnode,-513)

        endif

      endif
c     
c     check on all-pairs calculation request
      
      if(lnsq) then
        
        if(multt.eq.1) then
          
          kill = .true.
          call error(idnode,-422)
          
        endif
        
        if(keyfce/2.lt.2.or.keyfce/2.gt.3) then
          
          kill = .true.
          call error(idnode,-424)
          
        endif
        
      endif
c
c     check thole smearing function

      if (lthole) then
        if (nthole.ne.1.and.nthole.ne.3.and.nthole.ne.4) then
           write(nrite,*)
           write(nrite,*)'Error! need specify Thole smearing function!'
           kill = .true.
        endif
      endif
c     
c     check on steps before temperature scaling
      
      if(nstbts.eq.0) nstbts = nstrun+1
      
      if (kill) call error(idnode,0)
#ifdef VAMPIR
      call VTEND(3, ierr)
#endif

      return
      end

      

