      subroutine ewald_spme
     x (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,fplan,bplan,
     x  engcpe,vircpe,alpha,volm,epsq,nauxfft,key1,key2,key3,cell,chge,
     x  xxx,yyy,zzz,txx,tyy,tzz,fxx,fyy,fzz,stress,buffer,csp,ww1,ww2,
     x  ww3,bspx,bspy,bspz,bsdx,bsdy,bsdz,qqc,qqq,bscx,bscy,bscz,
     x  ffttable)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using the smoothed particle mesh ewald method
c     due to Essmann et al J. Chem. Phys. 103 (1995) 8577.
c     
c     parallel replicated data version (part 1)
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith july 1998
c     additional FFT code - j. geronowicz sept 1999
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c     wl
c     2002/05/31 13:59:31
c     1.6
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"
#include "fftw3.f"      
 
      logical newjob,lconsw
      
      dimension key1(kmaxd),key2(kmaxe),key3(kmaxf)
      dimension qqc(kmaxd,kmaxe,kmaxf),celprp(10),ffttable(mxftab)
      dimension chge(mxatms),cell(9),rcell(9),nauxfft(4)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension csp(mxspl),buffer(mxbuff),stress(9),omg(9)
      dimension bspx(mxspme,mxspl),bspy(mxspme,mxspl),bspz(mxspme,mxspl)
      dimension bsdx(mxspme,mxspl),bsdy(mxspme,mxspl),bsdz(mxspme,mxspl)
      complex*16 qqq(kmaxd,kmaxe,kmaxf),cpetot,vterm
      complex*16 bscx(kmaxd),bscy(kmaxe),bscz(kmaxf)
      complex*16 ww1(kmaxd),ww2(kmaxe),ww3(kmaxf)

#ifdef FFTW
      FFTW_PLAN_TYPE fplan, bplan
#else
      integer fplan, bplan
#endif

      real*8 bb1,bb2,bb3

      save newjob,engsic
      
      data newjob/.true./

#ifdef FFTW
c done in simdef.f
c      call dfftw_plan_dft_3d(fplan,kmaxd,kmaxe,kmaxf,
c     x     qqq,qqq,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
c      call dfftw_plan_dft_3d(bplan,kmaxd,kmaxe,kmaxf,
c     x     qqq,qqq,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
#endif

#ifdef VAMPIR
      call VTBEGIN(138, ierr)
#endif
      npass=1
      lconsw=.true.
      twopi=2.d0*pi

      if(newjob)then

        newjob=.false.

c     calculate self interaction correction
        
        engsic=0.d0
        
        do i=1,natms
          
          engsic=engsic+chge(i)**2
          
        enddo

        engsic=-r4pie0/epsq*alpha*engsic/sqrpi
#if defined ESSL || defined CRAY || defined FFTW || defined SGICRAY
c     initialise the complex exponential arrays

        call spl_cexp(kmax1,kmax2,kmax3,ww1,ww2,ww3)
#else
c     initialise the default fft routine

      call dlpfft3(1,1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)
#endif

c     calculate B-spline coefficients

        call bspcoe
     x    (nospl,kmax1,kmax2,kmax3,csp,bscx,bscy,bscz,ww1,ww2,ww3)

      endif

c     initialise coulombic potential energy
      
      engcpe=0.d0
      vircpe=0.d0

c     initalize stress tensor working arrays

      do i = 1,9
        omg(i) = 0.d0
      enddo

c     set working parameters

      rvolm=twopi/volm
      ralph=-0.25d0/alpha**2

c     set switch for TO, RD and HP boundary conditions

      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) then

        npass=2
        lconsw=.false.
        rvolm=0.5d0*rvolm
        shiftx=0.5d0*dble(kmax1)
        shifty=0.5d0*dble(kmax2)
        shiftz=0.5d0*dble(kmax3)
        if(imcon.eq.7)shiftz=0.d0

      endif

c     convert cell coordinates to fractional coordinates

      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)

      do i=1,natms
        
        txx(i)=dble(kmax1)*(rcell(1)*xxx(i)+rcell(4)*yyy(i)+
     x    rcell(7)*zzz(i)+0.5d0)
        tyy(i)=dble(kmax2)*(rcell(2)*xxx(i)+rcell(5)*yyy(i)+
     x    rcell(8)*zzz(i)+0.5d0)
        tzz(i)=dble(kmax3)*(rcell(3)*xxx(i)+rcell(6)*yyy(i)+
     x    rcell(9)*zzz(i)+0.5d0)

      enddo
      
c     construct B-splines for atoms
      
      call bspgen
     x  (nospl,natms,txx,tyy,tzz,bspx,bspy,bspz,bsdx,bsdy,bsdz)
      
c     zero 3D charge array
      
      nnn=kmaxd*kmaxe*kmaxf
      call set_block(nnn,0.d0,qqc)
      
c     construct 3D charge array
      
      do ipass=1,npass

        do i=1,natms
          
          do l=1,nospl
            
            ll=int(tzz(i))-l+2
            if(ll.gt.kmax3)ll=1
            if(ll.lt.1)ll=ll+kmax3
            do k=1,nospl
              
              kk=int(tyy(i))-k+2
              if(kk.gt.kmax2)kk=1
              if(kk.lt.1)kk=kk+kmax2
              
              
              do j=1,nospl
                
                jj=int(txx(i))-j+2
                if(jj.gt.kmax1)jj=1
                if(jj.lt.1)jj=jj+kmax1
                
                qqc(jj,kk,ll)=qqc(jj,kk,ll)+
     x            chge(i)*bspx(i,j)*bspy(i,k)*bspz(i,l)
                
              enddo
              
            enddo
            
          enddo
          
        enddo

        if(.not.lconsw)then

          do i=1,natms

            tx=txx(i)-shiftx
            ty=tyy(i)-shifty
            tz=tzz(i)-shiftz
            txx(i)=txx(i)-sign(shiftx,tx)
            tyy(i)=tyy(i)-sign(shifty,ty)
            tzz(i)=tzz(i)-sign(shiftz,tz)

          enddo

        endif

      enddo

c     load charge array into complex array for FFT

      call cpy_rtc(nnn,qqc,qqq)

c     calculate inverse 3D FFT of charge array (in place).

#if defined FFTW
      call dfftw_execute(fplan)
#elif defined ESSL
      inc2=kmaxd
      inc3=kmaxd*kmaxe

      call dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmax1,kmax2,kmax3,
     x  -1,1.d0,buffer,mxbuff)
#elif defined SGICRAY
      call zzfft3d( -1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )
#elif defined CRAY
      call ccfft3d( -1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )
#else
      call dlpfft3(0,1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)
#endif
c     set reciprocal space cutoff

      call dcell(rcell,celprp)

      rcpcut=0.5d0*min(dble(kmax1)*celprp(7),dble(kmax2)*celprp(8),
     x  dble(kmax3)*celprp(9))
      rcpcut=rcpcut*1.05d0*twopi
      rcpct2=rcpcut**2

c     calculate convolution of charge array with gaussian function

      do l=1,kmax3

        ll=l-1
        if(l.gt.kmax3/2)ll=l-kmax3-1
        tmp=twopi*dble(ll)
        rkx1=tmp*rcell(3)
        rky1=tmp*rcell(6)
        rkz1=tmp*rcell(9)
        bb3 = bscz(l)*conjg(bscz(l))

        do k=1,kmax2

          kk=k-1
          if(k.gt.kmax2/2)kk=k-kmax2-1
          tmp=twopi*dble(kk)
          rkx2=rkx1+tmp*rcell(2)
          rky2=rky1+tmp*rcell(5)
          rkz2=rkz1+tmp*rcell(8)
          bb2 = bb3*bscy(k)*conjg(bscy(k))
          
          do j=1,kmax1

            jj=j-1
            if(j.gt.kmax1/2)jj=j-kmax1-1
            tmp=twopi*dble(jj)
            rkx3=rkx2+tmp*rcell(1)
            rky3=rky2+tmp*rcell(4)
            rkz3=rkz2+tmp*rcell(7)
            bb1 = bb2*bscx(j)*conjg(bscx(j))
                
            rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

            if(rksq.gt.1.d-6.and.rksq.le.rcpct2)then

              vterm=bb1*exp(ralph*rksq)/rksq*qqq(j,k,l)
              akv=2.d0*(1.d0/rksq-ralph)*real(vterm*conjg(qqq(j,k,l)))
              omg(1)=omg(1)-rkx3*rkx3*akv
              omg(5)=omg(5)-rky3*rky3*akv
              omg(9)=omg(9)-rkz3*rkz3*akv
              omg(2)=omg(2)-rkx3*rky3*akv
              omg(3)=omg(3)-rkx3*rkz3*akv
              omg(6)=omg(6)-rky3*rkz3*akv
              qqq(j,k,l)=vterm

            else

              qqq(j,k,l)=(0.d0,0.d0)

            endif

          enddo

        enddo

      enddo

#if defined FFTW
      call dfftw_execute(bplan)
#elif defined ESSL
      call dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmax1,kmax2,kmax3,
     x  1,1.d0,buffer,mxbuff)
#elif defined SGICRAY
      call zzfft3d( 1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )
#elif defined CRAY
      call ccfft3d( 1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )
#else
      call dlpfft3(0,-1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)
#endif

c     calculate atomic forces

      call spme_for
     x  (idnode,mxnode,nospl,natms,kmax1,kmax2,kmax3,rvolm,
     x  epsq,chge,txx,tyy,tzz,fxx,fyy,fzz,bspx,bspy,bspz,
     x  bsdx,bsdy,bsdz,rcell,qqq,buffer)

c     complete product of charge array and its gaussian convolution

      call ele_prd(nnn,qqq,qqc,qqq)

c     calculate total energy

      call scl_csum(nnn,cpetot,qqq)

      eng1=real(cpetot)
      den=1.d0/dble(npass)
      engcpe=(den*rvolm*r4pie0*eng1/epsq+engsic)/dble(mxnode)
c     
c     calculate stress tensor (symmetrical)

      scal1=den*rvolm*r4pie0/(epsq*dble(mxnode))
      stress(1) = stress(1)+scal1*(omg(1)+eng1)
      stress(2) = stress(2)+scal1*omg(2)
      stress(3) = stress(3)+scal1*omg(3)
      stress(4) = stress(4)+scal1*omg(2)
      stress(5) = stress(5)+scal1*(omg(5)+eng1)
      stress(6) = stress(6)+scal1*omg(6)
      stress(7) = stress(7)+scal1*omg(3)
      stress(8) = stress(8)+scal1*omg(6)
      stress(9) = stress(9)+scal1*(omg(9)+eng1)
c     
c     virial term

      vircpe=-scal1*(omg(1)+omg(5)+omg(9)+3.d0*eng1)

#ifdef VAMPIR
      call VTEND(138, ierr)
#endif
      return
      end
