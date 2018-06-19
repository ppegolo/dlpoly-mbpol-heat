      subroutine ewald1
     x  (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x  engcpe,vircpe,alpha,volm,epsq,cell,chge,xxx,
     x  yyy,zzz,fxx,fyy,fzz,elc,emc,enc,els,ems,ens,
     x  ckc,cks,clm,slm,stress,buffer,ewlbuf)

c
c***********************************************************************
c
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c
c     parallel replicated data version (part 1)
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c
c     version 2
c     author    - t. forester april 1993
c     T3D optimizaion . t.forester july 1994
c
c     part 1 - reciprocal space terms (fourier part)
c
c     note - in loop over all k vectors k=2pi(ll/cl,mm/cl,nn/cl)
c     the values of ll,mm and nn are selected so that the symmetry of
c     reciprocal lattice is taken into account i.e. the following
c     rules apply.
c
c     ll ranges over the values 0 to kmax1 only.
c
c     mm ranges over 0 to kmax2 when ll=0 and over
c     -kmax2 to kmax2 otherwise.
c     nn ranges over 1 to kmax3 when ll=mm=0 and over
c     -kmax3 to kmax3 otherwise.
c
c     hence the result of the summation must be doubled at the end.
c
c     stress tensor added t.forester may 1994
c
c     wl
c     2001/06/12 12:46:53
c     1.8
c     Exp
c
c***********************************************************************
c
#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_ew1p, update_energy_ew1p,
     x                       update_force_ew1p, update_forces
#endif /* HEAT_CURRENT */

#include "dl_params.inc"

      logical newjob,lconsw,safe,leven

      dimension chge(mxatms),cell(9),rcell(9)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld)
      dimension elc(mxewld,0:1),els(mxewld,0:1)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc)
      dimension buffer(mxbuff),ewlbuf(mxebuf)
      dimension stress(9),omg(9)

#ifdef HEAT_CURRENT
      real(8) :: ckcs, ckss, dipdotk, stressfactor, scalefactor
      real(8) :: strs1k, strs2k, strs3k, strs4k, strs5k, strs6k
      real(8) :: strs7k, strs8k, strs9k
      real(8) :: dipvec(1:3), kvec(1:3)

      real(8) :: Btensor(9), akv2
      real(8) :: Cqi, Sqi, CMi(3), SMi(3)
      real(8) :: Qc, Qs, Mc(3), Ms(3)
      real(8) :: stressCC(9), stressDD(9)
      real(8) :: stressCD(9), stressDC(9)
      real(8) :: dpkCMi, dpkMc, dpkMs, dpkSMi
      real(8) :: force_tmp(3), rij(3), qforce2
      integer :: j2,i2
#endif /* HEAT_CURRENT */

      save newjob,engsic

      data newjob/.true./,lconsw/.true./,safe/.true./,leven/.true./

#ifdef VAMPIR
      call VTBEGIN(97, ierr)
#endif
      twopi=2.d0*pi
c
c     initialise coulombic potential energy

      engcpe=0.d0
      vircpe=0.d0

      if(alpha.lt.1.d-6)return
      if(mxewld.ne.msatms) call error(idnode,330)
c
c     set up atoms numbers for nodes

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c
c     initalize stress tensor working arrays

      do i = 1,9
        omg(i) = 0.d0
      enddo
c
c     set working parameters

      rvolm=twopi/volm
      ralph=-0.25d0/alpha**2
c
c     set switch for TO, RD and HP boundary conditions

      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) then

        lconsw=.false.
        rvolm=0.5d0*rvolm

      endif
c
c     construct reciprocal lattice vectors and set k vector range
      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      call dcell(rcell,buffer)

      rcpcut=min(dble(kmax1)*buffer(7),dble(kmax2)*buffer(8),
     x  dble(kmax3)*buffer(9))
      rcpcut=rcpcut*1.05d0*twopi
      rcpct2=rcpcut**2

      if(newjob)then
c
c     calculate self interaction correction

        engsic=0.d0

        do i=iatm0,iatm1

          engsic=engsic+chge(i)**2

        enddo

        engsic=-r4pie0/epsq*alpha*engsic/sqrpi
        newjob=.false.

      endif

      do i=iatm0,iatm1
#ifdef HEAT_CURRENT
        call update_energy_ew1p(i,-r4pie0/epsq*
     x    alpha*chge(i)**2/sqrpi)
#endif /* HEAT_CURRENT */
      end do

c
c     calculate and store exponential factors
c     convert real to reciprocal space coordinates

      i=0

      do j=iatm0,iatm1

        i=i+1
        elc(i,0)=1.d0
        emc(i,0)=1.d0
        enc(i,0)=1.d0
        els(i,0)=0.d0
        ems(i,0)=0.d0
        ens(i,0)=0.d0
        ssx=rcell(1)*xxx(j)+rcell(4)*yyy(j)+rcell(7)*zzz(j)
        ssy=rcell(2)*xxx(j)+rcell(5)*yyy(j)+rcell(8)*zzz(j)
        ssz=rcell(3)*xxx(j)+rcell(6)*yyy(j)+rcell(9)*zzz(j)
        elc(i,1)=cos(twopi*ssx)
        emc(i,1)=cos(twopi*ssy)
        enc(i,1)=cos(twopi*ssz)
        els(i,1)=sin(twopi*ssx)
        ems(i,1)=sin(twopi*ssy)
        ens(i,1)=sin(twopi*ssz)

      enddo

      limit=i

      do l=2,kmax2

        do i=1,limit

          emc(i,l)=emc(i,l-1)*emc(i,1)-ems(i,l-1)*ems(i,1)
          ems(i,l)=ems(i,l-1)*emc(i,1)+emc(i,l-1)*ems(i,1)

        enddo

      enddo

      do l=2,kmax3

        do i=1,limit

          enc(i,l)=enc(i,l-1)*enc(i,1)-ens(i,l-1)*ens(i,1)
          ens(i,l)=ens(i,l-1)*enc(i,1)+enc(i,l-1)*ens(i,1)

        enddo

      enddo
c
c     start of main loop over k vectors

      npass=1
      if((mxnode.gt.16).or.(mxebuf.gt.5000))npass=2

      do ipass=1,npass

        kkk=0
        mmin=0
        nmin=1

        do ll=0,kmax1

          l=ll
          tmp = twopi*dble(ll)
          rkx1=tmp*rcell(1)
          rky1=tmp*rcell(4)
          rkz1=tmp*rcell(7)

c
c     put cos(i,L) terms into cos(i,0) array

          if(l.eq.1) then

            do i=1,limit

              elc(i,0)=elc(i,1)
              els(i,0)=els(i,1)

            enddo

          elseif(l.gt.1) then

            do i=1,limit

              cs  = elc(i,0)
              elc(i,0)=cs      *elc(i,1)-els(i,0)*els(i,1)
              els(i,0)=els(i,0)*elc(i,1)+cs      *els(i,1)

            enddo

          endif

          do mm=mmin,kmax2

            m=iabs(mm)
            tmp = twopi*dble(mm)
            rkx2=rkx1+tmp*rcell(2)
            rky2=rky1+tmp*rcell(5)
            rkz2=rkz1+tmp*rcell(8)
c
c     set temporary products of exponential terms

            if(mm.ge.0)then

              do i=1,limit

                clm(i)=elc(i,0)*emc(i,m)-els(i,0)*ems(i,m)
                slm(i)=els(i,0)*emc(i,m)+ems(i,m)*elc(i,0)

              enddo

            else

              do i=1,limit

                clm(i)=elc(i,0)*emc(i,m)+els(i,0)*ems(i,m)
                slm(i)=els(i,0)*emc(i,m)-ems(i,m)*elc(i,0)

              enddo

            endif

            do nn=nmin,kmax3

              n=iabs(nn)

              if(.not.lconsw)then

                if(imcon.eq.7)then

                  leven=(mod(l+m,2).eq.0)

                else

                  leven=(mod(l+m+n,2).eq.0)

                endif

              endif

              if(lconsw.or.leven)then

                tmp = twopi*dble(nn)
                rkx3=rkx2+tmp*rcell(3)
                rky3=rky2+tmp*rcell(6)
                rkz3=rkz2+tmp*rcell(9)

c
c     test on magnitude of k vector

                rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

                if(rksq.le.rcpct2)then

c
c     calculate exp(ikr) terms and product with charges

                  i=0
                  if(nn.ge.0)then

                    do j=iatm0,iatm1

                      i=i+1

                      tchge=chge(j)
                      tclm=clm(i)
                      tenc=enc(i,n)
                      tslm=slm(i)
                      tens=ens(i,n)

                      ckc(i)=tchge*(tclm*tenc-tslm*tens)
                      cks(i)=tchge*(tslm*tenc+tclm*tens)

                    enddo

                  else

                    do j=iatm0,iatm1

                      i=i+1
                      tchge=chge(j)
                      tclm=clm(i)
                      tenc=enc(i,n)
                      tslm=slm(i)
                      tens=ens(i,n)

                      ckc(i)=tchge*(tclm*tenc+tslm*tens)
                      cks(i)=tchge*(tslm*tenc-tclm*tens)

                    enddo

                  endif

                  if(ipass.eq.1)then

c
c     calculate vector sums

                    ckcs=0.d0
                    ckss=0.d0

                    do i=1,limit

                      ckcs=ckcs+ckc(i)
                      ckss=ckss+cks(i)

                    enddo

c
c     perform global summation of exp(ikr) terms or store if npass=2

                    if(npass.eq.2)then

                      if(kkk+2.le.mxebuf)then

                        ewlbuf(kkk+1)=ckcs
                        ewlbuf(kkk+2)=ckss

                      else

                        safe=.false.

                      endif

                    elseif(mxnode.gt.1)then

                      buffer(3)=ckcs
                      buffer(4)=ckss
                      call gdsum(buffer(3),2,buffer(1))
                      ckcs=buffer(3)
                      ckss=buffer(4)

                    endif

                  endif

                  if(ipass.eq.npass)then

                    if(npass.eq.2)then

                      ckcs=ewlbuf(kkk+1)
                      ckss=ewlbuf(kkk+2)

                    endif
c
c     calculate akk coefficients

                    rrksq = 1.d0/rksq
                    if(lconsw)then
                      akk=exp(ralph*rksq)*rrksq
                    else
                      akk=4.0d0*exp(ralph*rksq)*rrksq
                    endif
                    bkk=akk
                    akv=2.d0*akk*(rrksq-ralph)

#ifdef HEAT_CURRENT
! PP_: begin new loop over atoms to store the atomic virial -----28-03-2018---------------------------------------------
                  scalefactor = 2.d0*rvolm*r4pie0/epsq

                  i = 0
                  kvec = (/rkx3,rky3,rkz3/)
                  do j = iatm0, iatm1
                    i = i + 1

      strs1k = -scalefactor*(ckc(i)*ckcs+cks(i)*ckss)*(akv*rkx3
     x  *rkx3-akk)
       strs2k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rkx3*rky3
       strs3k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rkx3*rkz3
       strs4k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rkx3*rky3
       strs5k = -scalefactor*(ckc(i)*ckcs+cks(i)*ckss)*(akv*rky3
     x   *rky3-akk)
       strs6k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rky3*rkz3
       strs7k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rkx3*rkz3
       strs8k = -scalefactor*akv*(ckc(i)*ckcs+cks(i)*ckss)*rky3*rkz3
       strs9k = -scalefactor*(ckc(i)*ckcs+cks(i)*ckss)*(akv*rkz3
     x   *rkz3-akk)

                    call update_stress_ew1p(j,1,1,-strs1k)
                    call update_stress_ew1p(j,1,2,-strs2k)
                    call update_stress_ew1p(j,1,3,-strs3k)
                    call update_stress_ew1p(j,2,1,-strs4k)
                    call update_stress_ew1p(j,2,2,-strs5k)
                    call update_stress_ew1p(j,2,3,-strs6k)
                    call update_stress_ew1p(j,3,1,-strs7k)
                    call update_stress_ew1p(j,3,2,-strs8k)
                    call update_stress_ew1p(j,3,3,-strs9k)
                    !
                    call update_energy_ew1p(j,scalefactor*akk*
     x                (ckc(i)*ckcs+cks(i)*ckss))



                   qforce=bkk*(cks(i)*ckcs-ckc(i)*ckss)
                   fxx(j)=fxx(j)+rkx3*qforce
                   fyy(j)=fyy(j)+rky3*qforce
                   fzz(j)=fzz(j)+rkz3*qforce

                    /*i2=0
                     do j2 = iatm0, iatm1
                      i2 = i2+1
                      !if (j2==j) cycle
                      !rij = (/xxx(j)-xxx(j2),yyy(j)-yyy(j2),
!    x                  zzz(j)-zzz(j2)/)
          !force_tmp = 2*scalefactor*bkk*dot_product(kvec,rij)*
!     x (cks(j)*ckc(j2)-ckc(j)*cks(j2))*rij/dot_product(rij,rij)
                  force_tmp = bkk*(cks(i)*
     x                  ckc(i2)-ckc(i)*cks(i2))*kvec
            force_tmp=force_tmp*4.d0*rvolm*r4pie0/epsq
                  call update_forces(j,j2,force_tmp)
                  call update_forces(j2,j,-force_tmp)
                    end do */
                  end do
! PP_: end loop over atoms to store the atomic virial --------------
!#else /* HEAT_CURRENT */

              i=0
              do j=iatm0,iatm1
                i=i+1
                qforce=bkk*(cks(i)*ckcs-ckc(i)*ckss)
                fxx(j)=fxx(j)+rkx3*qforce
                fyy(j)=fyy(j)+rky3*qforce
                fzz(j)=fzz(j)+rkz3*qforce
                force_tmp =(/rkx3*qforce,rky3*qforce,rkz3*qforce/)
                force_tmp = 2*4.d0*rvolm*r4pie0/epsq*force_tmp
                call update_force_ew1p(j,force_tmp)
                ! PP_:
                /* i2=0
                do j2=iatm0,iatm1
                  i2=i2+1
                  if (j2==j1) cycle
                  qforce2=bkk*(cks(i)*ckc(i2)-ckc(i)*cks(i2))
                  force_tmp =(/rkx3*qforce2,rky3*qforce2,rkz3*qforce2/)
                  force_tmp = 2*4.d0*rvolm*r4pie0/epsq*force_tmp
                  call update_forces(j,j2,force_tmp)
                  call update_forces(j2,j,-force_tmp)
                enddo */
              enddo
#endif /*HEAT_CURRENT*/

c
c     accumulate potential energy and virial terms

                    engcpe=engcpe+akk*(ckcs*ckcs+ckss*ckss)
                    virprs=akv*(ckcs*ckcs+ckss*ckss)

                    omg(1)=omg(1)-virprs*rkx3*rkx3
                    omg(5)=omg(5)-virprs*rky3*rky3
                    omg(9)=omg(9)-virprs*rkz3*rkz3
#ifdef STRESS
                    omg(2)=omg(2)-virprs*rkx3*rky3
                    omg(3)=omg(3)-virprs*rkx3*rkz3
                    omg(6)=omg(6)-virprs*rky3*rkz3
#endif


c
c     end vector loop

                  endif

                  kkk=kkk+2

                endif

              endif

            enddo

            nmin=-kmax3

          enddo

          mmin=-kmax2

        enddo

c
c     delayed global sum of exp(ikr) terms for npass=2 case

        if(ipass.eq.1.and.npass.eq.2)then

          if(safe)then

            call gdsum(ewlbuf,kkk,buffer)

            do i=1,limit

              elc(i,0)=1.d0
              els(i,0)=0.d0

            enddo

          else

            if(idnode.eq.0)then

              write(nrite,'(a,i10)')
     x          'dimension of ewlbuf array required ',kkk
              write(nrite,'(a,i10)')
     x          'dimension of current  ewlbuf array ',mxebuf

            endif

            call error(idnode,46)

          endif

        endif

      enddo

      engcpe=engcpe/dble(mxnode)
      do i = 1,9
        omg(i) = omg(i)/dble(mxnode)
      enddo
c
c     add self interaction correction to potential

      if(lconsw)then

        eng1 = engcpe
        engcpe=2.d0*rvolm*r4pie0*engcpe/epsq+engsic
        scal1 = 2.d0*rvolm*r4pie0/epsq
        scale = 4.d0*rvolm*r4pie0/epsq

      else

        eng1 = engcpe
        engcpe=rvolm*r4pie0*engcpe/epsq+engsic
        scal1 = rvolm*r4pie0/epsq
        scale = 2.d0*rvolm*r4pie0/epsq

      endif


c
c     calculate final forces

      do i=iatm0,iatm1

        fxx(i)=scale*fxx(i)
        fyy(i)=scale*fyy(i)
        fzz(i)=scale*fzz(i)

      enddo
#ifdef STRESS
c
c     calculate stress tensor (symmetrical)

      stress(1) = stress(1)+scal1*(omg(1)+eng1)
      stress(2) = stress(2)+scal1*omg(2)
      stress(3) = stress(3)+scal1*omg(3)
      stress(4) = stress(4)+scal1*omg(2)
      stress(5) = stress(5)+scal1*(omg(5)+eng1)
      stress(6) = stress(6)+scal1*omg(6)
      stress(7) = stress(7)+scal1*omg(3)
      stress(8) = stress(8)+scal1*omg(6)
      stress(9) = stress(9)+scal1*(omg(9)+eng1)
#endif
c
c     virial term

      vircpe = -scal1*(omg(1)+omg(5)+omg(9) + 3.d0*eng1)

#ifdef VAMPIR
      call VTEND(97, ierr)
#endif
      return
      end
