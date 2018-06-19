! 06 FEB 06 - IUCHI - ESP FOR DMS
! 14 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 14 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
!
      subroutine ewald1cp
     x     (lpolar,idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,
     x     engcpe,vircpe,alpha,volm,epsq,cell,chge,xxx,yyy,zzz,
     x     fxx,fyy,fzz,elc,emc,enc,els,ems,ens,ckc,cks,clm,slm,
     x     stress,buffer,ewlbuf,efdcrecx,efdcrecy,efdcrecz,
     x     efddmurecx,efddmurecy,efddmurecz,dipx,dipy,dipz,polr,polr2,
     x     ckr,skr)
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
c     version 3 - voth group
c     author    - c. j. burnham
c
c
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
!
!     Last updated: 6 Feb 2006 by S. Iuchi
c     
c***********************************************************************
c     
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces, only: fx_ttm_per, fy_ttm_per, fz_ttm_per, 
     $                      fx_ttm_ind, fy_ttm_ind, fz_ttm_ind
#endif /* TTM_FORCE_DECOMPOSITION */
      use ps_type_dms, only: vesp_k, vesp_s, vesp_dc_k, ldms
      
#include "dl_params.inc"
      
      logical newjob,lconsw,safe,leven
      logical lpolar
      
      dimension chge(mxatms),cell(9),rcell(9)
      dimension polr(mxatms),polr2(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ckc(mxewld),cks(mxewld),clm(mxewld),slm(mxewld)
      dimension ckr(mxewld),skr(mxewld)
      dimension elc(mxewld,0:1),els(mxewld,0:1)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc)
      dimension buffer(mxbuff),ewlbuf(mxebuf)
      dimension stress(9),omg(9)

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
      engdum=0.d0

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
!
!     with DMS

      if( ldms ) newjob = .true.

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
                  dmuckrsx=0.
                  dmuckrsy=0.
                  dmuckrsz=0.

                  dmuskrsx=0.
                  dmuskrsy=0.
                  dmuskrsz=0.

                  skdmuckrs=0.
                  skdmuskrs=0.
                  
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

                      ckr(i)=(tclm*tenc-tslm*tens)
                      skr(i)=(tslm*tenc+tclm*tens)

                      ckc(i)=tchge*ckr(i)
                      cks(i)=tchge*skr(i)
c
c     calculate dipole sums

                      if (polr2(j) .gt. 1.d-6) then

                         dmuckrsx=dmuckrsx+ckr(i)*dipx(j)
                         dmuckrsy=dmuckrsy+ckr(i)*dipy(j)
                         dmuckrsz=dmuckrsz+ckr(i)*dipz(j)

                         dmuskrsx=dmuskrsx+skr(i)*dipx(j)
                         dmuskrsy=dmuskrsy+skr(i)*dipy(j)
                         dmuskrsz=dmuskrsz+skr(i)*dipz(j)

                      endif

                    enddo
                    
                  else
                    
                    do j=iatm0,iatm1
                      
                      i=i+1
                      tchge=chge(j)
                      tclm=clm(i)
                      tenc=enc(i,n)
                      tslm=slm(i)
                      tens=ens(i,n)

                      ckr(i)=(tclm*tenc+tslm*tens)
                      skr(i)=(tslm*tenc-tclm*tens)

                      ckc(i)=tchge*ckr(i)
                      cks(i)=tchge*skr(i)
c
c     calculate dipole sums

                      if (polr2(j) .gt. 1.d-6) then

                         dmuckrsx=dmuckrsx+ckr(i)*dipx(j)
                         dmuckrsy=dmuckrsy+ckr(i)*dipy(j)
                         dmuckrsz=dmuckrsz+ckr(i)*dipz(j)

                         dmuskrsx=dmuskrsx+skr(i)*dipx(j)
                         dmuskrsy=dmuskrsy+skr(i)*dipy(j)
                         dmuskrsz=dmuskrsz+skr(i)*dipz(j)

                      endif

                    enddo
                    
                  endif

                  if(ipass.eq.1)then

c     
c     calculate vector sums
                    
                    ckcs=0.d0
                    ckss=0.d0
                    
                    do i=1,limit
c
c     calculate point charge sums
                      
                      ckcs=ckcs+ckc(i)
                      ckss=ckss+cks(i)

                    enddo
c
c      take dot product of k with dipole sums

                    skdmuckrs=rkx3*dmuckrsx+rky3*dmuckrsy
     x                        +rkz3*dmuckrsz
                    skdmuskrs=rkx3*dmuskrsx+rky3*dmuskrsy
     x                        +rkz3*dmuskrsz

                    summu1=(ckcs-skdmuskrs)
                    summu2=(ckss+skdmuckrs)
                    
c     
c     perform global summation of exp(ikr) terms or store if npass=2
                    
                    if(npass.eq.2)then

                      if(kkk+12.le.mxebuf)then
                        
                        ewlbuf(kkk+1)=ckcs
                        ewlbuf(kkk+2)=ckss
                        ewlbuf(kkk+3)=skdmuckrs
                        ewlbuf(kkk+4)=skdmuskrs
                        ewlbuf(kkk+5)=summu1
                        ewlbuf(kkk+6)=summu2
                        ewlbuf(kkk+7)=dmuckrsx
                        ewlbuf(kkk+8)=dmuckrsy
                        ewlbuf(kkk+9)=dmuckrsz
                        ewlbuf(kkk+10)=dmuskrsx
                        ewlbuf(kkk+11)=dmuskrsy
                        ewlbuf(kkk+12)=dmuskrsz

                      else

                        safe=.false.

                      endif

                    elseif(mxnode.gt.1)then
                      
                      buffer(13)=ckcs
                      buffer(14)=ckss
                      buffer(15)=skdmuckrs
                      buffer(16)=skdmuskrs
                      buffer(17)=summu1
                      buffer(18)=summu2
                      buffer(19)=dmuckrsx
                      buffer(20)=dmuckrsy
                      buffer(21)=dmuckrsz
                      buffer(22)=dmuskrsx
                      buffer(23)=dmuskrsy
                      buffer(24)=dmuskrsz
                      call gdsum(buffer(13),12,buffer(1))
                      ckcs=buffer(13)
                      ckss=buffer(14)
                      skdmuckrs=buffer(15)
                      skdmuskrs=buffer(16)
                      summu1=buffer(17)
                      summu2=buffer(18)
                      dmuckrsx=buffer(19)
                      dmuckrsy=buffer(20)
                      dmuckrsz=buffer(21)
                      dmuskrsx=buffer(22)
                      dmuskrsy=buffer(23)
                      dmuskrsz=buffer(24)

                    endif

                  endif

                  if(ipass.eq.npass)then

                    if(npass.eq.2)then

                      ckcs=ewlbuf(kkk+1)
                      ckss=ewlbuf(kkk+2)
                      skdmuckrs=ewlbuf(kkk+3)
                      skdmuskrs=ewlbuf(kkk+4)
                      summu1=ewlbuf(kkk+5)
                      summu2=ewlbuf(kkk+6)
                      dmuckrsx=ewlbuf(kkk+7)
                      dmuckrsy=ewlbuf(kkk+8)
                      dmuckrsz=ewlbuf(kkk+9)
                      dmuskrsx=ewlbuf(kkk+10)
                      dmuskrsy=ewlbuf(kkk+11)
                      dmuskrsz=ewlbuf(kkk+12)

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
c     
c     accumulate potential energy and virial terms

                 engcpe=engcpe+akk*(ckcs*ckcs+ckss*ckss)
                 engdum=engdum+akk*(summu1*summu1+summu2*summu2)
                 virprs=akv*(summu1*summu1+summu2*summu2)

                 omg(1)=omg(1)-virprs*rkx3*rkx3
     x                     -2.d0*akk*rkx3*dmuskrsx*summu1
     x                     +2.d0*akk*rkx3*dmuckrsx*summu2
                 omg(5)=omg(5)-virprs*rky3*rky3
     x                     -2.d0*akk*rky3*dmuskrsy*summu1
     x                     +2.d0*akk*rky3*dmuckrsy*summu2
                 omg(9)=omg(9)-virprs*rkz3*rkz3
     x                     -2.d0*akk*rkz3*dmuskrsz*summu1
     x                     +2.d0*akk*rkz3*dmuckrsz*summu2
#ifdef STRESS
                 omg(2)=omg(2)-virprs*rkx3*rky3
     x                     -2.d0*akk*rkx3*dmuskrsy*summu1
     x                     +2.d0*akk*rkx3*dmuckrsy*summu2
                 omg(3)=omg(3)-virprs*rkx3*rkz3
     x                     -2.d0*akk*rkx3*dmuskrsz*summu1
     x                     +2.d0*akk*rkx3*dmuckrsz*summu2
c                 omg(4)=omg(4)-virprs*rkx3*rky3
c     x                     -2.d0*akk*rky3*dmuskrsx*summu1
c     x                     +2.d0*akk*rky3*dmuckrsx*summu2
                 omg(6)=omg(6)-virprs*rky3*rkz3
     x                     -2.d0*akk*rky3*dmuskrsz*summu1
     x                     +2.d0*akk*rky3*dmuckrsz*summu2
c                 omg(7)=omg(7)-virprs*rkx3*rkz3
c     x                     -2.d0*akk*rkz3*dmuskrsx*summu1
c     x                     +2.d0*akk*rkz3*dmuckrsx*summu2
c                 omg(8)=omg(8)-virprs*rky3*rkz3
c     x                     -2.d0*akk*rkz3*dmuskrsy*summu1
c     x                     +2.d0*akk*rkz3*dmuckrsy*summu2
#endif
c     
c     calculate force on each site
                    
                    i=0
                    do j=iatm0,iatm1
                      
                      i=i+1

                      djmuk=0.d0
                      if (polr2(j).gt.1.d-6) djmuk=dipx(j)*rkx3+
     x                               dipy(j)*rky3+dipz(j)*rkz3

c
c     increment electric field in k-space

                      tt=-skr(i)*ckcs+ckr(i)*ckss

                      efdcrecx(j)=efdcrecx(j)-rkx3*akk*tt
                      efdcrecy(j)=efdcrecy(j)-rky3*akk*tt
                      efdcrecz(j)=efdcrecz(j)-rkz3*akk*tt
c
c     increment dipole field in k-space

                      if (polr2(j).gt.1.d-6 ) then

                         vvmu=skr(i)*skdmuskrs+ckr(i)*skdmuckrs

                         efddmurecx(j)=efddmurecx(j)-rkx3*akk*vvmu
                         efddmurecy(j)=efddmurecy(j)-rky3*akk*vvmu
                         efddmurecz(j)=efddmurecz(j)-rkz3*akk*vvmu

                      endif

                      ssmu=summu1*(-skr(i)*chge(j)-ckr(i)*djmuk)+
     x                     summu2*(-skr(i)*djmuk+ckr(i)*chge(j))

                      fxx(j)=fxx(j)-rkx3*akk*ssmu
                      fyy(j)=fyy(j)-rky3*akk*ssmu
                      fzz(j)=fzz(j)-rkz3*akk*ssmu

#ifdef TTM_FORCE_DECOMPOSITION
!     parmanent part
                      dumtmp = skr(i) * ckcs - ckr(i) * ckss

                      fx_ttm_per(j) = fx_ttm_per(j) +
     $                     rkx3 * akk * chge(j) * dumtmp
                      fy_ttm_per(j) = fy_ttm_per(j) +
     $                     rky3 * akk * chge(j) * dumtmp
                      fz_ttm_per(j) = fz_ttm_per(j) +
     $                     rkz3 * akk * chge(j) * dumtmp

!     induced part ( former 4 terms: D-C,  last two terms: D-D )

                      dumtmp = djmuk * skr(i) * ckss
     $                       + djmuk * ckr(i) * ckcs
     $                       - chge(j) * ckr(i) * skdmuckrs
     $                       - chge(j) * skr(i) * skdmuskrs
     $                       + djmuk * skr(i) * skdmuckrs
     $                       - djmuk * ckr(i) * skdmuskrs  

                      fx_ttm_ind(j) = fx_ttm_ind(j)
     $                     + rkx3 * akk * dumtmp
                      fy_ttm_ind(j) = fy_ttm_ind(j)
     $                     + rky3 * akk * dumtmp
                      fz_ttm_ind(j) = fz_ttm_ind(j)
     $                     + rkz3 * akk * dumtmp
#endif /* TTM_FORCE_DECOMPOSITION */

!      ESP for DMS
                      if( ldms ) then
                         vesp_k(j) = vesp_k(j) + 
     $                        akk * ( ckr(i) * ckcs + skr(i) * ckss )
                         vesp_dc_k(j) = vesp_dc_k(j) + akk *
     $                        ( skr(i)*skdmuckrs - ckr(i)*skdmuskrs )
                      endif

                    enddo
c     
c     end vector loop

                  endif

                  kkk=kkk+12

                endif
                
              endif

            enddo
            
            nmin=-kmax3
            
          enddo
          
          mmin=-kmax2
          
        enddo
c
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
      engdum=engdum/dble(mxnode)
      do i = 1,9
         omg(i) = omg(i)/dble(mxnode)
      enddo

c     
c     add self interaction correction to potential
      
      if(lconsw)then

        eng1 = engdum
        engcpe=2.d0*rvolm*r4pie0*engcpe/epsq +engsic
        scal1 = 2.d0*rvolm*r4pie0/epsq
        scale = 4.d0*rvolm*r4pie0/epsq

      else

        eng1 = engdum
        engcpe=rvolm*r4pie0*engcpe/epsq+engsic
        scal1 = rvolm*r4pie0/epsq
        scale = 2.d0*rvolm*r4pie0/epsq

      endif

      
c     
c     calculate final forces

      do i=iatm0,iatm1
        
        efdcrecx(i)=scale*efdcrecx(i)
        efdcrecy(i)=scale*efdcrecy(i)
        efdcrecz(i)=scale*efdcrecz(i)

        efddmurecx(i)=scale*efddmurecx(i)
        efddmurecy(i)=scale*efddmurecy(i)
        efddmurecz(i)=scale*efddmurecz(i)

        fxx(i)=scale*fxx(i)
        fyy(i)=scale*fyy(i)
        fzz(i)=scale*fzz(i)

!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
        fx_ttm_per(i) = fx_ttm_per(i) * scale
        fy_ttm_per(i) = fy_ttm_per(i) * scale
        fz_ttm_per(i) = fz_ttm_per(i) * scale
        
        fx_ttm_ind(i) = fx_ttm_ind(i) * scale
        fy_ttm_ind(i) = fy_ttm_ind(i) * scale
        fz_ttm_ind(i) = fz_ttm_ind(i) * scale
#endif /* TTM_FORCE_DECOMPOSITION */
!
!     ESP for DMS
          
        if( ldms ) then 
           vesp_k(i)    = vesp_k(i)    * scale 
           vesp_dc_k(i) = vesp_dc_k(i) * scale
        endif

      enddo
c
c     correct self-interactions

      dum=4.d0*r4pie0*alpha**3/sqrpi/3.d0
      do i=iatm0,iatm1

         efddmurecx(i)=efddmurecx(i)+dum*dipx(i)
         efddmurecy(i)=efddmurecy(i)+dum*dipy(i)
         efddmurecz(i)=efddmurecz(i)+dum*dipz(i)

      enddo
c
c     global exchange of electric field

! VB: no need to merge them; the fields are valid for iatm0:iatm1
!      if(mxnode.gt.1) then
!
!
!        call merge(idnode,mxnode,natms,mxbuff,efdcrecx,
!     x             efdcrecy,efdcrecz,buffer)
!
!        call merge(idnode,mxnode,natms,mxbuff,efddmurecx,
!     x             efddmurecy,efddmurecz,buffer)
!
!      endif
!
!     for ESP (self term)

      if (ldms) then
         do i=iatm0,iatm1
            vesp_s(i)=-2.0d0*r4pie0*chge(i)*alpha/sqrpi/epsq
         enddo
      endif ! ldms

#ifdef STRESS
c     
c     calculate stress tensor (symmetrical)

      stress(1) = stress(1)+scal1*(omg(1)+eng1)
      stress(2) = stress(2)+scal1*omg(2)
      stress(3) = stress(3)+scal1*omg(3)
      stress(4) = stress(4)+scal1*omg(2)
c      stress(4) = stress(4)+scal1*omg(4)
      stress(5) = stress(5)+scal1*(omg(5)+eng1)
      stress(6) = stress(6)+scal1*omg(6)
      stress(7) = stress(7)+scal1*omg(3)
c      stress(7) = stress(7)+scal1*omg(7)
      stress(8) = stress(8)+scal1*omg(6)
c      stress(8) = stress(8)+scal1*omg(8)
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
