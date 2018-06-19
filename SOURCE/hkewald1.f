      subroutine hkewald1
     x  (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x  engcpe,vircpe,alpha,epsq,cell,ahk,chge,xxx,
     x  yyy,zzz,fxx,fyy,fzz,elc,emc,els,ems,ckc,cks,
     x  stress,crn,pp,znp,zgc,zgs,buffer)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using Hautman Klein Ewald method
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith february 2000
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c     note - in loop over all k vectors k=2pi(ll/cl,mm/cl)
c     the values of ll and mm are selected so that the symmetry of
c     reciprocal lattice is taken into account i.e. the following
c     rules apply.
c     
c     ll ranges over the values 0 to kmax1 only.
c     
c     mm ranges over 1 to kmax2 when ll=0 and over
c     -kmax2 to kmax2 otherwise.
c     
c     wl
c     2001/06/12 13:10:08
c     1.1
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"
      
      logical newjob,safe
      
      dimension elc(mxewld,0:1),els(mxewld,0:1)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb)
      dimension buffer(mxbuff),cprop(10),ahk(0:mxhko)
      dimension chge(mxatms),cell(9),rcell(9),stress(9),omg(9)
      dimension ckc(mxewld),cks(mxewld),crn(0:mxhko,0:mxhko)
      dimension pp(2*mxhko),znp(mxhke,0:2*mxhko)
      dimension zgc(0:2*mxhko),zgs(0:2*mxhko)

      save newjob,engsic

      data newjob/.true./,safe/.true./

#ifdef VAMPIR
      call VTBEGIN(171, ierr)
#endif
c     
c     initialise coulombic potential energy

      engcpe=0.d0
      vircpe=0.d0
      if(alpha.lt.1.d-8)return
c     
c     set working parameters

      twopi=2.d0*pi
      ralph=0.5d0/alpha
      call dcell(cell,cprop)
      area=cprop(1)*cprop(2)*sqrt(1.d0-cprop(4)**2)
      rarea=pi/area
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
c     construct reciprocal lattice vectors and set k vector range

      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      call dcell(rcell,cprop)
      rcpcut=min(dble(kmax1)*cprop(7),dble(kmax2)*cprop(8))*
     x  1.05d0*twopi
      rcpct2=rcpcut**2
c     
c     compute quantities for first entry

      if(newjob)then

        newjob=.false.
c     
c     pbc check and array bound checks

        if(imcon.ne.6)call error(idnode,66)
        if(mxhke.ne.msatms) call error(idnode,331)
        if(mxewld.ne.msatms) call error(idnode,330)
c
c     check hk screening function at cutoff

c$$$        aaa=cerfr(ralph,rcpcut)
        aaa=erfc(ralph*rcpcut)/rcpcut
        if(aaa.gt.1.d-4)then

          call warning(idnode,105,aaa,0.d0,0.d0)
c          call error(idnode,487)
          
        endif
c     
c     calculate self interaction correction
        
        engsic=0.d0
        
        do i=iatm0,iatm1
          
          engsic=engsic+chge(i)**2
          
        enddo

        engsic=-r4pie0*alpha*engsic/(sqrpi*epsq)
c     
c     binomial coefficients

        k=0
        crn(0,0)=0.5d0
        do i=1,2*nhko

          pp(i)=1.d0
          pm1=pp(1)

          do j=2,i

            pm2=pp(j)
            pp(j)=pm2+pm1
            pm1=pm2

          enddo

          if(mod(i,2).eq.0)then

            k=k+1
            do j=0,k

              term=pp(j+1)*(-1.d0)**j
              crn(j,k)=term
              crn(k,j)=term

            enddo

            crn(k,k)=0.5d0*crn(k,k)

          endif

        enddo

      endif
c     
c     calculate and store powers of z_i
      
      i=0

      do j=iatm0,iatm1

        i=i+1
        znp(i,0)=1.d0
        znp(i,1)=zzz(j)

      enddo

      limit=i

      do k=2,2*nhko

        do i=1,limit

          znp(i,k)=znp(i,k-1)*znp(i,1)

        enddo

      enddo
c     
c     calculate and store exponential factors

      i=0

      do j=iatm0,iatm1
        
        i=i+1
        elc(i,0)=1.d0
        emc(i,0)=1.d0
        els(i,0)=0.d0
        ems(i,0)=0.d0
        ssx=rcell(1)*xxx(j)+rcell(4)*yyy(j)
        ssy=rcell(2)*xxx(j)+rcell(5)*yyy(j)
        elc(i,1)=cos(twopi*ssx)
        emc(i,1)=cos(twopi*ssy)
        els(i,1)=sin(twopi*ssx)
        ems(i,1)=sin(twopi*ssy)
        
      enddo

      do l=2,kmax2
        
        do i=1,limit
          
          emc(i,l)=emc(i,l-1)*emc(i,1)-ems(i,l-1)*ems(i,1)
          ems(i,l)=ems(i,l-1)*emc(i,1)+emc(i,l-1)*ems(i,1)
          
        enddo
        
      enddo
c     
c     start of main loop over k vectors

      mmin=1
      
      do ll=0,kmax1
        
        l=ll
        tmp = twopi*dble(ll)
        rkx1=tmp*rcell(1)
        rky1=tmp*rcell(4)
c     
c     put cos(i,L) terms into cos(i,0) array

        if(l.eq.1) then
          
          do i=1,limit

            elc(i,0)=elc(i,1)
            els(i,0)=els(i,1)
            
          enddo

        elseif(l.gt.1) then

          do i=1,limit

            cs=elc(i,0)
            elc(i,0)=cs*elc(i,1)-els(i,0)*els(i,1)
            els(i,0)=els(i,0)*elc(i,1)+cs*els(i,1)
            
          enddo
          
        endif

        do mm=mmin,kmax2
          
          m=iabs(mm)
          tmp = twopi*dble(mm)
          rkx2=rkx1+tmp*rcell(2)
          rky2=rky1+tmp*rcell(5)
c     
c     test on magnitude of k vector
          
          rksq=rkx2*rkx2+rky2*rky2

          if(rksq.le.rcpct2)then
c     
c     calculate exp(ikr) terms and product with charges
            
            i=0

            if(mm.ge.0)then
              
              do j=iatm0,iatm1
                
                i=i+1
                ckc(i)=chge(j)*(elc(i,0)*emc(i,m)-els(i,0)*ems(i,m))
                cks(i)=chge(j)*(els(i,0)*emc(i,m)+ems(i,m)*elc(i,0))
                
              enddo
              
            else
              
              do j=iatm0,iatm1
                
                i=i+1
                
                ckc(i)=chge(j)*(elc(i,0)*emc(i,m)+els(i,0)*ems(i,m))
                cks(i)=chge(j)*(els(i,0)*emc(i,m)-ems(i,m)*elc(i,0))
                
              enddo
              
            endif
c     
c     calculate sum of products of powers of z_i and q_i exp(ik.s_i)
            
            do k=0,2*nhko
              
              zgc(k)=0.d0
              zgs(k)=0.d0

              do i=1,limit

                zgc(k)=zgc(k)+ckc(i)*znp(i,k)
                zgs(k)=zgs(k)+cks(i)*znp(i,k)

              enddo

            enddo
c     
c     perform global summation of zgc and zgs arrays
            
            if(mxnode.gt.1)then
              
              call gdsum(zgc(0),2*nhko+1,buffer)
              call gdsum(zgs(0),2*nhko+1,buffer)
              
            endif
c     
c     calculate 0th order screening function

            rkk=sqrt(rksq)
c$$$            fn0=cerfr(ralph,rkk)
            fn0=erfc(ralph*rkk)/rkk
            gaus=exp(-(ralph*rkk)**2)/(alpha*sqrpi)
c     
c     sum terms for orders of the screening function

            fac=1.d0

            do k=0,nhko
c     
c     sum over z_i binomial contributions

              eterm=0.d0
              fng=fac*fn0
              do m=0,k

                n=2*k-m
c
c     sum energy terms

                eterm=eterm+crn(m,k)*(zgc(m)*zgc(n)+zgs(m)*zgs(n))
c     
c     calculate force contribution to each site
                
                i=0
                bkk=-fng*crn(m,k)

                do j=iatm0,iatm1
                  
                  i=i+1
                  force0=bkk*(znp(i,n)*(zgs(m)*ckc(i)-zgc(m)*cks(i))+
     x              znp(i,m)*(zgs(n)*ckc(i)-zgc(n)*cks(i)))
                  fxx(j)=fxx(j)+rkx2*force0
                  fyy(j)=fyy(j)+rky2*force0
#ifdef STRESS
                  omg(3)=omg(3)+rkx2*force0*zzz(j)
                  omg(6)=omg(6)+rky2*force0*zzz(j)
#endif
                  if(k.gt.0)then

                    if(m.eq.0)then
                      
                      forcez=bkk*dble(n)*znp(i,n-1)*(zgc(m)*ckc(i)+
     x                  zgs(m)*cks(i))
                      
                    else
                      
                      forcez=bkk*(dble(m)*znp(i,m-1)*(zgc(n)*ckc(i)+
     x                  zgs(n)*cks(i))+dble(n)*znp(i,n-1)*(zgc(m)*
     x                  ckc(i)+zgs(m)*cks(i)))
                    
                    endif
                    
                    omg(9)=omg(9)+forcez*zzz(j)
                    fzz(j)=fzz(j)+forcez

                  endif

                enddo

              enddo
c     
c     accumulate potential energy and stress tensor
              
              engcpe=engcpe+fng*eterm
              pterm=(dble(2*k-1)*fng-fac*gaus)/rksq
              omg(1)=omg(1)+eterm*(fng+pterm*rkx2*rkx2)
              omg(5)=omg(5)+eterm*(fng+pterm*rky2*rky2)
#ifdef STRESS
              omg(2)=omg(2)+eterm*pterm*rky2*rkx2
#endif
              fac=fac*rksq/(dble(2*(k+1))*dble(2*k+1))
c     
c     end of loop over orders of screening function

            enddo
c
c     end of if-block for  rksq <  rcpct2

          endif
c
c     end of inner loop over reciprocal lattice vectors

        enddo
        
        mmin=-kmax2
c
c     end of outer loop over reciprocal lattice vectors

      enddo
      
      engcpe=engcpe/dble(mxnode)
      do i = 1,9
        omg(i) = omg(i)/dble(mxnode)
      enddo
c     
c     add self interaction correction to potential

      scale=4.d0*rarea*r4pie0/epsq
      engcpe=scale*engcpe+engsic
c     
c     virial term

      vircpe=-scale*(omg(1)+omg(5)+omg(9))
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

      stress(1) = stress(1)+scale*omg(1)
      stress(2) = stress(2)+scale*omg(2)
      stress(3) = stress(3)+scale*omg(3)
      stress(4) = stress(4)+scale*omg(2)
      stress(5) = stress(5)+scale*omg(5)
      stress(6) = stress(6)+scale*omg(6)
      stress(7) = stress(7)+scale*omg(3)
      stress(8) = stress(8)+scale*omg(6)
      stress(9) = stress(9)+scale*omg(9)
#endif

#ifdef VAMPIR
      call VTEND(171, ierr)
#endif
      return
      end
