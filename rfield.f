#include "assert.h"
       subroutine rfield
     x  (idnode,imcon,mxnode,natms,lentry,list,ilist,rcut,epsq,
     x  nexatm,lexatm,jlist,chge,polr,polr2,xxx,yyy,zzz,alpha,cell,
     x  buffer,efieldkx,efieldky,efieldkz,ercp,drewd,eps2,
     x  dipx,dipy,dipz,emux,emuy,emuz,iloop,lexatm2,nexatm2,
     x  lthole,athole,athole12,athole13,
     x  athole_ion,ithole,n_ions,athole_ionwat,
     x  nthole,lttm,nttm2,listttm2,
     x  lads,ascd,n_water)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating electric field in a
c     periodic system using ewald's method
c     
c     parallel replicated data version (part 2)
c     
c     copyright - voth group
c     author    - c. j. burnham
c     parallel  - t. yan
c     
c     electric field and dipole tensor in r space
c     
c     2003/10/01 19:47:27
c     
c***********************************************************************
c     

#include "dl_params.inc"

      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension nexatm2(msatms),lexatm2(msatms,mxexcl)
      dimension ilist(mxxdf),jlist(mxxdf)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension chge(mxatms),polr(mxatms),polr2(mxatms),cell(9)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension ercp(mxegrd,0:3)
      dimension buffer(mxbuff)
      dimension listttm2(mxatms)

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)

      logical lthole,lttm,lads

!FP_fix_start
      real(8) athole_ion,athole_ionwat
      integer n_water, n_water_atoms, no_water_atoms, n_ions, ithole
!FP_fix_end

#ifdef VAMPIR
      call VTBEGIN(100, ierr)
#endif

c     
c     outer loop over atoms

      ii=0

      do i=idnode+1,natms,mxnode

        ii=ii+1
c
c     calculate interatomic distances

        do k=1,lentry(ii)

          j=list(ii,k)
          ilist(k) = j

          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)

        enddo
c
c     periodic boundary conditions

        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
c
c     square of distances

        do k=1,lentry(ii)

          r2 = xdf(k)**2+ydf(k)**2+zdf(k)**2
          rr = dsqrt(r2)
c
c     apply truncation of potential

          if (rcut.gt.rr) then
c
c     calculate the 'Smith' B (error) functions and gradients
c     for real space Ewald sum.  (I've modified
c     Smith's recursion to give -erf, not erfc functions.)
c     comment by c. j. burnham
c     modified to use 3-point interpolation by t. yan

             ll = int(rr/drewd)
             l1=ll+1
             l2=ll+2
             ppp = rr/drewd - dble(ll)
c
c     calculate interaction energy using 3-point interpolation

             do n=0,3

                vk = ercp(ll,n)
                vk1 = ercp(l1,n)
                vk2 = ercp(l2,n)

                t1 = vk + (vk1-vk)*ppp
                t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)

                bbb(n) = t1 + (t2-t1)*ppp*0.5d0

             enddo

             rri = 1.d0/rr
             rsqi = 1.d0/r2
             r3i = 1.d0/(rr**3)
             r5i = r3i*rsqi
c
c     no screening or switching functions for charge and dipole

             screencc=1.d0
             screendc=1.d0
             screendd=1.d0

c
c     assume same screening functions for chgs and dips

             do m=0,3

                chgchg(m)=bbb(m)
                chgdip(m)=bbb(m)
                dipdip(m)=bbb(m)

             enddo
            
             scc=screencc
             sdc=1.d0

             j=ilist(k)

             if (iloop.eq.1) then
c
c     charge field

                eij_ki = r4pie0*chge(j)*(chgdip(1)+sdc*r3i)
                eij_kj = r4pie0*chge(i)*(chgdip(1)+sdc*r3i)
c
c     here is the all dipole smeared option, i.e. lads=.true.

                if (lthole .and. lads) then

!                  write(5001,'(2i8,2f15.5)') i,j, polr(i), polr(j)

                  bb=(polr(i)*polr(j))**(1.d0/6.d0)

                  if (bb.lt.1.d-6) call fixme(__FILE__,__LINE__)

                  if (nthole.eq.3) then

                    arob3=ascd*(rr/bb)**3
                    dum1=dexp(-arob3)

                  else if (nthole.eq.4) then

                    arob4=ascd*(rr/bb)**4
                    dum1=dexp(-arob4)

                  else

                    arob=ascd*rr/bb
                    arob2=arob*arob
                    dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)

                  endif

                  eij_ki = eij_ki-r4pie0*chge(j)*dum1*r3i
                  eij_kj = eij_kj-r4pie0*chge(i)*dum1*r3i

                endif
 
                efieldkx(i) = efieldkx(i) + eij_ki*xdf(k)
                efieldky(i) = efieldky(i) + eij_ki*ydf(k)
                efieldkz(i) = efieldkz(i) + eij_ki*zdf(k)

                efieldkx(j) = efieldkx(j) - eij_kj*xdf(k)
                efieldky(j) = efieldky(j) - eij_kj*ydf(k)
                efieldkz(j) = efieldkz(j) - eij_kj*zdf(k)

             endif ! iloop.eq.1
c
c     dipole-dipole field
c
             if(polr2(i).gt.1.d-6 .and. polr2(j).gt.1.d-6) then

                sdd=1.d0

                rirj1 = xdf(k)*xdf(k)
                dtens1 = r4pie0*((3.d0*rirj1*r5i-r3i)*sdd
     x                           + rirj1*dipdip(2) - dipdip(1))
                rirj2 = xdf(k)*ydf(k)
                dtens2 = r4pie0*((3.d0*rirj2*r5i)*sdd
     x                           + rirj2*dipdip(2))
                rirj3 = xdf(k)*zdf(k)
                dtens3 = r4pie0*((3.d0*rirj3*r5i)*sdd
     x                           + rirj3*dipdip(2))
                rirj4 = ydf(k)*ydf(k)
                dtens4 = r4pie0*((3.d0*rirj4*r5i-r3i)*sdd
     x                           + rirj4*dipdip(2) - dipdip(1))
                rirj5 = ydf(k)*zdf(k)
                dtens5 = r4pie0*((3.d0*rirj5*r5i)*sdd
     x                           + rirj5*dipdip(2))
                rirj6 = zdf(k)*zdf(k)
                dtens6 = r4pie0*((3.d0*rirj6*r5i-r3i)*sdd
     x                           + rirj6*dipdip(2) - dipdip(1))
c
c      here is the all dipole smeared option, i.e. lads=.true.

                if (lthole .and. lads) then

                  bb=(polr(i)*polr(j))**(1.d0/6.d0)

                  if (nthole.eq.3) then
                    arob3=athole*(rr/bb)**3
                    dum1=dexp(-arob3)
                    dum2=(1.d0+arob3)*dum1
                  else if (nthole.eq.4) then
!FP_fix_start: different Thole damping and smearing for ion-ion.
                    if ((i.le.n_ions) .and. (j.le.n_ions)) then
                      if (ithole.eq.3) then
                        arob3=athole_ion*(rr/bb)**3
                        dum1=dexp(-arob3)
                        dum2=(1.d0+arob3)*dum1
c                       write(911,'(4i8,5f12.5)') 
c    x                        i, j, n_ions, ithole,
c    x                        polr(i), polr(j), athole_ion
                      else
                        arob4=athole_ion*(rr/bb)**4
                        dum1=dexp(-arob4)
                        dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
c                       write(912,'(4i8,5f12.5)') 
c    x                        i, j, n_ions, ithole,
c    x                        polr(i), polr(j), athole_ion
                      end if
c  for ion-water 
                    elseif (((i.le.n_ions) .and. (j.gt.n_ions)).or.
     x                       ((i.gt.n_ions) .and. (j.le.n_ions) )) then
                        arob4=athole_ionwat*(rr/bb)**4
                        dum1=dexp(-arob4)
                        dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
c                      write(2550,'(2i8,5f12.5)') 
c     x                      i, j, polr(i), polr(j), athole_ionwat
c                       write(2551,'(2i8,5f12.5)') 
c     x                      i, j, arob4,dum1,dum2,athole_ionwat

                    else
                     arob4=athole*(rr/bb)**4
                     dum1=dexp(-arob4)
                     dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
c                       write(2550,'(2i8,5f12.5)') 
c     x                      i, j, arob4,dum1,dum2
c                      write(2550,'(2i8,5f12.5)') 
c     x                      i, j, polr(i), polr(j), athole 
                   end if
!FP_fix_end: different Thole damping and smearing for ion-ion.
                  else
                    arob=athole*rr/bb
                    arob2=arob*arob
                    arob3=arob2*arob
                    dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)
                    dum2=dum1+dexp(-arob)*arob3/6.d0
                  endif

                  dtens1=dtens1-(3.d0*rirj1*r5i*dum2-dum1*r3i)*r4pie0
                  dtens2=dtens2-3.d0*rirj2*r5i*dum2*r4pie0
                  dtens3=dtens3-3.d0*rirj3*r5i*dum2*r4pie0
                  dtens4=dtens4-(3.d0*rirj4*r5i*dum2-dum1*r3i)*r4pie0
                  dtens5=dtens5-3.d0*rirj5*r5i*dum2*r4pie0
                  dtens6=dtens6-(3.d0*rirj6*r5i*dum2-dum1*r3i)*r4pie0

                endif

                dmui1=dipx(i)
                dmui2=dipy(i)
                dmui3=dipz(i)
                dmuj1=dipx(j)
                dmuj2=dipy(j)
                dmuj3=dipz(j)

                emux(i)=emux(i)+dtens1*dmuj1
     x                         +dtens2*dmuj2
     x                         +dtens3*dmuj3
                emuy(i)=emuy(i)+dtens2*dmuj1
     x                         +dtens4*dmuj2
     x                         +dtens5*dmuj3
                emuz(i)=emuz(i)+dtens3*dmuj1
     x                         +dtens5*dmuj2
     x                         +dtens6*dmuj3

                emux(j)=emux(j)+dtens1*dmui1
     x                         +dtens2*dmui2
     x                         +dtens3*dmui3
                emuy(j)=emuy(j)+dtens2*dmui1
     x                         +dtens4*dmui2
     x                         +dtens5*dmui3
                emuz(j)=emuz(j)+dtens3*dmui1
     x                         +dtens5*dmui2
     x                         +dtens6*dmui3

             endif

          endif

        enddo

      enddo
c
c     include intramolecular field
c
      ii=0

      do i=idnode+1,natms,mxnode

        ii=ii+1
c
c     calculate interatomic distances

        do k=1,nexatm(ii)

          j=lexatm(ii,k)
          jlist(k) = j

          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)

        enddo
c
c     periodic boundary conditions

        call images(imcon,0,1,nexatm(ii),cell,xdf,ydf,zdf)
c
c     square of distances

        do k=1,nexatm(ii)

          r2 = xdf(k)**2+ydf(k)**2+zdf(k)**2
          rr = dsqrt(r2)

          rri = 1.d0/rr
          rsqi = 1.d0/r2
          rr3 = rr*r2
          r3i = 1.d0/rr3
          r5i = r3i*rsqi
c
c     calculate the 'Smith' B (error) functions and gradients
c     for real space Ewald sum.  (I've modified
c     Smith's recursion to give -erf, not erfc functions.)
c     comment by c. j. burnham


          a=alpha*rr
          bbb(0)=(erfcc(a)-1.d0)*rri
          exp2a=exp(-a*a)

          do m=1,3

             fm=dble(m)
             bbb(m)=rsqi*((2.0*fm-1.d0)*bbb(m-1)
     x        +((2.0*alpha**2.0)**fm)*exp2a/alpha/sqrpi)

          enddo
c
c     assume same screening functions for chgs and dips

             screencc=1.d0
             screendc=1.d0
             screendd=1.d0

          do m=0,3

             chgchg(m)=bbb(m)
             chgdip(m)=bbb(m)
             dipdip(m)=bbb(m)

          enddo

          scc=screencc
          sdc=1.d0

          j=jlist(k)

          if (iloop.eq.1) then
c
c     charge field

             eij_ki = r4pie0*chge(j)*chgdip(1)/epsq
             eij_kj = r4pie0*chge(i)*chgdip(1)/epsq

             efieldkx(i) = efieldkx(i) + eij_ki*xdf(k)
             efieldky(i) = efieldky(i) + eij_ki*ydf(k)
             efieldkz(i) = efieldkz(i) + eij_ki*zdf(k)

             efieldkx(j) = efieldkx(j) - eij_kj*xdf(k)
             efieldky(j) = efieldky(j) - eij_kj*ydf(k)
             efieldkz(j) = efieldkz(j) - eij_kj*zdf(k)

          endif
c
c     add intramolecular dipole-dipole terms (not computed above)
c
          if (polr2(i).gt.1.d-6 .and. polr2(j).gt.1.d-6) then
             rirj1 = xdf(k)*xdf(k)
             dtens1 = rirj1*dipdip(2) - dipdip(1)

             rirj2 = xdf(k)*ydf(k)
             dtens2 = rirj2*dipdip(2)

             rirj3 = xdf(k)*zdf(k)
             dtens3 = rirj3*dipdip(2)

             rirj4 = ydf(k)*ydf(k)
             dtens4 = rirj4*dipdip(2) - dipdip(1)

             rirj5 = ydf(k)*zdf(k)
             dtens5 = rirj5*dipdip(2)

             rirj6 = zdf(k)*zdf(k)
             dtens6 = rirj6*dipdip(2) - dipdip(1)
c
c    Applequist's point dipole model

             dtens1=dtens1+3.d0*rirj1*r5i-r3i
             dtens2=dtens2+3.d0*rirj2*r5i
             dtens3=dtens3+3.d0*rirj3*r5i
             dtens4=dtens4+3.d0*rirj4*r5i-r3i
             dtens5=dtens5+3.d0*rirj5*r5i
             dtens6=dtens6+3.d0*rirj6*r5i-r3i

             if (lthole) then
c
c     Thole's modified model
!FP_fix_start
               n_water_atoms = 4 * n_water
               no_water_atoms = natms - n_water_atoms

               bb=(polr(i)*polr(j))**(1.d0/6.d0)

               athole_intra=athole13
               if (mod((i-no_water_atoms)-1,3) .eq. 0
     x          .or. mod((j-no_water_atoms)-1,3) .eq. 0)
     x             athole_intra=athole12
c original.
c              athole_intra=athole13
c              if(mod(i-1,3).eq.0.or.mod(j-1,3).eq.0)
c    x             athole_intra=athole12
c
c              write(811,'(2i8,5f12.5)') i, j, polr2(i), polr2(j),
c    x            athole_intra, athole12, athole13
!FP_fix_end

               if (nthole.eq.3) then
                 arob3=athole_intra*(rr/bb)**3
                 dum1=dexp(-arob3)
                 dum2=(1.d0+arob3)*dum1
               else if (nthole.eq.4) then
                 arob4=athole_intra*(rr/bb)**4
                 dum1=dexp(-arob4)
                 dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
               else
                 arob=athole_intra*rr/bb
                 arob2=arob*arob
                 arob3=arob2*arob
                 dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)
                 dum2=dum1+dexp(-arob)*arob3/6.d0
               endif

               dtens1=dtens1-3.d0*rirj1*r5i*dum2+dum1*r3i
               dtens2=dtens2-3.d0*rirj2*r5i*dum2
               dtens3=dtens3-3.d0*rirj3*r5i*dum2
               dtens4=dtens4-3.d0*rirj4*r5i*dum2+dum1*r3i
               dtens5=dtens5-3.d0*rirj5*r5i*dum2
               dtens6=dtens6-3.d0*rirj6*r5i*dum2+dum1*r3i

             endif

             dmui1=dipx(i)*r4pie0/epsq
             dmui2=dipy(i)*r4pie0/epsq
             dmui3=dipz(i)*r4pie0/epsq
             dmuj1=dipx(j)*r4pie0/epsq
             dmuj2=dipy(j)*r4pie0/epsq
             dmuj3=dipz(j)*r4pie0/epsq

             emux(i)=emux(i)+dtens1*dmuj1
     x                      +dtens2*dmuj2
     x                      +dtens3*dmuj3
             emuy(i)=emuy(i)+dtens2*dmuj1
     x                      +dtens4*dmuj2
     x                      +dtens5*dmuj3
             emuz(i)=emuz(i)+dtens3*dmuj1
     x                      +dtens5*dmuj2
     x                      +dtens6*dmuj3

             emux(j)=emux(j)+dtens1*dmui1
     x                      +dtens2*dmui2
     x                      +dtens3*dmui3
             emuy(j)=emuy(j)+dtens2*dmui1
     x                      +dtens4*dmui2
     x                      +dtens5*dmui3
             emuz(j)=emuz(j)+dtens3*dmui1
     x                      +dtens5*dmui2
     x                      +dtens6*dmui3

          endif ! dipole-dipole

        enddo

      enddo

c
c     exclude intramolecular charge field
c

      if (iloop.eq.1) then

        ii=0

        do i=idnode+1,natms,mxnode

          ii=ii+1
c
c     calculate interatomic distances

          do k=1,nexatm2(ii)

            j=lexatm2(ii,k)
            jlist(k) = j

            xdf(k)=xxx(i)-xxx(j)
            ydf(k)=yyy(i)-yyy(j)
            zdf(k)=zzz(i)-zzz(j)

          enddo
c
c     periodic boundary conditions

        call images(imcon,0,1,nexatm2(ii),cell,xdf,ydf,zdf)

c
c     square of distances

          do k=1,nexatm2(ii)

            r2 = xdf(k)**2+ydf(k)**2+zdf(k)**2
            rr = dsqrt(r2)

            r3i = 1.d0/(rr**3)

              j=jlist(k) 
c
c     exclude intramolecular charge field

              eij_ki = -r4pie0*chge(j)*r3i/eps2
              eij_kj = -r4pie0*chge(i)*r3i/eps2

              efieldkx(i) = efieldkx(i) + eij_ki*xdf(k)
              efieldky(i) = efieldky(i) + eij_ki*ydf(k)
              efieldkz(i) = efieldkz(i) + eij_ki*zdf(k)

              efieldkx(j) = efieldkx(j) - eij_kj*xdf(k)
              efieldky(j) = efieldky(j) - eij_kj*ydf(k)
              efieldkz(j) = efieldkz(j) - eij_kj*zdf(k)

          enddo

        enddo

      endif ! iloop.eq.1

c
c     sum contributions to electric field

      if(mxnode.gt.1) then

!VB: reduce only iatm0:iatm1 range (see pdsum.f)
#if 1
         if (iloop.eq.1) then
            call pdsum6(idnode,mxnode,natms,
     x         efieldkx,efieldky,efieldkz,emux,emuy,emuz,mxbuff,buffer)
         else
            call pdsum3(idnode,mxnode,natms,
     x         emux,emuy,emuz,mxbuff,buffer)
         endif ! iloop.eq.1
#else
        if(iloop.eq.1) then 

        j=0
          do i=1,natms

            buffer(j+1)=efieldkx(i)
            buffer(j+2)=efieldky(i)
            buffer(j+3)=efieldkz(i)
            j=j+3

          enddo

          call gdsum(buffer(1),3*natms,buffer(3*natms+1))

          j=0
          do i=1,natms

            efieldkx(i)=buffer(j+1)
            efieldky(i)=buffer(j+2)
            efieldkz(i)=buffer(j+3)
            j=j+3

          enddo

        endif

        j=0
        do i=1,natms

          buffer(j+1)=emux(i)
          buffer(j+2)=emuy(i)
          buffer(j+3)=emuz(i)
          j=j+3

        enddo

        call gdsum(buffer(1),3*natms,buffer(3*natms+1))

        j=0
        do i=1,natms

          emux(i)=buffer(j+1)
          emuy(i)=buffer(j+2)
          emuz(i)=buffer(j+3)
          j=j+3

        enddo
#endif
      endif

#ifdef VAMPIR
      call VTEND(100, ierr)
#endif
      return
      end
