      subroutine ddfield
     x    (idnode,mxnode,natms,imcon,lentry,list,ilist,rcut,
     x     xxx,yyy,zzz,dipx,dipy,dipz,
     x     emux,emuy,emuz,polr,alpha)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating induced-dipole field
c
c     parallel replicated data version
c
c     copyright - voth group 2003
c     author    - c. j. burnham  sept 2003
c
c***********************************************************************
c
#include "dl_params.inc"

      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension lentry(msatms),list(msatms,mxlist)
      dimension ilist(mxxdf)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension polr(mxatms)

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)
c
#ifdef VAMPIR
      call VTBEGIN(100, ierr)
#endif
c
c     evaluate the emu fields from the dipole tensor

      do i=1,natms

         emux(i)=0.d0
         emuy(i)=0.d0
         emuz(i)=0.d0

      enddo

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

	
	if (ii.eq.7603) then
	write(nrite,*)'7603  ',lentry(ii)
	do jj=1,lentry(ii)
	r2 = xdf(jj)**2+ydf(jj)**2+zdf(jj)**2
        rr=sqrt(r2)
	write(nrite,'2i8,3x,f12.6')jj,list(ii,jj),rr
	enddo
	endif

        enddo
c
c     periodic boundary conditions

        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)

        do k=1,lentry(ii)
c
c     square of distances

          r2 = xdf(k)**2+ydf(k)**2+zdf(k)**2
          rr = sqrt(r2)
	if (rr.lt.1.d-1) then
	write(nrite,*)'ddf',i,ilist(k),rr
	stop
	endif

          j=ilist(k)
c
c     apply truncation of potential

          if (rr.lt.rcut) then

             if(polr(i).gt.1.d-6 .and. polr(j).gt.1.d-6) then

                rri = 1./rr
                rsqi = 1./r2
                r3i = 1./(rr**3)
                r5i = r3i*rsqi

c
c     calculate the 'Smith' B (error) functions and gradients
c     for real space Ewald sum.  (I've modified
c     Smith's recursion to give -erf, not erfc functions.)
c     comment by c. j. burnham


                a=alpha*rr
                bbb(0)=(erfcc(a)-1.)*rri
                exp2a=exp(-a*a)
 
                do m=1,3

                   fm=dble(m)
                   bbb(m)=rsqi*((2.0*fm-1.0)*bbb(m-1)
     x              +((2.0*alpha**2.0)**fm)*exp2a/alpha/sqrpi)

                enddo
c
c     assume same screening functions for chgs and dips

                do m=0,3

                   chgchg(m)=bbb(m)
                   chgdip(m)=bbb(m)
                   dipdip(m)=bbb(m)

                enddo
c
c     dipole-dipole field

                sdd=1.

                rirj = xdf(k)*xdf(k)
                dtens1 = r4pie0*((3.*rirj*r5i-r3i)*sdd
     x                           + rirj*dipdip(2) - dipdip(1))
                rirj = xdf(k)*ydf(k)
                dtens2 = r4pie0*((3.*rirj*r5i)*sdd
     x                           + rirj*dipdip(2))
                rirj = xdf(k)*zdf(k)
                dtens3 = r4pie0*((3.*rirj*r5i)*sdd
     x                           + rirj*dipdip(2))
                rirj = ydf(k)*ydf(k)
                dtens4 = r4pie0*((3.*rirj*r5i-r3i)*sdd
     x                           + rirj*dipdip(2) - dipdip(1))
                rirj = ydf(k)*zdf(k)
                dtens5 = r4pie0*((3.*rirj*r5i)*sdd
     x                           + rirj*dipdip(2))
                rirj = zdf(k)*zdf(k)
                dtens6 = r4pie0*((3.*rirj*r5i-r3i)*sdd
     x                           + rirj*dipdip(2) - dipdip(1))

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

      return
      end

