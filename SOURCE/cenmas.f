      subroutine cenmas
     x  (idnode,imcon,mxnode,ntbond,listbnd,cell,
     x  natms,ntpmls,nummols,numsit,wgtsit,numbonds,nums,
     x  listbnd2,listbnd3,listin2,tstep,xmsd,ymsd,zmsd,
     x  xxx,yyy,zzz,vxx,vyy,vzz,xxx1,yyy1,zzz1,xdab,ydab,zdab,
     x  xdab2,ydab2,zdab2,xcm,ycm,zcm,vxcm,vycm,vzcm,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating c.o.m coordinates and velocities
c
c     t. yan 
c
c     the voth group
c
c     feb. 2004
c     
c     replicated data - blocked  data version
c
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      dimension nummols(mxtmls),numsit(mxtmls),numbonds(mxtmls)
      dimension wgtsit(mxsite)
      dimension listbnd(mxbond,3)
      dimension listin2(4*mxatms),listbnd2(mxatms),listbnd3(mxatms) 
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension xxx1(mxatms),yyy1(mxatms),zzz1(mxatms)
      dimension xcm(mxatms),ycm(mxatms),zcm(mxatms)
      dimension vxcm(mxatms),vycm(mxatms),vzcm(mxatms)
      dimension xmsd(mxatms),ymsd(mxatms),zmsd(mxatms)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension xdab2(mxatms),ydab2(mxatms),zdab2(mxatms)
      dimension buffer(mxbuff)
      
#ifdef VAMPIR
      call VTBEGIN(21, ierr)
#endif
c
c     check size of work arrays

      if((ntbond-mxnode+1)/mxnode.gt.msbad) call error(idnode,418)
c
c     initialize

      do i=1,natms

         xxx1(i)=xxx(i)
         yyy1(i)=yyy(i)
         zzz1(i)=zzz(i)

         xcm(i)=0.d0
         ycm(i)=0.d0
         zcm(i)=0.d0

         vxcm(i)=0.d0
         vycm(i)=0.d0
         vzcm(i)=0.d0

         xdab2(i)=0.d0
         ydab2(i)=0.d0
         zdab2(i)=0.d0

         listbnd2(i)=0
         listbnd3(i)=0

      enddo
c     
c     block indices

      ibnd1 = (idnode*ntbond)/mxnode + 1
      ibnd2 = ((idnode+1)*ntbond)/mxnode
c     
c     calculate atom separation vectors
      
      ii=0
      do i=ibnd1,ibnd2
        
        ii=ii+1
c     
c     indices of atoms involved
        
        ia=listbnd(ii,2)
        ib=listbnd(ii,3)

        listbnd2(i)=ia
        listbnd3(i)=ib
c     
c     components of bond vector
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

      enddo
c     
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)

      ii=0
      do i=ibnd1,ibnd2

        ii=ii+1
c
c     components of bond vector

        xdab2(i)=xdab(ii)
        ydab2(i)=ydab(ii)
        zdab2(i)=zdab(ii)

      enddo

c      if(mxnode.gt.1) then
c        call merge(idnode,mxnode,ntbond,mxbuff,xdab2,
c     x             ydab2,zdab2,buffer)
c        call merge(idnode,mxnode,ntbond,mxbuff,xdab2,
c     x             listbnd2,listbnd3,buffer)
c      endif

      if(mxnode.gt.1) then

        j=0
        k=0
        do i=1,ntbond

          buffer(j+1)=xdab2(i)
          buffer(j+2)=ydab2(i)
          buffer(j+3)=zdab2(i)
          listin2(k+1)=listbnd2(i)
          listin2(k+2)=listbnd3(i)
          j=j+3
          k=k+2

        enddo

        call gdsum(buffer(1),3*ntbond,buffer(3*ntbond+1))

        call gisum(listin2(1),2*ntbond,listin2(2*ntbond+1))

        j=0
        k=0
        do i=1,ntbond

          xdab2(i)=buffer(j+1)
          ydab2(i)=buffer(j+2)
          zdab2(i)=buffer(j+3)
          listbnd2(i)=listin2(k+1)
          listbnd3(i)=listin2(k+2)
          j=j+3
          k=k+2

        enddo

      endif

      do i=1,ntbond
c
c     indices of atoms involved

        ia=listbnd2(i)
        ib=listbnd3(i)
c
c     correct for pbc

        xxx1(ib)=xxx1(ia)-xdab2(i)
        yyy1(ib)=yyy1(ia)-ydab2(i)
        zzz1(ib)=zzz1(ia)-zdab2(i)

      enddo
c     
c     c.m. coordinates

      ndd=0
      nss=0
      nums=0

      do i=1,ntpmls

        im1=idnode*nummols(i)/mxnode+1
        im2=(idnode+1)*nummols(i)/mxnode

        do j=im1,im2

           jj=nums+j

           xcm(jj)=0.d0
           ycm(jj)=0.d0
           zcm(jj)=0.d0

           wtt=0.d0

           do k=1,numsit(i)

              idum=ndd+(j-1)*numsit(i)+k
              kdum=nss+k

              wtt=wtt+wgtsit(kdum)

              xcm(jj)=xcm(jj)+wgtsit(kdum)*xxx1(idum)
              ycm(jj)=ycm(jj)+wgtsit(kdum)*yyy1(idum)
              zcm(jj)=zcm(jj)+wgtsit(kdum)*zzz1(idum)

              vxcm(jj)=vxcm(jj)+wgtsit(kdum)*vxx(idum)
              vycm(jj)=vycm(jj)+wgtsit(kdum)*vyy(idum)
              vzcm(jj)=vzcm(jj)+wgtsit(kdum)*vzz(idum)

           enddo

           xcm(jj)=xcm(jj)/wtt
           ycm(jj)=ycm(jj)/wtt
           zcm(jj)=zcm(jj)/wtt

           vxcm(jj)=vxcm(jj)/wtt
           vycm(jj)=vycm(jj)/wtt
           vzcm(jj)=vzcm(jj)/wtt

        enddo 

        ndd=ndd+nummols(i)*numsit(i)
        nss=nss+numsit(i)
        nums=nums+nummols(i)

      enddo

c     global exchange of c.m. coordinates and velocities

      if(mxnode.gt.1) then

c        call merge(idnode,mxnode,nums,mxbuff,xcm,ycm,zcm,buffer)
c        call merge(idnode,mxnode,nums,mxbuff,vxcm,vycm,vzcm,buffer)

          j=0
          do i=1,nums

            buffer(j+1)=xcm(i)
            buffer(j+2)=ycm(i)
            buffer(j+3)=zcm(i)
            buffer(j+4)=vxcm(i)
            buffer(j+5)=vycm(i)
            buffer(j+6)=vzcm(i)
            j=j+6

          enddo

          call gdsum(buffer(1),6*nums,buffer(6*nums+1))

          j=0
          do i=1,nums

            xcm(i)=buffer(j+1)
            ycm(i)=buffer(j+2)
            zcm(i)=buffer(j+3)
            vxcm(i)=buffer(j+4)
            vycm(i)=buffer(j+5)
            vzcm(i)=buffer(j+6)
            j=j+6

          enddo

      endif
c
c     accumulate center-of-mass displacements

      inum0 = (idnode*nums)/mxnode + 1
      inum1 = ((idnode+1)*nums)/mxnode
      do i=inum0,inum1
        xmsd(i)=xmsd(i)+vxcm(i)*tstep
        ymsd(i)=ymsd(i)+vycm(i)*tstep
        zmsd(i)=zmsd(i)+vzcm(i)*tstep
      enddo
      if (mxnode.gt.1) then
        nbuff=mxbuff
        call merge(idnode,mxnode,nums,nbuff,xmsd,ymsd,zmsd,buffer)
      endif

#ifdef VAMPIR
      call VTEND(21, ierr)
#endif
      return
      end


