      subroutine ttm2list
     x  (nttm2,ntpmls,nummols,numsit,listttm2,weight)
c
c     construct list for ttm2 water model, c.j.burnham and s.s.xantheas
c     j. chem. phys. 116(2002)5124
c
c     t. yan  june 2004

#include "dl_params.inc"

      dimension nummols(mxtmls),numsit(mxtmls)
      dimension weight(mxatms)
      dimension listttm2(mxatms)
c
c     check

      if (numsit(ntpmls-1).ne.3) then
        write(nrite,*)'TTM2 should be the second to the last type!'
        stop
      endif
      if (numsit(ntpmls).ne.1) then
        write(nrite,*)'For TTM2 water, the last type is the M site!'
        stop
      endif
      if (nummols(ntpmls-1).ne.nummols(ntpmls)) then
        write(nrite,*)'For TTM2 water, the number of water = # M site!'
        stop
      endif
c
c     identify ttm2 water model

      nttm2=0
      nnn=0
      do i=1,ntpmls
        do j=1,nummols(i)
          do k=1,numsit(i)
            nnn=nnn+1
            if (i.eq.ntpmls-1) then
              nttm2=nttm2+1
              listttm2(nttm2)=nnn
            endif
          enddo
        enddo
      enddo
c
c     check

      if (mod(nttm2,3).ne.0) then
        write(nrite,*)'check subroutine ttm2list.f!'
        stop
      endif
c
c     identify oxygen atom

      dum=0.d0
      do i=1,3
         j=listttm2(i)
         if (weight(j).gt.dum) then
            dum=weight(j)
            ioxygen=i
         endif
      enddo
c
c     sort oxygen to the first place of h2o

      if (ioxygen.ne.1) then
         do i=1,nttm2,3
            nh=listttm2(i)
            ndo=ioxygen+i-1
            listttm2(i)=listttm2(ndo)
            listttm2(ndum2)=nh
         enddo
      endif
      
      return
      end
      
