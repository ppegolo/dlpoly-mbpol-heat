      subroutine xscale
     x  (idnode,mxnode,natms,keyens,imcon,tstep,
     x  cell,xxx,yyy,zzz,weight,eta,buffer)
c     
c***********************************************************************
c     
c     dl_poly routine to scale positions with change in box shape
c     
c     parallel replicated data version
c     
c     copyright daresbury laboratory 1995
c     author t.forester      october 1995
c     
c     wl
c     2000/01/18 14:06:00
c     1.3
c     Exp
c     
c***********************************************************************
c     
#include "dl_params.inc"

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension weight(mxatms)
      dimension eta(9),cell(9)
      dimension buffer(mxbuff)
#ifdef VAMPIR
      call VTBEGIN(65, ierr)
#endif
c     
c     assign block of atoms to processor

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

      if((keyens.eq.4).or.(keyens.eq.6)) then
c     
c     berendsen npt/nst

        do i = iatm0,iatm1

          xa = eta(1)*xxx(i) + eta(2)*yyy(i) + eta(3)*zzz(i)
          ya = eta(4)*xxx(i) + eta(5)*yyy(i) + eta(6)*zzz(i)
          za = eta(7)*xxx(i) + eta(8)*yyy(i) + eta(9)*zzz(i)

          xxx(i) = xa
          yyy(i) = ya
          zzz(i) = za

        enddo

      elseif(keyens.eq.5.or.keyens.eq.7) then
c     
c     hoover npt/nst
        
        totmas = 0.d0
        do i = 1,natms
          totmas = totmas + weight(i)
        enddo
        
        xcmo = 0.d0
        ycmo = 0.d0
        zcmo = 0.d0

        do i = 1,natms
          xcmo = xcmo + weight(i)*xxx(i)
          ycmo = ycmo + weight(i)*yyy(i)
          zcmo = zcmo + weight(i)*zzz(i)
        enddo
        xcmo = xcmo/totmas
        ycmo = ycmo/totmas
        zcmo = zcmo/totmas

        do i = iatm0,iatm1

          xa = xxx(i) - xcmo
          ya = yyy(i) - ycmo
          za = zzz(i) - zcmo

          xxx(i)=xxx(i)+tstep*(eta(1)*xa+eta(2)*ya+eta(3)*za)
          yyy(i)=yyy(i)+tstep*(eta(2)*xa+eta(5)*ya+eta(6)*za)
          zzz(i)=zzz(i)+tstep*(eta(3)*xa+eta(6)*ya+eta(9)*za)

        enddo

        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

      endif

      nbuff = mxbuff
      if(mxnode.gt.1)
     x  call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)

#ifdef VAMPIR
      call VTEND(65, ierr)
#endif
      return 
      end
