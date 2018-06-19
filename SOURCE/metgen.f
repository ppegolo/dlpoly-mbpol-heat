      subroutine metgen
     x  (idnode,ntpvdw,ntpatm,dlrpot,rcut,lstvdw,ltpvdw,prmvdw,vvv,ggg)
c     
c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for metal potentials
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith june 1995
c     
c     wl
c     2003/05/08 08:45:11
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension lstvdw(mxvdw),ltpvdw(mxvdw),prmvdw(mxvdw,mxpvdw)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
c     
c     sutton-chen potentials
      
      vv1(r,a,b,c)=a*(b/r)**c
      gg1(r,a,b,c)=a*c*(b/r)**c
      
#ifdef VAMPIR
      call VTBEGIN(151, ierr)
#endif
c     
c     define grid resolution for potential arrays
      
      dlrpot=rcut/dble(mxgrid-4)
c     
c     construct arrays for metal potentials
      
      kvdw=0
      do katm1=1,ntpatm

        do katm2=1,katm1

          kvdw=kvdw+1

          if(lstvdw(kvdw).gt.0)then

            ivdw=lstvdw(kvdw)
            
            if(ltpvdw(ivdw).eq.100)then
              
              jvdw=ivdw+ntpvdw+1

              do i=1,mxgrid
                
                rrr=dble(i)*dlrpot
                vvv(i,ivdw)=vv1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
                vvv(i,jvdw)=vv1(rrr,1.d0,prmvdw(ivdw,2),
     x            prmvdw(ivdw,4))
                ggg(i,ivdw)=gg1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
                ggg(i,jvdw)=gg1(rrr,0.5d0,prmvdw(ivdw,2),
     x            prmvdw(ivdw,4))

              enddo
              
              if(katm1.eq.katm2)then

                vvv(1,jvdw)=(prmvdw(ivdw,1)*prmvdw(ivdw,5))**2
                vvv(2,jvdw)=(prmvdw(ivdw,1)*prmvdw(ivdw,5))**2

              else

                nvdw=lstvdw((katm1*(katm1+1))/2)
                mvdw=lstvdw((katm2*(katm2+1))/2)
                vvv(1,jvdw)=(prmvdw(nvdw,1)*prmvdw(nvdw,5))**2
                vvv(2,jvdw)=(prmvdw(mvdw,1)*prmvdw(mvdw,5))**2

              endif

            else if(ltpvdw(ivdw).gt.100)then
              
              call error(idnode,151)
              
            endif

          endif
          
        enddo

      enddo

#ifdef VAMPIR
      call VTEND(151, ierr)
#endif
      return
      end
