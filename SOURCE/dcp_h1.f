      subroutine dcp_h0
     x  (dipmas,polr,dipx,dipy,dipz,dspx,dspy,dspz,vdxx,vdyy,vdzz,
     x  udxx,udyy,udzz,fdxx,fdyy,fdzz,chitd,engdke,emux,emuy,emuz,
     x  efieldkx,efieldky,efieldkz,efdcrecx,efdcrecy,efdcrecz,
     x  efddmurecx,efddmurecy,efddmurecz,sigmad,tautd,consvd,conintd,
     x  idnode,mxnode,natms,imcon,tstep,buffer)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     modified for dipole degrees of freedom - t. yan dec 2003
c     
c     wl
c     2000/01/18 14:05:47
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"

      dimension polr(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension dspx(mxatms),dspy(mxatms),dspz(mxatms)
      dimension vdxx(mxatms),vdyy(mxatms),vdzz(mxatms)
      dimension udxx(mxatms),udyy(mxatms),udzz(mxatms)
      dimension fdxx(mxatms),fdyy(mxatms),fdzz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      dimension buffer(mxbuff)

#ifdef VAMPIR
      call VTBEGIN(32, ierr)
#endif
c     
c     block indices
      
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c
c     store initial values of dipole and velocity

      do i=iatm0,iatm1

        dspx(i) = dipx(i)
        dspy(i) = dipy(i)
        dspz(i) = dipz(i)

        udxx(i) = vdxx(i)
        udyy(i) = vdyy(i)
        udzz(i) = vdzz(i)

      enddo
c
c     estimate kinetic energy at current timestep

      engdke = 0.d0
      do i = iatm0,iatm1

        fdxx(i)=-dipx(i)/polr(i)*r4pie0+efieldkx(i)+efdcrecx(i)+
     x       emux(i)+efddmurecx(i)
        fdyy(i)=-dipy(i)/polr(i)*r4pie0+efieldky(i)+efdcrecy(i)+
     x       emuy(i)+efddmurecy(i)
        fdzz(i)=-dipz(i)/polr(i)*r4pie0+efieldkz(i)+efdcrecz(i)+
     x       emuz(i)+efddmurecz(i)

        vdxt = vdxx(i) + 0.5d0*tstep*(fdxx(i)/dipmas)
        vdyt = vdyy(i) + 0.5d0*tstep*(fdyy(i)/dipmas)
        vdzt = vdzz(i) + 0.5d0*tstep*(fdzz(i)/dipmas)
c
c     kinetic energy at current timestep

        engdke = engdke + dipmas*(vdxt**2+vdyt**2+vdzt**2)
c
c     first estimate of new velocities - no thermostat

        vdxx(i) = udxx(i) + tstep*(fdxx(i)/dipmas)
        vdyy(i) = udyy(i) + tstep*(fdyy(i)/dipmas)
        vdzz(i) = udzz(i) + tstep*(fdzz(i)/dipmas)
c
c     advance dipoles using leapfrog

        dipx(i) = dspx(i) + tstep*vdxx(i)
        dipy(i) = dspy(i) + tstep*vdyy(i)
        dipz(i) = dspz(i) + tstep*vdzz(i)

      enddo
c
c     kinetic energy

      engdke =engdke*0.5d0
      if(mxnode.gt.1) call gdsum(engdke,1,buffer)
c
c     propagate chitd

      chitp = (engdke/sigmad - 1.d0)/tautd**2
      chitnew = chitd + tstep*chitp
      chit0 = 0.5d0*(chitd+chitnew)
c
c     begin thermostat iteration

      maxit = 3

      do iter = 1,maxit

        engdke=0.d0
c
c     find better estimate to kinetic energy at half timestep

        do i=iatm0,iatm1
c
c     calculate new velocities at half step

          vdxt=0.5d0*(udxx(i)+vdxx(i))
          vdyt=0.5d0*(udyy(i)+vdyy(i))
          vdzt=0.5d0*(udzz(i)+vdzz(i))

          vdxx(i)=udxx(i)+tstep*(fdxx(i)/dipmas-chit0*vdxt)
          vdyy(i)=udyy(i)+tstep*(fdyy(i)/dipmas-chit0*vdyt)
          vdzz(i)=udzz(i)+tstep*(fdzz(i)/dipmas-chit0*vdzt)
c
c     advance dipoles using leapfrog

          dipx(i) = dspx(i) + tstep*vdxx(i)
          dipy(i) = dspy(i) + tstep*vdyy(i)
          dipz(i) = dspz(i) + tstep*vdzz(i)
c
c     kinetic energy at current timestep

          vdxt=0.5d0*(udxx(i)+vdxx(i))
          vdyt=0.5d0*(udyy(i)+vdyy(i))
          vdzt=0.5d0*(udzz(i)+vdzz(i))

          engdke = engdke + dipmas*(vdxt**2+vdyt**2+vdzt**2)

        enddo

        engdke = engdke*0.5d0
        if(mxnode.gt.1) call gdsum(engdke,1,buffer)

c        if(idnode.eq.0)
c     x	write(nrite,*)dipx(1),tstep,vdxx(i)
c
c     improved prediction of chitd

        chitp = (engdke/sigmad - 1.d0)/tautd**2
        chitnew = chitd + tstep*chitp
        chit0 = 0.5d0*(chitd+chitnew)
c
c     end of thermostat iterations

      enddo
c
c     update thermostat variable

      chitd = chitnew
c
c     conserved quantity less kinetic and potential energy terms

      conintd = conintd + 2.d0*sigmad*tstep*chit0
      cons0 = conintd
      cons1 = sigmad*(tautd*chit0)**2
      consvd = cons0 + cons1
c
c     global exchange of dipoles and dipole velocities
      
      if(mxnode.gt.1)then
        
        nbuff=mxbuff
        call merge(idnode,mxnode,natms,nbuff,dipx,dipy,dipz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vdxx,vdyy,vdzz,buffer)
        
      endif
c      if(idnode.eq.0)
c     x write(nrite,*)'in dcp_h0, engdke, consvd',engdke,consvd

#ifdef VAMPIR
      call VTEND(32, ierr)
#endif
      return
      end
