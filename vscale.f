      subroutine vscale
     x  (idnode,mxnode,imcon,natms,degfre,temp,lstfrz,
     x  buffer,cell,weight,vxx,vyy,vzz,xxx,yyy,zzz)
      
c     
c*********************************************************************
c     
c     dl_poly subroutine for scaling the velocity arrays to the
c     desired temperature
c     Only for systems without rigid body units
c     
c     zeroes angular momentum in non-periodic system.
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1992.
c     author - w.smith july 1992
c     amended - t.forester oct 1993
c     amended - t.forester dec 1994 : block data
c     
c     wl
c     2000/01/18 14:06:00
c     1.4
c     Exp
c     
c*********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension buffer(mxbuff)
      dimension lstfrz(mxatms)
      dimension rot(9),rotinv(9),cell(9)
      dimension weight(natms),vxx(natms),vyy(natms),vzz(natms)
      dimension xxx(natms),yyy(natms),zzz(natms)
#ifdef VAMPIR
      call VTBEGIN(63, ierr)
#endif
c     
c     block indices
      
      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode
c     
c     calculate centre of mass position and motion from the system

  100 continue
      cmx=0.d0
      cmy=0.d0
      cmz=0.d0
      cmvx=0.d0
      cmvy=0.d0
      cmvz=0.d0
      sysmas =0.d0
      
      nnatm=0
      do i=iatm1,iatm2
        
        if(lstfrz(i).eq.0) then
          
          nnatm = nnatm + 1
          cmx = cmx + weight(i)*xxx(i)
          cmy = cmy + weight(i)*yyy(i)
          cmz = cmz + weight(i)*zzz(i)
          sysmas = sysmas + weight(i)
          
          cmvx=cmvx+vxx(i)*weight(i)
          cmvy=cmvy+vyy(i)*weight(i)
          cmvz=cmvz+vzz(i)*weight(i)
          
        endif
        
      enddo
      
      if(mxnode.gt.1) then
        buffer(8) = sysmas
        buffer(9) = cmx
        buffer(10) = cmy
        buffer(11) = cmz
        buffer(12) = cmvx
        buffer(13) = cmvy
        buffer(14) = cmvz
        
        call gdsum(buffer(8),7,buffer(1))
        
        sysmas =  buffer(8) 
        cmx    =  buffer(9) 
        cmy    =  buffer(10) 
        cmz    =  buffer(11) 
        cmvx   =  buffer(12) 
        cmvy   =  buffer(13) 
        cmvz   =  buffer(14) 
      endif
      
      cmx = cmx/sysmas
      cmy = cmy/sysmas
      cmz = cmz/sysmas
      
      cmvx=cmvx/sysmas
      cmvy=cmvy/sysmas
      cmvz=cmvz/sysmas
      
c     
c     remove centre of mass motion  

      do i=1,natms
        
        if(lstfrz(i).eq.0) then
          
          vxx(i)=vxx(i)-cmvx
          vyy(i)=vyy(i)-cmvy
          vzz(i)=vzz(i)-cmvz
          
        elseif(lstfrz(i).ne.0) then
          
          vxx(i)=0.d0
          vyy(i)=0.d0
          vzz(i)=0.d0
          
        endif
        
      enddo
      
c     
c     zero angular momentum about centre of mass - non-periodic system
      
      if(imcon.eq.0) then
c     
c     move to centre of mass origin
        
        do i = 1,natms
          
          xxx(i) = xxx(i) - cmx
          yyy(i) = yyy(i) - cmy
          zzz(i) = zzz(i) - cmz
          
        enddo
        
c     
c     angular momentum accumulators
        
        amx = 0.d0
        amy = 0.d0
        amz = 0.d0
c     
c     rotational inertia accumulators
        
        do i = 1,9
          rot(i) = 0.d0
        enddo
        
        do i = iatm1,iatm2
          
          amx = amx + weight(i)*(yyy(i)*vzz(i) - zzz(i)*vyy(i))
          amy = amy + weight(i)*(zzz(i)*vxx(i) - xxx(i)*vzz(i))
          amz = amz + weight(i)*(xxx(i)*vyy(i) - yyy(i)*vxx(i))
          
          rsq = xxx(i)**2 + yyy(i)**2 + zzz(i)**2
          rot(1) = rot(1) + weight(i)*(xxx(i)*xxx(i) - rsq)
          rot(2) = rot(2) + weight(i)* xxx(i)*yyy(i)
          rot(3) = rot(3) + weight(i)* xxx(i)*zzz(i)
          rot(5) = rot(5) + weight(i)*(yyy(i)*yyy(i) - rsq)
          rot(6) = rot(6) + weight(i)* yyy(i)*zzz(i)
          rot(9) = rot(9) + weight(i)*(zzz(i)*zzz(i) - rsq)
          
        enddo
c     
c     complete rotational inertia matrix
        
        rot(4) = rot(2)
        rot(7) = rot(3)
        rot(8) = rot(6)
c     
c     globally sum
        
        if(mxnode.gt.1) then
          buffer(13) = amx
          buffer(14) = amy
          buffer(15) = amz
          do i = 1,9
            buffer(15+i) = rot(i)
          enddo
          
          call gdsum(buffer(13),12,buffer(1))
          
          amx =  buffer(13) 
          amy =  buffer(14) 
          amz =  buffer(15) 
          do i = 1,9
            rot(i) = buffer(15+i)
          enddo
        endif
c     
c     invert rotational inertia matrix
        
        call invert (rot,rotinv,det)
c     
c     correction to angular velocity
        
        wxx = rotinv(1)*amx + rotinv(2)*amy + rotinv(3)*amz
        wyy = rotinv(4)*amx + rotinv(5)*amy + rotinv(6)*amz
        wzz = rotinv(7)*amx + rotinv(8)*amy + rotinv(9)*amz
c     
c     correction to linear velocity
        
        do i = 1,natms
          
          if(lstfrz(i).eq.0) then
            
            vxx(i) = vxx(i) + (wyy*zzz(i) - wzz*yyy(i))
            vyy(i) = vyy(i) + (wzz*xxx(i) - wxx*zzz(i))
            vzz(i) = vzz(i) + (wxx*yyy(i) - wyy*xxx(i))
            
          endif
          
        enddo
c     
c     reset positions to original reference frame
        
        do i = 1,natms
          
          xxx(i) = xxx(i) + cmx
          yyy(i) = yyy(i) + cmy
          zzz(i) = zzz(i) + cmz
          
        enddo
        
      endif
c     
c     evaluate translational temperature
c     and scale velocities

      engke = 0.d0
      do i = iatm1,iatm2
        engke = engke + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
      enddo
      engke = engke *0.5d0
      
      if(mxnode.gt.1) call gdsum (engke,1,buffer)
      
      scale=1.d0
      sigma = 0.5d0*temp*boltz*degfre
      if(engke.gt.1.d-6)scale=sqrt(sigma/engke)

      do i = iatm1,iatm2

        vxx(i) = vxx(i)*scale
        vyy(i) = vyy(i)*scale
        vzz(i) = vzz(i)*scale

      enddo
      
      if(mxnode.gt.1) then
c     
c     merge new atomic coordinates and velocities

        nbuff = mxbuff
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
      endif
      
#ifdef VAMPIR
      call VTEND(63, ierr)
#endif
      return
      end
