! 19 DEC 05 - IUCHI - output tabulated potential and forces for check
!
      subroutine fortab
     x  (idnode,ntpvdw,ntpatm,mxnode,dlrpot,rcut,engunit,
     x  unqatm,lstvdw,ltpvdw,prmvdw,vvv,ggg,buffer)
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith march 1994
c     
c     wl
c     2001/05/30 12:40:06
c     1.6
c     Exp
!
!     Last updated: 20 Dec 2005 by S. Iuchi
c     
c***********************************************************************
c     
      
#include "dl_params.inc"

      logical safe
      character*256 record
      character*8 unqatm(mxsite),atom1,atom2
      dimension lstvdw(mxvdw),ltpvdw(mxvdw),prmvdw(mxvdw,mxpvdw)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw),buffer(mxbuff)
      
#ifdef VAMPIR
      call VTBEGIN(142, ierr)
#endif
      if(idnode.eq.0)open (ntable,file='TABLE')
c     
c     skip header record
      
      call getrec(safe,idnode,mxnode,ntable,record)
      if(.not.safe)go to 100
c     
c     read mesh resolution
      
      call getrec(safe,idnode,mxnode,ntable,record)
      if(.not.safe)go to 100
      delpot=dblstr(record(1:1),40,idum)
      cutpot=dblstr(record(idum:idum),40,idum1)
      idum1 = idum+idum1
      ngrid=intstr(record(idum1:idum1),40,idum)
      dlrpot=rcut/dble(mxgrid-4)
      if(abs(delpot-dlrpot).le.1.d-8)delpot=dlrpot
      if((delpot.gt.dlrpot).or.(ngrid-4.ne.nint(cutpot/delpot)))then
        if(idnode.eq.0) then
          write(nrite,"('expected radial increment : ',1p,e15.7,/,
     x                '   TABLE radial increment : ',1p,e15.7,/,/,    
     x                'expected number of grid points : ',0p,i10,/,
     x                'grid points in TABLE           : ',i10)")
     x      dlrpot, delpot, mxgrid, ngrid
        endif
        
        call error(idnode,22)

      endif

      if(cutpot.lt.rcut) call error(idnode,504)
      if(idnode.eq.0) then
        if(abs(1.d0-(delpot/dlrpot)).gt.1d-7) write(nrite,
     x    "(/,' TABLE arrays resized for mxgrid = ',i10)") mxgrid

      endif
c     
c     read potential arrays for all pairs
      
      do ivdw=1,ntpvdw

c     
c     read potential arrays if potential not already defined
        
        if(ltpvdw(ivdw).eq.0)then
          
c     
c     read pair potential labels and long range corrections
          
          call getrec(safe,idnode,mxnode,ntable,record)
          if(.not.safe)go to 100

          atom1=record(1:8)
          atom2=record(9:16)
          call strip(atom1,8)
          call strip(atom2,8)
          prmvdw(ivdw,1)=dblstr(record(17:17),15,idum)
          prmvdw(ivdw,2)=dblstr(record(32:32),15,idum)
          
          katom1=0
          katom2=0
          
          do jtpatm=1,ntpatm
            
            if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
            if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
            
          enddo
          
          if(katom1.eq.0.or.katom2.eq.0) then
            if(idnode.eq.0) 
     x        write(nrite,'(a)') '****',atom1,'***',atom2,'****'
            call  error(idnode,81)
          endif

          
          keyvdw=(max(katom1,katom2)*(max(katom1,katom2)-1))/2+
     x      min(katom1,katom2)
          
          if(lstvdw(keyvdw).ne.ivdw) call error(idnode,23)
          
c     
c     read potential arrays
          
          
          if(idnode.eq.0)then

            if(mxbuff.lt.ngrid)  then
              
              write(nrite,*) 'mxbuff must be >= ',ngrid,' in fortab'
              call error(idnode,48)
              
            endif


            read(ntable,'(4e15.8)',end=100)(buffer(i),i=1,ngrid)
c
c     reconstruct arrays using 3pt interpolation

            rdr = 1.d0/delpot
            do i= 1,mxgrid
              rrr = dble(i)*dlrpot
              l = int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              vk = buffer(l)
              vk1 = buffer(l+1)
              vk2 = buffer(l+2)
            
              t1 = vk + (vk1-vk)*ppp
              t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)
              vvv(i,ivdw) = t1 + (t2-t1)*ppp*0.5d0

            enddo

            read(ntable,'(4e15.8)',end=100)(buffer(i),i=1,ngrid)
c
c     reconstruct ggg arrays using 3pt interpolation

            do i= 1,mxgrid
              rrr = dble(i)*dlrpot
              l = int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              vk = buffer(l)
              vk1 = buffer(l+1)
              vk2 = buffer(l+2)
            
              t1 = vk + (vk1-vk)*ppp
              t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)
            
              ggg(i,ivdw) = t1 + (t2-t1)*ppp*0.5d0

            enddo

            call gdsum(vvv(1,ivdw),mxgrid,buffer)
            call gdsum(ggg(1,ivdw),mxgrid,buffer)

          else
            
            if(mxbuff.lt.mxgrid) call error(idnode,48)

            do i=1,mxgrid

              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            enddo

            call gdsum(vvv(1,ivdw),mxgrid,buffer)
            call gdsum(ggg(1,ivdw),mxgrid,buffer)

          endif

        endif
        
      enddo
c     
c     convert to internal units
      
      do k=1,ntpvdw
        
        if(ltpvdw(k).eq.0)then

          do i=1,mxgrid
            
            vvv(i,k)=vvv(i,k)*engunit
            ggg(i,k)=ggg(i,k)*engunit
            
          enddo
          
        endif
        
      enddo
      
      if(idnode.eq.0)close (ntable)
      
      if(idnode.eq.0)write(nrite,'(/,/,1x,a)')
     x  'potential tables read from TABLE file'

!              
!     output tabulated potential and forces for check
!     (input unit is used)

      do k=1,ntpvdw

         write(6,'(a2,i2,a30)') '# ',k,'Tabulated potential and forces'
         
         do i=1,mxgrid
            
            rrr = dble(i)*dlrpot
            write(6,'(5d20.10)') rrr, vvv(i,k),
     $           ggg(i,k) / rrr , vvv(i,k) / engunit,
     $           ggg(i,k) / rrr / engunit
            
         enddo
      enddo
     
#ifdef VAMPIR
      call VTEND(142, ierr)
#endif
      return
      
c     
c     end of file error exit
      
  100 continue

      if(idnode.eq.0)close (ntable)
      
      call error(idnode,24)
      
      end

