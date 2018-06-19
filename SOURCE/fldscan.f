      subroutine fldscan
     x   (idnode,mxnode,nfield,ntable,mxtbnd,mxtang,mxtdih,
     x   mxngp,mxgrp,mxatms,mxgatm,mxneut,mxteth,mxtet1,mxsvdw,
     x   mxvdw,mxn1,mxtbp,mxexcl,mxsite,mxbond,mxcons,mxangl,mxungp,
     x   mxdihd,mxtshl,mxshl,mxpmf,mxspmf,mxtinv,mxinv,mxfbp,ngrid,
     x   mxtmls,mxtcon,rctbp,rcfbp)
c     
c***********************************************************************
c     
c     dl_poly routine for scanning the field file to determine the
c     required parameters
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c     
c     wl
c     2003/05/08 08:45:11
c     1.9
c     Exp
c
c***********************************************************************
c     
      implicit real*8(a-h,o-z)
      
      parameter (mega=10000,mmk=1000)
      
      character*80 header
      character*256 record
      character*8 name,chr(mmk)
      logical check,ltable,safe,lneut
c sandeep

      integer,parameter :: nofields=100
      integer           :: countfields,fieldlength,tmpint
      integer           :: starts(nofields),ends(nofields)


      
#ifdef VAMPIR
      call VTBEGIN(139, ierr)
#endif
      mxtmls=0
      mxatms=0
      mxgrp=0
      mxtcon=0
      mxtbnd=0
      mxtang=0
      mxtdih=0
      mxtinv=0
      mxpmf=0
      mxspmf=0
      mxungp=0
      mxngp=0
      mxneut=0
      nmetal=0
      mxn1=0
      nxn1=0
      nold = -1
      mxgatm=0
      mxteth=0
      mxtet1=0
      mxsvdw=0
      mxvdw=0
      mxtbp=0
      mxexcl=0
      mxsite=0
      mxbond=0
      mxcons=0
      mxangl=0
      mxdihd=0
      mxinv=0
      mxshl=0
      mxtshl=0
      mxfbp=0
      ngrid=0
      rctbp=0.d0
      rcfbp=0.d0
      safe=.true.
      lneut=.false.
      ltable=.false.

c     
c     open force field data file
      
      if(idnode.eq.0)open (nfield,file='FIELD')
      
      call getrec(safe,idnode,mxnode,nfield,record)
      if(.not.safe)go to 3000
      
c     
c     read and process directives from field file
      
      do nrecs=1,mega
        
        call getrec(safe,idnode,mxnode,nfield,record)
        if(.not.safe)go to 3000

        call cal_field(record,nofields,countfields,starts,ends)
        call lowcase(record,40)
c        call strip(record,256)
        if(record(starts(1):ends(1)).eq.'neut')then

          lneut=.true.

        elseif(record(starts(1):ends(1)).eq.'molecular')then
          
          mxtmls=intstr(record(starts(2):ends(2)),80,idum)
          
          do itmols=1,mxtmls
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 3000
            
            call cal_field(record,nofields,countfields,starts,ends)
 
            do itrec=1,mega
              call getrec(safe,idnode,mxnode,nfield,record)
              if(.not.safe)go to 3000
              call cal_field(record,nofields,countfields,starts,ends)

              call lowcase(record,40)
c              call strip(record,256)
              
              ksite=0
              
              if(record(starts(1):ends(1)).eq.'nummols')then
                
                nummols=intstr(record(starts(2):ends(2)),80,idum)
                
              elseif(record(starts(1):ends(1)).eq.'atoms')then
                
                numsit=intstr(record(starts(2):ends(2)),80,idum)
                mxatms = mxatms+numsit*nummols
                mxsite=mxsite+numsit
                ksite=0
                do isite=1,numsit
                  
                  if(ksite.lt.numsit)then
                    
                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 3000
                    call cal_field
     x                   (record,nofields,countfields,starts,ends)
                    
                    name=record(starts(1):ends(1))
                    call strip(name,8)
c sandeep: set variables to default if they are undefined
                    tmpint=scan(record(starts(4):ends(4)),'.')
                    if(tmpint==0) then 
                      nrept=intstr(record(starts(4):ends(4)),80,idum)
                      ifrz =intstr(record(starts(5):ends(5)),80,idum)
                      nneu=intstr(record(starts(6):ends(6)),80,idum)
                     else
                      nrept=intstr(record(starts(6):ends(6)),80,idum)
                      ifrz =intstr(record(starts(7):ends(7)),80,idum)
                      nneu=intstr(record(starts(8):ends(8)),80,idum)
                      if(starts(6)==0) nrept=1
                      if(starts(7)==0) ifrz=0
                      if(starts(8)==0) neugp=0
                    endif
                    if(nrept.eq.0)nrept=1
                    
                    if(lneut)then
                      if(nneu.ne.nold) nxn1 =0
                      nxn1 = nxn1+nrept
                      mxn1 = max(mxn1,nxn1)
                      nold = nneu
                    endif

                    if(mxsvdw.eq.0)then

                      mxsvdw=1
                      chr(1)=name

                    else

                      check=.true.
                      do j=1,mxsvdw

                        if(name.eq.chr(j))check=.false.

                      enddo
                      if(check)then

                        mxsvdw=mxsvdw+1
                        if(mxsvdw.le.mmk)chr(mxsvdw)=name

                      endif

                    endif
                    if(nrept.eq.0)nrept=1
                    ksite=ksite+nrept
                    
                  endif
                  
                enddo

                if(mmk.lt.mxsvdw)call error(idnode,34)

                if(lneut)mxneut = mxneut+nneu*nummols
                
              elseif(record(starts(1):ends(1)).eq.'shell')then
                numshl=intstr(record(starts(2):ends(2)),80,idum)
                mxtshl=mxtshl+numshl
                mxshl=mxshl+nummols*numshl

                do ishls=1,numshl

                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'bonds')then
                
                numbonds=intstr(record(starts(2):ends(2)),80,idum)
                mxtbnd=mxtbnd+numbonds
                mxbond=mxbond+nummols*numbonds
                
                do ibonds=1,numbonds
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'constraints')then
                
                numcon=intstr(record(starts(2):ends(2)),80,idum)
                mxtcon=mxtcon+numcon
                mxcons=mxcons+nummols*numcon
                
                do icon=1,numcon
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'angles')then
                
                numang=intstr(record(starts(2):ends(2)),80,idum)
                mxtang=mxtang+numang
                mxangl=mxangl+nummols*numang
                
                do iang=1,numang
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'dihedrals')then
                
                numdih=intstr(record(starts(2):ends(2)),80,idum)
                mxtdih=mxtdih+numdih
                mxdihd=mxdihd+nummols*numdih
                
                do idih=1,numdih
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'inversions')then
                
                numinv=intstr(record(starts(2):ends(2)),80,idum)
                mxtinv=mxtinv+numinv
                mxinv=mxinv+nummols*numinv
                
                do iinv=1,numinv
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'rigid')then
                
                numgrp=intstr(record(starts(2):ends(2)),80,idum)
                mxungp=mxungp+numgrp
                mxgrp = mxgrp+numgrp*nummols
                
                do kgrp=1,numgrp
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  numgsit=intstr(record(starts(1):ends(1)),5,idum)
                  mxgatm = mxgatm+numgsit*nummols
                  mxngp=max(mxngp,numgsit)
                  do j=16,numgsit,16

                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 3000

                  enddo
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'teth')then
                
                numteth=intstr(record(starts(2):ends(2)),80,idum)
                mxteth=mxteth+numteth
                mxtet1 = mxtet1+numteth*nummols
                
                do iteth=1,numteth
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  
                enddo
                
              elseif(record(starts(1):ends(1)).eq.'pmf')then
                
                do ipmf = 1,2
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 3000
                  call cal_field
     x                 (record,nofields,countfields,starts,ends)
c                  call strip(record,40)
                  call lowcase(record,40)
! pmf unit are two units, use 3 
                  npmf=intstr(record(starts(3):ends(3)),80,idum)       
                  mxspmf=mxspmf+npmf

                  do jpmf=1,npmf
                    
                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 3000
                    
                  enddo
                  
                enddo

                mxpmf=mxpmf+nummols
                
              elseif(record(1:6).eq.'finish')then
                
                go to 1000
                
              endif
              
            enddo
            
 1000       continue
            
          enddo
          
        elseif((record(starts(1):ends(1)).eq.'vdw').or.
     x         (record(starts(1):ends(1)).eq.'vdwtable')) then
c        elseif(record(1:3).eq.'vdw') then
c          call strip(record(4:4),37)
          if(record(starts(1):ends(1)).eq.'vdwtable')ltable=.true.
c          if(record(4:6).eq.'tab')ltable=.true.
          ntpvdw=intstr(record(starts(2):ends(2)),80,idum)
          mxvdw=max(ntpvdw,(mxsvdw*(mxsvdw+1))/2)
          do itpvdw=1,ntpvdw
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 3000
            call cal_field(record,nofields,countfields,starts,ends)
            call lowcase(record,40)
            if(record(starts(3):ends(3)).eq.'stch')nmetal=nmetal+1
            if(record(starts(3):ends(3)).eq.'tab')ltable=.true.
c  sandeep:    warning:  Empty space won't mean anything now.  should be changed later.
            if(record(starts(3):ends(3)).eq.'   ')ltable=.true.
            
          enddo
          if(nmetal.gt.0)mxvdw=2*mxvdw
          if(ltable)then

            if(idnode.eq.0)open(ntable,file='TABLE')

            call getrec(safe,idnode,mxnode,ntable,record)
            if(.not.safe)go to 4000
            call getrec(safe,idnode,mxnode,ntable,record)
            if(.not.safe)go to 4000
            ngrid=intstr(record(31:31),10,idum)

            close (ntable)

          endif

          
        elseif(record(starts(1):ends(1)).eq.'tbp') then
          
          ntptbp=intstr(record(starts(2):ends(2)),80,idum)
          mxtbp=ntptbp
          
          do itptbp=1,ntptbp
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 3000
            call cal_field(record,nofields,countfields,starts,ends)
            rct=dblstr(record(starts(9):ends(9)),12,idum)
            rctbp=max(rctbp,rct)

          enddo
          
        elseif(record(starts(1):ends(1)).eq.'fbp') then
          
          ntpfbp=intstr(record(starts(2):ends(2)),80,idum)
          mxfbp=ntpfbp
          do itpfbp=1,ntpfbp
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 3000
            call cal_field(record,nofields,countfields,starts,ends)
            rct=dblstr(record(starts(8):ends(8)),80,idum)
            rcfbp=max(rcfbp,rct)

          enddo
          
        elseif(record(starts(1):ends(1)).eq.'extern') then
          
          call getrec(safe,idnode,mxnode,nfield,record)
          if(.not.safe)go to 3000
          call cal_field(record,nofields,countfields,starts,ends)
          nfld=intstr(record(starts(2):ends(2)),80,idum)
          if(nfld.eq.0)nfld=5
          
          do k=1,nfld,5

            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000

          enddo

        elseif(record(starts(1):ends(1)).eq.'close')then
          
          go to 2000
          
        endif
        
      enddo
      
 2000 continue
      if(idnode.eq.0)close (nfield)

      if(mxpmf.gt.0)mxpmf=mxatms
      if(mxtcon.gt.0)mxexcl=max(mxexcl,6)
      if(mxtbnd.gt.0)mxexcl=max(mxexcl,6)
      if(mxtang.gt.0)mxexcl=max(mxexcl,16)
      if(mxtdih.gt.0)mxexcl=max(mxexcl,50)
      if(mxtinv.gt.0)mxexcl=max(mxexcl,50)
      if(mxneut.gt.0)mxexcl=max(mxexcl,10*mxn1*mxn1)
      if(mxgrp.gt.0)mxexcl=max(mxexcl,mxngp)

#ifdef VAMPIR
      call VTEND(139, ierr)
#endif
      return

 3000 continue
      if(idnode.eq.0) close (nfield)
      call error(idnode,52)
      return

 4000 continue
      if(idnode.eq.0) close (ntable)
      call error(idnode,24)
      return

      end
