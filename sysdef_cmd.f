! PAESANI GROUP: 04/18/2015.
      subroutine sysdef_cmd
     x  (lneut,lmetal,lnsq,molnam,mbnrg_index,sitnam,
     x  unqatm,idnode,mxnode,keyfce,
     x  keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpfbp,nshels,
     x  nhko,nlatt,alpha,dlrpot,drewd,engunit,prmpmf,rcut,rvdw,rcuttb,
     x  rcutfb,indpmf,keyang,keybnd,keydih,keyinv,keytet,lfzsit,listyp,
     x  lstang,lstbnd,lstcon,lstdih,lstinv,lstgst,lstvdw,lsttbp,lstfbp,
     x  lstshl,ltpsit,lsttet,ltpvdw,ltptbp,ltpfbp,npmf,nugrp,numang,
     x  numbonds,numcon,numdih,numinv,numgrp,numgsit,nummols,numsit,
     x  numpmf,numteth,numshl,chgsit,lpolar,polarsit,polarsit2,
     x  erc,fer,ercp,pmfwght,ggg,prmang,prmbnd,
     x  prmcon,prmdih,prminv,prmfld,prmtet,prmvdw,prmtbp,prmfbp,prmshl,
     x  vvv,wgtsit,rcut3b,rcut4b,buffer,ahk,hon,dhn,fon,
     x  keyumb,prmumb,lmbnrg)

c     
c***********************************************************************
c     
c     dl_poly subroutine for reading in the molecular specifications
c     of the system to be simulated
c     version for rigid unit data and neutral groups
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     amended   - w.smith march 1994 
c     amended   - t.forester april 1994
c     amended   - w.smith  dec 1994 - get rec etc
c     amended   - a.smondyrev may 2000 - keydih=5 for 
c                 ryckaert-bellemans potential in dihedrals
c     amended   - a.smondyrev may 2000 - keydih=6 for 
c                 fluorinated ryckaert-bellemans potential in dihedrals
c     
c     wl
c     2001/06/12 13:00:03
c     1.19
c     Exp
!
!     Last updated: 23 Feb 2006 by S. Iuchi
c     
c***********************************************************************
c     
! from module

      use ps_type_pes, only: lpspes, init_ps_type_pes
      
#include "dl_params.inc"
      
      parameter (mega=10000)
      
      character*4  keyword
      character*8  atom0,atom1,atom2,atom3
      character*8  sitnam(mxsite)
      character*8  unqatm(mxsite)
      character*40 molnam(mxtmls)
      character*256 record
      character*256 ctmp
      logical atmchk,lunits,lmols,lneut,ltable,lmetal,lnsq,lshl,safe
      logical lpmf,novdw,lpolar

      dimension ahk(0:mxhko)
      dimension hon(mxegrd,0:mxhko),dhn(mxegrd,0:mxhko),fon(mxegrd,0:7)
      dimension nummols(mxtmls),numsit(mxtmls),numcon(mxtmls)
      dimension numpmf(mxtmls),numinv(mxtmls)
      dimension numdih(mxtmls),numang(mxtmls),numbonds(mxtmls)
      dimension ltpsit(mxsite),lstcon(mxtcon,2),prmcon(mxtcon)
      dimension wgtsit(mxsite),chgsit(mxsite),parpot(30)
      dimension polarsit(mxsite),polarsit2(mxsite)
      dimension lfzsit(mxsite),nugrp(mxsite),buffer(mxbuff)
      dimension keybnd(mxtbnd),lstbnd(mxtbnd,2),prmbnd(mxtbnd,mxpbnd)
      dimension keyang(mxtang),lstang(mxtang,3),prmang(mxtang,mxpang)
      dimension keydih(mxtdih),lstdih(mxtdih,4),prmdih(mxtdih,mxpdih)
      dimension keyinv(mxtinv),lstinv(mxtinv,4),prminv(mxtinv,mxpinv)
      dimension lstvdw(mxvdw),ltpvdw(mxvdw),prmvdw(mxvdw,mxpvdw)
      dimension lsttbp(mxtbp),ltptbp(mxtbp),prmtbp(mxtbp,mxptbp)
      dimension lstfbp(mxfbp),ltpfbp(mxfbp),prmfbp(mxfbp,mxpfbp)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw),prmfld(mxfld)
      dimension erc(mxegrd),fer(mxegrd),rcut3b(mxtbp),rcut4b(mxfbp)
      dimension ercp(mxegrd,0:3)
      dimension numgrp(mxtmls),numgsit(mxungp)
      dimension listyp(mxungp),lstgst(mxungp,mxngp)
      dimension lsttet(mxteth),numteth(mxtmls),keytet(mxteth)
      dimension prmtet(mxteth,mxpbnd)
      dimension prmshl(mxtshl),lstshl(mxtshl,2),numshl(mxtmls)
      dimension npmf(2),indpmf(mxspmf)
      dimension pmfwght(mxspmf)
c      character(len=8) :: mbnrg_atname(mxtmls,mxsite)


! sandeep 

      integer,parameter :: nofields=100
      integer           :: countfields,fieldlength
      integer           :: starts(nofields),ends(nofields)
      integer           :: mbnrg_index(mxtmls)
  
      logical           :: lmbnrg,ltmp

      character(len=3)  :: tmpchr1
      character(len=4)  :: tmpchr2
      character(len=8)  :: tmpchr3
      

! PIMD/CMD
      integer :: keyumb
      real(8) :: prmumb(5)
#ifdef VAMPIR
      call VTBEGIN(2, ierr)
#endif
c     
c     initialise system counters: atomic site index, number of 
c     constraints, bond angles, dihedrals, inversions, chemical bonds,
c     unique atom types, total number of atoms,
c     total number of rigid groups, number of tethered atoms,
c     number of three body potentials
      
      nsite =0
      nconst=0
      nangle=0
      ndihed=0
      ninver=0
      nbonds=0
      ntpatm=0
      natms =0
      ngrp  =0
      nteth =0
      ntptbp=0
      ntpfbp=0
      nshels=0
      nspmf = 0
      keyfld = 0
      mbnrg_index=0
      
      lunits=.false.
      lmols =.false.
      lneut =.false.
      ltable=.false.
      lmetal=.false.
      lshl  =.false.
      lpmf  =.false.
      novdw =.true.
      lpspes = .false.   ! Partridge and Schwenke type PES
      engunit = 1.d0
      ltmp=.true.

      do i = 1,mxtmls
        numbonds(i) =0
        numpmf(i) =0
        numcon(i) =0
        numang(i) =0
        numdih(i) =0
        numinv(i) =0
        numgrp(i) =0
        numsit(i) =0
        numteth(i) =0
        numshl(i) =0
      enddo

      do ipmf = 1,2
        npmf(ipmf) = 0
      enddo
      do jpmf = 1,mxspmf
        indpmf(jpmf) = 0
      enddo

c     
c     open force field data file
      
      if(idnode.eq.0)open (nfield,file='FIELD',status='old')
      
      if(idnode.eq.0) 
     x  write(nrite,"(/,/,'SYSTEM SPECIFICATION')")
      
      call getrec(safe,idnode,mxnode,nfield,record)
      if(.not.safe)go to 2000

c     
c     read and process directives from field file
      
      do nrecs=1,mega
        
        call getrec(safe,idnode,mxnode,nfield,record)
        if(.not.safe)go to 2000
c     
c     convert to lowercase
        
        call lowcase(record,40)
c     
c     strip out blank
        
c        call strip(record,40)

c sandeep
        call cal_field(record,nofields,countfields,starts,ends)

        if(countfields==0)  cycle           ! takes care of blank lines
        
        if(record(1:1).eq.'#') then
c        if(record(1:1).eq.'#'.or.record(1:1).eq.' ') then
c     
c     record is commented out
          
          
        elseif(record(starts(1):ends(1)).eq.'units') then
c     
c     identify energy unit for input/output
          
          lunits=.true.
c          record(1:35) = record(6:40)
c          call strip(record,35)
          
          
          if(record(starts(2):ends(2)).eq.'ev') then 
            engunit = 9648.530821d0
            if(idnode.eq.0) 
     x        write(nrite,"(/,' energy units = electron volts ')")
            
          elseif(record(starts(2):ends(2)).eq.'kcal')then
            engunit = 418.4d0
            if(idnode.eq.0) 
     x        write(nrite,"(/,' energy units = kcal/ mol ')")
            
          elseif(record(starts(2):ends(2)).eq.'kj')then
            engunit = 1.d2
            if(idnode.eq.0)
     x        write(nrite,"(/,' energy units = kjoule/mol ')")
            
          elseif(record(starts(2):ends(2)).eq.'internal') then
            if(idnode.eq.0) 
     x        write(nrite,"(/,' energy units = dl_poly internal',
     x        ' units ')")

! sandeep warning: changing the syntax            
          elseif(starts(2).eq.0) then
c          elseif(record(starts(2):ends(2)).eq.' ') then
            if(idnode.eq.0) 
     x        write(nrite,"(/,' energy units = dl_poly internal ',
     x        'units ')")
            
          else
            
            if(idnode.eq.0) write(nrite,'(a)') record
            call error(idnode,5)
            
          endif

c     
c     neutral group control option
          
        elseif(record(starts(1):ends(1)).eq.'neut') then
          
          lneut =.true.
          if(idnode.eq.0) 
     x      write(nrite,"(/,' neutral group implementation in use')")

c     
c     can't have neutral groups with all-pairs

          if(lnsq) call error(idnode,426)
          
c     
c     specify molecular species
          
        elseif(record(starts(1):ends(1)).eq.'molecular')then
          
c     
c     number of molecular types
          
          if(lmols) call error(idnode,11)
          lmols=.true.
c  integer is at 3rd field, not at 2nd
          ntpmls=intstr(record(starts(3):ends(3)),80,idum)     
          
          if(idnode.eq.0) 
     x      write(nrite,"(/,/,1x,'number of molecular types',6x,i10)") 
     x      ntpmls
          
          if(ntpmls.gt.mxtmls) call error(idnode,10)


         
          
c
c     initialise total system charge
          
          sumchg=0.d0
          
c     
c     read in molecular characteristics
          
          do itmols=1,ntpmls
            
            if(idnode.eq.0) 
     x        write(nrite,"(/,1x,'molecular species type',9x,i10)") 
     x        itmols

c     
c     name of molecular species
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000

            call cal_field(record,nofields,countfields,starts,ends)
            
! mbnrg -- default is to use polynomials 

            molnam(itmols)=record(starts(1):ends(1))
            mbnrg_index(itmols)=0 ; ctmp=''

! expects that ion is the first molecule

! MRR & DZ - Add ion ion mb-nrg & multiple MB-nrg potentials for ions
!            By default, is going to use TTM/MBPOL.
!            if(ltmp.and.lmbnrg) then 
            if(lmbnrg) then 
              mbnrg_index(itmols)=0
!              mbnrg_index(itmols)=1 
              if(countfields>1) ctmp=record(starts(2):ends(2))
              if(trim(adjustl(ctmp)) == 'ttm-nrg') mbnrg_index(itmols)=0
              if(trim(adjustl(ctmp)) == 'mb-nrg') mbnrg_index(itmols)=1
!              ltmp=.false. ! Commented on 20072017
! END MRR & DZ

            endif 

!            molnam(itmols)=record(1:40)

            if(idnode.eq.0) 
     x        write(nrite,"(/,/,1x,'name of species:',13x,a40)") 
     x        molnam(itmols)
            if((idnode.eq.0).and.(mbnrg_index(itmols)>0)) 
     x        write(nrite,"(/,/,1x,'Polynomials active for species:'
     x        ,13x,a40)") molnam(itmols)
            
c     
c     stop processing if energy unit has not been specified
            
            if(.not.lunits)call error(idnode,6)
            
c     
c     read molecular data
            
            do itrec=1,mega
              
              call getrec(safe,idnode,mxnode,nfield,record)
              if(.not.safe)go to 2000

              call cal_field(record,nofields,countfields,starts,ends)
        
              call lowcase(record,40)
c              call strip(record,40)
              
              ksite=0
              
              if(record(starts(1):ends(1)).eq.'nummols')then
                
                nummols(itmols)=
     x            intstr(record(starts(2):ends(2)),80,idum)
                if(idnode.eq.0)
     x            write(nrite,"(/,1x,'number of molecules  ',
     x            10x,i10)") nummols(itmols)
                
              elseif(record(starts(1):ends(1)).eq.'atoms')then
c     
c     read in atomic details
                
                numsit(itmols)=intstr(record(starts(2):ends(2)),40,idum)
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of atoms/sites',
     x              10x,i10)") numsit(itmols)
                  if(.not.lneut .and. .not.lpolar)
     x              write(nrite,"(/,/,1x,'atomic characteristics:',
     x              /,/,21x,' site',5x,'name',10x,'mass',8x,
     x              'charge',4x,'repeat',4x,'freeze'/)")
                  if(.not.lneut .and. lpolar)
     x              write(nrite,"(/,/,1x,'atomic characteristics:',
     x              /,/,21x,' site',5x,'name',10x,'mass',8x,
     x              'charge',6x,'polar',7x,'polar2',7x,'repeat',
     x              4x,'freeze'/)")
                  if(lneut)
     x              write(nrite,"(/,/,1x,'atomic characteristics:',/
     x              /,21x,' site',5x,'name',10x,'mass',8x,'charge',
     x              4x,'repeat',4x,'freeze',3x,'chg grp')")
                  
                endif

                
                do isite=1,numsit(itmols)
                  
                  if(ksite.lt.numsit(itmols))then
c     
c     read atom name, site number, mass, charge, freeze option
                    
                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 2000

                    call cal_field
     x                   (record,nofields,countfields,starts,ends)
                    atom1=record(starts(1):ends(1))
                    weight=dblstr(record(starts(2):ends(2)),80,idum)
                    charge=dblstr(record(starts(3):ends(3)),80,idum)
                    if (lpolar) then
                       polar=dblstr(record(starts(4):ends(4)),80,idum)
                      polar2=dblstr(record(starts(5):ends(5)),80,idum)
                       nrept=intstr(record(starts(6):ends(6)),80,idum)
                       ifrz =intstr(record(starts(7):ends(7)),80,idum)
                       neugp=intstr(record(starts(8):ends(8)),80,idum)
c sandeep: if any variable is undefined, it's value is set to zero or one 
                       if(starts(4)==0) polar=0.d0
                       if(starts(5)==0) polar2=0.d0
                       if(starts(6)==0) nrept=1
                       if(starts(7)==0) ifrz=0
                       if(starts(8)==0) neugp=0
                    else
                       nrept=intstr(record(starts(4):ends(4)),80,idum)
                       ifrz =intstr(record(starts(5):ends(5)),80,idum)
                       neugp=intstr(record(starts(6):ends(6)),80,idum)
                       if(starts(4)==0) nrept=1
                       if(starts(5)==0) ifrz=0
                       if(starts(6)==0) neugp=0
                    endif
                    
c                    call strip(atom1,8)
                    if(nrept.eq.0)nrept=1
                    ksite=ksite+nrept
                    
                    if(idnode.eq.0) then
                      
                      if(.not.lneut .and. .not.lpolar) then

                        write(nrite,
     x                    "(21x,i5,5x,a8,2f12.5,2i10)")
     x                    nsite+1,atom1,weight,charge,nrept,
     x                    ifrz

                      elseif(.not.lneut .and. lpolar) then

                        write(nrite,
     x                    "(21x,i5,5x,a8,4f12.5,2i10)")
     x                    nsite+1,atom1,weight,charge,polar,polar2,
     x                    nrept,ifrz

                      else

                        write(nrite,
     x                    "(21x,i5,5x,a8,2f12.5,3i10)")
     x                    nsite+1,atom1,weight,charge,nrept,
     x                    ifrz,neugp

                      endif

                    endif
                    
                    
                    do irept=1,nrept
                      
                      nsite=nsite+1
                      if(nsite.gt.mxsite) call error(idnode,20)
                      
c mbnrg
c                    if(mbnrg_index(itmols)>0)
c     x                  mbnrg_atname(itmols,nsite)=atom1 
c                    endif 
                      sitnam(nsite)=atom1
                      wgtsit(nsite)=weight
                      chgsit(nsite)=charge
                      if(lpolar) then
                         polarsit(nsite)=polar
                         polarsit2(nsite)=polar2
                      end if
                      lfzsit(nsite)=ifrz
                      nugrp(nsite)=neugp
                      
                    enddo
                    
c     
c     establish list of unique atom types
                    
                    atmchk=.true.
                    
                    do jsite=1,ntpatm
                      
                      if(atom1.eq.unqatm(jsite)) then
                        
                        atmchk=.false.
                        do irept=nsite,nsite-nrept+1,-1
                          
                          ltpsit(irept)=jsite
                          
                        enddo
                        
                      endif
                      
                    enddo
                    
                    if(atmchk)then
                      
                      ntpatm=ntpatm+1
                      if(ntpatm.gt.mxsvdw)call error(idnode,14)
                      unqatm(ntpatm)=atom1
                      
                      do irept=nsite,nsite-nrept+1,-1
                        
                        ltpsit(irept)=ntpatm
                        
                      enddo
                      
                    endif
                    
                  endif
                  
                enddo
                
c     read core - shell spring parameters
                
              elseif(record(starts(1):ends(1)).eq.'shell')then
                
                lshl=.true.
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numshl(itmols)=numshl(itmols)+ntmp
                if(idnode.eq.0) then
                  
                  write(nrite,
     x              "(/,1x,'number of core-shell units',5x,i10)")
     x              ntmp
                  write(nrite,
     x              "(/,/,1x,'core-shell details:',/,/,21x,
     x              5x,'index',5x,'index',6x,'parameter')")
                  
                endif
                
                do ishls=1,numshl(itmols)
                  
                  nshels=nshels+1
                  
                  if(nshels.gt.mxtshl) call error(idnode,57)
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1=intstr(record(starts(1):ends(1)),80,idum)
                  iatm2=intstr(record(starts(2):ends(2)),80,idum)
                  lstshl(nshels,1)=iatm1
                  lstshl(nshels,2)=iatm2
                  prmshl(nshels)=
     x              dblstr(record(starts(3):ends(3)),80,idum)
                  if(idnode.eq.0) write(nrite,
     x              "(21x,2i10,f15.4)")
     x              lstshl(nshels,1),lstshl(nshels,2),
     x              prmshl(nshels)

c     
c     test for frozen cores or shells
                  
                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2
                  if(lfzsit(isite1)*lfzsit(isite2).ne.0)
     x              call error(idnode,49)
                  
c     
c     convert energy units to internal units
                  
                  prmshl(nshels)=prmshl(nshels)*engunit
                  
                enddo
c     
c     read chemical bond force constant and bondlength
                
              elseif(record(starts(1):ends(1)).eq.'bonds')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numbonds(itmols)=numbonds(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of chemical bonds',
     x              7x,i10)") ntmp
                  write(nrite,"(/,/,1x,'chemical bond details:',
     x              /,/,21x,7x,'key',5x,'index',5x,'index',28x,
     x              'parameters', /)")
                endif
                
                ibond1 = numbonds(itmols)
                  
                do ibond=1,ibond1
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(2):ends(2)),80,idum)
                  iatm2 = intstr(record(starts(3):ends(3)),80,idum)
c     
c     test for frozen atom pairs

                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2
                  
                  if(lfzsit(isite1)*lfzsit(isite2).ne.0) then

                    numbonds(itmols) = numbonds(itmols) -1
                    if(idnode.eq.0) write(nrite,'(12x,2a)')
     x                '*** frozen *** ',record(1:40)

                  else

                    nbonds=nbonds+1
                    if(nbonds.gt.mxtbnd) call error(idnode,30)
                    
                    keyword=record(starts(1):ends(1))
                    call lowcase(keyword,4)
c                    call strip(keyword,4)
                    if    (keyword.eq.'harm') then
                      keybnd(nbonds)=1
                    elseif(keyword.eq.'-hrm') then
                      keybnd(nbonds)=-1
                    elseif(keyword.eq.'mors') then
                      keybnd(nbonds)=2
                    elseif(keyword.eq.'-mrs') then
                      keybnd(nbonds)=-2
                    elseif(keyword.eq.'12-6') then
                      keybnd(nbonds)=3
                    elseif(keyword.eq.'-126') then
                      keybnd(nbonds)=-3
                    elseif(keyword.eq.'rhrm') then
                      keybnd(nbonds)=4
                    elseif(keyword.eq.'-rhm') then
                      keybnd(nbonds)=-4
                    elseif(keyword.eq.'quar') then
                      keybnd(nbonds)=5
                    elseif(keyword.eq.'-qur') then
                      keybnd(nbonds)=-5
                    elseif(keyword.eq.'expo') then
                      keybnd(nbonds)=7
                    elseif(keyword.eq.'-exp') then
                      keybnd(nbonds)=-7
                    elseif(keyword.eq.'qmrs') then
                      keybnd(nbonds)=8
                    else
                      if(idnode.eq.0) write(nrite,*) record
                      call error(idnode,444)
                    endif
                    
                    lstbnd(nbonds,1)= iatm1
                    lstbnd(nbonds,2)= iatm2
                    prmbnd(nbonds,1)=
     x                dblstr(record(starts(4):ends(4)),80,idum)
                    prmbnd(nbonds,2)=
     x                dblstr(record(starts(5):ends(5)),80,idum)
                    prmbnd(nbonds,3)=
     x                dblstr(record(starts(6):ends(6)),80,idum)
                    prmbnd(nbonds,4)=
     x                dblstr(record(starts(7):ends(7)),80,idum)
                    
                    do i=1,mega
                      if(starts(i)==0) then
                        prmbnd(nbonds,i-3:)=0.d0
                        exit
                      endif
                    enddo
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(27x,a4,2i10,2x,1p,10e15.6)")
     x                keyword,lstbnd(nbonds,1),
     x                lstbnd(nbonds,2),(prmbnd(nbonds,j),j=1,mxpbnd)
                    
c     
c     convert energy units to internal units
                    
                    if(abs(keybnd(nbonds)).eq.3) then
                      prmbnd(nbonds,2)=prmbnd(nbonds,2)*engunit
                    endif
                    if(abs(keybnd(nbonds)).eq.5) then
                      prmbnd(nbonds,3)=prmbnd(nbonds,3)*engunit
                      prmbnd(nbonds,4)=prmbnd(nbonds,4)*engunit
                    endif
                    
                    prmbnd(nbonds,1)=prmbnd(nbonds,1)*engunit
                    
                  endif

                enddo
                
c     
c     read bond atom indices and constraint bondlength
                
              elseif(record(starts(1):ends(1)).eq.'constraints')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numcon(itmols)=numcon(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of bond constraints',
     x              5x,i10)") ntmp
                  write(nrite,"(/,/,1x,'constraint bond details:',
     x              /,/,21x,5x,'index',5x,'index',2x,'bondlength',/)
     x              ")
                endif

                icnst1 = numcon(itmols)
                do icnst=1,icnst1
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(1):ends(1)),80,idum)
                  iatm2 = intstr(record(starts(2):ends(2)),80,idum)
c     
c     test for frozen atom pairs

                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2

                  if(lfzsit(isite1)*lfzsit(isite2).ne.0) then
                    
                    numcon(itmols) = numcon(itmols) -1
                    if(idnode.eq.0) write(nrite,'(14x,3a)')
     x                'frozen ',record(1:10),' ****'

                  else

                    nconst=nconst+1
                    
                    if(nconst.gt.mxtcon) call error(idnode,40)

                    lstcon(nconst,1)= iatm1
                    lstcon(nconst,2)= iatm2
                    prmcon(nconst)=
     x                dblstr(record(starts(3):ends(3)),80,idum)
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(21x,2i10,f12.6)")
     x                lstcon(nconst,1),lstcon(nconst,2),
     x                prmcon(nconst)
                    
                  endif
                  
                enddo
                
c     
c     read pmf bond atom indices, weights and constraint bondlength
                
              elseif(record(starts(1):ends(1)).eq.'pmf')then
                
                if(lpmf) call error(idnode,484)
                lpmf = .true.
                numpmf(itmols)=1

                prmpmf = dblstr(record(starts(2):ends(2)),80,idum)

                if(idnode.eq.0) then
                  write(nrite,"(/,1x,' PMF      bondlength :',
     x              5x,f20.10)") prmpmf
                  write(nrite,
     x              "(/,/,12x,'unit, site and weight details:'
     x              ,/,/,16x,'unit',6x,'index',5x,'weight')")
                endif

                do ipmf = 1,2

                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000
c                  call strip(record,40)
                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  call lowcase(record,40)

                  tmpchr1=record(starts(1):ends(1))
                  tmpchr2=record(starts(2):ends(2))
                  tmpchr3=tmpchr1//' '//tmpchr2
                   

                  if(tmpchr3.ne.'pmf unit') call error(idnode,462)
                  npmf(ipmf) =intstr(record(starts(3):ends(3)),80,idum)

                  do jpmf=1,npmf(ipmf)
                    
                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 2000

                    call cal_field
     x                   (record,nofields,countfields,starts,ends)

                    iatm1 = intstr(record(starts(1):ends(1)),80,idum)
                    wght  = dblstr(record(starts(2):ends(2)),80,idum)
                    if(wght.eq.0.d0) wght = 1.d0
                    
                    nspmf=nspmf+1
                    
                    if(nspmf.gt.mxspmf) call error(idnode,460)

                    indpmf(nspmf) = iatm1
                    pmfwght(nspmf) = wght
                    
                    if(idnode.eq.0) then

                      if(jpmf.eq.1) then
                        write(nrite,"(16x,i5,i10,f12.6)")
     x                  ipmf,indpmf(nspmf),pmfwght(nspmf)
                      else
                        write(nrite,"(21x,i10,f12.6)")
     x                  indpmf(nspmf),pmfwght(nspmf)
                      endif

                    endif

                  enddo
                  
                enddo
                
c     
c     read intramolecular angular potential parameters
                
              elseif(record(starts(1):ends(1)).eq.'angles')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numang(itmols)=numang(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of bond angles',
     x              10x,i10)") ntmp
                  write(nrite,"(/,/,1x,'bond angle details:',
     x              /,/,21x,7x,'key',5x,'index',5x,'index',5x,
     x              'index',5x,'f-const',7x,'angle',/)")
                endif
                
                iang1 = numang(itmols)
                do iang=1,iang1
                  
c     
c     read bond angle potential parameters
                  
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(2):ends(2)),80,idum)
                  iatm2 = intstr(record(starts(3):ends(3)),80,idum)
                  iatm3 = intstr(record(starts(4):ends(4)),80,idum)
c     
c     test for frozen atom pairs

                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2
                  isite3 = nsite - numsit(itmols) + iatm3

                  if(lfzsit(isite1)*lfzsit(isite2)*
     x              lfzsit(isite3).ne.0) then
                    
                    numang(itmols) = numang(itmols) -1
                    if(idnode.eq.0) write(nrite,'(14x,3a)')
     x                'frozen ',record(1:20),' ****'

                  else

                    nangle=nangle+1
                    
                    if(nangle.gt.mxtang) call error(idnode,50)
                    
                    keyword=record(starts(1):ends(1))
                    call lowcase(keyword,4)
c                    call strip(keyword,4)
                    if    (keyword.eq.'harm') then
                      keyang(nangle)=1
                    elseif(keyword.eq.'-hrm') then
                      keyang(nangle)=-1
                    elseif(keyword.eq.'quar') then
                      keyang(nangle)=2
                    elseif(keyword.eq.'-qur') then
                      keyang(nangle)=-2
                    elseif(keyword.eq.'thrm') then
                      keyang(nangle)=3
                    elseif(keyword.eq.'-thm') then
                      keyang(nangle)=-3
                    elseif(keyword.eq.'shrm') then
                      keyang(nangle)=4
                    elseif(keyword.eq.'-shm') then
                      keyang(nangle)=-4
                    elseif(keyword.eq.'bvs1') then
                      keyang(nangle)=5
                    elseif(keyword.eq.'-bv1') then
                      keyang(nangle)=-5
                    elseif(keyword.eq.'bvs2') then
                      keyang(nangle)=6
                    elseif(keyword.eq.'-bv2') then
                      keyang(nangle)=-6
                    elseif(keyword.eq.'hcos') then
                      keyang(nangle)=7
                    elseif(keyword.eq.'-hcs') then
                      keyang(nangle)=-7
                    elseif(keyword.eq.'cos ') then
                      keyang(nangle)=8
                    elseif(keyword.eq.'-cos') then
                      keyang(nangle)=-8
                    elseif(keyword.eq.'schw') then
                      keyang(nangle)=9
                      if(.not.lpspes) call init_ps_type_pes
                      lpspes=.true.
                    elseif(keyword.eq.'-sch') then
                      keyang(nangle)=-9
                      lpspes=.true.
                      if(.not.lpspes) call init_ps_type_pes
                    else
                      if(idnode.eq.0) write(nrite,*) record
                      call error(idnode,440)
                    endif

                    lstang(nangle,1)=iatm1
                    lstang(nangle,2)=iatm2
                    lstang(nangle,3)=iatm3
                    prmang(nangle,1)=
     x                dblstr(record(starts(5):ends(5)),80,idum)
                    prmang(nangle,2)=
     x                dblstr(record(starts(6):ends(6)),80,idum)
                    prmang(nangle,3)=
     x                dblstr(record(starts(7):ends(7)),80,idum)
                    prmang(nangle,4)=
     x                dblstr(record(starts(8):ends(8)),80,idum)
                    
                    do i=1,mega
                      if(starts(i)==0) then
                        prmang(nangle,i-4:)=0.d0
                        exit
                      endif
                    enddo
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(27x,a4,3i10,1p,e12.4,0p,9f12.6)")
     x                keyword,(lstang(nangle,ia),ia=1,3),
     x                (prmang(nangle,ja),ja=1,mxpang)
c     
c     convert energies to internal units
                    
                    prmang(nangle,1) = prmang(nangle,1)*engunit
                    if(abs(keyang(nangle)).eq.2) then
                      prmang(nangle,3) = prmang(nangle,3)*engunit
                      prmang(nangle,4) = prmang(nangle,4)*engunit
                    endif
c     
c     convert angles to radians
                    
                    prmang(nangle,2)=prmang(nangle,2)*(pi/180.d0)
                    

                  endif

                enddo
                
c     
c     read intramolecular dihedral potential parameters
                
              elseif(record(starts(1):ends(1)).eq.'dihedrals')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numdih(itmols)=numdih(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of dihedral angles',
     x              6x,i10)") ntmp
                  write(nrite,"(/,/,1x,'dihedral angle details:',
     x              /,/,21x,7x,'key',5x,'index',5x,'index',5x,
     x              'index',5x,'index',5x,'f-const',7x,'angle',
     x              8x,'trig',4x,'1-4 elec',5x,'1-4 vdw',/)")
                endif
                
                idih1 = numdih(itmols)
                do idih=1,idih1
                  
c     
c     read dihedral bond angle potential parameters
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(2):ends(2)),80,idum)
                  iatm2 = intstr(record(starts(3):ends(3)),80,idum)
                  iatm3 = intstr(record(starts(4):ends(4)),80,idum)
                  iatm4 = intstr(record(starts(5):ends(5)),80,idum)
c     
c     test for frozen atom pairs

                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2
                  isite3 = nsite - numsit(itmols) + iatm3
                  isite4 = nsite - numsit(itmols) + iatm4

                  if(lfzsit(isite1)*lfzsit(isite2)*
     x              lfzsit(isite3)*lfzsit(isite4).ne.0) then
                    
                    numdih(itmols) = numdih(itmols) -1
                    if(idnode.eq.0) write(nrite,'(14x,3a)')
     x                'frozen ',record(1:25),' ****'

                  else

                    ndihed=ndihed+1
                    
                    if(ndihed.gt.mxtdih) call error(idnode,60)

                    keyword=record(starts(1):ends(1))
                    call lowcase(keyword,4)
c                    call strip(keyword,4)

                    if(keyword.eq.'cos ') then
                      keydih(ndihed)=1
                    elseif(keyword.eq.'harm') then
                      keydih(ndihed)=2
                    elseif(keyword.eq.'hcos') then
                      keydih(ndihed)=3
                    elseif(keyword.eq.'cos3') then
                      keydih(ndihed)=4
                    elseif(keyword.eq.'ryck') then
                      keydih(ndihed)=5
                    elseif(keyword.eq.'rbf') then 
                      keydih(ndihed)=6
                    else
                      if(idnode.eq.0) write(nrite,*) record
                      call error(idnode,448)
                    endif

                    lstdih(ndihed,1)=iatm1
                    lstdih(ndihed,2)=iatm2
                    lstdih(ndihed,3)=iatm3
                    lstdih(ndihed,4)=iatm4
                    prmdih(ndihed,1)=
     x                dblstr(record(starts(6):ends(6)),80,idum)
                    prmdih(ndihed,2)=
     x                dblstr(record(starts(7):ends(7)),80,idum)
                    prmdih(ndihed,3)=
     x                dblstr(record(starts(8):ends(8)),80,idum)
                    prmdih(ndihed,4)=
     x                dblstr(record(starts(9):ends(9)),80,idum)
                    prmdih(ndihed,5)=
     x                dblstr(record(starts(10):ends(10)),80,idum)

                    do i=1,mega
                      if(starts(i)==0) then
                        prmdih(ndihed,i-5:)=0.d0
                        exit
                      endif
                    enddo
                    
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(27x,a4,4i10,1p,e12.4,0p,9f12.6)")
     x                keyword,(lstdih(ndihed,ia),ia=1,4),
     x                (prmdih(ndihed,ja),ja=1,mxpdih)
c     
c     convert energies to internal units and angles to radians
                    
                    prmdih(ndihed,1)=prmdih(ndihed,1)*engunit

                    if(keydih(ndihed).eq.4)then

                      prmdih(ndihed,2)=prmdih(ndihed,2)*engunit
                      prmdih(ndihed,3)=prmdih(ndihed,3)*engunit

                    else

                      prmdih(ndihed,2)=prmdih(ndihed,2)*(pi/180.d0)

                    endif
                    
                  endif

                enddo
                
c     
c     read intramolecular inversion potential parameters
                
              elseif(record(starts(1):ends(1)).eq.'inversions')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numinv(itmols)=numinv(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of inversion terms',
     x              6x,i10)") ntmp
                  write(nrite,"(/,/,1x,'inversion potential details:',
     x              /,/,21x,7x,'key',5x,'index',5x,'index',5x,
     x              'index',5x,'index',5x,'f-const',7x,'angle',/)")
                endif
                
                inv1 = numinv(itmols)
                do inv=1,inv1
                  
c     
c     read inversion potential parameters
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(2):ends(2)),80,idum)
                  iatm2 = intstr(record(starts(3):ends(3)),80,idum)
                  iatm3 = intstr(record(starts(4):ends(4)),80,idum)
                  iatm4 = intstr(record(starts(5):ends(5)),80,idum)
c     
c     test for frozen atom pairs

                  isite1 = nsite - numsit(itmols) + iatm1
                  isite2 = nsite - numsit(itmols) + iatm2
                  isite3 = nsite - numsit(itmols) + iatm3
                  isite4 = nsite - numsit(itmols) + iatm4

                  if(lfzsit(isite1)*lfzsit(isite2)*
     x              lfzsit(isite3)*lfzsit(isite4).ne.0) then
                    
                    numinv(itmols) = numinv(itmols) -1
                    if(idnode.eq.0) write(nrite,'(14x,3a)')
     x                'frozen ',record(1:25),' ****'

                  else

                    ninver=ninver+1
                    
                    if(ninver.gt.mxtinv) call error(idnode,73)

                    keyword=record(starts(1):ends(1))
                    call lowcase(keyword,4)
c                    call strip(keyword,4)

                    if(keyword.eq.'harm') then
                      keyinv(ninver)=1
                    elseif(keyword.eq.'hcos') then
                      keyinv(ninver)=2
                    elseif(keyword.eq.'plan') then
                      keyinv(ninver)=3
                    else
                      if(idnode.eq.0) write(nrite,*) record
                      call error(idnode,449)
                    endif

                    lstinv(ninver,1)=iatm1
                    lstinv(ninver,2)=iatm2
                    lstinv(ninver,3)=iatm3
                    lstinv(ninver,4)=iatm4
                    prminv(ninver,1)=
     x                dblstr(record(starts(6):ends(6)),80,idum)
                    prminv(ninver,2)=
     x                dblstr(record(starts(7):ends(7)),80,idum)
                    do i=1,mega
                      if(starts(i)==0) then
                        prminv(ninver,i-5:)=0.d0
                        exit
                      endif
                    enddo
                    
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(27x,a4,4i10,1p,e12.4,0p,9f12.6)")
     x                keyword,(lstinv(ninver,ia),ia=1,4),
     x                (prminv(ninver,ja),ja=1,mxpinv)
c     
c     convert energies to internal units and angles to radians
                    
                    prminv(ninver,1)=prminv(ninver,1)*engunit

                    if(keyinv(ninver).eq.2)then

                      prminv(ninver,2)=cos(prminv(ninver,2)*(pi/180.d0))

                    endif
                    
                  endif

                enddo
                
c     
c     read rigid body data
              elseif(record(starts(1):ends(1)).eq.'rigid')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numgrp(itmols)=numgrp(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of rigid units    ',
     x              6x,i10)") ntmp
                  write(nrite,"(/,' rigid body details:',/,/,21x,
     x              6x,'unit',3x,'indices',/) ")
                endif
                
                do igrp=1,numgrp(itmols)
                  
                  ngrp=ngrp+1
                  
                  if(ngrp.gt.mxungp) call error(idnode,301)
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000
                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  numgsit(ngrp)=intstr(record(starts(1):ends(1)),5,idum)
                  
                  if(numgsit(ngrp).gt.mxngp) 
     x              call error (idnode,302)
                  
                  listyp(ngrp) = ngrp
                  do j=1,min(15,numgsit(ngrp))
                    
                    lstgst(ngrp,j)=
     x                intstr(record(starts(j+1):ends(j+1)),80,idum)
                    
                  enddo
                  do j=16,numgsit(ngrp),16
                    call getrec(safe,idnode,mxnode,nfield,record)
                    if(.not.safe)go to 2000
                    call cal_field
     x                   (record,nofields,countfields,starts,ends)

                    k=1
                    do m=j,min(numgsit(ngrp),j+15)
                      
                      lstgst(ngrp,m)=
     x                  intstr(record(starts(k):ends(k)),80,idum)
                      k=k+1
                      
                    enddo
                  enddo
                  
                  if(idnode.eq.0) 
     x              write(nrite,"(21x,10i10,100(/,21x,10i10))")
     x              listyp(ngrp),(lstgst(ngrp,j),j=1,
     x              numgsit(ngrp))
                  
                enddo
c     
c     read tethered atom indices and tethering parameters
                
              elseif(record(starts(1):ends(1)).eq.'teth')then
                
                ntmp=intstr(record(starts(2):ends(2)),80,idum)
                numteth(itmols)=numteth(itmols)+ntmp
                if(idnode.eq.0) then
                  write(nrite,"(/,1x,'number of tethered atoms ',
     x              6x,i10)") ntmp
                  write(nrite,"(/,' tethered atom details:',/,/,
     x              21x,7x,'key',6x,'atom',19x,'parameters',/) ")
                endif

                iteth1 = numteth(itmols)
                do iteth=1,iteth1
                  
                  call getrec(safe,idnode,mxnode,nfield,record)
                  if(.not.safe)go to 2000
c                  call strip(record,40)

                  call cal_field
     x                 (record,nofields,countfields,starts,ends)

                  iatm1 = intstr(record(starts(2):ends(2)),80,idum)
c     
c     test for frozen atom 

                  isite1 = nsite - numsit(itmols) + iatm1

                  
                  if(lfzsit(isite1).ne.0) then

                    numteth(itmols) = numteth(itmols) -1
                    if(idnode.eq.0) write(nrite,'(12x,2a)')
     x                '*** frozen *** ',record(1:40)

                  else

                    nteth=nteth+1
                    if(nteth.gt.mxteth) call error(idnode,62)

                    keyword=record(starts(1):ends(1))
                    call lowcase(keyword,4)
c                    call strip(keyword,4)
                    
                    if(keyword.eq.'harm') then
                      keytet(nteth) = 1
                    elseif(keyword.eq.'rhrm') then
                      keytet(nteth) = 2
                    elseif(keyword.eq.'quar') then
                      keytet(nteth) = 3
                    else
                      if(idnode.eq.0) write(nrite,*) record
                      call error(idnode,450)
                    endif

                    lsttet(nteth)=iatm1
                    prmtet(nteth,1)=
     x                dblstr(record(starts(3):ends(3)),80,idum)
                    prmtet(nteth,2)=
     x                dblstr(record(starts(4):ends(4)),80,idum)
                    prmtet(nteth,3)=
     x                dblstr(record(starts(5):ends(5)),80,idum)

                    do i=1,mega
                      if(starts(i-2)==0) then
                        prmtet(nteth,i:)=0.d0
                        exit
                      endif
                    enddo
                    
                    if(idnode.eq.0) 
     x                write(nrite,"(27x,a4,i10,1p,9e12.4)")
     x                keyword,lsttet(nteth),
     x                (prmtet(nteth,j),j=1,mxpbnd)
c
c     convert energy units to internal units
                    
                    if(abs(keytet(nteth)).eq.1) then
                      prmtet(nteth,1)=prmtet(nteth,1)*engunit
                    elseif(abs(keytet(nteth)).eq.2) then
                      prmtet(nteth,1)=prmtet(nteth,1)*engunit
                    elseif(abs(keytet(nteth)).eq.3) then
                      prmtet(nteth,1)=prmtet(nteth,1)*engunit
                      prmtet(nteth,2)=prmtet(nteth,2)*engunit
                      prmtet(nteth,3)=prmtet(nteth,3)*engunit
                    endif

                  endif

                enddo
                
c     
c     finish of data for one molecular type
                
              elseif(record(starts(1):ends(1)).eq.'finish')then
                
c     
c     running total of number of atoms in system
                
                natms=natms+nummols(itmols)*numsit(itmols)
                if(natms.gt.mxatms) call error(idnode,75)
                
c     
c     check core-shell units are not both in same rigid body unit

                nshels=nshels-numshl(itmols)
                do k1 = 1,numshl(itmols)

                  nshels=nshels+1
                  ia = lstshl(nshels,1)
                  ib = lstshl(nshels,2)

                  ngrp = ngrp - numgrp(itmols)

                  do kk = 1,numgrp(itmols)
                    
                    ngrp = ngrp + 1
                    id = listyp(ngrp)
                    
                    do jj = 1,numgsit(id)-1
                      
                      ia1 = lstgst(ngrp,jj)
                      if(ia1.eq.ia) then

                        do jk = jj,numgsit(id)
                          
                          ib1 = lstgst(ngrp,jk)
                          if(ib1.eq.ib) then 
                            
                            if(idnode.eq.0)write(nrite,'(/,13x,a,2i10)')
     x                        'error: sites ',ia,ib
                            call error(idnode,456)

                          endif

                        enddo

                      elseif(ia1.eq.ib) then

                        do jk = jj,numgsit(id)
                          
                          ib1 = lstgst(ngrp,jk)
                          if(ib1.eq.ia) then 
                            
                            if(idnode.eq.0)write(nrite,'(/,13x,a,2i10)')
     x                        'error: sites ',ia,ib
                            call error(idnode,456)

                          endif

                        enddo

                      endif

                    enddo
                  enddo
                enddo

                go to 1000
                
              else
c     
c     error exit for unidentified directive in molecular data
                
                if(idnode.eq.0) write(nrite,'(12x,a)') record
                call error(idnode,12)
                
              endif
              
            enddo
            
 1000       continue
            
          enddo
          
c     
c     calculate system charge
          
          jsite=0
          do itmols=1,ntpmls
            
            do lsite=1,numsit(itmols)
              
              jsite=jsite+1
              sumchg=sumchg+dble(nummols(itmols))*chgsit(jsite)
              
            enddo
            
          enddo
          
          if(abs(sumchg).gt.1.0d-6) then
            
c           call error(idnode,90)
            call warning(idnode,60,sumchg,0.d0,0.d0)
            
          endif
          
          
c     
c     read in the nonbonded potential energy parameters
c   taking care of 'vdw' and also, 'vdwtable'
          
        elseif((record(starts(1):ends(1)).eq.'vdw').or.
     x    (record(starts(1):ends(1)).eq.'vdwtable')) then
          
          novdw=.false.
          ntpvdw=intstr(record(starts(2):ends(2)),80,idum)
          ltable=(record(starts(1):ends(1)).eq.'vdwtable')
c          ltable=(record(4:8).eq.'table')
          
          if(idnode.eq.0) then
            
            write(nrite,"(/,/,1x,'number of specified pair ',
     x        'potentials',i10)") ntpvdw
            write(nrite,"(/,/,16x,'atom 1  ','atom 2  ',3x,
     x        ' key',30x,'parameters'/,/)")
            
          endif      

          if(ntpvdw.gt.mxvdw) call error(idnode,80)
          if(.not.lunits) call error(idnode,6)
          if(.not.lmols) call error(idnode,13)
          
          do ivdw=1,mxvdw
            lstvdw(ivdw)=0
          enddo
          
          do itpvdw=1,ntpvdw
            
            do i=1,mxpvdw
              parpot(i)=0.d0
            enddo
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000
            call cal_field(record,nofields,countfields,starts,ends)
            atom1=record(starts(1):ends(1))
            atom2=record(starts(2):ends(2))

            keyword=record(starts(3):ends(3))
            call lowcase(keyword,4)
c            call strip(keyword,4)

            if(keyword.eq.'12-6') then
              keypot = 1
            elseif(keyword.eq.'lj  ') then
              keypot = 2
            elseif(keyword.eq.'nm  ') then
              keypot = 3
            elseif(keyword.eq.'buck') then
              keypot = 4
            elseif(keyword.eq.'bhm ') then
              keypot = 5
            elseif(keyword.eq.'hbnd') then
              keypot = 6
            elseif(keyword.eq.'snm ') then
              keypot = 7
            elseif(keyword.eq.'hcnm') then
              keypot = 7
            elseif(keyword.eq.'mors') then
              keypot = 8
            elseif(keyword.eq.'lj-8') then
              keypot = 9
            elseif(keyword.eq.'ljex') then
              keypot = 10
            elseif(keyword.eq.'ttm4') then
              keypot = 11
            elseif(keyword.eq.'tt68') then
              keypot = 12
            elseif(keyword.eq.'btt6') then
              keypot = 13
            elseif(keyword.eq.'tab ') then
              keypot = 0
            elseif(keyword.eq.'stch') then
              keypot = 100
            else
              if(idnode.eq.0) write(nrite,*) record
              call error(idnode,452)
            endif

            do i=4,countfields
           
               parpot(i-3)=dblstr(record(starts(i):ends(i)),80,idum)

            enddo

            do i=1,mega
              if(starts(i)==0) then
                parpot(i-3:)=0.d0
                exit
              endif
            enddo
                    
c            call strip(atom1,8)
c            call strip(atom2,8)
            if(idnode.eq.0) 
     x        write(nrite,"(16x,2a8,2x,a4,3x,1p,13e13.5)") 
     x        atom1,atom2,keyword,(parpot(j),j=1,mxpvdw)

            katom1=0
            katom2=0
            
            do jtpatm=1,ntpatm
              
              if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
              if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
              
            enddo
            
            if(katom1.eq.0.or.katom2.eq.0) then
              call  error(idnode,81)
            endif
            
            keyvdw=(max(katom1,katom2)*(max(katom1,katom2)-1))/2+
     x        min(katom1,katom2)
c     
c     convert energies to internal unit

            if(keyvdw.gt.mxvdw) call error(idnode,82)
            
            parpot(1) = parpot(1)*engunit
            
            if(keypot.eq.1) then
              
              parpot(2) = parpot(2)*engunit
              
            else if(keypot.eq.4) then
              
              parpot(3) = parpot(3)*engunit
              
            else if(keypot.eq.5) then
              
              parpot(4) = parpot(4)*engunit
              parpot(5) = parpot(5)*engunit
              
            else if(keypot.eq.6) then
              
              parpot(2) = parpot(2)*engunit

            else if(keypot.eq.9) then

              parpot(2) = parpot(2)*engunit   
              parpot(3) = parpot(3)*engunit
              parpot(4) = parpot(4)*engunit 

            else if(keypot.eq.10) then

              parpot(3) = parpot(3)*engunit  ! C
              parpot(4) = parpot(4)*engunit  ! D (see above for A)

            else if(keypot.eq.11) then

              parpot(2) = parpot(2)*engunit
              parpot(3) = parpot(3)*engunit
              parpot(4) = parpot(4)*engunit
              parpot(5) = parpot(5)*engunit
              parpot(6) = parpot(6)*engunit
              
            else if(keypot.eq.12) then      !TT-damped 6/8 dispersion
              !(a,b,c,d) => -a/r^6*TT(b*r)-c/r^8*TT(d*r)
              parpot(3) = parpot(3)*engunit

            else if(keypot.eq.13) then      !buck + RT-damped 6 dispersion
              !(a,b,c,d) => a*exp(-r*b) - c/r^6*TT(d*r)
              parpot(3) = parpot(3)*engunit

            else if(keypot.eq.100) then
              
              lmetal=.true.
              
            endif

            if(lmetal.and.(2*ntpvdw.ge.mxvdw)) call error(idnode,71)
            
            ltable=(ltable.or.(keypot.eq.0))

            if(lstvdw(keyvdw).ne.0) call error(idnode,15)
            lstvdw(keyvdw)=itpvdw
            ltpvdw(itpvdw)=keypot
            
            do i=1,mxpvdw
              
              prmvdw(itpvdw,i)=parpot(i)
              
            enddo
            
          enddo

c     
c     generate metal force arrays

          if(lmetal)then

            call metgen
     x        (idnode,ntpvdw,ntpatm,dlrpot,rvdw,lstvdw,ltpvdw,
     x        prmvdw,vvv,ggg)

          endif
c     
c     generate nonbonded force arrays

          if((ntpvdw.gt.0.and.mod(keyfce,2).eq.1).or.(keyfce.eq.2))
     x      then
            
            call forgen
     x        (ltable,idnode,ntpvdw,dlrpot,rvdw,ltpvdw,
     x        prmvdw,vvv,ggg)
            
            if(ltable)then
              
              call fortab
     x          (idnode,ntpvdw,ntpatm,mxnode,dlrpot,rvdw,engunit,
     x          unqatm,lstvdw,ltpvdw,prmvdw,vvv,ggg,buffer)

            endif
            
          endif
          
c     
c     check for unspecified atom-atom potentials
          
          ntab =(ntpatm*(ntpatm+1))/2
          
          if(ntpvdw.lt.ntab) then
            
            call warning(idnode,110,0.d0,0.d0,0.d0)

            if(mxvdw.le.ntpvdw) call error(idnode,82)

            do i = 1,ntab

! bug fix              
              if(lstvdw(i).eq.0) lstvdw(i)=ntpvdw+1

!              if(lstvdw(i).eq.0) then 

!                 lstvdw(i)=ntpvdw+1
!                 ltpvdw(i)=-1

!              endif              
               
            enddo
c     
c     define zero potential for undefined interactions
            
            do i = 1,mxgrid
              
              ggg(i,ntpvdw+1) = 0.d0
              vvv(i,ntpvdw+1) = 0.d0
              
            enddo
            
          endif
c     
c     generate error function complement tables for ewald sum

          if((keyfce/2.eq.1.or.keyfce/2.eq.6) .and. .not.lpolar) 
     x      call erfcgen(keyfce,alpha,drewd,rcut,erc,fer)

          if((keyfce/2.eq.1.or.keyfce/2.eq.6) .and. lpolar)
     x      call erfcgenp(keyfce,alpha,drewd,rcut,ercp)
c
c     generate screening function tables for hautman-klein-ewald

          if(keyfce/2.eq.7) call hkgen
     x      (idnode,nhko,nlatt,alpha,drewd,rcut,ahk,hon,fon,dhn)
c     
c     read in the three body potential energy parameters
          
        elseif(record(starts(1):ends(1)).eq.'tbp') then
          
          ntptbp=intstr(record(starts(2):ends(2)),80,idum)
          
          if(idnode.eq.0) then
            
            write(nrite,"(/,/,1x,'number of specified three ',
     x        'body potentials',i10)") ntptbp
            write(nrite,"(/,/,16x,'atom 1  ','atom 2  ','atom 3  ',
     x        3x,' key',30x,'parameters'/,/)")
            
          endif      
          if(ntptbp.gt.mxtbp) call error(idnode,83)
          if(.not.lunits) call error(idnode,6)
          if(.not.lmols) call error(idnode,13)
          
          do itbp=1,mxtbp
            lsttbp(itbp)=0
          enddo
          
          do itbp=1,mxtbp,mx2tbp
            lsttbp(itbp)=-1
          enddo
          
          rcuttb=0.d0
          
          do itptbp=1,ntptbp
            
            do i=1,mxptbp
              parpot(i)=0.d0
            enddo
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000
            call cal_field(record,nofields,countfields,starts,ends)
c     
c     Note the order!! atom0 is the central atom
            atom1=record(starts(1):ends(1))
            atom0=record(starts(2):ends(2))
            atom2=record(starts(3):ends(3))
            keyword=record(starts(4):ends(4))
            call lowcase(keyword,4)
c            call strip(keyword,4)
            if(keyword.eq.'harm') then
              keypot=0
            elseif(keyword.eq.'thrm') then
              keypot=1
            elseif(keyword.eq.'shrm') then
              keypot=2
            elseif(keyword.eq.'bvs1') then
              keypot=3
            elseif(keyword.eq.'bvs2') then
              keypot=4
            elseif(keyword.eq.'hbnd') then
              keypot=5
            else
              if(idnode.eq.0) write(nrite,*) record
              call error(idnode,442)
            endif

            parpot(1)=dblstr(record(starts(5):ends(5)),80,idum)
            parpot(2)=dblstr(record(starts(6):ends(6)),80,idum)
            parpot(3)=dblstr(record(starts(7):ends(7)),80,idum)
            parpot(4)=dblstr(record(starts(8):ends(8)),80,idum)
            parpot(5)=dblstr(record(starts(9):ends(9)),80,idum)
            
            do i=1,mega
              if(starts(i)==0) then
                parpot(i-4:)=0.d0
                exit
              endif
            enddo
                    
c            call strip(atom0,8)
c            call strip(atom1,8)
c            call strip(atom2,8)
            if(idnode.eq.0) 
     x        write(nrite,"(16x,3a8,4x,a4,1x,1p,9e13.5)") 
     x        atom1,atom0,atom2,keyword,(parpot(j),j=1,mxptbp)
            
            katom0=0
            katom1=0
            katom2=0
            
            do jtpatm=1,ntpatm
              
              if(atom0.eq.unqatm(jtpatm))katom0=jtpatm
              if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
              if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
              
            enddo
            
            if(katom0.eq.0.or.katom1.eq.0.or.katom2.eq.0) 
     x        call error(idnode,84)
            
            keytbp=(max(katom1,katom2)*(max(katom1,katom2)-1))/2+
     x        min(katom1,katom2)+(katom0-1)*mx2tbp
            
            if(keytbp.gt.mxtbp) call error(idnode,86)
            
c     
c     convert parameters to internal units
            
            parpot(1) = parpot(1)*engunit
            if(keypot.ne.5)parpot(2) = parpot(2)*(pi/180.d0)
            
            if(lsttbp(keytbp).gt.0) call error(idnode,18)
            lsttbp(keytbp)=itptbp
            ltptbp(itptbp)=keypot
            ktbp=mx2tbp*((keytbp-1)/mx2tbp)+1
            if(lsttbp(ktbp).lt.0)lsttbp(ktbp)=0
            
c     
c     calculate max three body cutoff
            
            
            rcuttb=max(rcuttb,parpot(5))
            rcut3b(itptbp)=parpot(5)
            
c     
c     store three body potential parameters
            
            do i=1,4
              prmtbp(itptbp,i)=parpot(i)
            enddo
            if(mxptbp.ge.6) then
              do i=6,mxptbp
                prmtbp(itptbp,i-1)=parpot(i-1)
              enddo
            endif
          enddo

          if(rcuttb.lt.1.d-6)call error(idnode,451)
          
c     
c     read in the four body potential energy parameters
          
        elseif(record(starts(1):ends(1)).eq.'fbp') then
          
          ntpfbp=intstr(record(starts(2):ends(2)),80,idum)
          
          if(idnode.eq.0) then
            
            write(nrite,"(/,/,1x,'number of specified four ',
     x        'body potentials',i10)") ntpfbp
            write(nrite,"(/,/,16x,'atom 1  ','atom 2  ','atom 3  ',
     x        'atom 4  ',3x,' key',30x,'parameters'/,/)")
            
          endif      
          if(ntpfbp.gt.mxfbp) call error(idnode,89)
          if(.not.lunits) call error(idnode,6)
          if(.not.lmols) call error(idnode,13)
          
          do ifbp=1,mxfbp
            lstfbp(ifbp)=0
          enddo
          
          do ifbp=1,mxfbp,mx3fbp
            lstfbp(ifbp)=-1
          enddo
          
          rcutfb=0.d0
          
          do itpfbp=1,ntpfbp
            
            do i=1,mxpfbp
              parpot(i)=0.d0
            enddo
            
            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000
            call cal_field(record,nofields,countfields,starts,ends)
c     
c     Note the order!! atom0 is the central atom
            atom0=record(starts(1):ends(1))
            atom1=record(starts(2):ends(2))
            atom2=record(starts(3):ends(3))
            atom3=record(starts(4):ends(4))
            keyword=record(starts(5):ends(5))
            call lowcase(keyword,4)
c            call strip(keyword,4)
            if(keyword.eq.'harm') then
              keypot=1
            elseif(keyword.eq.'hcos') then
              keypot=2
            elseif(keyword.eq.'plan') then
              keypot=3
            else
              if(idnode.eq.0) write(nrite,*) record
              call error(idnode,443)
            endif

            parpot(1)=dblstr(record(starts(6):ends(6)),80,idum)
            parpot(2)=dblstr(record(starts(7):ends(7)),80,idum)
            parpot(3)=dblstr(record(starts(8):ends(8)),80,idum)
            do i=1,mega
              if(starts(i)==0) then
                parpot(i-5:)=0.d0
                exit
              endif
            enddo
                    
            
c            call strip(atom0,8)
c            call strip(atom1,8)
c            call strip(atom2,8)
c            call strip(atom3,8)
            if(idnode.eq.0) 
     x        write(nrite,"(16x,4a8,4x,a4,1x,1p,9e13.5)") 
     x        atom0,atom1,atom2,atom3,keyword,(parpot(j),j=1,mxpfbp)
            
            katom0=0
            katom1=0
            katom2=0
            katom3=0
            
            do jtpatm=1,ntpatm
              
              if(atom0.eq.unqatm(jtpatm))katom0=jtpatm
              if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
              if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
              if(atom3.eq.unqatm(jtpatm))katom3=jtpatm
              
            enddo
            
            if(katom0.eq.0.or.katom1.eq.0.or.katom2.eq.0.or.
     x         katom3.eq.0) call error(idnode,91)
            
            ka1=max(katom1,katom2,katom3)
            ka3=min(katom1,katom2,katom3)
            ka2=katom1+katom2+katom3-ka1-ka3
            keyfbp=ka3+(ka2*(ka2-1))/2+(ka1*(ka1**2-1))/6+
     x        (katom0-1)*mx3fbp

            if(keyfbp.gt.mxfbp) call error(idnode,101)

c     
c     convert parameters to internal units
            
            parpot(1) = parpot(1)*engunit
            parpot(2) = parpot(2)*(pi/180.d0)

            if(keypot.eq.2)then

              parpot(2)=cos(parpot(2))

            endif

            if(lstfbp(keyfbp).gt.0) call error(idnode,19)
            lstfbp(keyfbp)=itpfbp
            ltpfbp(itpfbp)=keypot
            kfbp=mx3fbp*((keyfbp-1)/mx3fbp)+1
            if(lstfbp(kfbp).lt.0)lstfbp(kfbp)=0
            
c     
c     calculate max four body cutoff
            
            
            rcutfb=max(rcutfb,parpot(3))
            rcut4b(itpfbp)=parpot(3)
            
c     
c     store four body potential parameters
            
            do i=1,mxpfbp
              prmfbp(itpfbp,i)=parpot(i)
            enddo

          enddo

          if(rcutfb.lt.1.d-6)call error(idnode,453)
          
c     
c     read external field data
          
        elseif(record(starts(1):ends(1)).eq.'extern') then
          
          call getrec(safe,idnode,mxnode,nfield,record)
          if(.not.safe)go to 2000
          call cal_field(record,nofields,countfields,starts,ends)
          call lowcase(record,40)
c          call strip(record,40)
          keyword = record(starts(1):ends(1))

          if(keyword.eq.'elec') then
            keyfld =1 
          elseif(keyword.eq.'oshm') then
            keyfld=2
          elseif(keyword.eq.'shrx') then
            keyfld=3
          elseif(keyword.eq.'grav') then
            keyfld=4
          elseif(keyword.eq.'magn') then
            keyfld=5
          elseif(keyword.eq.'sphr') then
            keyfld=6
          elseif(keyword.eq.'zbnd') then
            keyfld=7
          else
            if(idnode.eq.0) write(nrite,*) record
            call error(idnode,454)
          endif

          do i = 1,mxfld
            prmfld(i)=0.d0
          enddo
          
          nfld=intstr(record(starts(2):ends(2)),40,idum)
          if(nfld.eq.0)nfld=5
          do k=1,nfld,5

            call getrec(safe,idnode,mxnode,nfield,record)
            if(.not.safe)go to 2000
            call cal_field(record,nofields,countfields,starts,ends)
            do m=1,5
              
              prmfld(k+m-1)=dblstr(record(starts(m):ends(m)),80,idum)
              
            enddo

          enddo

          
          if(idnode.eq.0) then
            
            write(nrite,"(/,/,1x,'external field key ',13x,a4,
     x        /,/,30x,'external field parameters')") keyword
            write(nrite,"(2(/,1x,1p,5e15.5))") prmfld
            
          endif      
c     
c     convert to internal units
          
          if(keyfld.eq.1.or.keyfld.eq.4.or.keyfld.eq.5) then
            
            if(.not.lunits)call error(idnode,6)
            
            do i = 1,3
              prmfld(i) = prmfld(i)*engunit
            enddo
            
          elseif(keyfld.eq.2.or.keyfld.eq.6.or.keyfld.eq.7) then
            
            prmfld(1) = prmfld(1)*engunit
            
          endif
          
! PIMD/CMD
        elseif(record(starts(1):ends(1)).eq.'umbrel') then

          keyumb=1
          do i=1,5
            prmumb(i)=0.d0
          end do

          call getrec(safe,idnode,mxnode,nfield,record)
          if (.not.safe) go to 2000
          call cal_field(record,nofields,countfields,starts,ends)

          do m=1,5
           prmumb(m)=dblstr(record(starts(m):ends(m)),80,idum)
          end do

          prmumb(1)=prmumb(1)*engunit
c     
c     close force field file
          
        elseif(record(starts(1):ends(1)).eq.'close')then
          
          if(idnode.eq.0) close (nfield)
          if(novdw.and.mod(keyfce,2).eq.1)call error(idnode,145)

#ifdef VAMPIR
          call VTEND(2, ierr)
#endif
          return
          
c     
c     error exit for unidentified directive
          
        else
          
          if(idnode.eq.0) write(nrite,'(a80)') record(1:80)
          call error(idnode,4)
          
        endif
        
      enddo
      
c     
c     uncontrolled error exit from field file procesing
      
      if(idnode.eq.0) close (nfield)
      call error(idnode,16)
      
c     
c     end of field file error exit
      
 2000 continue
      if(idnode.eq.0) close (nfield)
      call error(idnode,52)

      return
      end
