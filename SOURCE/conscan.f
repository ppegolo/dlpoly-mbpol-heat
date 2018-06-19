      subroutine conscan
     x  (lewald,lspme,lhke,idnode,mxnode,nread,imcon,mxstak,kmaxa,
     x  kmaxb,kmaxc,kmaxd,kmaxe,kmaxf,nhko,rcut,rvdw,delr,cell)
c     
c***********************************************************************
c     
c     dl_poly subroutine for scanning the contents of the control file
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith  june       1997
c     
c     wl
c     2001/06/12 12:39:22
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
      implicit real*8(a-h,o-z)

      parameter (mega=10000,pi =3.141592653589793d0)
      
      logical safe,lewald,lspme,lhke
      character*256 record
      dimension celprp(10),cell(9)

#ifdef VAMPIR
      call VTBEGIN(132, ierr)
#endif
      nhko=0
      mxstak=1
      kmaxa=0
      kmaxb=1
      kmaxc=1
      kmaxd=1
      kmaxe=1
      kmaxf=1
      rcut=0.0
      rvdw=0.0
      delr=0.0
      lhke=.false.
      lspme=.false.
c     
c     open the simulation input file
      
      if(idnode.eq.0)open (nread,file='CONTROL')
      
      call getrec(safe,idnode,mxnode,nread,record)
      if(.not.safe)go to 2000
      
      do nrecs= 1,mega
        
        call getrec(safe,idnode,mxnode,nread,record)
        if(.not.safe)go to 2000
        call lowcase(record,256)
        call strip(record,256)

        if(record(1:5).eq.'stack') then
          
          mxstak= intstr(record,256,idum)

        elseif(record(1:5).eq.'ewald'.or. 
     x      record(1:4).eq.'spme'.or.
     x      record(1:3).eq.'hke') then
c     
c     read Ewald or HK-Ewald or SPM-Ewald sum parameters
          
          lhke=(record(1:3).eq.'hke')
          lspme=(record(1:4).eq.'spme')
          lewald=(record(1:5).eq.'ewald')
          
          if(record(7:15).eq.'precision' .or.
     x      record(6:14).eq.'precision' .or.
     x      record(5:13).eq.'precision') then
            
            record(1:87)=record(14:100)
            eps = dblstr(record(1:1),87,idum)
            if(lhke) then

              ilen=87-idum 
              record(1:ilen)=record(idum:87)
              nhko=intstr(record(1:1),ilen,idum)
              jlen=ilen-idum
              record(1:jlen)=record(idum:ilen)
              nlatt=intstr(record(1:1),jlen,idum)
              nlatt=min(nlatt,2)
              
            endif
            
            if(rcut.lt.1.d-6)rcut=10.d0
c     
c     compute alpha and the kmax

            if(lewald.or.lspme)then
              
              call dcell(cell,celprp)
              eps = min(abs(eps),0.5d0)
              tol = sqrt(abs(log(eps*rcut)))
              alpha = sqrt(abs(log(eps*rcut*tol)))/rcut
              tol1 = sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
              fac = 1.d0
              if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) 
     x          fac = 2.d0**(1.d0/3.d0)
              kmax1 = nint(0.25d0 + fac*celprp(1)*alpha*tol1/pi)
              kmax2 = nint(0.25d0 + fac*celprp(2)*alpha*tol1/pi)
              kmax3 = nint(0.25d0 + fac*celprp(3)*alpha*tol1/pi)
              
            elseif(lhke)then
              
              if(nhko.eq.0)then
                if(eps.le.1.d-6)then
                  alpha=3.46d0/rcut
                elseif(eps.le.1.d-5)then
                  alpha=3.14d0/rcut
                else
                  alpha=2.76d0/rcut
                endif
              elseif(nhko.eq.1)then
                if(eps.le.1.d-6)then
                  alpha=4.37d0/rcut
                elseif(eps.le.1.d-5)then
                  alpha=4.08d0/rcut
                else
                  alpha=3.75d0/rcut
                endif                
              elseif(nhko.eq.2)then
                if(eps.le.1.d-6)then
                  alpha=5.01d0/rcut
                elseif(eps.le.1.d-5)then
                  alpha=4.74d0/rcut
                else
                  alpha=4.44d0/rcut
                endif
              elseif(nhko.eq.3)then
                if(eps.le.1.d-6)then
                  alpha=5.55d0/rcut
                elseif(eps.le.1.d-5)then
                  alpha=5.28d0/rcut
                else
                  alpha=5.00d0/rcut
                endif
              endif
              alpha=alpha/dble(2*nlatt+1)
              if(abs(cell(9)).lt.1.d-8)cell(9)=1.d0
              call dcell(cell,celprp)
              tol=2.d0*alpha*sqrt(abs(log(eps*alpha)))
              tol1=2.d0*alpha*sqrt(abs(log(eps*alpha*tol)))
              kmax1 = nint(0.25d0 + 0.5d0*celprp(1)*tol1/pi)
              kmax2 = nint(0.25d0 + 0.5d0*celprp(2)*tol1/pi)
              kmax3=1
              
            endif
            
          else
            
            alpha= dblstr(record,100,idum)
            ilen= 100-idum 
            record(1:ilen)= record(idum:100)
            
            kmax1=intstr(record,ilen,idum)
            ilen1= ilen-idum 
            record(1:ilen1)= record(idum:ilen)
            
            kmax2=intstr(record,ilen1,idum)
            ilen= ilen1-idum 
            record(1:ilen)= record(idum:ilen1)
            
            if(lhke)then

              kmax3=1
              nhko=intstr(record,ilen,idum)

            else

              kmax3=intstr(record,ilen,idum)

            endif
            
          endif
c     
c     for spme double kmax and set to next power of 2, with current
c     upper limit of 512

          if(lspme)then

            kmaxpow2 = 1
            do while (kmax1.gt.kmaxpow2.and.kmaxpow2.lt.256)
              kmaxpow2 = kmaxpow2 * 2
            end do
            kmaxd = 2 * kmaxpow2

            kmaxpow2 = 1
            do while (kmax2.gt.kmaxpow2.and.kmaxpow2.lt.256)
              kmaxpow2 = kmaxpow2 * 2
            end do
            kmaxe = 2 * kmaxpow2

            kmaxpow2 = 1
            do while (kmax3.gt.kmaxpow2.and.kmaxpow2.lt.256)
              kmaxpow2 = kmaxpow2 * 2
            end do
            kmaxf = 2 * kmaxpow2

          elseif(lhke) then

            kmaxa=kmax1
            kmaxb=kmax2
            kmaxc=1

          else

            kmaxa=kmax1
            kmaxb=kmax2
            kmaxc=kmax3

          endif

        elseif(record(1:3).eq.'cut') then

          rcut= dblstr(record,100,idum)

        elseif(record(1:4).eq.'rvdw') then

          rvdw= dblstr(record,100,idum)
          
        elseif(record(1:4).eq.'delr')then
          
          delr= dblstr(record,100,idum)
          
        elseif(record(1:6).eq.'finish')then

          go to 1000
          
        endif
        
      enddo
      
 1000 continue
      if(idnode.eq.0)close (nread)
      if(rvdw.eq.0.0)rvdw=rcut
#ifdef VAMPIR
      call VTEND(132, ierr)
#endif
      return

 2000 continue
      if(idnode.eq.0)close(nread)
      call error(idnode,17)
      return
      end
