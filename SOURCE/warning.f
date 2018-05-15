      subroutine warning(idnode,kode,a,b,c)
c     
c***********************************************************************
c     
c     dl_poly subroutine for printing warning messages and returning
c     control back to the main program
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    april 1994
c     
c     wl
c     2001/06/12 13:05:16
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      if(idnode.eq.0) then

        if (kode.eq. 10) then
          write(nrite,'(/,1x,a)')
     x      ' *** warning - no pair forces in use ***'
          
        elseif (kode.eq. 20) then
          
          ia = nint(a)
          ib = nint(b)
          ic = nint(c)
          write(nrite,'(/,1x,a50,i5,5x,a6,2i10)')
     x      ' *** warning - : 1..4 scale factors reset for molecule ',
     x      ia,'sites ',ib,ic
          
        elseif (kode.eq. 30) then
          write(nrite,'(/,1x,a)')
     x      ' *** warning - atomistic cutoff with electrostatics ***'

        elseif(kode.eq. 40) then
          write(nrite,'(/,1x,a,/,1x,a,f12.6)')
     x      ' *** warning - radial cutoff reset ***',
     x      'new potential cutoff radius    ',a

        elseif(kode.eq. 50) then
          write(nrite,'(/,1x,a,/,1x,a,f12.6)')
     x      ' *** warning - short range cutoff reset ***',
     x      'new cutoff radius (rvdw)       ',a

        elseif(kode.eq. 60) then
          write(nrite,'(/,1x,a,f12.6,a)')
     x      ' *** warning - total system charge:',a,' ***'

        elseif(kode.eq. 70) then
          write(nrite,'(/,1x,a,f12.6,a)')
     x      ' *** warning - switching length reset to: ',a,' ***'

        elseif(kode.eq. 80) then
          write(nrite,'(/,1x,a,f12.6,a)')
     x      ' *** warning -  requested thermostat unavailable ***'

        elseif(kode.eq. 90) then

          ia=nint(a)
          ib=nint(b)
          write(nrite,'(/,1x,a)')
     x      ' *** warning -  cannot activate link cell option ***'
          write(nrite,'(/,1x,a,i6,a,i6)')
     x      ' *** you must change parameter mxcell from ',ib,' to ',ia
          write(nrite,'(/,1x,a)')
     x      ' *** using standard Brode-Ahlrichs list builder  ***'

        elseif(kode.eq.100) then

          write(nrite,'(/,1x,a,a,1p,e12.4,a)')
     x      ' *** warning - HK real space screening function problem',
     x      ' at cut off:',a,' ***'

        elseif(kode.eq.105) then

          write(nrite,'(/,1x,a,a,1p,e12.4,a)')
     x      ' *** warning - HK recip space screening function problem',
     x      ' at cut off:',a,' ***'

        elseif(kode.eq.110) then

          write(nrite,'(/,1x,a)')
     x    ' *** warning - undefined atom-atom interactions set to zero'
     x    //' ***'

        else

          write(nrite,'(/,1x,a)')
     x      ' *** unspecified warning encountered ***'
        endif
        
      endif
      
      return
      end
