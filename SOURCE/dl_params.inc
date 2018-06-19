! 09 FEB 06 - IUCHI - COMMENT OUT ndipole=901
! 01 DEC 05 - IUCHI - noutopt=913
! 17 NOV 05 - IUCHI - ndip3=912
! 16 NOV 05 - IUCHI - DIPOLE MOMENTS FOR TTM MODEL: ndip1, 2=910,911
! 10 NOV 05 - IUCHI - DECOMPOSED FORCES FOR TTM MODEL: nforce4,5=908,909
! 09 NOV 05 - IUCHI - DECOMPOSED FORCES FOR TTM MODEL: nforce1,2,3=905,906,907
! 01 NOV 05 - IUCHI - INTRAPOTENTIAL FOR TTM MODEL: nintra=904
!
!*********************************************************************
!
!     dl_poly insert file specifying fundamental parameters for the
!     entire package
!
!     copyright - daresbury laboratory 1997
!     authors - w. smith
!
!
!     note the following internal units apply everywhere
!
!     unit of time      (to)    =          1 x 10**(-12) seconds
!     unit of length    (lo)    =          1 x 10**(-10) metres
!     unit of mass      (mo)    = 1.6605402  x 10**(-27) kilograms
!     unit of charge    (qo)    = 1.60217733 x 10**(-19) coulombs
!     unit of energy    (eo)    = 1.6605402  x 10**(-23) joules
!     unit of pressure  (po)    = 1.6605402  x 10**(  7) pascals
!                               = 163.882576             atmospheres
!
!     wl
!     2001/06/12 12:41:05
!     1.11
!     Exp
!
!     Last updated: 1 Dec 2005 by S. iuchi
!
!*********************************************************************
!
      implicit real*8 (a-h,o-z)

!     array allocation parameters (set by subroutine parset)

      common/params1/kmaxa,kmaxb,kmaxc,minnode,msatms,msbad,msgrp
      common/params2/mspmf,msteth,mxangl,mxatms,mxbond,mxbuff,mxcell
      common/params3/mxcons,mxdihd,mxewld,mxexcl,mxfbp,mxfld,mxgatm
      common/params4/mxgrid,mxgrp,mxinv,mxlist,mxlshp,mxneut,mxngp
      common/params5/mxnstk,mxpang,mxpbnd,mxpdih,mxpfbp,mxpinv,mxpmf
      common/params6/mxproc,mxptbp,mxpvdw,mxrdf,mxshl,mxsite,mxspmf
      common/params7/mxstak,mxsvdw,mxtang,mxtbnd,mxtbp,mxtcon,mxtdih
      common/params8/mxteth,mxtinv,mxtmls,mxtshl,mxungp,mxvdw,mxxdf
      common/params9/mx2tbp,mx3fbp,mxebuf,mxquat,mxshak,mxspl,kmaxd
      common/params10/kmaxe,kmaxf,mxspme,mxftab,mxhko,mxegrd,mxhke

!
!     standard pi values
 
      parameter (pi =3.141592653589793d0,sqrpi =1.7724538509055159d0)
!
!     conversion factor for coulombic terms in internal units
!     i.e. (unit(charge)**2/(4 pi eps0 unit(length))/unit(energy)
 
      parameter (r4pie0 =138935.4835d0)
!
!     boltzmann constant in internal units
 
      parameter (boltz =8.31451115d-1)
!
!     planck's constant in internal units

      parameter (hbar=6.350780668d0)
!
!     conversion factor for pressure from internal units to katms

      parameter (prsunt = 1.63882576d-1)
!
!     main input channel
 
      parameter (nread = 7)
!
!     main output channel
 
      parameter (nrite = 8)
!
!     force field input channel
 
      parameter (nfield = 9)
!
!     configuration file input channel
 
      parameter (nconf = 10)
!
!     statistical data file output channel
 
      parameter (nstats = 20)
!
!     trajectory history file channel
 
      parameter (nhist = 21)
!
!     acummulators restart dump file
 
      parameter (nrest = 22)
!
!     tabulated potential file channel
 
      parameter (ntable = 23)
!
!     rdf file channel number
 
      parameter (nrdfdt = 24)
!
!     z density file channel number
 
      parameter (nzdndt = 25)
!
!     data dumping interval in event of system crash
 
      parameter (ndump = 10000)
!
!     data for dipole moments
!      
!      parameter (ndipole=901)
!
!     data for COM

      parameter (ncomwat=902)
!
!     data for MOVIE

      parameter (nmovi=903)
!
!     data for intramolecular potential

      parameter (nintra=904)
!
!     data for ttm forces (CC part)

      parameter (nforce1=905)
!
!     data for ttm forces (DC and DD parts)

      parameter (nforce2=906)
!
!     data for ttm forces (CC+DC+DD+INTRA for check)

      parameter (nforce3=907)
!
!     data for ttm forces (INTRA part)

      parameter (nforce4=908)
!
!     data for ttm forces (LJ part)

      parameter (nforce5=909)
!
!     data for ttm induced dipole

      parameter (ndip1=910)
!
!     data for ttm molecular dipole

      parameter (ndip2=911)
!
!     data for ttm molecular dipole

      parameter (ndip3=912)
!
!     data for steepest-descent like optimization

      parameter (noutopt=913)




