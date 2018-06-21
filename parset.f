! 28 SEP 06 - IUCHI - FROM MXGRID TO MXGRID+2 FOR MXBUFF
!
      subroutine parset(idnode,mxnode,buffer)

c     
c***********************************************************************
c     
c     dl_poly subroutine to determine required array sizes for 
c     allocation of memory manager
c
c     author - w.smith june 1997
c     copyright daresbury laboratory 1997
c     
c     wl
c     2001/06/12 12:55:54
c     1.11
c     Exp
c
!     Last modifed: 28 Sep 2006 by S. Iuchi
!     
c***********************************************************************
c     
      
#include "dl_params.inc"
      logical loglnk,lewald,lspme,lhke
      dimension cell(9),celprp(10),buffer(10)
#ifdef VAMPIR
      call VTBEGIN(1, ierr)
#endif
      lhke=.false.
      lspme=.false.
      lewald=.false.
c
c     specify maximum and minimum nodes

      mxproc =mxnode
      minnode=mxnode

c     
c     scan the FIELD file data
      
      call fldscan
     x   (idnode,mxnode,nfield,ntable,mxtbnd,mxtang,mxtdih,
     x   mxngp,mxgrp,mxatms,mxgatm,mxneut0,mxteth,mxtet1,mxsvdw,
     x   mxvdw,mxn1,mxtbp0,mxexcl,mxsite,mxbond,mxcons,mxangl,mxungp,
     x   mxdihd,mxtshl,mxshl,mxpmf,mxspmf,mxtinv,mxinv,mxfbp0,ngrid,
     x   mxtmls,mxtcon,rctbp,rcfbp)
c     
c     scan CONFIG file data
      
      call cfgscan
     x   (.true.,idnode,mxnode,nconf,imcon,volm,xhi,yhi,zhi,cell,buffer)
c
c     scan CONTROL file data

      call conscan
     x     (lewald,lspme,lhke,idnode,mxnode,nread,imcon,mxstak,kmaxa,
     x     kmaxb,kmaxc,kmaxd,kmaxe,kmaxf,nhko,rcut,rvdw,delr,cell)
c
c     set dimension of working coordinate arrays

      msatms=max(1,(mxatms+minnode-1)/minnode)
c
c     maximum number of neutral groups

      mxneut =max(mxneut0+1,1)
c
c     maximum number of molecule types

      mxtmls=max(mxtmls,1)
c
c     maximum number of specified bondlength constraints

      mxtcon=max(mxtcon,1)
c
c     maximum number of chemical bond potentials

      mxtbnd=max(mxtbnd,1)
c
c     maximum number of different bond angle potentials

      mxtang=max(mxtang,1)
c
c     maximum number of different torsional potentials

      mxtdih=max(mxtdih,1)
c
c     maximum number of different inversion potentials

      mxtinv=max(mxtinv,1)
c
c     maximum number of unique rigid body units

      mxungp=max(mxungp,1)
c
c     maximum number of tethered atom potentials

      mxteth=max(mxteth,1)
c
c     set maximum number of unique atom types

      mxsvdw=max(1,mxsvdw)
c
c     maximum number of different pair potentials

      mxvdw=max(mxvdw,1)+1
c
c     maximum number of three body potentials

      if(mxtbp0.eq.0)then

         mxtbp=1
         mx2tbp=1

      else

         mx2tbp=(mxsvdw*(mxsvdw+1))/2
         mxtbp=mx2tbp*mxsvdw

      endif
c
c     maximum number of four body potentials

      if(mxfbp0.eq.0)then

         mxfbp=1
         mx3fbp=1

      else

         mx3fbp=(mxsvdw*(mxsvdw+1)*(mxsvdw+2))/6
         mxfbp=mxsvdw*mx3fbp

      endif
c
c     maximum number of angular potential parameters

      mxpang = 4
c
c     maximum number of three body potential parameters

      mxptbp = mxpang+1
c
c     maximum number of four body potential parameters

      mxpfbp = 3
c
c     maximum number of parameters for dihedrals

      mxpdih = 5
c
c     maximum number of parameters for inversion potentials

      mxpinv = 2
c
c     maximum number of parameters for bond potentials

      mxpbnd = 4
c
c     maximum number of parameters for vdw potentials

      mxpvdw = 13 

c
c     maximum number of external field parameters

      mxfld = 10

c
c     maximum number of excluded atoms per atom

      mxexcl=max(mxexcl,1)
      mxexcl=4048 !VB for MOF
c
c     maximum number of different sites in system

      mxsite=max(mxsite,1)
c
c     maximum number of chemical bonds per node

      mxbond =max(1,(mxbond+minnode-1)/minnode)
c
c     maximum number of bond angles per node

      mxangl =max(1,(mxangl+minnode-1)/minnode)
c
c     maximum number of torsion angles per node

      mxdihd =max(1,(mxdihd+minnode-1)/minnode)
c
c     maximum number of inversion potentials per node

      mxinv =max(1,(mxinv+minnode-1)/minnode)
c
c     maximum number of constraints per node

      mxcons = max(1,2*((mxcons+minnode-1)/minnode))
c
c     maximum number of tethered atoms per node

      msteth = max(1,(mxtet1+minnode-1)/minnode)
c
c     maximum size for working arrays for bonds, angles, dihedrals
c     inversion potentials, tethers and core-shell units

      msbad = max(mxbond,mxangl,mxdihd,mxinv,msteth,mxshl)
c
c     maximum number of grid points in potentials arrays

      if(ngrid.eq.0)then

        ngrid = max(1000,int(rcut/0.0001d0+0.5d0)+4)

      endif
      mxgrid = ngrid
      mxegrd = 1
      if(lewald.or.lspme.or.lhke)mxegrd=ngrid
c
c     maximum dimension of rdf arrays

      mxrdf =max(128,int(rvdw/0.05d0+0.5d0))
c
c     maximum number of rigid groups in system

      mxgrp=max(mxgrp,1)
c
c     maximum number of rigid groups per node

      msgrp=max(1,(mxgrp+minnode-1)/minnode)
c
c     maximum number of sites per rigid unit

      mxngp=max(mxngp,3)
c
c     maximum number of sites in rigid units

      mxgatm = max(1,mxgatm)
c
c     maximum number of timesteps in stack arrays

      mxstak=max(100,mxstak)
c
c     maximum number of variables in stack arrays

      mxnstk =45+mxsvdw
c
c     dimension of shake shared atoms array

      mxlshp = max(mxcons*2,1)
c
c     set dimension of working arrays in ewald sum

      mxewld =1
      mxebuf =1
      if(lewald)then

        mxftab = 1
        mxewld = msatms
        mxebuf = (12*kmaxa+1)*(12*kmaxb+1)*(12*kmaxc+1)-1
        if(mxnode.le.16.and.mxebuf.le.5000)mxebuf=1

      endif
c
c     set dimension of working arrays in spme

      mxspme=1
      if(lspme)then

        mxspme=mxatms
        mxftab=2*(kmaxd+kmaxe+kmaxf)

      endif
c
c     set dimension of working arrays for HK ewald
      
      mxhko=1
      mxhke=1
      if(lhke)then

        mxhko=2
        mxewld=msatms
        mxhke=msatms
        if(nhko.gt.0)mxhko=max(2,nhko)
        mxebuf = (2*kmaxa+1)*(2*kmaxb+1)-1
        if(mxnode.le.16.and.mxebuf.le.5000)mxebuf=1

      endif
c
c     maximum dimension of principal transfer buffer

      mxbuff=max(6*mxatms,8*(mxcons+1),8*(mxgrp+1),mxnstk*mxstak,
     x     mxebuf,mxgrid+2,2*kmaxa*kmaxb*kmaxc,2*kmaxd*kmaxe*kmaxf,
     x     10000)
!      mxbuff=max(6*mxatms,8*(mxcons+1),8*(mxgrp+1),mxnstk*mxstak,
!     x     mxebuf,mxgrid,2*kmaxa*kmaxb*kmaxc,2*kmaxd*kmaxe*kmaxf,
!     x     10000)
c
c     maximum size of verlet neighbour/link cell list for each atom

c     decide if link-cells in use or not

      loglnk=.true.
      cut=rcut+delr
      dens=dble(mxatms)/volm
      ratio=1.5d0*dens*(4.d0*pi/3.d0)*cut**3
      mxlist=min(nint(ratio),(mxatms+1)/2)
      if(imcon.eq.0) then
        
        cell(1) = max(xhi+2.d0*cut,3.d0*cut)
        cell(5) = max(yhi+2.d0*cut,3.d0*cut)
        cell(9) = max(zhi+2.d0*cut,3.d0*cut)
        
      endif
      if(imcon.eq.6)then

        cell(9) = max(zhi+2.d0*cut,3.d0*cut,cell(9))

      endif

      call dcell(cell,celprp)
      ilx = max(3,int(celprp(7)/cut))
      ily = max(3,int(celprp(8)/cut))
      ilz = max(3,int(celprp(9)/cut))
      ncells = ilx*ily*ilz

      if(ncells.eq.27) loglnk=.false.
      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) loglnk = .false.
      if(mxneut0.gt.0.and.ncells.le.36) loglnk=.false.
      
      mxcell=1
      if(loglnk)then

        mxlist=64*nint(1.5d0*dens*celprp(10)/dble(ncells))
        mxcell=(ilx+2)*(ily+2)*(ilz+2)

      endif

      if(mxneut0.gt.0) mxlist = (mxneut+1)/2
      mxlist=2*mxlist
      mxlist=max(256,mxlist) !VB: for trimer in huge box
      if(mxtbp0.gt.0.or.mxfbp0.gt.0)then

        if(mxtbp0.gt.0)cut=min(cut,rctbp)
        if(mxfbp0.gt.0)cut=min(cut,rcfbp)
        ilx = max(3,int(celprp(7)/cut))
        ily = max(3,int(celprp(8)/cut))
        ilz = max(3,int(celprp(9)/cut))
        mxcell=max(mxcell,(ilx+2)*(ily+2)*(ilz+2))

      endif
c
c     maximum size for coordinate difference arrays

      mxxdf=max(mxlist,mxatms,mxcons,mxn1*mxn1*(mxneut+1)/2)
c
c     maximum number of core-shell unit types

      mxtshl=max(mxtshl,1)
c
c     maximum number of core-shell units

      mxshl=max(mxshl,1)
c
c     potential of mean force array parameter

      mxpmf=max(mxpmf,1)
c
c     number of pmf constraints on a processor

      mspmf=max(1,(mxpmf+minnode-1)/minnode)
c
c     maximum number of sites to define pmf units

      mxspmf=max(mxspmf,1)
c
c     maximum iterations in quaternion integration
 
      mxquat = 100
c
c     maximum number of shake cycles
 
      mxshak = 100
c
c     maximum b-spline interpolation order

      mxspl = 12
#ifdef VAMPIR
      call VTEND(1, ierr)
#endif
      return
      end
