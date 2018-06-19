#ifndef SERIAL
c*********************************************************************
c
c     dl_poly include file for MPI, PVM and SHMEM
c
c     trf
c     1996/05/01 12:03:26
c     1.2
c     Exp
c
c*********************************************************************
#ifdef MPI
#ifdef MPIU       
#define mpipriv mpipriv_
#endif
 
#include "mpif.h"
ckt IF CRAYT3D AND MPI DEFINED, DO NOT DECLARE Dlen HERE
#ifdef CRAYT3D

#else
      integer Dlen
      parameter (Dlen = 8)
#endif
#endif

#ifdef INTEL

      integer Dlen, Ilen
      parameter (Dlen = 8, Ilen = 4)

#endif
#if PVM

      include 'fpvm3.h'

#endif
#if SHMEM || SGISHMEM

      include 'mpif.h'

#endif
#ifdef CRAYT3D
c     data type and byte length specification

      integer Ilen,Dlen,Llen
      parameter (Ilen = 8, Dlen = 8, Llen = 8)

#endif
#ifdef PVM

      integer Itype,Dtype,Ilen,Dlen,Llen
      parameter (Itype = INTEGER4, Dtype = REAL8)
      parameter (Ilen = 4, Dlen = 8, Llen = 4)

c     global goup -- should be commented out for Cray T3D

      character*8 pvmall
      data pvmall/'global'/

c     PVM specifications and common blocks

c     number of nodes

      integer Nodes
      parameter (Nodes = 256)

c     parameters for Communication routines

c     Nbranch is the number of branches at each node in the gsum 
c     routines
c     Encoding is the default encoding used for PVM messages

      integer Nbranch, Encoding
      parameter( Nbranch=2, Encoding=PVMRAW )

c     Main PVM common block

c     tidtable is a list of all the tids
c     my_tid is my tid
c     root_tid is the tid of the root processor

      integer tidtable, my_tid, root_tid
      common/COMMSTIDS/tidtable(0:Nodes-1), my_tid, root_tid

      integer no_of_nodes, hcdim, me
      common /COMMSNODES/ no_of_nodes, hcdim, me

#endif
c     parameters for message tags

      integer Xgcopy_tag, Merge1_tag
      integer Dgsum_tag1, Dgsum_tag2
      integer Igmax_tag1, Igmax_tag2
      integer Igsum_tag1, Igsum_tag2
      integer Shmove_tag, Merge_tag
      integer Merge4_tag, Shell_tag
      integer Gstate_tag1, Gstate_tag2
      integer Passcon_tag


c     MPI tagsizes

      parameter(
     $     Xgcopy_tag = 105,
     $     Dgsum_tag1  = 6001,
     $     Dgsum_tag2  = 9001,
     $     Igsum_tag1  = 12001,
     $     Igsum_tag2  = 15001,
     $     Shmove_tag  = 18001,
     $     Merge_tag   = 6001,
     $     Merge1_tag  = 9001,
     $     Merge4_tag  = 12001,
     $     Gstate_tag1 = 15001,
     $     Gstate_tag2 = 18001,
     $     Shell_tag   = 21001,
     $     Passcon_tag =  3001,
     $     Igmax_tag1  = 24001,
     $     Igmax_tag2  = 27001)


#endif

