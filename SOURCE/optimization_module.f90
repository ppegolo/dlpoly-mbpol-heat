!
!  ==  optimization subroutine ==
!

module optimization

  implicit none 

  integer,private,save        ::  maxstep
  real(kind=8),private,save   ::  maxforce,rmsforce,maxdis,rmsdis

!l-bfgs parameters 

      integer,private,save                   :: totvar
      integer,private,save                   :: ndim, nwork,msave
!      real(kind=8),allocatable :: x(:),g(:),diag(:),w(:)
      real(kind=8),allocatable,private,save  :: diag(:),w(:)
      real(kind=8),private,save              :: eps,xtol
      integer,private,save                   :: iprint(2),iflag,n,m,j
      logical,private,save                   :: diagco

      integer,private,save                   :: output_index,output_rank,output_freq



 ! subroutines

  public :: opt_init,lbfgs_init
  public :: opt_force_energy
  public :: opt_finish

! external lb2
! common /lb3/mp,lp,gtol,stpmin,stpmax

 ! values needed for linesearch. Default values will work for almost all cases. 
 !  uncomment above line in case you want to change alpha in linesearch

  contains

 
  subroutine opt_init(natms_r,nrite,idnode,mxnode,nread)
 
    implicit none 

    integer,intent(in) :: natms_r,nrite

    integer :: mega=10000,nrecs

    logical :: safe 

    integer :: idnode,mxnode,nread,i

    integer,parameter :: nofields=100
    integer           :: countfields,fieldlength
    integer           :: starts(nofields),ends(nofields)
    character(len=256):: record
    character(len=256):: tmprecord


    output_index=nrite
    output_rank=idnode
  
    !default values 

    maxstep=10

    maxforce=1.500d2         ! in internal units  equivalent to .000302228 au
    rmsforce=1.000d2
! right now, distance are not checked for convergence -- if needed, modify it later 
    maxdis=0.001d0
    rmsdis=0.002d0 
    output_freq=1 

    totvar=3*natms_r

    if(idnode.eq.0) open(nread,file='CONTROL',status='old')

    if(idnode.eq.0) rewind(nread) 

     do nrecs=1,mega

        call getrec(safe,idnode,mxnode,nread,record)

        call lowcase(record,40)

        call cal_field(record,nofields,countfields,starts,ends)

        if(countfields==0)  cycle           ! takes care of blank lines

        if(record(1:1).eq.'#') then

         elseif(record(starts(1):ends(1)).eq.'geomopt') then

            tmprecord=record(starts(2):ends(2))
            if(starts(2)/=0) read(tmprecord,*) maxstep
        
            tmprecord=record(starts(3):ends(3))
            if(starts(3)/=0)  read(tmprecord,*) output_freq

            tmprecord=record(starts(4):ends(4))
            if(starts(4)/=0)  read(tmprecord,*) maxforce
        
!            tmprecord=record(starts(4):ends(4))
!            if(starts(4)/=0)  read(tmprecord,*) rmsforce

         elseif(record(starts(1):ends(1)).eq.'finish')then

            goto 2000

       endif

     enddo

     2000 continue 
    
     if(idnode.eq.0) then
       close(nread) 

       write(nrite,*)' ==============================================    '
       write(nrite,*)'    GEOMETRY OPTIMIZATION' 
       write(nrite,*)' ==============================================    '
       write(nrite,*)
       write(nrite,*)'CONVERGENCE CRITERIA SET  '
       write(nrite,*)
       write(nrite,*)'MAX STEPS      =   ',maxstep
       write(nrite,*)'MAX FORCE      =   ',maxforce,'j/mol/A'
!       write(nrite,*)'RMS FORCE      =   ',rmsforce,'j/mol/A'
!       write(nrite,*)'max distance   =   ',maxdis
!       write(nrite,*)'rms distance   =   ',rmsdis
       write(nrite,*)
       write(nrite,*)
       if(maxforce<100.d0) then
          write(nrite,*) 'Note: Maximum force cutoff is low --- '
       endif
       write(nrite,*)
     endif

   end subroutine opt_init


   subroutine lbfgs_init()

!
      m=20           ! no of l-bfgs vectors -- increase it in case of failure.. it should work in most cases
      iprint(1)= 1
      iprint(2)= 0
      diagco= .false.
      eps= maxforce
!      eps= 1.0d-5
      xtol= 1.0d-16
      iflag=0

      ndim=3*totvar
      msave=7
      nwork=ndim*(2*msave +1)+2*msave

      allocate(diag(ndim),w(nwork))

   end subroutine lbfgs_init 



!selects otpimization algorithm  and line search 


   subroutine opt_force_energy(atmnam,r,f,pot,opt_error)
      implicit none

      real(kind=8) :: r(*),f(*),pot
      character(len=8) :: atmnam(*)
      integer      :: opt_error


      call lbfgs(atmnam,output_freq,output_index,output_rank,maxstep,totvar,m,r,pot,f,diagco,diag,iprint,eps,xtol,w,iflag)

      opt_error=iflag
 
   end subroutine opt_force_energy


   subroutine opt_finish()
     implicit none 

     deallocate(diag,w)

   end subroutine opt_finish

end module 
