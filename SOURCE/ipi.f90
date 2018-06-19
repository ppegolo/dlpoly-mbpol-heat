
      MODULE ipi
        USE f90sockets, ONLY : open_socket, writebuffer, readbuffer

      IMPLICIT NONE 

      
! I-PI INTERFACE
      integer, parameter :: dp = selected_real_kind(15, 307)
      INTEGER, PARAMETER :: MSGLEN=12
      INTEGER        :: nwat 
      LOGICAL :: isinit=.false., imginit, hasdata=.false., ionode, imglock
      CHARACTER*12 :: header
      CHARACTER*1024 :: parbuffer
      INTEGER socket, nat, natwrap, ccmd, nimg, parbufflen, oimg, cbuf,i
      CHARACTER*1024 :: host  
      CHARACTER*300000  :: initbuffer              ! length can store around 3500 water dipoles
      REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot, mtxbuf(9)
      REAL*8, ALLOCATABLE :: combuf(:)
      REAL(KIND=dp) :: sleeptime


      !added by VK
      INTERFACE
      SUBROUTINE uwait(sec) BIND(C, NAME="uwait")
         USE ISO_C_BINDING, ONLY: C_DOUBLE
      REAL(C_DOUBLE)                                     :: sec

      END SUBROUTINE
      END INTERFACE
      CONTAINS
 
      SUBROUTINE IPI_INIT(serveraddr, idnode)
       CHARACTER*1024  serveraddr  
       INTEGER inet, port, idnode
       
       if (.not.(isinit)) then
         isinit = .true.
         ionode=.false.
         oimg = -1
         imginit = .false.
         imglock = .false.
         sleeptime=0.001d0
         if (idnode.eq.0) ionode=.true.
         inet=1
         host=serveraddr(1:INDEX(serveraddr,":",.true.)-1)//achar(0)
         read(serveraddr(INDEX(serveraddr,':',.true.)+1:),*) port   
       
         IF (adjustl(trim(serveraddr(1:INDEX(serveraddr,':')-1))).eq."unix") &
          inet=0
!         if (ionode) write(*,*) adjustl(trim(serveraddr(1:INDEX(serveraddr,':')-1))),  &
!           adjustl(trim(serveraddr(1:INDEX(serveraddr,':')-1))).eq."unix"
         host=serveraddr(INDEX(serveraddr,':')+1: &
         INDEX(serveraddr,':',.true.)-1)//achar(0)    
!         IF (ionode) write(*,*)  &
!       " @ DRIVER MODE: Connecting to host:port:inet ",  &
!       trim(host), port, inet
       
         if (ionode) call open_socket(socket, inet, port, host)   
       END IF

      
      END SUBROUTINE
      
      SUBROUTINE IPI_STEP(volm, pot, stress, cell, com3n, dipole)
          IMPLICIT NONE
          REAL*8 volm, pot, stress(:), cell(:), ipipot, celldata(10)
          REAL*8, ALLOCATABLE :: com3n(:,:), dipole(:,:)
          INTEGER ierror, nprocs, mpiindex, fwait, iwait, mpistatus(128) !oversized MPISTATUS
          INTEGER idnode, wait_msg, wait_req

          REAL,PARAMETER :: conv_length_a_au=0.52917721067d0
          REAL,PARAMETER :: conv_energy_dlpoly_au=3.8087989d-06
          REAL,PARAMETER :: cons_ratio=418.6799938/418.4d0
!SR: this ratio (418.67999/418.4) is used because (1 kcal/mol)/(1 kjoule/mol) is 
! not exactly 4.184 which is used in dl_poly 
! to be consistent, taken this ratio from ipi so that output energies are same across dlpoly and ipi 

#ifdef MPI
          include "mpif.h"
          call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
          call MPI_COMM_RANK(MPI_COMM_WORLD, idnode, ierror)
#endif
      do       
#ifdef MPI       
          IF (ionode) THEN
            call readbuffer(socket, header, MSGLEN)
            wait_msg = 0
            DO iwait = 0, nprocs-1
               IF (iwait /= 0) THEN
                  call MPI_Send(wait_msg, 1, MPI_INT, iwait, 666, MPI_COMM_WORLD, ierror)
               ENDIF
            ENDDO
          ELSE
            !figure out array of requests
            !<pre>INCLUDE ’mpif.h’ MPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST,ERROR
            CALL MPI_Irecv(wait_msg, 1, MPI_INTEGER, 0, 666, MPI_COMM_WORLD, wait_req, ierror)
            CALL MPI_Test(wait_req, fwait, mpistatus, ierror)
            DO WHILE (fwait == 0)
               CALL MPI_Test(wait_req, fwait, mpistatus, ierror)
               CALL uwait(sleeptime)
            ENDDO
          ENDIF
          call MPI_Bcast(header,MSGLEN,MPI_CHAR,0,MPI_COMM_WORLD,ierror) !broadcast MSGLEN chars to all nodes
#else
         if (ionode) call readbuffer(socket, header, MSGLEN)
#endif
!          if (ionode) write(*,*) "Read header ", header
          if (trim(header) == "STATUS") then
            if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
                if (hasdata) then
                   call writebuffer(socket,"HAVEDATA    ",MSGLEN)
                else if (imginit) then
                   call writebuffer(socket,"READY       ",MSGLEN)
                else 
                   call writebuffer(socket,"NEEDINIT    ",MSGLEN)             
                endif
            endif
         elseif (trim(header) == "INIT") then
            if (ionode) call readbuffer(socket, nimg) ! actually this reads the replica id   
            if (ionode) call readbuffer(socket, parbufflen) ! length of parameter string -- ignored at present!
            if (ionode) call readbuffer(socket, parbuffer, parbufflen)   
#ifdef MPI            
            call MPI_Bcast(nimg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
            
            imginit = .true.
            if (oimg == nimg) then 
               imglock = .true.
            else
               imglock = .false.
            endif
            oimg = nimg
         elseif (trim(header)=="GETFORCE") then
         if (hasdata) then 
            combuf=RESHAPE(com3n,(/3*natwrap/))*conv_energy_dlpoly_au*conv_length_a_au*cons_ratio
            vir=transpose(reshape(stress,(/3,3/)))*conv_energy_dlpoly_au*cons_ratio

            ipipot=pot*conv_energy_dlpoly_au*cons_ratio
!         if (ionode) write(*,*) 
!     x     " @ DRIVER MODE: Returning v,forces,stress "
            if (ionode) then      
               call writebuffer(socket,"FORCEREADY  ",MSGLEN)            
               call writebuffer(socket,ipipot)
               call writebuffer(socket,natwrap)            
               call writebuffer(socket,combuf,3*natwrap)
               call writebuffer(socket,reshape(vir,(/9/)),9)
               initbuffer = " "
               WRITE(initbuffer,*) (dipole(i,:),i=1,nwat)
               cbuf = LEN_TRIM(initbuffer)
               CALL writebuffer(socket,cbuf) ! Writes back the water dipole all individual molecules
               CALL writebuffer(socket,initbuffer,cbuf)
            endif
            hasdata=.false.
            imginit = .false.
         endif
      elseif  (trim(header)=="POSDATA") then
         if (ionode) then        
            call readbuffer(socket, mtxbuf, 9)
            cellh=RESHAPE(mtxbuf,(/3,3/))
            call readbuffer(socket, mtxbuf, 9)
            cellih=RESHAPE(mtxbuf,(/3,3/))
            call readbuffer(socket, natwrap)
!            cellh=transpose(cellh)
!            cellih=transpose(cellih)
         endif
#ifdef MPI
         call MPI_Bcast(cellh,9,MPI_DOUBLE,0,MPI_COMM_WORLD,ierror)
         call MPI_Bcast(cellih,9,MPI_DOUBLE,0,MPI_COMM_WORLD,ierror)
         call MPI_Bcast(natwrap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
         nat = natwrap
         nwat = nat/3
         if (.not.allocated(combuf)) then
            allocate(combuf(3*natwrap))
            allocate(com3n(3,natwrap))
         endif
         if (.not.allocated(dipole)) then
            allocate(dipole(nwat,3))
         endif
    
         if (ionode) call readbuffer(socket, combuf, natwrap*3)
#ifdef MPI
         call MPI_Bcast(combuf,3*natwrap,MPI_DOUBLE,0,MPI_COMM_WORLD,ierror)
#endif
         cell=reshape(transpose(cellh), (/9/))*conv_length_a_au !cell to Angstrom
         call dcell(cell,celldata)
         volm = celldata(10)
           
         com3n = RESHAPE(combuf, (/ 3 , natwrap /) )*conv_length_a_au ! positions to Angstrom         
         exit
      else
         write(*,*) "Unexpected message from wrapper", trim(header)
         stop "ended"

      endif   
      enddo
      END SUBROUTINE

      END MODULE ipi
