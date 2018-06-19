
  MODULE Langevin_Thermostat

    IMPLICIT NONE

    SAVE 

    ! taken from  centroid_module.f90 
    real(8), parameter :: pi  = 3.141592653589793d0
    real(8), parameter :: hbar = 6.350780668d0
    real(8), parameter :: kB = 8.31451115d-1

    INTEGER :: thermoType,numS_gle,TAYLORN,MATSQP 
    INTEGER :: thermo_centroid,npslices,rank,nproc
    real*8 :: dt,wn,omega0,beta,c1,c2,wk,coswk,wksinwk,sinoverwk
    real*8 :: c1_cell,c2_cell, cell_mass, press_ext, tau_cell, omega_cell
    real*8,ALLOCATABLE,DIMENSION(:) :: one_over_sqrt_nmmass
    real*8,ALLOCATABLE,DIMENSION(:) :: one_over_nmmass
    real*8,ALLOCATABLE,DIMENSION(:) :: Sx,Sy,Sz,vcm_gle,vcm_gle_ 
    real*8,ALLOCATABLE,DIMENSION(:) :: Rx,Ry,Rz
    real*8,ALLOCATABLE,DIMENSION(:,:) :: Amat_gle,Cmat_gle 
    real*8,ALLOCATABLE,DIMENSION(:,:) :: Tmat_gle,Smat_gle 
    CHARACTER (LEN=256) AMAT_FILE 
    CHARACTER (LEN=256) CMAT_FILE 

    ! scaling constants for normal mode transformations, masses, freq, etc
    ! the scaling constants take the values from nose-hoover to langevin 
    !logical :: langevin,evolve_in_nm,add_spring_to_forces 
    !real*8 :: freq_scale, force_scale 
    real(8) :: freq_scale = 1.0d0
    real(8) :: force_scale = 1.0d0
    logical :: langevin = .false.
    logical :: evolve_in_nm = .true.
    logical :: add_spring_to_forces = .true.

    PUBLIC :: langevin, evolve_in_nm, add_spring_to_forces,  &
              freq_scale, force_scale,  &
              Langevin_init, Langevin_thermostat_integrate,  &
              Langevin_evolve_first,Langevin_evolve_second,  &
              Langevin_Barostat_integrate

  CONTAINS

    SUBROUTINE Langevin_init(lhead,pimd_head,ibead,nproc_,rank_,idnode,dt_,np_,omega0_,natoms,  &
                        thermoType_,numS_,mass,file_out,temperature,Aunit,Cunit, &
                        Afile,Cfile,press,tau_cell_, omega_cell_ , cell_mass_)

      use multibead, only:comm_bead,comm_mb

      IMPLICIT NONE

#ifdef MPI
#  include "mpif.h"
#endif

      integer, intent(in) :: thermoType_,np_,numS_,file_out,idnode,natoms,ibead
      integer, intent(in) :: nproc_,rank_ 
      real*8, intent(in) :: dt_,omega0_,mass(:),temperature 
      real*8, intent(in) :: press, tau_cell_, omega_cell_
      real*8, intent(out) :: cell_mass_ 
      logical, intent(in) :: lhead,pimd_head 

      INTEGER N,NDIM,IERR,i,j,k,l,readCAR,p
      INTEGER np,ni,ns,k_seed,glemat_comm
      CHARACTER*256 Aunit,Cunit
      CHARACTER*256 Afile,Cfile 
      LOGICAL readCp,ife,calculate_matrix
      real*8 vtmp(3),vnew(3),facA,facC,dummy
      !real*8 vtmp(3),vnew(3),facA,facC,tp,fct,azero,aone,aminus,dummy
      !real*8,ALLOCATABLE,DIMENSION(:,:) :: AI,BI
      integer, allocatable :: seeds(:)

      if(thermoType_ < 1 .or. thermoType_ > 6 ) then
        langevin = .false.
        return 
      endif
      
      ! defaults
      CMAT_FILE = 'GLE_CMAT'
      AMAT_FILE = 'GLE_AMAT'
      !CMAT_FILE = Afile
      !AMAT_FILE = Cfile 
      numS_gle = 0 
      TAYLORN=10
      MATSQP=10
      readCp=.false.
      readCAR=0
      omega0 = 100.0d0
      thermo_centroid = 1
      langevin = .true.
      evolve_in_nm = .false.

      rank = rank_
      nproc = nproc_
      dt = dt_ 
      npslices = np_
      omega0 = omega0_
      thermoType = thermoType_
      numS_gle = numS_
      tau_cell = tau_cell_
      omega_cell = omega_cell_
      press_ext = press

      ! by default using the harmonic propagator
      add_spring_to_forces = .false.
      
      allocate(seeds(2*nproc))
      !call generate_seed(rank,nproc,seeds)
      call generate_seed(seeds)
      dummy = random(seeds(rank+1),seeds(rank+1+nproc))
      deallocate(seeds)


      ! scaling functions 
      force_scale = dble(npslices) 
      freq_scale = dsqrt(force_scale)  

      !force_scale = 1.d0
      !freq_scale = 1.d0 

      ! these are the GLE-based algorithms 
      if(thermoType .eq.2 .OR.   & 
         thermoType .eq. 3 .OR.  & 
         thermoType .eq. 5 .OR.  &
         thermoType .eq. 6) then
        ! needed in 1/ps
        if( Aunit(1:2) .eq. 'ps' ) then
          facA=1.0d0    
        elseif( Aunit(1:2) .eq. 'fs' ) then
          facA=1000.0d0    
        elseif( Aunit(1:1) .eq. 's' ) then
          facA=1.0d-12    
        elseif( Aunit(1:2) .eq. 'au' ) then
          facA=1.0d0/(2.41888432650500e-05)    
        else
          call error(idnode,709)
        endif
        if(thermoType .gt. 2) then
          ! needed in 10 J/mol
          if( Cunit(1:1) .eq. 'k' ) then
           facC=0.831446222087654d0
          elseif( Cunit(1:2) .eq. 'ev' ) then
            facC=9648.53363969351d0
          elseif( Cunit(1:1) .eq. 'au' ) then
            facC=262549.964037578d0
          else
            call error(idnode,710)
          endif
        endif
      endif 

      if(pimd_head) then
        open(file_out,file='OUTPUT_PIMD',position='append')
      endif

      if(thermoType .eq. 1) then
        if (pimd_head) then
          write(file_out,*)'Using Langevin thermostat on the normal modes.'
          write(file_out,*)'Thermostat natural frequency, w0 = ',omega0
        endif
      else if(thermoType .eq. 2) then

        IF(numS_gle .lt. 1) then
          CALL error(idnode,704) 
        endif
        if(pimd_head) then
          write(file_out,*)'Using GLE method to thermostat the PI calculation.'  
          write(file_out,*)'The fluctuation-dissipation theorem will be obeyed,'
          write(file_out,*)'so a canonical distribution will be achieved in the'
          write(file_out,*)'associated classical problem.' 
          write(file_out,*)'GLE - Number of additional momenta: ',numS_gle
          write(file_out,*)'GLE - Assuming that the A matrix is in GLE_AMAT' 
          write(file_out,*)'      in units of ',TRIM(Aunit) 
          write(file_out,*)'facA: ',facA   
        endif
      else if(thermoType .eq. 3) then
        thermoType = 2
        readCp=.true.

        IF(numS_gle .lt. 1) then
          CALL error(idnode,704) 
        endif
        if(pimd_head) then
          write(file_out,*)'Using the GLE+PI method to accelerate PI convergence.'
          write(file_out,*)'GLE - Number of additional momenta: ',numS_gle
          write(file_out,*)'GLE - Assuming that the A matrix is in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Aunit) 
          write(file_out,*)'facA: ',facA   
          write(file_out,*)'GLE - Assuming that the Cp matrix is also in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Cunit) 
          write(file_out,*)'facC: ',facC 
        endif
      ! PILE_NPT in NM 
      elseif(thermoType.eq.4 ) then
        if (pimd_head) then
          write(file_out,*)'Using PILE_NPT thermostat.' 
          write(file_out,*)'A Langevin thermostat on the normal modes.'
          write(file_out,*)'Thermostat natural frequency, w0 = ',omega0
          write(file_out,*)'Barostat natural frequency: ',omega_cell
          write(file_out,*)'Cell relaxation time: ',tau_cell
        endif
      elseif(thermoType.eq.5 ) then
        readCp=.true.
        IF(numS_gle .lt. 1) then
          CALL error(idnode,704)
        endif
        if (pimd_head) then
          write(file_out,*)'Using PIGLET-NVT thermostat.'
          write(file_out,*)'An accelerated Langevin thermostat on the normal modes.'
          write(file_out,*)'GLE - Number of additional momenta: ',numS_gle
          write(file_out,*)'GLE - Assuming that the A matrix is in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Aunit)
          write(file_out,*)'facA: ',facA
          write(file_out,*)'GLE - Assuming that the Cp matrix is also in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Cunit)
          write(file_out,*)'facC: ',facC
        endif
      ! PIGLET_NPT in NM 
      elseif(thermoType.eq.6 ) then
        readCp=.true.
        IF(numS_gle .lt. 1) then
          CALL error(idnode,704)
        endif
        if (pimd_head) then
          write(file_out,*)'Using PIGLET_NPT thermostat.'
          write(file_out,*)'An accelerated Langevin thermostat on the normal modes.'
          write(file_out,*)'GLE - Number of additional momenta: ',numS_gle
          write(file_out,*)'GLE - Assuming that the A matrix is in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Aunit)
          write(file_out,*)'facA: ',facA
          write(file_out,*)'GLE - Assuming that the Cp matrix is also in GLE_AMAT'
          write(file_out,*)'      in units of ',TRIM(Cunit)
          write(file_out,*)'facC: ',facC
          write(file_out,*)'Thermostat natural frequency, w0 = ',omega0
          write(file_out,*)'Barostat natural frequency: ',omega_cell
          write(file_out,*)'Cell relaxation time: ',tau_cell
        endif
      else
        CALL error(idnode,705) 
      endif

      allocate(one_over_sqrt_nmmass(natoms),one_over_nmmass(natoms))

        ! assuming a simulation of the classical system at the physical
        ! temperature, for a RPMD-like mass-tensor, change these below
        beta = 1.0d0/(temperature*kB)  ! 1.0/kB*T in 10 J/mol 
        wn = sqrt( dble(npslices))/beta/hbar*freq_scale       ! n*kB*T/hbar in ps^-1 
        do i = 1,natoms
         one_over_nmmass(i) = 1.0d0/mass(i)
         one_over_sqrt_nmmass(i) = 1.0d0/dsqrt(mass(i))
        enddo
        

         ! coswk, wksinwk, and sinoverwk used if integrating harmonic hamiltonian
         ! nslices-->[1:NP], k-->[0:NP-1] 
        if(ibead .eq. 1) then
          c1 = dexp(-dt*0.5d0*omega0)
          c2 = dsqrt((1.0d0 - c1*c1)*force_scale/beta)
          wk = 0.0d0
          coswk = 1.0d0
          wksinwk = 0.0d0
          sinoverwk = dt
        else
          wk = 2.0d0*wn*dsin((ibead-1.0d0)*pi/npslices)
          c1 = dexp(-dt*wk)
          c2 = dsqrt((1.0d0 - c1*c1)*force_scale/beta)
          coswk = dcos(wk*dt)
          wksinwk = wk*dsin(wk*dt)
          sinoverwk = dsin(wk*dt)/wk
        endif
        
         if(thermoType.eq.4 .or. thermoType.eq.6) then
           evolve_in_nm=.true.
           cell_mass = 3.0d0*natoms*tau_cell*tau_cell/beta 
           cell_mass_ = cell_mass 
           c1_cell = dexp(-dt*0.5d0*omega_cell)
           c2_cell = dsqrt((1.0d0 - c1_cell*c1_cell)*force_scale*cell_mass/beta)
         endif

        !debug
        !write(100+ibead+1,*)'wn,beta,one_over_sqrt_nmmass(1),one_over_sqrt_nmmass(2),c1,c2:', &
        ! wn,beta,one_over_sqrt_nmmass(1),one_over_sqrt_nmmass(2),c1,c2
        !call flush(100+ibead+1)

        if(thermoType.eq.2 .OR.      &   ! GLE or PI+GLE
               thermoType.eq.5 .OR.      &   ! PIGLET-NVT 
               thermoType.eq.6 ) then    ! PIGLET-NPT 

          glemat_comm = comm_mb
          calculate_matrix = .false. 
          if(thermoType.eq.5 .or. thermoType.eq.6) then
            glemat_comm = comm_bead
            if(lhead) calculate_matrix = .true.
          else
            if(pimd_head) calculate_matrix = .true.
          endif

          ALLOCATE(Sx(numS_gle+1),Sy(numS_gle+1),Sz(numS_gle+1),vcm_gle(3*(numS_gle+1)+1),   &
                  vcm_gle_(3*(numS_gle+1)+1),Tmat_gle(numS_gle+1,numS_gle+1), & 
                  Smat_gle(numS_gle+1,numS_gle+1),Rx(numS_gle+1),Ry(numS_gle+1),Rz(numS_gle+1))

          if(calculate_matrix) then
         
            ALLOCATE(Amat_gle(numS_gle+1,numS_gle+1), &
             Cmat_gle(numS_gle+1,numS_gle+1))
           
            inquire(file=TRIM(AMAT_FILE),exist=ife)
            if(.not.ife) call error(idnode,706)
    
!          write(file_out,*) 'about to reading matrices'
!          call flush(file_out)
            OPEN(UNIT=551,FILE=TRIM(AMAT_FILE),STATUS='OLD')
            if(thermoType.eq.5 .or. thermoType.eq.6) then
              call read_GLE_MAT(551, ibead, npslices,readCp)
            else
              call read_GLE_MAT(551, 1, 1, readCp)
            endif
!          write(file_out,*) 'finished reading matrices'
!          call flush(file_out)
            Amat_gle = Amat_gle*facA*freq_scale/sqrt(dble(npslices))
            if(readCp) then
             Cmat_gle = Cmat_gle*facC*force_scale/dble(npslices)
            else ! Cp = nkT
             Cmat_gle = 0.d0 
             do i=1,numS_gle+1
              Cmat_gle(i,i) = force_scale/beta
             enddo
            endif

            CLOSE(551)

            call build_GLE_matrices(IERR)
          !write(*,*) 'finished with GLE matrices' 
 
            j=numS_gle+1 
            if(IERR .ne. 0) then 
              WRITE(file_out,*)'Error in Cholesky decomposition of gle matrix, info:',IERR
              WRITE(file_out,*) ' Amat:'
              do i=1,j
               do k=1,j
                write(file_out,77,ADVANCE="No") Amat_gle(i,k)
               enddo
               write(file_out,*) ''
              enddo
              WRITE(file_out,*) ' Cmat:'
              do i=1,j
               do k=1,j
                write(file_out,77,ADVANCE="No") Cmat_gle(i,k)
               enddo
               write(file_out,*) ''
              enddo
              WRITE(file_out,*) ' Tmat:'
              do i=1,j
               do k=1,j
                write(file_out,77,ADVANCE="No") Tmat_gle(i,k)
               enddo
               write(file_out,*) ''
              enddo
              WRITE(file_out,*) ' SSt:'
              do i=1,j
               do k=1,j
                write(file_out,77,ADVANCE="No") Smat_gle(i,k)
               enddo
               write(file_out,*) ''
              enddo
              call error(idnode,708)
            endif

            ! output matrices for safety
            WRITE(file_out,*) '************************************************ '
            WRITE(file_out,*) ' GLE matrices:'
            WRITE(file_out,*) ' Amat:'
            do i=1,j
             do k=1,j
              write(file_out,77,ADVANCE="No") Amat_gle(i,k)
             enddo
             write(file_out,*) '' 
            enddo
            WRITE(file_out,*) ' Cmat:'
            do i=1,j
             do k=1,j
              write(file_out,77,ADVANCE="No") Cmat_gle(i,k)
             enddo
             write(file_out,*) ''
            enddo
            WRITE(file_out,*) ' Tmat:'
            do i=1,j
             do k=1,j
              write(file_out,77,ADVANCE="No") Tmat_gle(i,k)
             enddo
             write(file_out,*) ''
            enddo
            WRITE(file_out,*) ' Smat:'
            do i=1,j
             do k=1,j
              write(file_out,77,ADVANCE="No") Smat_gle(i,k)
             enddo
             write(file_out,*) ''
            enddo
            WRITE(file_out,*) '************************************************ '

            j=(numS_gle+1)*(numS_gle+1)  
            call MPI_Bcast( Tmat_gle(1,1), j, MPI_DOUBLE_PRECISION, 0,   & 
              glemat_comm , ierr )
            call MPI_Bcast( Smat_gle(1,1), j, MPI_DOUBLE_PRECISION, 0,   & 
              glemat_comm , ierr )
     
          else

            j=(numS_gle+1)*(numS_gle+1)  
            call MPI_Bcast( Tmat_gle(1,1), j, MPI_DOUBLE_PRECISION, 0,   & 
              glemat_comm , ierr )
            call MPI_Bcast( Smat_gle(1,1), j, MPI_DOUBLE_PRECISION, 0,   & 
              glemat_comm , ierr )

          endif  ! pimd_head 

        endif

        !write(*,*) 'finished init: ',rank

        if(pimd_head) then
          close(file_out)
        endif

 77     FORMAT((E13.5,A2))

    END SUBROUTINE


    SUBROUTINE Langevin_thermostat_integrate(idnode,ibead,iatm1,iatm2,nvx,nvy,nvz,vxx,vyy,vzz, & 
                  gle_vx,gle_vy,gle_vz,fix_com,fict_mass_nmode,alpha_cell)

      use multibead,    only: comm_bead

      IMPLICIT NONE

#ifdef MPI
#  include "mpif.h"
#endif

      integer, intent(in) :: ibead,iatm1,iatm2,idnode  
      real*8 , intent(inout) :: nvx(:),nvy(:),nvz(:)
      real*8 , intent(inout) :: vxx(:),vyy(:),vzz(:)
      real*8 , intent(inout) :: gle_vx(:),gle_vy(:),gle_vz(:) 
      real*8 , intent(inout) :: alpha_cell 
      real*8, intent(in) :: fict_mass_nmode(:,:)
      logical, intent(in) :: fix_com

      integer :: i,j,l,ns,nk,ierr
      real*8 :: c2oversqrtm,rnd(3),vcm(4),vcm_(4),mass,q1,q2

      ! Langevin dynamics on normal modes 
      if(thermoType.eq.1) then
        call transform_velocity_from_cart_to_nmode(vxx,vyy,vzz,nvx,nvy,nvz)
        if(ibead.eq.1) then
          if(thermo_centroid .eq. 1) then
          ! PILE-L  
            do i=iatm1,iatm2
              c2oversqrtm = c2*one_over_sqrt_nmmass(i)
              CALL getRandom(c2oversqrtm,rnd)
              nvx(i) = c1*nvx(i)+rnd(1)
              nvy(i) = c1*nvy(i)+rnd(2)
              nvz(i) = c1*nvz(i)+rnd(3)
            enddo
          else
          ! PILE-G  
          ! NVE (CMD)  
            ! other thermostats for centroid, e.g. CMD
          endif
        else ! ibead != 1
          do i=iatm1,iatm2
            c2oversqrtm = c2*one_over_sqrt_nmmass(i)
            CALL getRandom(c2oversqrtm,rnd)
            nvx(i) = c1*nvx(i)+rnd(1)
            nvy(i) = c1*nvy(i)+rnd(2)
            nvz(i) = c1*nvz(i)+rnd(3)
          enddo
        endif
        ! transform to cartesian 
        call transform_velocity_from_nmode_to_cart(nvx,nvy,nvz,vxx,vyy,vzz)

        if(fix_com) then
          vcm(1:4)=0.0
          do i=iatm1,iatm2
            mass = fict_mass_nmode(i,ibead)
            vcm(4) = vcm(4) + mass
            vcm(1) = vcm(1) + vxx(i)*mass
            vcm(2) = vcm(2) + vyy(i)*mass
            vcm(3) = vcm(3) + vzz(i)*mass
          enddo
#ifdef MPI
          call MPI_ALLREDUCE(vcm(1),vcm_(1),4,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,MPI_COMM_WORLD ,ierr)
#endif

          vcm(1:3) = vcm_(1:3)/vcm_(4)
          do i=iatm1,iatm2
            vxx(i) = vxx(i) - vcm(1)
            vyy(i) = vyy(i) - vcm(2)
            vzz(i) = vzz(i) - vcm(3)
          enddo
        endif

 ! **************************  PI+GLE  ******************************** 
      else if(thermoType.eq.2) then

        ! I should be here in cartesian 
      ! in PI+GLE, gle_vi is always kept in cartesian space, where the
      ! thermostat is applied (not in NM coords, PIGLET is applied in NM coords)   
        ns = (iatm1-1)*(numS_gle+1)+1
        do i=iatm1,iatm2
          gle_vx(ns) = vxx(i)
          gle_vy(ns) = vyy(i)
          gle_vz(ns) = vzz(i)
          ns = ns + (numS_gle+1)
        enddo

        do i=iatm1,iatm2 
         ! thermostat half step
         Sx(:) = 0.d0
         Sy(:) = 0.d0
         Sz(:) = 0.d0
         ! silly, but doint it like this for now
         do l=1,numS_gle+1
           CALL getRandom(1.d0,rnd)
           Rx(l) = rnd(1)
           Ry(l) = rnd(2)
           Rz(l) = rnd(3)
         enddo 
         do l=1,numS_gle+1
          ns=(i-1)*(numS_gle+1)+1
          do j=1,numS_gle+1
           q1 = Tmat_gle(l,j)
           q2 = Smat_gle(l,j)*one_over_sqrt_nmmass(i)
           Sx(l) = Sx(l) + q1*gle_vx(ns) + q2*Rx(j)
           Sy(l) = Sy(l) + q1*gle_vy(ns) + q2*Ry(j)
           Sz(l) = Sz(l) + q1*gle_vz(ns) + q2*Rz(j)
           ns=ns+1
          enddo
         enddo
         gle_vx((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sx(1:numS_gle+1)
         gle_vy((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sy(1:numS_gle+1)
         gle_vz((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sz(1:numS_gle+1)
        enddo


        if(fix_com) then
          vcm(:)=0.0
          nk = (iatm1-1)*(numS_gle+1)+1
          do i=iatm1,iatm2
            mass = fict_mass_nmode(i,1)
            ns=1
            do l=1,numS_gle+1 
              vcm_gle(ns) = vcm_gle(ns) + gle_vx(nk)*mass
              vcm_gle(ns+1) = vcm_gle(ns+1) + gle_vy(nk)*mass
              vcm_gle(ns+2) = vcm_gle(ns+2) + gle_vz(nk)*mass
              ns = ns + 3
              nk = nk + 1
            enddo
            vcm_gle(ns) = vcm_gle(ns) + mass
          enddo
#ifdef MPI
          call MPI_ALLREDUCE(vcm_gle(1),vcm_gle_(1),3*(numS_gle+1)+1,     &
                      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD ,ierr)   
#endif

          vcm_gle(1:3*(numS_gle+1)) =                                         &
           vcm_gle_(1:3*(numS_gle+1))/vcm_gle_(3*(numS_gle+1)+1)
          nk = (iatm1-1)*(numS_gle+1)+1
          do i=iatm1,iatm2
            ns = 1
            do l=1,numS_gle+1
              gle_vx(nk) = gle_vx(nk) - vcm_gle(ns)
              gle_vy(nk) = gle_vy(nk) - vcm_gle(ns+1)
              gle_vz(nk) = gle_vz(nk) - vcm_gle(ns+2)
              ns = ns + 3
              nk = nk + 1 
            enddo
          enddo
        endif

        ns = (iatm1-1)*(numS_gle+1)+1
        do i=iatm1,iatm2
          vxx(i) = gle_vx(ns)
          vyy(i) = gle_vy(ns)
          vzz(i) = gle_vz(ns)
          ns = ns + (numS_gle+1)
        enddo

      elseif(thermoType .eq. 4) then

        if(ibead.eq.1) then

          if(rank.eq.0) then
            CALL getRandom(c2_cell,rnd)
            alpha_cell = c1_cell*alpha_cell+rnd(1)
          endif

          call MPI_BCAST(alpha_cell,1,MPI_DOUBLE_PRECISION,  &
                       0,comm_bead,ierr)

          if(thermo_centroid .eq. 1) then
            do i=iatm1,iatm2
              c2oversqrtm = c2*one_over_sqrt_nmmass(i)
              CALL getRandom(c2oversqrtm,rnd)
              nvx(i) = c1*nvx(i)+rnd(1)
              nvy(i) = c1*nvy(i)+rnd(2)
              nvz(i) = c1*nvz(i)+rnd(3)
            enddo
          else
            ! other thermostats for centroid, e.g. CMD
          endif
        else ! ibead != 1
          do i=iatm1,iatm2
            c2oversqrtm = c2*one_over_sqrt_nmmass(i)
            CALL getRandom(c2oversqrtm,rnd)
            nvx(i) = c1*nvx(i)+rnd(1)
            nvy(i) = c1*nvy(i)+rnd(2)
            nvz(i) = c1*nvz(i)+rnd(3)
          enddo
        endif
        ! transform to cartesian 

        if(fix_com) then
          vcm(1:4)=0.0
          do i=iatm1,iatm2
            mass = fict_mass_nmode(i,ibead)
            vcm(4) = vcm(4) + mass
            vcm(1) = vcm(1) + nvx(i)*mass
            vcm(2) = vcm(2) + nvy(i)*mass
            vcm(3) = vcm(3) + nvz(i)*mass
          enddo
#ifdef MPI
          call MPI_ALLREDUCE(vcm(1),vcm_(1),4,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,MPI_COMM_WORLD ,ierr)
#endif

          vcm(1:3) = vcm_(1:3)/vcm_(4)
          do i=iatm1,iatm2
            nvx(i) = nvx(i) - vcm(1)
            nvy(i) = nvy(i) - vcm(2)
            nvz(i) = nvz(i) - vcm(3)
          enddo
        endif

 ! *******************  PIGLET-NVT anf NPT  ******************************** 
      else if(thermoType.eq.5 .or. thermoType.eq.6) then


        if(ibead.eq.1.and.thermoType.eq.6) then

          if(rank.eq.0) then
            CALL getRandom(c2_cell,rnd)
            alpha_cell = c1_cell*alpha_cell+rnd(1)
          endif

          call MPI_BCAST(alpha_cell,1,MPI_DOUBLE_PRECISION,  &
                       0,comm_bead,ierr)
        endif

        ns = (iatm1-1)*(numS_gle+1)+1
        do i=iatm1,iatm2
          gle_vx(ns) = nvx(i)
          gle_vy(ns) = nvy(i)
          gle_vz(ns) = nvz(i)
          ns = ns + (numS_gle+1)
        enddo

        do i=iatm1,iatm2 
         ! thermostat half step
         Sx(:) = 0.d0
         Sy(:) = 0.d0
         Sz(:) = 0.d0
         ! silly, but doint it like this for now
         do l=1,numS_gle+1
           CALL getRandom(1.d0,rnd)
           Rx(l) = rnd(1)
           Ry(l) = rnd(2)
           Rz(l) = rnd(3)
         enddo 
         do l=1,numS_gle+1
          ns=(i-1)*(numS_gle+1)+1
          do j=1,numS_gle+1
           q1 = Tmat_gle(l,j)
           q2 = Smat_gle(l,j)*one_over_sqrt_nmmass(i)
           Sx(l) = Sx(l) + q1*gle_vx(ns) + q2*Rx(j)
           Sy(l) = Sy(l) + q1*gle_vy(ns) + q2*Ry(j)
           Sz(l) = Sz(l) + q1*gle_vz(ns) + q2*Rz(j)
           ns=ns+1
          enddo
         enddo
         gle_vx((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sx(1:numS_gle+1)
         gle_vy((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sy(1:numS_gle+1)
         gle_vz((i-1)*(numS_gle+1)+1:i*(numS_gle+1)) = Sz(1:numS_gle+1)
        enddo


        if(fix_com) then
          vcm(:)=0.0
          nk = (iatm1-1)*(numS_gle+1)+1
          do i=iatm1,iatm2
            mass = fict_mass_nmode(i,1)
            ns=1
            do l=1,numS_gle+1 
              vcm_gle(ns) = vcm_gle(ns) + gle_vx(nk)*mass
              vcm_gle(ns+1) = vcm_gle(ns+1) + gle_vy(nk)*mass
              vcm_gle(ns+2) = vcm_gle(ns+2) + gle_vz(nk)*mass
              ns = ns + 3
              nk = nk + 1
            enddo
            vcm_gle(ns) = vcm_gle(ns) + mass
          enddo
#ifdef MPI
          call MPI_ALLREDUCE(vcm_gle(1),vcm_gle_(1),3*(numS_gle+1)+1,     &
                      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD ,ierr)   
#endif

          vcm_gle(1:3*(numS_gle+1)) =                                         &
           vcm_gle_(1:3*(numS_gle+1))/vcm_gle_(3*(numS_gle+1)+1)
          nk = (iatm1-1)*(numS_gle+1)+1
          do i=iatm1,iatm2
            ns = 1
            do l=1,numS_gle+1
              gle_vx(nk) = gle_vx(nk) - vcm_gle(ns)
              gle_vy(nk) = gle_vy(nk) - vcm_gle(ns+1)
              gle_vz(nk) = gle_vz(nk) - vcm_gle(ns+2)
              ns = ns + 3
              nk = nk + 1 
            enddo
          enddo
        endif

        ns = (iatm1-1)*(numS_gle+1)+1
        do i=iatm1,iatm2
          nvx(i) = gle_vx(ns)
          nvy(i) = gle_vy(ns)
          nvz(i) = gle_vz(ns)
          ns = ns + (numS_gle+1)
        enddo

      else
         call error(idnode,705)
      endif

    END SUBROUTINE

    SUBROUTINE Langevin_Barostat_integrate() 
      IMPLICIT NONE

    END SUBROUTINE

    SUBROUTINE  Langevin_evolve_first(idnode,ibead,iatm1,iatm2,   &
                  xxx,yyy,zzz,npx,npy,npz,nvx,nvy,nvz,vxx,vyy,vzz,   &
                  fxx,fyy,fzz,nfx,nfy,nfz,alpha_cell,volume,Pint  &
                  ,Pvir,Epot_deriv ) 
  
      use multibead,    only: comm_bead

      IMPLICIT NONE

#ifdef MPI
#  include "mpif.h"
#endif

      integer, intent(in) :: ibead,iatm1,iatm2,idnode
      real*8 , intent(inout) :: xxx(:),yyy(:),zzz(:)
      real*8 , intent(inout) :: npx(:),npy(:),npz(:)
      real*8 , intent(inout) :: nvx(:),nvy(:),nvz(:)
      real*8 , intent(inout) :: vxx(:),vyy(:),vzz(:)
      real*8 , intent(in) :: fxx(:),fyy(:),fzz(:)
      real*8 , intent(in) :: nfx(:),nfy(:),nfz(:),Pvir,Epot_deriv
      real*8, intent(inout) :: alpha_cell,volume
      real*8, intent(out) :: Pint

      integer :: iatom,j,l,ns,nk,ierr
      real*8 :: dtx2,wfac,wfac2,vforce,p2,com1(10),com2(10),dtx6
      real*8 :: vnew,pnew

      dtx2 = dt*0.5d0
      ! Langevin dynamics on normal modes 
      if(thermoType.eq.1 .OR.    &  
         thermoType.eq.2 .OR.    &
         thermoType.eq.3) then

        ! propagate velocities 
        ! not using harmonic propagator for now
        do iatom = iatm1, iatm2
          wfac = dtx2 * one_over_nmmass(iatom)
          vxx(iatom) = vxx(iatom) + fxx(iatom)*wfac
          vyy(iatom) = vyy(iatom) + fyy(iatom)*wfac
          vzz(iatom) = vzz(iatom) + fzz(iatom)*wfac
        enddo

        call transform_position_from_cart_to_nmode(xxx,yyy,zzz,npx,npy,npz)
        call transform_velocity_from_cart_to_nmode(vxx,vyy,vzz,nvx,nvy,nvz)

        do iatom = iatm1, iatm2
          vnew = coswk*nvx(iatom) - wksinwk*npx(iatom)
          pnew = coswk*npx(iatom) + sinoverwk*nvx(iatom)
          nvx(iatom) = vnew
          npx(iatom) = pnew

          vnew = coswk*nvy(iatom) - wksinwk*npy(iatom)
          pnew = coswk*npy(iatom) + sinoverwk*nvy(iatom)
          nvy(iatom) = vnew
          npy(iatom) = pnew

          vnew = coswk*nvz(iatom) - wksinwk*npz(iatom)
          pnew = coswk*npz(iatom) + sinoverwk*nvz(iatom)
          nvz(iatom) = vnew
          npz(iatom) = pnew
        enddo

        call transform_position_from_nmode_to_cart(npx,npy,npz,xxx,yyy,zzz)
        call transform_velocity_from_nmode_to_cart(nvx,nvy,nvz,vxx,vyy,vzz)
      

      ! PILE_NPT in NM 
      elseif(thermoType.eq.4 .or. thermoType.eq.6 ) then

        if(ibead.eq.1) then

          dtx6 = dt/6.0d0
          vforce = 0.d0
          p2=0.d0
          do iatom=iatm1,iatm2
            vforce = vforce &
              + nfx(iatom)*(nvx(iatom)+dtx6*nfx(iatom)*one_over_nmmass(iatom))  & 
              + nfy(iatom)*(nvy(iatom)+dtx6*nfy(iatom)*one_over_nmmass(iatom))  & 
              + nfz(iatom)*(nvz(iatom)+dtx6*nfz(iatom)*one_over_nmmass(iatom))  
            p2 = p2 + (nvx(iatom)**2 + nvy(iatom)**2 +        &
                    nvz(iatom)**2) / one_over_nmmass(iatom)   
          enddo
          vforce = vforce*dtx2*dtx2
          com1(1) = vforce
          com1(2) = p2

#ifdef MPI
          call MPI_ALLREDUCE(com1(1),com2(1),2,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,comm_bead ,ierr)
#endif

          vforce = com2(1)
          Pint = (com2(2)/npslices + 2.0d0*Epot_deriv - Pvir) / 3.0d0 / volume

          alpha_cell = alpha_cell                     & 
            + 1.5d0*npslices*dt*(volume*(Pint-press_ext)+1.0d0/beta) + vforce

          ! forces do not contain spring term
          do iatom = iatm1, iatm2
            wfac = dtx2 * one_over_nmmass(iatom)
            nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
            nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
            nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
          enddo

          wfac = exp(dt*alpha_cell/cell_mass) 
          wfac2 = sinh( dt*alpha_cell/cell_mass  ) *cell_mass/alpha_cell
          do iatom = iatm1, iatm2
            npx(iatom) = npx(iatom) * wfac + nvx(iatom)*wfac2
            nvx(iatom) = nvx(iatom) / wfac

            npy(iatom) = npy(iatom) * wfac + nvy(iatom)*wfac2
            nvy(iatom) = nvy(iatom) / wfac

            npz(iatom) = npz(iatom) * wfac + nvz(iatom)*wfac2
            nvz(iatom) = nvz(iatom) / wfac
          enddo          
          volume = volume*wfac*wfac*wfac


        else

          do iatom = iatm1, iatm2
            wfac = dtx2 * one_over_nmmass(iatom)
            nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
            vnew = coswk*nvx(iatom) - wksinwk*npx(iatom)
            pnew = coswk*npx(iatom) + sinoverwk*nvx(iatom) 
            nvx(iatom) = vnew
            npx(iatom) = pnew 
            nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
            vnew = coswk*nvy(iatom) - wksinwk*npy(iatom)
            pnew = coswk*npy(iatom) + sinoverwk*nvy(iatom) 
            nvy(iatom) = vnew
            npy(iatom) = pnew 
            nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
            vnew = coswk*nvz(iatom) - wksinwk*npz(iatom)
            pnew = coswk*npz(iatom) + sinoverwk*nvz(iatom) 
            nvz(iatom) = vnew
            npz(iatom) = pnew 
          enddo

        endif

        com1(1) = volume
        com1(2) = Pint 
        com1(3) = alpha_cell 
        call MPI_Bcast( com1(1), 3, MPI_DOUBLE_PRECISION, 0,   &
               MPI_COMM_WORLD,ierr )
        volume = com1(1)
        Pint = com1(2)
        alpha_cell = com1(3)  

      ! PIGLET - NVT 
      elseif(thermoType.eq.5 ) then

        call transform_velocity_from_nmode_to_cart(nvx,nvy,nvz,vxx,vyy,vzz)

        ! propagate velocities 
        do iatom = iatm1, iatm2
          wfac = dtx2 * one_over_nmmass(iatom)
          vxx(iatom) = vxx(iatom) + fxx(iatom)*wfac
          vyy(iatom) = vyy(iatom) + fyy(iatom)*wfac
          vzz(iatom) = vzz(iatom) + fzz(iatom)*wfac
        enddo

        call transform_velocity_from_cart_to_nmode(vxx,vyy,vzz,nvx,nvy,nvz)

        do iatom = iatm1, iatm2
          vnew = coswk*nvx(iatom) - wksinwk*npx(iatom)
          pnew = coswk*npx(iatom) + sinoverwk*nvx(iatom)
          nvx(iatom) = vnew
          npx(iatom) = pnew

          vnew = coswk*nvy(iatom) - wksinwk*npy(iatom)
          pnew = coswk*npy(iatom) + sinoverwk*nvy(iatom)
          nvy(iatom) = vnew
          npy(iatom) = pnew

          vnew = coswk*nvz(iatom) - wksinwk*npz(iatom)
          pnew = coswk*npz(iatom) + sinoverwk*nvz(iatom)
          nvz(iatom) = vnew
          npz(iatom) = pnew
        enddo


      else
         call error(idnode,705)
      endif


    END SUBROUTINE

    SUBROUTINE  Langevin_evolve_second(idnode,ibead,iatm1,iatm2,   &
                  xxx,yyy,zzz,npx,npy,npz,nvx,nvy,nvz,vxx,vyy,vzz,   &
                  fxx,fyy,fzz,nfx,nfy,nfz,alpha_cell,volume,Pint  &
                  ,Pvir,Epot_deriv ) 
  
      use multibead,    only: comm_bead

      IMPLICIT NONE

#ifdef MPI
#  include "mpif.h"
#endif

      integer, intent(in) :: ibead,iatm1,iatm2,idnode
      real*8 , intent(inout) :: xxx(:),yyy(:),zzz(:)
      real*8 , intent(inout) :: npx(:),npy(:),npz(:)
      real*8 , intent(inout) :: nvx(:),nvy(:),nvz(:)
      real*8 , intent(inout) :: vxx(:),vyy(:),vzz(:)
      real*8 , intent(in) :: fxx(:),fyy(:),fzz(:)
      real*8 , intent(in) :: nfx(:),nfy(:),nfz(:),Pvir,Epot_deriv
      real*8, intent(inout) :: alpha_cell,volume
      real*8, intent(out) :: Pint

      integer :: iatom,j,l,ns,nk,ierr
      real*8 :: dtx2,wfac,wfac2,vforce,p2,com1(10),com2(10),dtx6
      real*8 :: vnew,pnew

        !write(200+rank,*) 'ev2 ',npx(iatm1),nvx(iatm1)
        !write(200+rank,*) '       ',alpha_cell


      dtx2 = dt*0.5d0
      ! Langevin dynamics on normal modes 
      if(thermoType.eq.1 .OR.    &
         thermoType.eq.2 .OR.    &
         thermoType.eq.3) then

        ! propagate velocities 
        do iatom = iatm1, iatm2
          wfac = dtx2 * one_over_nmmass(iatom) 
          vxx(iatom) = vxx(iatom) + fxx(iatom)*wfac
          vyy(iatom) = vyy(iatom) + fyy(iatom)*wfac
          vzz(iatom) = vzz(iatom) + fzz(iatom)*wfac
        enddo

      ! PILE_NPT in NM 
      elseif(thermoType.eq.4 .or. thermoType.eq.6 ) then

        if(ibead.eq.1) then

          dtx6 = dt/6.0d0
          vforce = 0.d0
          p2=0.d0
          do iatom=iatm1,iatm2
            vforce = vforce &
              + nfx(iatom)*(nvx(iatom)+dtx6*nfx(iatom)*one_over_nmmass(iatom))  & 
              + nfy(iatom)*(nvy(iatom)+dtx6*nfy(iatom)*one_over_nmmass(iatom))  & 
              + nfz(iatom)*(nvz(iatom)+dtx6*nfz(iatom)*one_over_nmmass(iatom))  
            p2 = p2 + (nvx(iatom)**2 + nvy(iatom)**2 +        &
                    nvz(iatom)**2) / one_over_nmmass(iatom)   
          enddo
          vforce = vforce*dtx2*dtx2
          com1(1) = vforce
          com1(2) = p2

#ifdef MPI
          call MPI_ALLREDUCE(com1(1),com2(1),2,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,comm_bead ,ierr)
#endif

          vforce = com2(1)
          Pint = (com2(2)/npslices + 2.0d0*Epot_deriv - Pvir) / 3.0d0 / volume

          alpha_cell = alpha_cell                     & 
            + 1.5d0*npslices*dt*(volume*(Pint-press_ext)+1.0d0/beta) + vforce

          ! forces do not contain spring term
          do iatom = iatm1, iatm2
            wfac = dtx2 * one_over_nmmass(iatom)
            nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
            nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
            nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
          enddo

        else

          do iatom = iatm1, iatm2
            wfac = dtx2 * one_over_nmmass(iatom)
            nvx(iatom) = nvx(iatom) + nfx(iatom)*wfac
            nvy(iatom) = nvy(iatom) + nfy(iatom)*wfac
            nvz(iatom) = nvz(iatom) + nfz(iatom)*wfac
          enddo

        endif

#ifdef MPI
        com1(1) = Pint 
        com1(2) = alpha_cell 
        call MPI_Bcast( com1(1), 2, MPI_DOUBLE_PRECISION, 0,   &
               MPI_COMM_WORLD,ierr )
        Pint = com1(1)
        alpha_cell = com1(2)  
#endif

        !write(200+rank,'(7g)') alpha_cell,Pint,vforce,1.0d0/beta,volume*(Pint-press_ext) &
        ! ,npslices*dt*(volume*(Pint-press_ext)+1.0d0/beta),volume

      ! PIGLET - NVT 
      elseif(thermoType.eq.5 ) then

        call transform_velocity_from_nmode_to_cart(nvx,nvy,nvz,vxx,vyy,vzz)

        ! propagate velocities 
        do iatom = iatm1, iatm2
          wfac = dtx2 * one_over_nmmass(iatom)
          vxx(iatom) = vxx(iatom) + fxx(iatom)*wfac
          vyy(iatom) = vyy(iatom) + fyy(iatom)*wfac
          vzz(iatom) = vzz(iatom) + fzz(iatom)*wfac
        enddo

        call transform_velocity_from_cart_to_nmode(vxx,vyy,vzz,nvx,nvy,nvz)

      else
         call error(idnode,705)
      endif

    END SUBROUTINE

    function random(seed1,seed2)
      implicit none

      integer, intent(in),optional :: seed1,seed2 
      integer s1,s2
      real*8 random
      logical new
      real*4 u(97)
      integer i,j,k,l,ir,jr,m,ii,jj 
      real*8 s,t,c,cd,cm,uni  
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

      if(new)then
     
!     initial values of i,j,k must be in range 1 to 178 (not all 1)
!     initial value of l must be in range 0 to 168.
!        i=12
!        j=34
!        k=56
!        l=78

        if(present(seed1)) then
         s1 = MOD(seed1, 31328) 
         i = MOD(s1/177, 177) + 2
         j = MOD(s1    , 177) + 2
        else
         i=12
         j=34
        endif 

        if(present(seed2)) then
         s2 = MOD(seed2, 30081) 
         k = MOD(s2/169, 178) + 1
         l = MOD(s2,     169)
        else
         k=56
         l=78
        endif
     
        ir=97
        jr=33
        new=.false.
        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
        random=0.d0
      else
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        random=dble(uni)
      endif
      return
    end function


    ! marsaglia alg taken from duni()
    ! reimplementing here to add random initial seed
    SUBROUTINE getRandom(w,r) 
      IMPLICIT NONE

      real*8, intent(in) :: w
      real*8, intent(out) :: r(3)
      integer i 
      real*8 twopi,r1 !,duni
      parameter ( twopi = 6.283185307179586d0 )
      !external duni 

      do i=1,3
        r1=random()  ! duni() 
! mmorales: can duni return zero? 
        do while(r1.lt.1.0d-10)
         r1=random()   !duni() 
        enddo 
        !r(i) = w* COS( twopi*duni() ) * SQRT( -2.d0*LOG(r1) )
        r(i) = w* COS( twopi*random() ) * SQRT( -2.d0*LOG(r1) )
      enddo

    END SUBROUTINE

    !SUBROUTINE generate_seed(rank,nproc,seeds) 
    SUBROUTINE generate_seed(seeds) 
      IMPLICIT NONE

#ifdef MPI
#  include "mpif.h"
#endif

!      integer, intent(in) :: rank, nproc
      integer, intent(out) :: seeds(:)
      integer i,ierr 
      integer,parameter :: s = 86456
      !integer, external :: irand
      integer :: irand

      ! develop better alternative later
      if(rank.eq.0) then 

        call srand(s)
        do i=1,2*nproc
          seeds(i) = irand(0) 
        enddo

        call MPI_Bcast( seeds, 2*nproc, MPI_INTEGER, 0,   &
               MPI_COMM_WORLD,ierr )

      else

        call MPI_Bcast( seeds, 2*nproc, MPI_INTEGER, 0,   & 
               MPI_COMM_WORLD,ierr )

      endif


    END SUBROUTINE

    SUBROUTINE build_GLE_matrices(ierr) 
      IMPLICIT NONE

      integer, intent(out) :: ierr
      integer :: i,j,k,l
      real*8 :: tp,azero,aone,aminus,fct
      real*8, allocatable :: AI(:,:),BI(:,:) 

      ALLOCATE(AI(numS_gle+1,numS_gle+1), &
       BI(numS_gle+1,numS_gle+1))

      ! build GLE matrices
      ! 1. Use Taylor expansion + product to construct exp(-dt*A/2) 
      tp = dt*0.5d0/DBLE(MATSQP)
      azero=0.d0
      aone=1.d0
      aminus=-1.d0
      AI=0.d0
      BI=0.d0
      do i=1,numS_gle+1
       AI(i,i) = 1.d0
       BI(i,i) = 1.d0
      enddo

      j=numS_gle+1
      do i=1,TAYLORN
      !(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       fct=-tp/DBLE(i)
#ifdef BGQ
#else
       CALL DGEMM('N','N',j,j,j,fct, &
           AI,j,Amat_gle,j,azero,Tmat_gle,j)
#endif
       AI = Tmat_gle
       BI = BI+AI
      enddo

      ! make sure that elements of AI are small before finishing
      ! BI should now be converged to exp(-dt*0.5*A/MATSQP)  
      Tmat_gle = BI
      do i=1,MATSQP-1
#ifdef BGQ
#else
        CALL DGEMM('N','N',j,j,j,aone,  &
            Tmat_gle,j,BI,j,azero,AI,j)
#endif
        Tmat_gle = AI
      enddo

      ! now Tmat_gle is exp(-dt*0.5*A)  
      ! build SSt = Cp - TCpTt
#ifdef BGQ
#else
      CALL DGEMM('N','T',j,j,j,aone,  &
         Cmat_gle,j,Tmat_gle,j,azero,AI,j)
      Smat_gle = Cmat_gle
      CALL DGEMM('N','N',j,j,j,aminus,  &
         Tmat_gle,j,AI,j,aone,Smat_gle,j)
      CALL  DPOTRF('L',j,Smat_gle,j,IERR);
#endif

      ! Smat_gle is strictly lower triangular
      do i=1,j
       do k=i+1,j
        Smat_gle(i,k) = 0.d0
       enddo
      enddo

      DEALLOCATE(AI,BI)

    END SUBROUTINE

    SUBROUTINE read_GLE_MAT(Funit, n, ntot, readCp) 
      IMPLICIT NONE

      integer, intent(in) :: Funit, n, ntot
      logical, intent(in) :: readCp
      integer i,j,k,l,ierr,mn
      parameter(mn=100) 
      character :: p(mn)*(100) 
      !real*8 :: rlread 
      !integer :: ipickoff


      ! read A matrices
      do i=1,ntot
        do j=1,numS_gle+1
          ierr=ipickoff(Funit,p,k,mn)
          if(ierr.ne.0) call error(713,rank) 
          if(k.ne.numS_gle+1) call error(714,rank) 
          if(i.eq.n) then
            do l=1,numS_gle+1    
              Amat_gle(j,l) = rlread(p(l))
            enddo
          endif 
        enddo  
      enddo    

      if(.NOT.readCp) return 

      ! read C matrices
      do i=1,ntot
        do j=1,numS_gle+1
          ierr=ipickoff(Funit,p,k,mn)
          if(ierr.ne.0) call error(713,rank)
          if(k.ne.numS_gle+1) call error(714,rank)
          if(i.eq.n) then
            do l=1,numS_gle+1
              Cmat_gle(j,l) = rlread(p(l))
            enddo
          endif
        enddo
      enddo

    END SUBROUTINE

! taken from bopimc for now, not sure what's the easiest way to do this 
      function ipickoff(iu,p,n,mn)
! on entry: iu is unit number for read, mn is the size of array p
!c on exit pickoff.eq.1 means eof was encountered
! p contains as characters the parameters, there are n of them
!  blanks or commas separate characters
! line length is maximally 256 and characters after a # or ! are ignored 
      implicit none
      integer line,ipickoff,iu,n, mn,ifact,i,k,ld,ist,ln,ls
      parameter (line=256)
      character data*257, p(mn)*(*)

100   continue
      ifact=0
      n=0
      i=0
      read (iu,6,END=2) data
      ipickoff=0
6     format(a256)
! starting from end find how long the line is
      do k=line,1,-1
        if(data(k:k).ne.' ') go to 401
      enddo
401   ld=k

!     Additional safety check to stop XLF problem
!     If test above fails for all k (due to empty line) then 
!     ld will be initialized to 0, and the next write will crash.
      if(ld<1)return

! now strip off comment part
      ls=0
      do k=1,ld
      if(data(k:k).eq.'!'.or.data(k:k).eq.'#') goto 402
       ls=k
      enddo
402   ls=ls+1
      data(ls:ls)=' '

1     i=i+1
      if(i.gt.ls)then
        if(n.gt.0)then
           return
         endif
         go to 100
       endif

      if(data(i:i).eq.' '.or.data(i:i).eq.',')then

       if(ifact.ne.0) then
! finish the word
        ln=i-ist
        if(ln.ge.len(p(n))) then
          ipickoff=-2
          return
        endif

        p(n)(1:ln)=data(ist:i-1)
! terminate with  blanks
        do 30 k=ln+1,len(p(n))
30      p(n)(k:k)=' '
        ifact=0
       endif

      else

       if(ifact.eq.0) then
! found a new parameter
         ifact=1
         ist=i
         n=n+1
         if(n.gt.mn) then
          ipickoff=-3
          return
         endif
       endif

      endif
      go to 1

2     ipickoff=1
      return
      end function

      function rlread(p)
      implicit none
      integer ln,pt
      real*8 rlread
      character p*(*),blank*25,both*25
      blank='                            '
      ln=index(p,' ')-1
      !If no decimal point, add one to end...
      pt=index(p,'.')
      if(pt.lt.1)then
         !Add to end...
         both=blank(1:25-ln-1)//p(1:ln)//'.'
      else
         both=blank(1:25-ln)//p(1:ln)
      endif
      read(both,'(e25.10)') rlread
      return
      end function

  END MODULE Langevin_Thermostat 

