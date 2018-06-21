!SR
!  ==  vibrational analysis subroutine ==
! following the vibrational analysis notes given here
!    http://depa.fquim.unam.mx/amyd/archivero/ANALISIS_VIBRACIONAL_3420.pdf

module vibanal

  implicit none 


  real(kind=8),private,save :: avogadro
  real(kind=8),private,save :: clight
  real(kind=8),private,save :: amu
  real(kind=8),private,save :: angstom
  real(kind=8),private,save :: pi


  integer,private,save        ::  maxstep
  real(kind=8),private,save   ::  maxforce,rmsforce,maxdis,rmsdis,disp
  real(kind=8),allocatable,private,save   ::  r0(:),fprev(:),hessian(:,:)


  integer,private,save        :: output_index,output_rank,output_freq,totvar,count,lcalc
  integer,private,save        :: countl
  logical,private,save        :: firststep

  integer,private,save ::  nvib

 ! subroutines

  public :: vib_init
  public :: vibrational_analysis
  public :: vib_finish


  contains

 
  subroutine vib_init(natms_r,nrite,idnode,mxnode,nread)
 
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

    avogadro=6.022140857d23    
    clight=2.99792458d8        !m
    amu=1.660539040d-27    ! in kg
    angstom=1.d-10         ! m 
    pi=atan(1.d0)*4.d0           
    disp=0.0052918d0         ! in angstrom

    totvar=3*natms_r
    nvib=totvar-6        ! subtract translations and rotations
    firststep=.True.
    lcalc = -1 
    count=0

    allocate(r0(totvar),fprev(totvar),hessian(totvar,totvar))

    if(idnode.eq.0) open(nread,file='CONTROL',status='old')

    if(idnode.eq.0) rewind(nread) 

     do nrecs=1,mega

        call getrec(safe,idnode,mxnode,nread,record)

        call lowcase(record,40)

        call cal_field(record,nofields,countfields,starts,ends)

        if(countfields==0)  cycle           ! takes care of blank lines

        if(record(1:1).eq.'#') then

         elseif(record(starts(1):ends(1)).eq.'vibanal') then

            tmprecord=record(starts(2):ends(2))
            if(starts(2)/=0) read(tmprecord,*) disp
        
         elseif(record(starts(1):ends(1)).eq.'finish')then

            goto 2000

       endif

     enddo

     2000 continue 
    
     if(idnode.eq.0) then
       close(nread) 

       write(nrite,*)' ==============================================    '
       write(nrite,*)'    VIBRATIONAL ANALYSIS' 
       write(nrite,*)' ==============================================    '
       write(nrite,*)
       write(nrite,*)'Displace each atom by      =   ',disp
       write(nrite,*)
       write(nrite,*)
     endif

   end subroutine vib_init


   subroutine vibrational_analysis(atmnam,mass,r,f,pot,flag)
     implicit none  

     real(kind=8) :: r(*),pot,f(*),mass(*)
     character(len=8) :: atmnam(*)
     integer,intent(out) :: flag
     integer          ::  i,j,k,l,mm,ss
 

     real(kind=8),allocatable     ::  WW(:)
     real(kind=8)     ::  sign1,IX(3)
     character(len=4) ::  nn

     real(kind=8),allocatable :: Rot(:,:),D(:,:),shessian(:,:),TransRot(:,:),r0com(:)

     real(kind=8) :: sum1(3),summ,rcom(3),iner(3,3),rm(3)
     real(kind=8) :: norm,prx,pry,prz


! remove translations and rotations 
     

!displace each atom in equilibrium structure by disp in + and - directions

      flag=0

      lcalc=-1*lcalc

      if(.not.firststep.and.lcalc==1) then

         hessian(countl,:)=(f(1:totvar)-fprev(1:totvar))/(2.d0*disp) ! f is grad

      endif 

      fprev=f(1:totvar)

      if(.not.firststep) then
         if(output_rank==0) then

            if(mod(countl,3)==1) then
               if(lcalc==-1) then
                   nn='X-DX'
                 else
                   nn='X+DX'
               endif
             endif

            if(mod(countl,3)==2) then
               if(lcalc==-1) then
                   nn='Y-DY'
                 else
                   nn='Y+DY'
               endif
             endif

            if(mod(countl,3)==0) then
               if(lcalc==-1) then
                   nn='Z-DZ'
                 else
                   nn='Z+DZ'
               endif
             endif

            write(output_index,*) 'ATOM  ',(countl+2)/3,'  COORDINATE  ',nn

         endif
       endif


      if(firststep) then
        r0=r(1:totvar)
        firststep=.False.
      endif


      count=count+1 


      if(count>2*totvar) then

! multiply with masses 
! Right now, only water (H, D, O) masses are included

!        do l=1,totvar/3
         
!         if(trim(atmnam(l))=="OW".or.trim(atmnam(l))=="O") then 
!
!             mass(3*l-2:3*l)=15.9994d0  
!
!          elseif(trim(atmnam(l))=="HW".or.trim(atmnam(l))=="H") then 
!
!             mass(3*l-2:3*l)=1.0079d0
!
!          elseif(trim(atmnam(l))=="DW".or.trim(atmnam(l))=="D") then 
!
!             mass(3*l-2:3*l)=2.0141d0
!
!         endif

!        enddo

   
        do i=1,totvar
          do j=1,totvar

             mm=(i+2)/3; ss=(j+2)/3 
       
             hessian(i,j)=hessian(i,j)/sqrt(mass(mm)*mass(ss))  

          enddo
        enddo

! symmetrize 

     do i=1,totvar
       do j=i+1,totvar

         hessian(i,j)=hessian(j,i)

       enddo
     enddo


! save hessian 

     allocate(shessian(totvar,totvar))
     shessian=hessian

! diagonalize the hessian         it uses lapack library 
         
     allocate(WW(totvar))
     call  diamat_all(hessian,WW)


  ! eigen vectors -- displacement 
!    do i = 1, totvar
!      j=(i+2)/3 
!      shessian(i,i) = 1.0d0/mass(j)
!    enddo
!
!    allocate(D(totvar,totvar))
!    ! Cartesian displacements of the normal modes
!    D = matmul(shessian,matmul(D,hessian))
!
!     write(*,*) "D"
!     write(*,*) D 
!
!    do i = 1, totvar
!       norm     = 1.0d0/sum(D(:,i)*D(:,i))
!       ! Reduced Masess
!       ! Renormalize displacements and convert in Angstrom
!       D(:,i)   = DSQRT(norm)*D(:,i)
!    enddo
!


! discard eigen vectors 

     hessian=shessian

     deallocate(shessian)

! print eigen values 

! convert to m^-1  

       WW=WW*10.d0/(avogadro*amu*angstom*angstom)    ! forces are in 10J/mol


!     write(2845,*) 
!     write(2845,*) 
!     write(2845,*) 'writing frequencies',WW(1:totvar)
!     write(2845,*) 
!     write(2845,*) 
!     write(2845,*) 'constants  '
!     write(2845,*) 'clight',clight,'avogardo',avogadro,'pi',pi,'amu',amu,'angstom',angstom

       do i=1,totvar

          sign1=sign(1.d0,WW(i))

          WW(i)=abs(WW(i))
      
          WW(i)=sign1*dsqrt(WW(i)/(4.d0*pi*pi*clight*clight))

       enddo

        WW=WW/100.d0          ! to cm^-1
 
! frequencies before removing the translations and rotations 

!       write(2845,*) 'writing frequencies',WW(1:totvar)

!      write(output_index,*) 
!      write(output_index,*) 'writing frequencies'
!      write(output_index,*) 
!      write(output_index,*) WW(1:totvar)


!     call print_normalmodes(totvar,totvar/3,atmnam,r0,D,WW)


       deallocate(WW)


! remove translational and rotational 

! calculate moment of inertia matrix 

       sum1=0.d0
       summ=0.d0
       do i=1,totvar,3
         mm=(i+2)/3
         sum1(1)=sum1(1)+mass(mm)*r0(i)
         sum1(2)=sum1(2)+mass(mm)*r0(i+1)
         sum1(3)=sum1(3)+mass(mm)*r0(i+2)
         summ=summ+mass(mm)

       enddo

       rcom=sum1/summ

       iner=0.d0
       do i=1,totvar,3

        mm=(i+2)/3

        rm(1)=r0(i)-rcom(1)
        rm(2)=r0(i+1)-rcom(2)
        rm(3)=r0(i+2)-rcom(3)

         iner(1,1)=iner(1,1)+mass(mm)*(rm(2)*rm(2)+rm(3)*rm(3))
         iner(2,2)=iner(2,2)+mass(mm)*(rm(1)*rm(1)+rm(3)*rm(3))
         iner(3,3)=iner(3,3)+mass(mm)*(rm(1)*rm(1)+rm(2)*rm(2))

         iner(1,2)=iner(1,2)-mass(mm)*rm(1)*rm(2)
         iner(1,3)=iner(1,3)-mass(mm)*rm(1)*rm(3)
         iner(2,3)=iner(2,3)-mass(mm)*rm(2)*rm(3)
       enddo

       iner(2,1)=iner(1,2)
       iner(3,1)=iner(1,3)
       iner(3,2)=iner(2,3)



       call  diamat_all(iner,IX)

!translational
       allocate(D(totvar,3))

       D=0.d0

       do i=1,totvar,3
        
         mm=(i+2)/3

         D(i,1)=sqrt(mass(mm))
         D(i+1,2)=sqrt(mass(mm))
         D(i+2,3)=sqrt(mass(mm))
      
       enddo

!normalization 

       do i=1,3

         norm=dsqrt(dot_product(D(:,i),D(:,i)))
         D(:,i)=D(:,i)/norm

       enddo


!rotational 

     allocate(Rot(totvar,3))

      do i=1,totvar,3

        mm=(i+2)/3

        rm(1)=r0(i)-rcom(1)
        rm(2)=r0(i+1)-rcom(2)
        rm(3)=r0(i+2)-rcom(3)

        prx=rm(1)*iner(1,1)+rm(2)*iner(2,1)+rm(3)*iner(3,1)
        pry=rm(1)*iner(1,2)+rm(2)*iner(2,2)+rm(3)*iner(3,2)
        prz=rm(1)*iner(1,3)+rm(2)*iner(2,3)+rm(3)*iner(3,3)

       Rot(i,1)=(pry*iner(1,3)-prz*iner(1,2))*sqrt(mass(mm))
       Rot(i+1,1)=(pry*iner(2,3)-prz*iner(2,2))*sqrt(mass(mm))
       Rot(i+2,1)=(pry*iner(3,3)-prz*iner(3,2))*sqrt(mass(mm))

       Rot(i,2)=(prz*iner(1,1)-prx*iner(1,3))*sqrt(mass(mm))
       Rot(i+1,2)=(prz*iner(2,1)-prx*iner(2,3))*sqrt(mass(mm))
       Rot(i+2,2)=(prz*iner(3,1)-prx*iner(3,3))*sqrt(mass(mm))

       Rot(i,3)=(prx*iner(1,2)-pry*iner(1,1))*sqrt(mass(mm))
       Rot(i+1,3)=(prx*iner(2,2)-pry*iner(2,1))*sqrt(mass(mm))
       Rot(i+2,3)=(prx*iner(3,2)-pry*iner(3,1))*sqrt(mass(mm))

      enddo


!normalization 

      do i=1,3 

         norm=dsqrt(dot_product(Rot(:,i),Rot(:,i)))
         Rot(:,i)=Rot(:,i)/norm

      enddo

      allocate(TransRot(totvar,6))

      TransRot=0.d0
      do i=1,6    ! trans+rot
        if(i<=3) then
          TransRot(:,i)=D(:,i)
         else
           l=i-3
           TransRot(:,i)=Rot(:,l)
        endif
      enddo

      deallocate(D,Rot)

! Schmidt orthoganalization

      allocate(D(totvar,nvib))      ! subracted translations and rotations

     call build_D_matrix(TransRot,6,D,full=.FALSE.,natoms=totvar/3)

      allocate(shessian(nvib,nvib))
     
      shessian=matmul(transpose(D),matmul(hessian,D))
      
!      deallocate(hessian,D)

      allocate(WW(nvib))
     call diamat_all(shessian,WW)

  ! eigen vectors -- displacement 
    hessian=0.d0 
     do i = 1, totvar
       mm=(i+2)/3
       hessian(i,i) = 1.0d0/sqrt(mass(mm))
     enddo

!     write(*,*) "displacement",nvib
!     do i=1,nvib
!        write(*,*) D(:,i) 
!     enddo

     ! Cartesian displacements
     D = matmul(hessian,matmul(D,shessian))

     do i = 1, nvib
        norm = dsqrt(dot_product(D(:,i),D(:,i)))
        ! Reduced Masess
        ! Renormalize displacements and convert in Angstrom
        D(:,i)   = D(:,i)/norm
     enddo
    
     ! now D contains eigen vectors in columns 


! convert to m^-1  

       WW=WW*10.d0/(avogadro*amu*angstom*angstom)    ! forces are in 10J/mol


!     write(2845,*) 
!     write(2845,*) 
!     write(2845,*) 'writing frequencies',WW(1:nvib)
!     write(2845,*) 
!     write(2845,*) 

       do i=1,nvib

          sign1=sign(1.d0,WW(i))

          WW(i)=abs(WW(i))
      
          WW(i)=sign1*dsqrt(WW(i)/(4.d0*pi*pi*clight*clight))

       enddo

        WW=WW/100.d0          ! to cm^-1
     
!       write(2845,*) 'writing frequencies',WW(1:nvib)

       open(unit=2845,file='vibfreq.dat',action='write')

       if(output_rank==0) then
          write(2845,*) 'zpe in kcal/mol',sum(WW(1:nvib))*0.00285911d0*0.5d0 
          write(2845,*) 'writing frequencies in cm-1'
          write(2845,"(3F10.2)") WW(1:nvib)

          write(output_index,*) 
          write(output_index,*) 'writing frequencies cm-1'
          write(output_index,*) 
          write(output_index,*) WW(1:nvib)

! print eigen vectors so as to read using JMOL 
       call print_normalmodes(nvib,totvar/3,atmnam,r0,D,WW)

       endif 

       deallocate(shessian,D,WW)


       flag=-1

       goto 140
      endif 



      r(1:totvar)=r0
      countl=(count+1)/2
      r(countl)=r0(countl)+disp*(-1)**(count)


!      write(2845,*) 'pot',pot
         
      140 continue   

   end subroutine vibrational_analysis


   subroutine vib_finish()
     implicit none 

     deallocate(r0,fprev)

   end subroutine vib_finish


  subroutine print_normalmodes(nvib,tot_atoms,atmnam,coor,nm_coor,freq )
     implicit none 
     character(len=300) :: file1,file2 
     character(len=200) :: cha(2)
 
     integer,intent(in) :: nvib,tot_atoms
     real(kind=8),intent(in)  :: freq(nvib),nm_coor(3*tot_atoms,nvib),coor(tot_atoms*3)
     character(len=8),intent(in) :: atmnam(tot_atoms)

     real(kind=8) :: intens
     integer      :: x1,x2,x3

     integer :: i,j,nmol,k,atoms_mol,dummy1,dummy2,natoms_mol,l,x,y,z,n,count
     integer,allocatable :: AN(:)
     character(len=17) :: rm(5),fr=' Frequencies --  '
     character(len=7) :: ck(3) 
     character(len=8) :: a  
     character(len=2) :: typ 
     real(kind=8) :: gz,val(3)  


     open(unit=2846,file='VIBJMOL.log',action='write')
     
     rm(1)=' Red. masses --  ' 
     rm(2)=' Frc consts  --  ' 
     rm(3)=' IR Inten    --  ' 
     rm(4)=' Raman Activ --  ' 
     rm(5)=' Depolar     --  '
     gz=0.d0 
     
     typ='?A'
     a=' Atom AN'
     ck(1)='      X'
     ck(2)='      Y'
     ck(3)='      Z'
     
     
     allocate(AN(tot_atoms))
     
     intens=0.d0     ! we are not calculating the intensities 

     write (2846,*)'Entering Gaussian System'
     write(2846,*)'this file is generated to visualize', &
                   ' the normal modes in JMOL software'
     write(2846,*)'Please note, that this is a "faked" output;'
     write(2846,*)'there are no intensities computed in here.'
     write(2846,*)'Standard orientation:'
     write(2846,*)'---------------------------------------', &
      '------------------------------'
     write(2846,'(A,2(5X,A),14X,A)') & 
               'Center','Atomic','Atomic','Coordinates (Angstroms)'
     write(2846,'(2(A,5X),1X,A,3X,3(11X,A))') & 
                   'Number','Number','Type','X','Y','Z'
     write(2846,*)'---------------------------------------', & 
                       '------------------------------' 

! it takes care of printing atomic numbers; only works for water at the moment.

     dummy2=0 

     do i=1,tot_atoms
       if(trim(atmnam(i))=="OW".or.trim(atmnam(i))=="O") AN(i)=8 
       if(trim(atmnam(i))=="HW".or.trim(atmnam(i))=="H") AN(i)=1 
       write(2846,22) i,AN(i),dummy2,coor(3*i-2),coor(3*i-1),coor(3*i)
     enddo

! reading unnecessary 8 lines in the input file 

      write(2846,*)'--------------------------------------------', & 
                   '-------------------------'
      write(2846,*)'      basis functions          primitive ', & 
                   'gaussians'
      write(2846,*)'      alpha electrons          beta electrons'
      write(2846,*)'********************************************', &
                   '**************************'
      write(2846,*)
      write(2846,*)'Harmonic frequencies (cm**-1), IR intensities ',&
                   '(KM/Mole),'
      write(2846,*)'Raman scattering activities (A**4/AMU), Raman ',&
                   'depolarization ratios,'
      write(2846,*)'reduced masses (AMU), force constants ',&
                   '(mDyne/A) and normal coordinates:'

 
      l=0 
      do i=1,nvib/3

        x1=3*i-2; x2=3*i-1; x3=3*i 

        write(2846,23) x1,x2,x3
        write(2846,24) typ, typ, typ 
        write(2846,25) fr,freq(x1),freq(x2),freq(x3)  
        write(2846,25) rm(1),gz,gz,gz 
        write(2846,25) rm(2),gz,gz,gz 
        write(2846,25) rm(3),gz,gz,gz 
        write(2846,25) rm(4),gz,gz,gz        ! intensities 
        write(2846,25) rm(5),gz,gz,gz 
        write(2846,26) a,(ck(n),n=1,3),(ck(n),n=1,3),(ck(n),n=1,3)
        do k=1,tot_atoms

          write(2846,27) k,AN(k),nm_coor(3*k-2,x1),nm_coor(3*k-1,x1),nm_coor(3*k,x1), &
        &                      nm_coor(3*k-2,x2),nm_coor(3*k-1,x2),nm_coor(3*k,x2), &
        &                      nm_coor(3*k-2,x3),nm_coor(3*k-1,x3),nm_coor(3*k,x3)   
        enddo 
      enddo  
      
!      write(2,*) 'Normal termination of Gaussian 98.'
      
      22   FORMAT(i5,i11,i14,4x,3(3x,f11.6))
      23   FORMAT(i22,2i23)
      24   FORMAT(20x,a2,2(21x,a2))
      25   FORMAT(a17,f9.4,2f23.4)
      26   FORMAT(a8,3a7,2(2x,3a7))
      27   FORMAT(2i4,3(f9.2,2f7.2))

! print xyz files of normal modes

     open(unit=2847,file='VIBxyz.xyz',action='write')

     write(2847,28) tot_atoms
     write(2847,29) ' Config: initial'
     do i=1,tot_atoms
       write(2847,30) atmnam(i),coor(3*i-2),coor(3*i-1),coor(3*i)
     enddo     

     do i=1,nvib

       write(2847,28) tot_atoms
       write(2847,31) ' Config:        ',i,'  Freq = ',freq(i)
       do k=1,tot_atoms
         write(2847,30) atmnam(k),nm_coor(3*k-2,i),nm_coor(3*k-1,i),nm_coor(3*k,i)
       enddo
     enddo

      28   FORMAT(i13)
      29   FORMAT(a16)
      30   FORMAT(a4,3(2f12.5))
      31   FORMAT(a16,i4,a9,2f10.4)

  end subroutine print_normalmodes



! Copied from cp2k --- diagonalization library 
! *****************************************************************************
!> \brief Diagonalize the symmetric n by n matrix a using the LAPACK
!>        library. Only the upper triangle of matrix a is used.
!>        Externals (LAPACK 3.0)
!> \param a ...
!> \param eigval ...
!> \param dac ...
!> \date    29.03.1999
!> \par Variables
!>      - a       : Symmetric matrix to be diagonalized (input; upper triangle) ->
!>      -           eigenvectors of the matrix a (output).
!>      - dac     : If true, then the divide-and-conquer algorithm is applied.
!>      - eigval  : Eigenvalues of the matrix a (output).
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE diamat_all(a,eigval,dac)
    INTEGER,PARAMETER                        :: dp=kind(0.d0)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: a
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigval
    LOGICAL, INTENT(IN), OPTIONAL            :: dac

    INTEGER                                  :: handle, info, liwork, lwork, &
                                                n, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwork
    INTEGER, EXTERNAL                        :: ilaenv
    LOGICAL                                  :: divide_and_conquer
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work

    EXTERNAL dsyev,dsyevd


    ! Get the size of the matrix a
    n = SIZE(a,1)

    ! Check the size of matrix a
    IF (SIZE(a,2) /= n) THEN
       write(*,*) "Check the size of matrix a (parameter #1)"
       stop
    END IF

    ! Check the size of vector eigval
    IF (SIZE(eigval) /= n) THEN
       write(*,*) "The dimension of vector eigval is too small"
       stop
    END IF

    ! Check, if the divide-and-conquer algorithm is requested

    IF (PRESENT(dac)) THEN
       divide_and_conquer = dac
    ELSE
       divide_and_conquer = .FALSE.
    END IF

    ! Get the optimal work storage size

    IF (divide_and_conquer) THEN
       lwork = 2*n**2 + 6*n + 1
       liwork = 5*n + 3
    ELSE
       nb = ilaenv(1,"DSYTRD","U",n,-1,-1,-1)
       lwork = (nb + 2)*n
    END IF

    ! Allocate work storage

    ALLOCATE (work(lwork))
    IF (divide_and_conquer) THEN
       ALLOCATE (iwork(liwork))
    END IF

    ! Diagonalize the matrix a

    IF (divide_and_conquer) THEN
       CALL dsyevd("V","U",n,a,n,eigval,work,lwork,iwork,liwork,info)
    ELSE
       CALL dsyev("V","U",n,a,n,eigval,work,lwork,info)
    END IF

    IF (info /= 0) THEN
       IF (divide_and_conquer) THEN
          write(*,*) "The matrix diagonalization with dsyevd failed"
          stop
       ELSE
          write(*,*) "The matrix diagonalization with dsyev failed"
          stop
       END IF
    END IF

    ! Release work storage

    DEALLOCATE (work)

    IF (divide_and_conquer) THEN
       DEALLOCATE (iwork)
    END IF

  END SUBROUTINE diamat_all

   ! *****************************************************************************
!> \brief Generates the transformation matrix from hessian in cartesian into
!>      internal coordinates (based on Gram-Schmidt orthogonalization)
!> \param mat ...
!> \param dof ...
!> \param Dout ...
!> \param full ...
!> \param natoms ...
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE build_D_matrix(mat,dof,Dout,full,natoms)
    INTEGER,PARAMETER                        :: dp=kind(0.d0)
    REAL(KIND=dp), DIMENSION(:, :),intent(in) :: mat
!    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat
    INTEGER, INTENT(IN)                      :: dof
    REAL(KIND=dp), DIMENSION(:, :),INTENT(out) :: Dout
!    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Dout
    LOGICAL, OPTIONAL                        :: full
    INTEGER, INTENT(IN)                      :: natoms

    INTEGER                                  :: handle, i, ifound, iseq, j, &
                                                nvib
    LOGICAL                                  :: my_full
    REAL(KIND=dp)                            :: norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: D

    ! Generate the missing vectors of the orthogonal basis set
    nvib = 3*natoms-dof
    ALLOCATE(work(3*natoms))
    ALLOCATE(D(3*natoms,3*natoms))
    ! Check First orthogonality in the first element of the basis set
    DO i = 1, dof
       D(:,i) = mat(:,i)
       DO j = i+1, dof
          norm = DOT_PRODUCT(mat(:,i),mat(:,j))
!          CPASSERT(ABS(norm)<thrs_motion)
       END DO
    END DO
    ! Generate the nvib orthogonal vectors
    iseq   = 0
    ifound = 0
    DO WHILE (ifound /= nvib)
       iseq = iseq + 1
!       CPASSERT(iseq<=3*natoms)
       work       = 0.0_dp
       work(iseq) = 1.0_dp
       ! Gram Schmidt orthogonalization
       DO i = 1, dof+ifound
          norm = DOT_PRODUCT(work,D(:,i))
          work(:) = work - norm * D(:,i)
       END DO
       ! Check norm of the new generated vector
       norm = SQRT(DOT_PRODUCT(work,work))
!       IF (norm>=10E4_dp*thrs_motion) THEN
        IF (norm>=1.d-6) then
          ! Accept new vector
          ifound = ifound + 1
          D(:,dof+ifound) = work / norm
       END IF
    END DO
!    CPASSERT(dof+ifound==3*natoms)
!    IF (my_full) THEN
!       ALLOCATE(Dout(3*natoms,3*natoms))
!       Dout = D
!    ELSE
!       ALLOCATE(Dout(3*natoms,nvib))
       Dout = D(:,dof+1:)
!    END IF
    DEALLOCATE(work)
    DEALLOCATE(D)
!    DEALLOCATE(mat)
  END SUBROUTINE build_D_matrix



end module vibanal
