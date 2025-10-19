! Tsunami code based on research of S. Koshimura et al 2004,02,06
! Modified for parallel computing by:
! - Carlos Davila
! - Julian Palacios
! - Fernando Garcia
! To compile and run this code use:
! mpif90 1_7_tun_mp.f90 1_7_tun_ut.f90 -O3 -o 1_7_tun.exe
! mpirun -n 48 ./1_7_tun.exe
!
program tunami_1d
  implicit none
  include "mpif.h"
  integer XOR1,YOR1
  integer II1,JJ1,JJ1tot,NP
  real(4) DX1,DY1
  integer nei,nbuf,ntot
  integer ierr,rank,size
  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
  ! Read values
  call read_params(XOR1,YOR1,II1,JJ1,JJ1tot,DX1,DY1,NP,rank,size,MPI_COMM_WORLD,ierr)
  call read_mpisize(nei,nbuf,ntot,rank,size,MPI_COMM_WORLD,ierr)
  ! Main program
  call main_allocate(II1,JJ1,JJ1tot,DX1,DY1,NP,nei,nbuf,ntot,rank,size,MPI_COMM_WORLD,ierr)
  ! Finalize MPI
  call MPI_Finalize(ierr)
endprogram
!
subroutine main_allocate(II1,JJ1,JJ1tot,DX1,DY1,NP,nei,nbuf,ntot,rank,size,comm,ierr)
  !
  implicit none
  include "mpif.h"
  ! Start parameters
  integer K6,KP,KL,KSNZ,KSTRT,KEND
  real(4) FM,DT,SH,GG
  parameter (K6=100,KP=20)
  parameter (FM=0.025,SH=-50.0)
  parameter (KL=3000,DT=0.20)!3600
  parameter (KSNZ=300,KSTRT=0,KEND=3000)!3600
  parameter (GG=9.81)
  ! Position parameters
  integer II1,JJ1,JJ1tot,NP
  real(4) DX1,DY1
  ! Parameters for REG 1
  real(4) HZ1(II1,JJ1tot)
  real(4) Z1(II1,JJ1tot,2),M1(II1,JJ1tot,2),N1(II1,JJ1tot,2)
  real(4) DZ1(II1,JJ1tot,2),ZOUT1(II1,JJ1tot)
  real(4) DM1(II1,JJ1tot,2),DN1(II1,JJ1tot,2)
  real(4) HM1(II1,JJ1tot),HN1(II1,JJ1tot),ZM1(II1,JJ1tot)
  real(4) VEL1(II1,JJ1tot,2),ZI1(II1,JJ1tot),VM1(II1,JJ1tot)
  real(4) N01(II1,JJ1tot),TH1(II1,JJ1tot)
  ! Output points
  integer IP(NP),JP(NP)
  ! integer IP(NP),JP(NP),PZ(NP)
  character(40) PNAME(NP)
  ! General variables
  integer KSNZ1,KSTRT1,KEND1
  integer IS1,IE1,JS1,JE1
  integer K,KPE,KK,NC1
  real(4) RX1,RY1
  ! MPI variables
  real(4) VAL(ntot)
  integer nei,nbuf,ntot,idnei(nei),index(0:nei)
  integer imp_I(nbuf),imp_J(nbuf),exp_I(nbuf),exp_J(nbuf)
  integer rank,size,comm,ierr
  character(len=30) fname
  real(4) t(17),et(KL,17)
  real(8) result
  integer I,J,tt
  ! Output point data file
  write(fname,'(A19,I2.2,A4)') './0_files/mar/outp_',rank,'.txt'
  if(NP.NE.0) open(10,file=fname,status='old',action='read')
  ! write(fname,'(A18,I2.2,A4)') './3_results/point_',rank+1,'.dat'
  ! if(NP.NE.0) open(11,file=fname,status='UNKNOWN',action='write')
  ! Bathymetry data files
  write(fname,'(A19,I2.2,A4)') './0_files/ascii/d1_',rank,'.dat'
  open(21,file=fname,status='old',action='read')
  ! Roughness coefficient data files
  ! open(28,file='./BATHY/n0.dat',status='old',action='read')
  ! open(29,file='./BATHY/theta.dat',status='old',action='read')
  ! Deform data files
  write(fname,'(A17,I2.2,A4)') './1_deforms/def1_',rank,'.dat'
  open(31,file=fname,status='old',action='read')
  !Result data files
  open(40,file='./2_results/rslt/max.dat',status='UNKNOWN',action='write')
  open(50,file='./2_results/rslt/vmax.dat',status='UNKNOWN',action='write')
  open(60,file='./2_results/rslt/inund.dat',status='UNKNOWN',action='write')
  !=================================================================
  ! Pre-processing
  !=================================================================
  ! Read index of MPI
  call read_index(nei,nbuf,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  ! Conditions
  KSNZ1=KSNZ; KSTRT1=KSTRT; KEND1=KEND
  IS1=1; IE1=II1
  JS1=merge(1,3,rank==0); JE1=JS1+JJ1-1
  KK=0
  if(rank.eq.0) write(*,'(A20,I2)') '-----Size:-----',size
  ! Read topobathy and outpoints, output: HZ
  call READ(II1,JJ1tot,HZ1,NP,IP,JP,DX1,PNAME,IS1,IE1,JS1,JE1,21,10,1,rank,size,comm,ierr)
  call sync_2d(HZ1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  ! Read deforms, output: HZ,Z,DZ
  call INTL(II1,JJ1tot,Z1,M1,N1,DZ1,HZ1,ZM1,ZOUT1,IS1,IE1,JS1,JE1,31,1,VEL1,VM1,N01,TH1,rank,size,comm,ierr)
  call sync_2d(HZ1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  call sync_3d(Z1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  call sync_3d(DZ1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  ! Creating the mean water depth, ouput: HM,HN
  call HMN(HZ1,HM1,HN1,II1,JJ1tot,IS1,IE1,JS1,JE1,rank,size,comm,ierr)
  call sync_2d(HM1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  call sync_2d(HN1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
  ! Output of water level, output
  call SNAPZ(1,II1,JJ1tot,Z1,DZ1,KK,KSNZ1,ZOUT1,HZ1,IS1,IE1,JS1,JE1,1,KSTRT1,KEND1,rank,size,comm,ierr)
  !
  if(rank.eq.0) then
    write(*,'(A40)')'========================================'
    write(*,*)'PRE-PROCESSING IS COMPLETED.'
    write(*,'(A40)')'========================================'
  endif
  !=================================================================
  ! Processing
  !=================================================================
  RX1=DT/DX1
  RY1=DT/DY1
  !
  NC1=0
  !
  write(fname,'(A22,I2.2,A4)') './2_results/rslt/time_',rank,'.dat'
  open(70,file=fname,status='UNKNOWN',action='write')
  t(16) = MPI_Wtime()
  do K=1,KL
    !
    t(1) = MPI_Wtime()
    KK=K
    ! Write simulation steps
    if(rank.eq.0) then
      if(MOD(KK,K6).eq.0) then
        KPE=INT(100.0*REAL(KK)/REAL(KL))
        write(*,'(A6,I6,A2,I3,A12)')'KK at ',KK,'--',KPE,' % completed'
      elseif(MOD(KK,50).eq.0) then
        write(*,*) KK,' STEP'
      endif
    endif
    t(2) = MPI_Wtime()
    ! Mass conservation, output: Z,DZ
    call NLMASS(II1,JJ1tot,Z1,M1,N1,DZ1,HZ1,RX1,RY1,KK,NC1,1,rank)
    t(3) = MPI_Wtime()
    call sync_3d(Z1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(DZ1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    t(4) = MPI_Wtime()
    ! Open boundary condition, output: Z
    call OPENBOUND(II1,JJ1tot,Z1,M1,N1,HZ1)
    t(5) = MPI_Wtime()
    call sync_3d(Z1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    t(6) = MPI_Wtime()
    ! Momentum conservation, output: M,N,DM,DN
    call NLMNT2(SH,GG,II1,JJ1tot,Z1,M1,N1,DZ1,DM1,DN1,&
                  HZ1,HM1,HN1,RX1,RY1,DT,FM)
    t(7) = MPI_Wtime()
    call sync_3d(M1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(N1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(DM1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(DN1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    t(8) = MPI_Wtime()
    ! Updating maximum value
    call ZMAX (II1,JJ1tot,IS1,IE1,JS1,JE1,Z1,ZM1)
    t(9) = MPI_Wtime()
    !Output of water discharge
    call SNAPV (0,II1,JJ1tot,VEL1,M1,N1,DM1,DN1,HZ1,KK,KSNZ1,IS1,IE1,JS1,JE1,1,KSTRT1,KEND1,rank,size,comm,ierr)
    t(10) = MPI_Wtime()
    !Exchange for last step data to next step data, output: Z,M,N,DZ
    call CHANGE (II1,JJ1tot,Z1,M1,N1,DZ1)
    t(11) = MPI_Wtime()
    call sync_3d(Z1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(M1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(N1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    call sync_3d(DZ1,II1,JJ1tot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,rank,size,comm,ierr)
    t(12) = MPI_Wtime()
    ! Output of water level
    call SNAPZ(1,II1,JJ1tot,Z1,DZ1,KK,KSNZ1,ZOUT1,HZ1,IS1,IE1,JS1,JE1,1,KSTRT1,KEND1,rank,size,comm,ierr)
    t(13) = MPI_Wtime()
    ! Updating maximum value
    call VMAX (II1,JJ1,IS1,IE1,JS1,JE1,VEL1,VM1)
    t(14) = MPI_Wtime()
    !
    if(NC1.eq.1) go to 999
    !
    do tt=1,13
      et(K,tt) = t(tt+1)-t(tt)
      !write(*,*) 'rango: ',rank,'time',tt,': ',et(K,tt)
    enddo
    et(K,14) = t(14)-t(1)
    write(70,'(14E12.5)') (et(K,tt),tt=1,14)
    !
  enddo
  !
  close(70)
  !
  t(17) = MPI_Wtime()
  ! Output of the run-up height and velocity
  if(NC1.EQ.0) call OUTDATA(II1,JJ1tot,IS1,IE1,JS1,JE1,ZM1,40,ZI1,HZ1,VM1,50,60)
  !
  999 continue
  !
  if(rank.eq.0) then
  write(*,'(A40)')'========================================'
    write(*,*)'PROCESSING IS COMPLETED.'
    write(*,*)'Total time:',t(17)-t(16)
    write(*,'(A40)')'========================================'
  endif
  !=================================================================

endsubroutine
