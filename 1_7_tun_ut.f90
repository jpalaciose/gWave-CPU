!=================================================================
! GENERAL ROUTINES
!=================================================================
! Read parameters
subroutine read_params(XOR1,YOR1,II1,JJ1,JJ1tot,DX1,DY1,NP,rank,size,comm,ierr)
  implicit none
  include "mpif.h"
  integer XOR1,YOR1
  integer II1,JJ1,JJ1tot,NP
  real(4) DX1,DY1
  integer ierr,rank,comm,size
  character(len=30) fname
  !Read values of position
  write(fname,'(A20,I2.2,A4)') './0_files/pos/posit_',rank,'.txt'
  open(1,file=fname,status='old',action='read')
  read(1,*); read(1,*) XOR1,YOR1
  read(1,*); read(1,*) II1,JJ1,JJ1tot
  read(1,*); read(1,*) DX1,DY1
  close(1)
  !
  write(fname,'(A20,I2.2,A4)') './0_files/mar/noutp_',rank,'.txt'
  open(2,file=fname,status='old',action='read')
  read(2,*)NP
  close(2)
endsubroutine
! Read MPI parameters
subroutine read_mpisize(nei,nbuf,ntot,rank,size,comm,ierr)
  implicit none
  integer nei,nbuf,ntot
  integer rank,size,comm,ierr
  character(len=30) fname
  !
  write(fname,'(A23,I2.2,A4)') './0_files/nei/neitable_',rank,'.txt'
  open(1,file=fname, status='old',action='read')
  read(1,*); read(1,*) nei
  read(1,*); read(1,*) nbuf
  read(1,*); read(1,*) ntot
  close(1)
endsubroutine
! Read index values
subroutine read_index(nei,nbuf,idnei,index,imp_I,imp_J,exp_I,exp_J,&
                        rank,size,comm,ierr)
  implicit none
  integer i_nei,item
  integer nei,nbuf,idnei(nei),index(0:nei)
  integer imp_I(nbuf),imp_J(nbuf),exp_I(nbuf),exp_J(nbuf)
  integer rank,size,comm,ierr
  character(len=30) fname
  !
  index = 0
  write(fname,'(A20,I2.2,A4)') './0_files/idx/index_',rank,'.txt'
  open(1,file=fname, status='old',action='read')
  read(1,*); read(1,*) (idnei(i_nei),i_nei= 1,nei)
  read(1,*); read(1,*) (index(i_nei),i_nei= 1,nei)
  read(1,*); read(1,*) (imp_I(item),item= 1,nbuf)
  read(1,*); read(1,*) (imp_J(item),item= 1,nbuf)
  read(1,*); read(1,*) (exp_I(item),item= 1,nbuf)
  read(1,*); read(1,*) (exp_J(item),item= 1,nbuf)
  close(1)
endsubroutine
! Synchronize 2D array
subroutine sync_2d(MTX,II,JJtot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,&
                  rank,size,comm,ierr)
  implicit none
  include "mpif.h"
  integer I,J,II,JJtot,i_nei,item,cont
  real(4) MTX(II,JJtot)
  integer nei,nbuf,ntot,idnei(nei)
  integer index(0:nei),imp_I(nbuf),imp_J(nbuf),exp_I(nbuf),exp_J(nbuf)
  integer iS,iE,BUFlength,tag
  integer irecvreq(nei),statsend(MPI_STATUS_SIZE,nei)
  integer isendreq(nei),statrecv(MPI_STATUS_SIZE,nei)
  real(4) SENDbuf(nbuf),RECVbuf(nbuf)
  integer rank,size,comm,ierr
  !
  if(size.eq.1) return
  !
  SENDbuf=0
  RECVbuf=0
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    do item=iS,iE
      SENDbuf(item)=MTX(exp_I(item),exp_J(item))
    enddo
  enddo
  !
  tag = 1
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_ISEND(SENDbuf(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,isendreq(i_nei),ierr)
  enddo
  !
  do i_nei= 1, nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_IRECV(RECVbuf(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,irecvreq(i_nei),ierr)
  enddo
  !
  if(nei.ne.0)then
    call MPI_WAITALL(nei,isendreq(1),statrecv,ierr)
    call MPI_WAITALL(nei,irecvreq(1),statsend,ierr)
  endif
  !
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    do item=iS,iE
      MTX(imp_I(item),imp_J(item)) = RECVbuf(item)
    enddo
  enddo
endsubroutine
! Synchronize 3D array
subroutine sync_3d(MTX,II,JJtot,nei,nbuf,ntot,idnei,index,imp_I,imp_J,exp_I,exp_J,&
                  rank,size,comm,ierr)
  implicit none
  include "mpif.h"
  integer I,J,II,JJtot,i_nei,item,cont
  real(4) MTX(II,JJtot,2)
  integer nei,nbuf,ntot,idnei(nei)
  integer index(0:nei),imp_I(nbuf),imp_J(nbuf),exp_I(nbuf),exp_J(nbuf)
  integer iS,iE,BUFlength,tag
  integer irecvreq(nei),statsend(MPI_STATUS_SIZE,nei)
  integer isendreq(nei),statrecv(MPI_STATUS_SIZE,nei)
  real(4) SENDbuf1(nbuf),SENDbuf2(nbuf)
  real(4) RECVbuf1(nbuf),RECVbuf2(nbuf)
  integer rank,size,comm,ierr
  !
  if(size.eq.1) return
  !Send
  SENDbuf1=0;SENDbuf2=0
  RECVbuf1=0;RECVbuf2=0
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    do item=iS,iE
      SENDbuf1(item)=MTX(exp_I(item),exp_J(item),1)
      SENDbuf2(item)=MTX(exp_I(item),exp_J(item),2)
    enddo
  enddo
  !
  tag = 1
  !First array
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_ISEND(SENDbuf1(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,isendreq(i_nei),ierr)
  enddo
  !
  do i_nei= 1, nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_IRECV(RECVbuf1(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,irecvreq(i_nei),ierr)
  enddo
  !
  if(nei.ne.0)then
    call MPI_WAITALL(nei,isendreq(1),statrecv,ierr)
    call MPI_WAITALL(nei,irecvreq(1),statsend,ierr)
  endif
  !Second array
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_ISEND(SENDbuf2(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,isendreq(i_nei),ierr)
  enddo
  !
  do i_nei= 1, nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    BUFlength= iE + 1 - iS
    call MPI_IRECV(RECVbuf2(iS),BUFlength,MPI_REAL4,idnei(i_nei),tag,comm,irecvreq(i_nei),ierr)
  enddo
  !
  if(nei.ne.0)then
    call MPI_WAITALL(nei,isendreq(1),statrecv,ierr)
    call MPI_WAITALL(nei,irecvreq(1),statsend,ierr)
  endif
  !Receive
  do i_nei=1,nei
    iS = index(i_nei-1)+1
    iE = index(i_nei)
    do item=iS,iE
      MTX(imp_I(item),imp_J(item),1) = RECVbuf1(item)
      MTX(imp_I(item),imp_J(item),2) = RECVbuf2(item)
    enddo
  enddo
endsubroutine
!
subroutine norm(MTX,II,JJtot,IZS,IZE,JZS,JZE,result,rank,size,comm,ierr)
  implicit none
  include "mpif.h"
  integer I,J,II,JJtot,IZS,IZE,JZS,JZE
  real(4) MTX(II,JJtot)
  real(8) norm2,norm2g,result
  integer rank,size,comm,ierr
  !
  norm2 = 0d0
  do I=IZS,IZE
    do J=JZS,JZE
      norm2 = norm2 + MTX(I,J)*MTX(I,J)
    enddo
  enddo
  call MPI_Allreduce(norm2,norm2g,1,MPI_REAL8,MPI_SUM,comm,ierr)
  result = norm2g**0.5
endsubroutine
!=================================================================
! READ DATA FROM ASCII AND OKADA
!=================================================================
! Read topobathy and outpoint
subroutine READ(II,JJtot,HZ,NP,IP,JP,DX,PNAME,&
                  IZS,IZE,JZS,JZE,NBAT,NOUT,NREG,rank,size,comm,ierr)
  !
  implicit none
  include "mpif.h"
  integer II,JJtot,NP,IP(NP)
  integer JP(NP),NBAT,NOUT,NREG
  integer I,J,IZS,IZE,JZS,JZE,N
  real(4) DX,HZ(II,JJtot)
  real(4) HMIN,HMAX
  character(40) PNAME(NP)
  integer rank,size,comm,ierr
  !Read topobathymetry
  read(NBAT,'(10F10.2)') ((HZ(I,J),I=IZS,IZE),J=JZE,JZS,-1)
  close(NBAT)
  !Min and max values
  HMIN = MINVAL(HZ(IZS:IZE, JZS:JZE))
  HMAX = MAXVAL(HZ(IZS:IZE, JZS:JZE))
  ! write(6,'(A4,I1,A8,I2,A12,F7.1,A12,F7.1,A12,F7.1)') 'REG(',NREG,') rank: ',rank,'HMIN(m)= ',HMIN,'HMAX(m)= ',HMAX,'Sqrt(2gh)=',DX/SQRT(2.0*9.8*HMAX)
  !Read output point data
  do N=1,NP
    read(NOUT,*)IP(N),JP(N),PNAME(N)
  enddo
  close(NOUT)
endsubroutine
! Read deforms
subroutine INTL(II,JJtot,Z,M,N,DZ,HZ,ZM,ZOUT,&
                  IZS,IZE,JZS,JZE,NDEF,NREG,VEL,VM,N0,TH,rank,size,comm,ierr)
  !
  implicit none
  include "mpif.h"
  integer I,J,JJtot,IZS,IZE,JZS,JZE
  integer II,NDEF,NREG
  real(4) Z(II,JJtot,2),VEL(II,JJtot,2),VM(II,JJtot)
  real(4) M(II,JJtot,2),N(II,JJtot,2),DZ(II,JJtot,2)
  real(4) HZ(II,JJtot),ZM(II,JJtot),ZOUT(II,JJtot)
  real(4) N0(II,JJtot),TH(II,JJtot)
  integer rank,size,comm,ierr
  !Initial values
  VM(:,:)    = 0.0
  VEL(:,:,:) = 0.0
  ZM(:,:)    = 0.0
  ZOUT(:,:)  = 0.0
  M(:,:,:)   = 0.0
  N(:,:,:)   = 0.0
  N0(:,:)    = 0.0
  TH(:,:)    = 0.0
  !Modify the value of the topography according to the surface deformation
  read(NDEF,'(10F7.3)') ((Z(I,J,1),I=IZS,IZE),J=JZE,JZS,-1)
  close(NDEF)
  !
  Z(:,:,2)  = Z(:,:,1)
  DZ(:,:,:) = 0.0
  !
  do J=JZS,JZE
    do I=IZS,IZE
      if(NREG.NE.5) then
        HZ(I,J) = HZ(I,J) - Z(I,J,1)
      endif
      if (HZ(I,J) .LE. 0.0) then
        Z(I,J,1) = 0.0
        Z(I,J,2) = Z(I,J,1)
      else
        DZ(I,J,1) = HZ(I,J) + Z(I,J,2)
        DZ(I,J,2) = HZ(I,J) + Z(I,J,2)
      endif
    enddo
  enddo
  !
  ! if(rank.eq.0) then
  !   write(*,'(A30)')'------------------------------'
  !   write(*,'(A66,I1,A12)') 'Initializing the Data Set and Reading Deformation Data in Region (',NREG,') terminated'
  ! endif
endsubroutine
!=================================================================
! TUNAMI ROUTINES
!=================================================================
! Creating the mean water depth
subroutine HMN(HZ,HM,HN,II,JJtot,IZS,IZE,JZS,JZE,&
                rank,size,comm,ierr)
  implicit none
  include "mpif.h"
  integer I,J,II,JJtot,IZS,IZE,JZS,JZE
  real(4) HM(II,JJtot),HN(II,JJtot),HZ(II,JJtot)
  integer rank,size,comm,ierr
  !Calculate HMN
  do I=1,II
    do J=1,JJtot
      if(I.LT.II)then
        HM(I,J) = 0.5*(HZ(I,J)+HZ(I+1,J))
      else
        HM(I,J) = HZ(I,J)
      endif
      if(J.LT.JJtot)then
        HN(I,J) = 0.5*(HZ(I,J)+HZ(I,J+1))
      else
        HN(I,J) = HZ(I,J)
      endif
    enddo
  enddo
endsubroutine
!Output of water level
subroutine SNAPZ(ISNP,II,JJtot,Z,DZ,KK,KSNZ,ZOUT,HZ,&
                  IZS,IZE,JZS,JZE,NREG,KSTRT,KEND,rank,size,comm,ierr)
  !
  implicit none
  include "mpif.h"
  integer KC,I,J
  real(4) GX
  integer II,JJtot,KK,KSNZ,IZS,IZE,JZS,JZE
  integer NREG,KSTRT,KEND,ISNP
  real(4) Z(II,JJtot,2),DZ(II,JJtot,2),ZOUT(II,JJtot),HZ(II,JJtot)
  character(len=30) fname
  parameter (GX=1.0E-5)
  integer rank,size,comm,ierr
  !
  if(ISNP.EQ.0) return
  if(KK.LT.KSTRT) return
  if(KK.GT.KEND) return
  if(MOD(KK,KSNZ).NE.0) return
  !
  KC=KK/KSNZ
  !
  write(fname,'(A17,I1,I4.4,A1,I2.2,A4)') '2_results/elr/elr',NREG,KC,'_',rank,'.dat'
  open(70,file=fname,status='UNKNOWN',action='write')
  !
  ZOUT(:,:) = -99.0
  do J=JZS, JZE
    do I=IZS, IZE
      if(DZ(I,J,2).GT.GX) ZOUT(I,J) = Z(I,J,2)
    enddo
  enddo
  !
  do J=JZS, JZE
    do I=IZS, IZE
      if(HZ(I,J).LE.0.0) then
        if(ZOUT(I,J)+HZ(I,J).GT.0.001) then
          ZOUT(I,J) = ZOUT(I,J) + HZ(I,J)
        else
          ZOUT(I,J) = -99.0
        endif
      endif
    enddo
  enddo
  !
  write(70,'(680F8.3)')((ZOUT(I,J),I=IZS,IZE),J=JZE,JZS,-1)
  close(70)
endsubroutine
!Mass conservation
subroutine NLMASS(II,JJtot,Z,M,N,DZ,HZ,RX,RY,KK,NC,NREG,rank)
  !
  implicit none
  integer II,JJtot,KK,NC,NREG,I,J
  real(4) Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2)
  real(4) DZ(II,JJtot,2),HZ(II,JJtot),RX,RY
  real(4) GX,ZZ,DD
  integer rank
  parameter (GX=1.0E-5)
  !
  do J=2,JJtot
    do I=2,II
    !If the vertical wall is placed at the shoreline,
    !set the following criterion as L.T. 0.0     
      if(HZ(I,J).GE.-100.0) then
        ZZ = Z(I,J,1) - RX*(M(I,J,1)-M(I-1,J,1)) &
             - RY*(N(I,J,1)-N(I,J-1,1))
        if(ABS(ZZ) .LT. GX) ZZ = 0.0
        DD = ZZ + HZ(I,J)
        if(DD .LT. GX) DD = 0.0
        DZ(I,J,2) = DD
        Z(I,J,2)  = DD - HZ(I,J)
        !*** Checking the Stability ***
        if(ABS(Z(I,J,2)).GT.100.0) then
          NC=1
          write(*,'(A24,3I6)')'Over flow Z at (K,I,J) :', KK,I,J
          write(*,*)'Z :', Z(I,J,2)
          write(*,*)'Point Number :', (JJtot-J)*II+I
          write(*,*)'Within Region :', NREG
          write(*,*)'Computation is unstable.'
          return
        endif
      endif
    enddo
  enddo
endsubroutine
! Open boundary condition
subroutine OPENBOUND(II,JJtot,Z,M,N,HZ)
  !
  implicit none
  real(4) HH,CC,UU,ZZ,HHH
  integer KK,I,J,JP,IP
  integer II,JJtot
  real(4) Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2),HZ(II,JJtot)
  !
  do KK=1,2
    if(KK.eq.1) then
      J  = 1
      JP = J
    elseif(KK.eq.2) then
      J = JJtot
      JP = J - 1
    endif
    !
    do I=2, II
      if(HZ(I,J) .GT. 0.0) then
        HH = HZ(I,J)
        CC = SQRT(9.8*HH)
        UU = 0.5*ABS(M(I,J,2)+M(I-1,J,2))
        UU = SQRT(UU**2+N(I,JP,2)**2)
        ZZ = UU/CC
        if((J.eq.1  .AND. N(I,JP,2).GT.0.0) .OR. &
            (J.eq.JJtot .AND. N(I,JP,2).LT.0.0)) then
          Z(I,J,2) = -ZZ
        else
          Z(I,J,2) = ZZ
        endif
      endif
    enddo
    !
    if(KK.eq.1) then
      I  = 1
      IP = I
    elseif(KK.eq.2) then
      I = II
      IP = I - 1
    endif
    !
    do J=2, JJtot
      if(HZ(I,J) .GT. 0.0) then
        HHH = HZ(I,J)
        CC = SQRT(9.8*HHH)
        UU = 0.5*ABS(N(I,J,2)+N(I,J-1,2))
        UU = SQRT(UU**2+M(IP,J,2)**2)
        ZZ = UU/CC
        if((I.eq.1  .AND. M(IP,J,2).GT.0.0) .OR. &
            (I.eq.II .AND. M(IP,J,2).LT.0.0)) then
          Z(I,J,2) = -ZZ
        else
          Z(I,J,2) = ZZ
        endif
      endif
    enddo
  enddo
  !
endsubroutine
! Momentum conservation
subroutine NLMNT2(SH,GG,II,JJtot,Z,M,N,DZ,DM,DN,HZ,HM,HN,RX,RY,DT,FM)
  !
  implicit none
  real(4) GX,DM1,DM2,DN1,DN2,FN
  integer II,JJtot,I,J
  real(4) SH,GG,Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2)
  real(4) DZ(II,JJtot,2),DM(II,JJtot,2),DN(II,JJtot,2)
  real(4) HZ(II,JJtot),HM(II,JJtot),HN(II,JJtot),RX,RY,DT,FM
  parameter (GX=1.0E-5)
  !
  DM(:,:,:) = 0.0
  DN(:,:,:) = 0.0
  !
  do J=1, JJtot
    do I=1, II
      if(I .LT. II) then
        DM1 = 0.25*(Z(I,J,1)+Z(I,J,2)+Z(I+1,J,1)+Z(I+1,J,2)) &
               + 0.5*(HZ(I,J)+HZ(I+1,J))
        DM2 = 0.5*(Z(I,J,2)+Z(I+1,J,2)+HZ(I,J)+HZ(I+1,J))
      else
        DM1 = Z(I,J,2) + HZ(I,J)
        DM2 = Z(I,J,2) + HZ(I,J)
      endif
      if(DM1 .GE. GX) DM(I,J,1) = DM1
      if(DM2 .GE. GX) DM(I,J,2) = DM2
      if(J .LT. JJtot) then
        DN1 = 0.25*(Z(I,J,1)+Z(I,J,2)+Z(I,J+1,1)+Z(I,J+1,2))+ 0.5*(HZ(I,J)+HZ(I,J+1))
        DN2 = 0.5*(Z(I,J,2)+Z(I,J+1,2)+HZ(I,J)+HZ(I,J+1))
      else
        DN1 = Z(I,J,2)+HZ(I,J)
        DN2 = Z(I,J,2)+HZ(I,J)
      endif
      if(DN1 .GE. GX) DN(I,J,1) = DN1
      if(DN2 .GE. GX) DN(I,J,2) = DN2
    enddo
  enddo
  !
  FN=0.5*DT*GG*FM**2
  !------- X-DIRECTION -------
  
  do J=1, JJtot
    do I=1, II
      if(HZ(I,J).GT.SH .AND. HM(I,J).GT.SH) then
        call XMMT(GG,I,J,II,JJtot,HZ,Z,DZ,DM,M,N,RX,RY,FN)
      endif
  !------- Y-DIRECTION -------
      if(HZ(I,J).GT.SH .AND. HN(I,J).GT.SH) then
        call YMMT(GG,I,J,II,JJtot,HZ,Z,DZ,DN,M,N,RX,RY,FN)
      endif
    enddo
  enddo
endsubroutine
! Momentum conservation (X)
subroutine XMMT(GG,I,J,II,JJtot,HZ,Z,DZ,DM,M,N,RX,RY,FN)
  !
  implicit none
  real(4) RX,RY,FN,GX,GX2
  real(4) DD,XNN,FF,XM,XDM,XNE,XMM0
  integer I,J,II,JJtot
  real(4) GG
  real(4) Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2)
  real(4) DZ(II,JJtot,2),DM(II,JJtot,2),HZ(II,JJtot)
  parameter (GX=1.0E-5,GX2=1.0E-5)
  !
  if(I.EQ.II) then
    M(I,J,2)=0.0
    return
  endif
  !
  if(DZ(I,J,2).GT.GX) then
    if(DZ(I+1,J,2).GT.GX) then
      DD=DM(I,J,2)
    else
      if(Z(I,J,2)+HZ(I+1,J).GT.GX) then
        DD=Z(I,J,2)+HZ(I+1,J)
      else
        M(I,J,2)=0.0
        return
      endif
    endif
  else
    if(DZ(I+1,J,2).GT.GX) then
      if(Z(I+1,J,2)+HZ(I,J).GT.GX) then
        DD=Z(I+1,J,2)+HZ(I,J)
      else
        M(I,J,2)=0.0
        return
      endif
    else
      M(I,J,2)=0.0
      return
    endif
  endif
  !------  LINEAR TERM  ------
  if(J.NE.1) then
     if(DD.GE.GX) then
        XNN=0.25*(N(I,J,1)+N(I+1,J,1)+N(I,J-1,1)+N(I+1,J-1,1))
        FF=FN*SQRT(M(I,J,1)**2+XNN**2)/DD**(7.0/3.0)
        XM=(1.0-FF)*M(I,J,1)-GG*RX*DD*(Z(I+1,J,2)-Z(I,J,2))
     else
        M(I,J,2)=0.0
        return
     endif
  else
     M(I,J,2)=0.0
     return
  endif
  !----  NON-LINEAR TERM  ----
  !THIS IF BLOCK (++++) IS ACTIVATED ONLY IN CASE THAT
  !THE LARGE REGION IS DESCRIBED WITH LINEAR THEORY
  !++++++++++++
  if(I.GT.3.AND.I.LT.II-3.AND.J.GT.3.AND.J.LT.JJtot-3) then
    !++++++++++++
    if(DM(I,J,1).GE.GX) then
      !++++++++++++
      if(M(I,J,1).GT.0.0) then
        if(I.NE.1) then
          if(DM(I-1,J,1).GE.GX) then
            XDM=M(I-1,J,1)**2/DM(I-1,J,1)
            if(DZ(I-1,J,2).LT.GX)XDM=0.0
            if(DZ(I,J,2).LT.GX)XDM=0.0
          else
            XDM=0.0
          endif
          XM=XM-RX*(M(I,J,1)**2/DM(I,J,1)-XDM) 
        endif
      else
        if(DM(I+1,J,1).GE.GX) then
          XDM=M(I+1,J,1)**2/DM(I+1,J,1)
          if(DZ(I+2,J,2).LT.GX)XDM=0.0
          if(DZ(I+1,J,2).LT.GX)XDM=0.0
        else
          XDM=0.0
        endif
        XM=XM-RX*(XDM-M(I,J,1)**2/DM(I,J,1))
      endif
      !
      if(XNN.GT.0.0) then
        if(J.NE.2) then
          XNE=0.25*(N(I,J-1,1)+N(I+1,J-1,1)+N(I,J-2,1)+N(I+1,J-2,1))
          if(DM(I,J-1,1).GE.GX) then
            XDM=M(I,J-1,1)*XNE/DM(I,J-1,1)
            if(DZ(I,J-2,2).LT.GX) then
              XDM=0.0
              XMM0=XM-RY*(M(I,J,1)*XNN/DM(I,J,1)-XDM)
            else
              XMM0=XM-RY*(M(I,J,1)*XNN/DM(I,J,1)-XDM)
              if(DZ(I,J-1,2).LT.GX)XMM0=XM
              if(DZ(I+1,J-1,2).LT.GX)XMM0=XM
              if(DZ(I+1,J-2,2).LT.GX)XMM0=XM
            endif
            XM=XMM0/(1.0+FF) 
          else
            XDM=0.0
            XM=XM-RY*(M(I,J,1)*XNN/DM(I,J,1)-XDM)
            XM=XM/(1.0+FF)
          endif
        endif
      else
        XNE=0.25*(N(I,J+1,1)+N(I+1,J+1,1)+N(I,J,1)+N(I+1,J,1))
        if(DM(I,J+1,1).GE.GX) then
          XDM=M(I,J+1,1)*XNE/DM(I,J+1,1)
          if(DZ(I,J+1,2).LT.GX) then
            XDM=0.0
            XMM0=XM-RY*(XDM-M(I,J,1)*XNN/DM(I,J,1))
          else
            XMM0=XM-RY*(XDM-M(I,J,1)*XNN/DM(I,J,1))
            if(DZ(I,J+2,2).LT.GX)XMM0=XM
            if(DZ(I+1,J+1,2).LT.GX)XMM0=XM
            if(DZ(I+1,J+2,2).LT.GX)XMM0=XM
          endif
          XM=XMM0/(1.0+FF) 
        else
          XDM=0.0
          XM=XM-RY*(XDM-M(I,J,1)*XNN/DM(I,J,1))
          XM=XM/(1.0+FF)
        endif
      endif
    else
      XM=XM/(1.0+FF)
    endif
  else
     XM=XM/(1.0+FF)
  endif
  !LIMITING OF DISCHARGE FLUX 
  if(ABS(XM).LT.GX)XM=0.0
  if(XM.GT.10.0*DD)XM=10.0*DD
  if(XM.LT.-10.0*DD)XM=-10.0*DD
  M(I,J,2)=XM
endsubroutine
! Momentum conservation (Y)
subroutine YMMT(GG,I,J,II,JJtot,HZ,Z,DZ,DN,M,N,RX,RY,FN)
  !
  implicit none
  real(4) RX,RY,FN,GX,GX2
  real(4) DD,XN,XMM,XME,XNN0,XDN,FF
  integer I,J,II,JJtot
  real(4) GG
  real(4) Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2)
  real(4) DZ(II,JJtot,2),DN(II,JJtot,2),HZ(II,JJtot)
  PARAMETER (GX=1.0E-5,GX2=1.0E-5)
  !
  if(J.EQ.JJtot) then
    N(I,J,2)=0.0
    return
  endif
  !
  if(DZ(I,J,2).GT.GX) then
    if(DZ(I,J+1,2).GT.GX) then
      DD=DN(I,J,2)
    else
      if(Z(I,J,2)+HZ(I,J+1).GT.GX) then
        DD=Z(I,J,2)+HZ(I,J+1)
      else
        N(I,J,2)=0.0
        return
      endif
    endif
  else
    if(DZ(I,J+1,2).GT.GX) then
      if(Z(I,J+1,2)+HZ(I,J).GT.GX) then
        DD=Z(I,J+1,2)+HZ(I,J)
      else
        N(I,J,2)=0.0
        return
      endif
    else
      N(I,J,2)=0.0
      return
    endif
  endif
  !------  LINEAR TERM  ------
  if(I.NE.1) then
    if(DD.GE.GX) then
      XMM=0.25*(M(I,J,1)+M(I,J+1,1)+M(I-1,J,1)+M(I-1,J+1,1))
      FF=FN*SQRT(N(I,J,1)**2+XMM**2)/DD**(7.0/3.0)
      XN=(1.0-FF)*N(I,J,1)-GG*RY*DD*(Z(I,J+1,2)-Z(I,J,2))
    else
      N(I,J,2)=0.0
      return
    endif
  else
    N(I,J,2)=0.0
    return
  endif
  !----  NON-LINEAR TERM  ----
  !THIS IF BLOCK (++++) IS ACTIVATED ONLY IN CASE THAT
  !THE LARGE REGION IS DESCRIBED WITH LINEAR THEORY
  !++++++++++++
  if(I.GT.3.AND.I.LT.II-3.AND.J.GT.3.AND.J.LT.JJtot-3) then
    !++++++++++++
    if(DN(I,J,1).GE.GX) then
      !++++++++++++
      if(N(I,J,1).GT.0.0) then
        if(J.NE.1) then
          if(DN(I,J-1,1).GE.GX) then
            XDN=N(I,J-1,1)**2/DN(I,J-1,1)
            if(DZ(I,J-1,2).LT.GX)XDN=0.0
            if(DZ(I,J,2).LT.GX)XDN=0.0
          else
            XDN=0.0
          endif
          XN=XN-RY*(N(I,J,1)**2/DN(I,J,1)-XDN) 
        endif
      else
        if(DN(I,J+1,1).GE.GX) then
          XDN=N(I,J+1,1)**2/DN(I,J+1,1)
          if(DZ(I,J+2,2).LT.GX)XDN=0.0
          if(DZ(I,J+1,2).LT.GX)XDN=0.0
        else
          XDN=0.0
        endif
        XN=XN-RY*(XDN-N(I,J,1)**2/DN(I,J,1))
      endif
      !
      if(XMM.GT.0.0) then
        if(I.NE.2) then
          XME=0.25*(M(I-1,J,1)+M(I-1,J+1,1)+M(I-2,J,1)+M(I-2,J+1,1))
          if(DN(I-1,J,1).GE.GX) then
            XDN=N(I-1,J,1)*XME/DN(I-1,J,1)
            if(DZ(I-2,J,2).LT.GX) then
              XDN=0.0
              XNN0=XN-RX*(N(I,J,1)*XMM/DN(I,J,1)-XDN)
            else
              XNN0=XN-RX*(N(I,J,1)*XMM/DN(I,J,1)-XDN)
              if(DZ(I-2,J+1,2).LT.GX)XNN0=XN
              if(DZ(I-1,J,2).LT.GX)XNN0=XN
              if(DZ(I-1,J+1,2).LT.GX)XNN0=XN
            endif
            XN=XNN0/(1.0+FF) 
          else
            XDN=0.0
            XN=XN-RX*(N(I,J,1)*XMM/DN(I,J,1)-XDN)
            XN=XN/(1.0+FF)
          endif
        endif
      else
        XME=0.25*(M(I+1,J,1)+M(I+1,J+1,1)+M(I,J,1)+M(I,J+1,1))
        if(DN(I+1,J,1).GE.GX) then
          XDN=N(I+1,J,1)*XME/DN(I+1,J,1)
          if(DZ(I+1,J,2).LT.GX) then
            XDN=0.0
            XNN0=XN-RX*(XDN-N(I,J,1)*XMM/DN(I,J,1))
          else
            XNN0=XN-RX*(XDN-N(I,J,1)*XMM/DN(I,J,1))
            if(DZ(I+2,J,2).LT.GX)XNN0=XN
            if(DZ(I+1,J+1,2).LT.GX)XNN0=XN
            if(DZ(I+2,J+1,2).LT.GX)XNN0=XN
          endif
          XN=XNN0/(1.0+FF) 
        else
          XDN=0.0
          XN=XN-RX*(XDN-N(I,J,1)*XMM/DN(I,J,1))
          XN=XN/(1.0+FF)
        endif
      endif
    else
      XN=XN/(1.0+FF)
    endif
  else
    XN=XN/(1.0+FF)
  endif
  !---- LIMITING OF DISCHARGE FLUX 
  if(ABS(XN).LT.GX)XN=0.0
  if(XN.GT.10.0*DD)XN=10.0*DD
  if(XN.LT.-10.0*DD)XN=-10.0*DD
  N(I,J,2)=XN
endsubroutine
! Updating maximum value
subroutine ZMAX(II,JJtot,IZS,IZE,JZS,JZE,Z,ZM)
  !
  implicit none
  integer I,J,II,JJtot,IZS,IZE,JZS,JZE
  real(4) Z(II,JJtot,2),ZM(II,JJtot)
  !
  do J=JZS,JZE
    do I=IZS,IZE
      if(ZM(I,J) .LT. Z(I,J,2)) ZM(I,J) = Z(I,J,2)
    enddo
  enddo
endsubroutine
! Output of water discharge
subroutine SNAPV(ISNV,II,JJtot,VEL,M,N,DM,DN,HZ,KK,KSNV,&
                  IVS,IVE,JVS,JVE,NREG,KSTRT,KEND,rank,size,comm,ierr)
  !
  implicit none
  integer I,J,KK,KC
  real(4) UU,VV,UUU,VVV,DDM1,DDN1,DDM2,DDN2,GX
  integer ISNV,KSNV,II,JJtot,IVS,IVE,JVS,JVE,NREG,KSTRT,KEND
  real(4) M(II,JJtot,2),N(II,JJtot,2),HZ(II,JJtot)
  real(4) DM(II,JJtot,2),DN(II,JJtot,2),VEL(II,JJtot,2)
  character(40) NM1,NM2
  parameter (GX=1.0E-2)
  integer rank,size,comm,ierr
  !
  do J=2,JJtot
    do I=2,II
      if(HZ(I,J).GE.-30.0) then
        !
        UU=0.0
        VV=0.0
        !
        if(DM(I,J,2).GE.GX) then
          DDM1=DM(I,J,2)
          UU=M(I,J,2)/DDM1
        endif
        if(DN(I,J,2).GE.GX) then
          DDN1=DN(I,J,2)
          VV=N(I,J,2)/DDN1
        endif
        !
        UUU=0.0
        VVV=0.0
        !
        if(DM(I-1,J,2).GE.GX) then
          DDM2=DM(I-1,J,2)
          UUU=M(I-1,J,2)/DDM2
        endif
        if(DN(I,J-1,2).GE.GX) then
          DDN2=DN(I,J-1,2)
          VVV=N(I,J-1,2)/DDN2
        endif
        !
        VEL(I,J,1)=0.5*(UU+UUU)
        VEL(I,J,2)=0.5*(VV+VVV)
        !
        if(ABS(VEL(I,J,1)).GT.20.0) then
          if(VEL(I,J,1).GT.0.0)VEL(I,J,1)=20.0
          if(VEL(I,J,1).LT.0.0)VEL(I,J,1)=-20.0
        endif 
        if(ABS(VEL(I,J,2)).GT.20.0) then
          if(VEL(I,J,2).GT.0.0)VEL(I,J,2)=20.0
          if(VEL(I,J,2).LT.0.0)VEL(I,J,2)=-20.0
        endif
      endif
    enddo
  enddo
  !
  if(ISNV.EQ.0) return
  if(KK.LT.KSTRT) return
  if(KK.GT.KEND) return
  if(MOD(KK,KSNV).NE.0) return
  !
  KC=KK/KSNV
  !
  write(NM1,'(A17,I1,I4.4,A1,I2.2,A4)')'2_results/vel/vxr',NREG,KC,'_',rank,'.dat'
  write(NM2,'(A17,I1,I4.4,A1,I2.2,A4)')'2_results/vel/vyr',NREG,KC,'_',rank,'.dat'
  !
  open(91,file=NM1,status='UNKNOWN',action='write')
  write(91,'(680F8.3)')((VEL(I,J,1),I=IVS,IVE),J=JVE,JVS,-1)
  close(91)
  !
  open(92,file=NM2,status='UNKNOWN',action='write')
  write(92,'(680F8.3)')((VEL(I,J,2),I=IVS,IVE),J=JVE,JVS,-1)
  close(92)
endsubroutine
! Exchange for last step data to next step data
subroutine CHANGE(II,JJtot,Z,M,N,DZ)
  !
  implicit none
  integer II,JJtot
  real(4) Z(II,JJtot,2),M(II,JJtot,2),N(II,JJtot,2),DZ(II,JJtot,2)
  !
  Z(:,:,1)  = Z(:,:,2)
  M(:,:,1)  = M(:,:,2)
  N(:,:,1)  = N(:,:,2)
  DZ(:,:,1) = DZ(:,:,2)
endsubroutine
! Updating maximum value
subroutine VMAX(II,JJtot,IZS,IZE,JZS,JZE,VEL,VM)
  !
  implicit none
  real(4) VVV
  integer I,J,II,JJtot,IZS,IZE,JZS,JZE
  real(4) VEL(II,JJtot,2),VM(II,JJtot)
  !
  do J=JZS,JZE
    do I=IZS,IZE
      VVV = MIN(SQRT(VEL(I,J,1)*VEL(I,J,1)+VEL(I,J,2)*VEL(I,J,2)),20.0)
      if(VM(I,J) .LT. VVV) VM(I,J) = VVV
    enddo
  enddo
endsubroutine
!Output of the run-up height and velocity
subroutine OUTDATA(II,JJtot,IZS,IZE,JZS,JZE,ZM,NZ,ZI,HZ,VM,NV,NR)
  !
  implicit none
  integer I,J,NV,NR
  integer II,JJtot,IZS,IZE,JZS,JZE,NZ
  real(4) ZM(II,JJtot),HZ(II,JJtot),VM(II,JJtot),ZI(II,JJtot)
  !
  do J=JZS,JZE
     do I=IZS,IZE
        if(ZM(I,J)+HZ(I,J).LE.0.001) then
           ZM(I,J)=0.0
        endif
     enddo
  enddo
  !
  write(NZ,'(680F8.3)')((ZM(I,J),I=IZS,IZE),J=JZE,JZS,-1)
  close(NZ)
  !
  write(NV,'(680F8.3)')((VM(I,J),I=IZS,IZE),J=JZE,JZS,-1)
  close(NV)
  !Output maximum inundation depth
  do J=JZE,JZS,-1
     do I=IZS,IZE
        if(HZ(I,J).LE.0.0) then
           if(ZM(I,J)+HZ(I,J).GT.0.001) then
              ZI(I,J)=ZM(I,J)+HZ(I,J)
           else
              ZI(I,J)=0.0
           endif
        else
           ZI(I,J)=-99.0
        endif
     enddo
  enddo
  !
  write(NR,'(680F8.3)')((ZI(I,J),I=IZS,IZE),J=JZE,JZS,-1)
  close(NR)
endsubroutine