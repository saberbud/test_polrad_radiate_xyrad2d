      program test_polrad_radiate_xyrad2d

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 E0,EP,TH_deg,THR,SIGMA_BORN_P,SIGMA_RAD_P,pi,low,high,W
c      DOUBLE PRECISION integral,integral1
      REAL*8 A,Z,EPS,PF
      REAL*8 SIGQFS,sigqf,SIGMOT,xsmot,inelas_cxsn,in_cxsn,xsin
      REAL*8 type_cxsn
      REAL*8 qf_cxsn,xsqf,xsqfs
      REAL*8 DELP,DELTHETA
      REAL*8 xsnr(25,10),xsr(25,10)
      REAL MOTT,xsmott
      REAL*8 xs_pol_perp,xs_pol_para
      real*8 g1sol,g2sol,g1_perp_para_sol,g2_perp_para_sol
      real*8 din1,din2
      INTEGER IFLAG,xystat,xstype,g12_type
      integer maid_on,he3g_on

      pi=3.141592653354
      E0 = 5892
      EP = 1000
      TH_deg=30.0

      din1=5
      din2=5

c      E0 = 11000
c      EP = 7700
c      TH_deg=15.0

      THR = TH_deg*PI/180.0
      W=E0-EP
      PF=130.0
      EPS=10.0
      Z=2
      A=3

      fEPS=10.0
      fEPSD=15.0
      fPF=130.0

      E0gev=E0/1000.
      EPgev=EP/1000.
      xynu=E0-EP

      call INIF1F209

      write(*,*)'E0,EP,THR,xynu',E0,EP,THR,xynu

c      go to 177

      xstype=2
      CALL POLSIG(E0,EP,THR,0,0,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)"POLRAD xs:SIGMA_BORN_P,SIGMA_RAD_P=",SIGMA_BORN_P,SIGMA_RAD_P

      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'RADIATE xs: SIGNR,SIGRAD',SIGNR,SIGRAD

      call xyrad2d(E0,THR,xynu,SIGNR,SIGRAD,xstype,din1,din2)
      write(*,*)'XYRAD2D xs: SIGNR,SIGRAD',SIGNR,SIGRAD

      stop

!     The code below are some detailed tests on different modules.
!     Each test can be run by putting go to XXX in line 50.

      !write(*,*)'E0,THR',E0,THR
      xsmot=SIGMOT(E0,THR)
      write(*,*)'xsmot=',xsmot
      xsmott=MOTT(E0gev,THR)
      write(*,*)'xsmott=',xsmott

      xsin=inelas_cxsn(E0gev,EPgev,THR)
      write(*,*)'xsin=',xsin

      xsin=in_cxsn(E0gev,EPgev,THR)
      write(*,*)'xsin no qf=',xsin

c      CALL RADIATE(E0,THR,xynu,xsin,SIGRAD,xF)
c      write(*,*)'SIGRAD',SIGRAD

      xsqf=qf_cxsn(E0gev,EPgev,THR)
      write(*,*)'xsqf=',xsqf

      xstype=0
      CALL POLSIG(E0,EP,THR,0,0,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)  !unpol
      write(*,*)"POLRAD:SIGMA_BORN_P,SIGMA_RAD_P=",SIGMA_BORN_P,SIGMA_RAD_P

      xstype=1
      CALL POLSIG(E0,EP,THR,0,0,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)"POLRAD qe:SIGMA_BORN_P,SIGMA_RAD_P=",SIGMA_BORN_P,SIGMA_RAD_P

      xstype=2
      CALL POLSIG(E0,EP,THR,0,0,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)  !unpol
      write(*,*)"POLRAD no qe:SIGMA_BORN_P,SIGMA_RAD_P=",SIGMA_BORN_P,SIGMA_RAD_P

      xstype=0
      xsin_0=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'xsin_0',xsin_0

      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD

      xstype=4
      xsin_4=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'xsin_4',xsin_4
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD

      xstype=1
      xsqf_1=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'xsqf_1',xsqf_1
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD

      xstype=5
      xsqf_5=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'xsqf_5',xsqf_5
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD

      xstype=2
      xsin_2=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'xsin_2',xsin_2
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD

c      call qfs_born(xsqf,Z,1,fEPS,fEPSD,fPF,E0,EP,THR,1)
c      write(*,*)'xsqf_born qf=',xsqf

      sigqf=SIGQFS(E0,TH_deg,W,Z,A,EPS,PF)
      write(*,*)'sigqf 3He=',sigqf

      xstype=3
      sigqf_2=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'sigqf_2=',sigqf_2
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'SIGNR,SIGRAD',SIGNR,SIGRAD
c      xyN=A-Z
c      DELP=0.1
c      DELTHETA=0.0

c      call qfs(xsqfs,Z,1,E0,EP,DELP,THR,DELTHETA)
c      write(*,*)'xsqfs=',xsqfs

c      call qfs_born(xsqfs,Z,1,fEPS,fEPSD,fPF,E0,EP,THR,0)
c      write(*,*)'xsqfs_born=',xsqfs

  177 write(*,*)'Polarized test:E0,EP,THR=',E0,EP,THR   

      x=0.17197679330499369
      Q2=1.5787566475098915

      call OPENDB(Ebeam)
      call OPENDBN(Ebeam)
      IFLAG=readg1g2dataf()
      IFLAG=readg1g2dssvf()
      IFLAG=readg1g2grsvf()

      xstype=0
      xynp=1
      xynn=0
      call g1g2types(x,Q2,xynp,xynn,xstype,g1p,g2p)
      write(*,*)'type 0:g1p,g2p=',g1p,g2p

      xynp=0
      xynn=1
      call g1g2types(x,Q2,xynp,xynn,xstype,g1n,g2n)
      write(*,*)'type 0:g1n,g2n=',g1n,g2n

      xstype=10
      xynp=1
      xynn=0
      call g1g2types(x,Q2,xynp,xynn,xstype,g1p,g2p)
      write(*,*)'type 10 3hecomb:g1p,g2p=',g1p,g2p

c      xstype=10
      xynp=0
      xynn=1
      call g1g2types(x,Q2,xynp,xynn,xstype,g1n,g2n)
      write(*,*)'type 10 3hecomb:g1n,g2n=',g1n,g2n

      xstype=0
      xynp=2
      xynn=1
      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
      write(*,*)'type 0 3hecomb:g1he3,g2he3=',g1he3,g2he3

      xstype=10
      xynp=2
      xynn=1
      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
      write(*,*)'type 10 3hecomb:g1he3,g2he3=',g1he3,g2he3

      sigma_perp=xs_pol_perp(E0gev,EPgev,THR,xstype)
      write(*,*)'sigma_perp=',sigma_perp

      sigma_para=xs_pol_para(E0gev,EPgev,THR,xstype)
      write(*,*)'sigma_para=',sigma_para

      g12_type=10
      xstype=6
      maid_on=1
      he3g_on=1
      call xsg1g2_init(g12_type)
      call g1g2types_init(maid_on,he3g_on)
      sigma_perp2=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'sigma_perp2=',sigma_perp2

      xstype=7
      sigma_para2=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'sigma_para2=',sigma_para2

      g1sol=g1_perp_para_sol(E0gev,EPgev,THR,sigma_perp,sigma_para)
      g2sol=g2_perp_para_sol(E0gev,EPgev,THR,sigma_perp,sigma_para)
      write(*,*)'g1sol,g2sol',g1sol,g2sol

c      stop

c      go to 197

c      xstype=2
c      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
c      write(*,*)'type 2 3hecomb:g1he3,g2he3=',g1he3,g2he3

c      xstype=12
c      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
c      write(*,*)'type 12 3hecomb:g1he3,g2he3=',g1he3,g2he3

c      sigma_perp=xs_pol_perp(E0gev,EPgev,THR,xstype)
c      write(*,*)'sigma_perp=',sigma_perp

c      xstype=3
c      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
c      write(*,*)'type 3 3hecomb:g1he3,g2he3=',g1he3,g2he3

c      xstype=13
c      call g1g2types(x,Q2,xynp,xynn,xstype,g1he3,g2he3)
c      write(*,*)'type 13 3hecomb:g1he3,g2he3=',g1he3,g2he3

c      sigma_perp=xs_pol_perp(E0gev,EPgev,THR,xstype)
c      write(*,*)'sigma_perp=',sigma_perp

c      xstype=10   !JAM model
c      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
c      write(*,*)'para, xstype=',xstype
c      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
c      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
c      write(*,*)'xstype=',xstype

      xstype=10   !JAM model
      CALL POLSIG(E0,EP,THR,0,3,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'perp, xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
c      write(*,*)'xstype=',xstype
      write(*,*)'E0,EP,THR,xynu',E0,EP,THR,xynu

      g12_type=10
      xstype=6
      call xsg1g2_init(g12_type)
      CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
      write(*,*)'RADIATE:SIGNR,SIGRAD',SIGNR,SIGRAD

      xstype=6
      call xyrad2d(E0,THR,xynu,SIGNR,SIGRAD,xstype,din1,din2)
      write(*,*)'xyrad2d:SIGNR,SIGRAD',SIGNR,SIGRAD
      write(*,*)'xstype,g12_type',xstype,g12_type

c      call xsg1g2_init(g12_type)
      sigma_perp2=type_cxsn(E0gev,EPgev,THR,xstype)
      write(*,*)'sigma_perp2=',sigma_perp2


c      go to 197
      go to 217

      stop

      xstype=0
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      xstype=12   !DSSV
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      xstype=2
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      xstype=13   !GRSV
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      xstype=3
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      xstype=11   !MAID
      CALL POLSIG(E0,EP,THR,0,1,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
      write(*,*)'xstype=',xstype
      write(*,*)"POLRAD:TH_RAD,SIGMA_BORN_P,SIGMA_RAD_P=",THR,SIGMA_BORN_P,SIGMA_RAD_P
      write(*,*)"POLRAD:Ratio=",SIGMA_RAD_P/SIGMA_BORN_P
      write(*,*)'xstype=',xstype

      stop

  197 write(*,*)'Loop start Polrad'

      maid_on=0
      he3g_on=1
      call g1g2types_init(maid_on,he3g_on)

      open(unit=34,file='bc_out.txt',IOSTAT=xystat,status='OLD')
      if(xystat.ne.0)then
        write(*,*)'cannot open file'
        stop
      endif

      open(22,file="result.txt")

      do l=1,10
        read(34,*)xytemp,EPgev,THR,xyxs
        EP=EPgev*1000.0d0
        xynu=E0-EP

        write(*,*)'read',l,E0gev,EPgev,THR

        xstype=16
        CALL POLSIG(E0,EP,THR,0,3,SB,SP,xfp,xstype)
        sigma_perp=xs_pol_perp(E0gev,EPgev,THR,xstype)

        write(*,*)'  ',SB,sigma_perp,(SB/sigma_perp)

        CALL POLSIG(E0,EP,THR,0,1,SB1,SP1,xfp,xstype)
        sigma_para=xs_pol_para(E0gev,EPgev,THR,xstype)

        write(*,*)'  ',SB1,sigma_para,(SB1/sigma_para)

        g1sol=g1_perp_para_sol(E0gev,EPgev,THR,sigma_perp,sigma_para)
        g2sol=g2_perp_para_sol(E0gev,EPgev,THR,sigma_perp,sigma_para)

        write(*,*)'g1sol,g2sol',g1sol,g2sol

        write(22,*)g1sol,g2sol
        write(22,*)l,SB,sigma_perp,SP,(SB/sigma_perp),(SP/SB)
c        write(22,*)g1sol,g2sol
        write(22,*)l,SB1,sigma_para,SP1,(SB1/sigma_para),(SP1/SB1)

c        write(22,*)l,SB,sigma_perp,SP,(SB/sigma_perp),(SP/SB)
c        write(22,*)l,SB,sigma_perp,SP

c        do lt=1,3
c          xstype=lt-1
c          write(*,*)'xstype=',xstype
c          CALL RADIATE(E0,THR,xynu,SIGNR,SIGRAD,xF,xstype)
c          CALL POLSIG(E0,EP,THR,0,0,SIGMA_BORN_P,SIGMA_RAD_P,xfp,xstype)
c          write(*,*)'SIGMA_BORN_P,SIGMA_RAD_P',SIGMA_BORN_P,SIGMA_RAD_P
c          write(22,*)l,xstype,EPgev,THR,SIGMA_BORN_P,SIGMA_RAD_P
c        enddo

      enddo
      close(34)
      close(22)

      stop

  217 write(*,*)'Loop start: radiate xyrad2d'

      g12_type=15
      xstype=6
      maid_on=1
      he3g_on=1
      call xsg1g2_init(g12_type)
      call g1g2types_init(maid_on,he3g_on)

      open(unit=34,file='bc_out.txt',IOSTAT=xystat,status='OLD')
      if(xystat.ne.0)then
        write(*,*)'cannot open file'
        stop
      endif

      open(22,file="result.txt")

      do l=1,10
        read(34,*)xytemp,EPgev,THR,xyxs
        EP=EPgev*1000.0d0
        xynu=E0-EP

        write(*,*)'read',l,E0gev,EPgev,THR

        xstype=6  !xs perp
        CALL RADIATE(E0,THR,xynu,SB,SP,xF,xstype)
        write(*,*)'RADIATE,finished'
        call xyrad2d(E0,THR,xynu,SB2,SP2,xstype,din1,din2)
        write(*,*)'xyrad2d, finished'
        xstype=g12_type
        sigma_perp=xs_pol_perp(E0gev,EPgev,THR,xstype)

        write(*,*)'  ',sigma_perp,(SB/sigma_perp),(SB2/sigma_perp)

        xstype=7  !xs para
c        sigma_para=xs_pol_para(E0gev,EPgev,THR,g12_type)
c        write(*,*)'sigma_para',sigma_para
        CALL RADIATE(E0,THR,xynu,SB1,SP1,xF,xstype)
        write(*,*)'RADIATE,finished'
        call xyrad2d(E0,THR,xynu,SB12,SP12,xstype,din1,din2)
        write(*,*)'xyrad2d, finished'
        xstype=g12_type
        sigma_para=xs_pol_para(E0gev,EPgev,THR,xstype)

        write(*,*)'  ',sigma_para,(SB1/sigma_para),(SB12/sigma_para)

        write(22,*)l,sigma_perp,SB,SP,SB2,SP2
        write(22,*)l,sigma_para,SB1,SP1,SB12,SP12

      enddo

      close(34)
      close(22)

      stop

      end program test_polrad_radiate_xyrad2d

