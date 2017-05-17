      !> SIGX = DIS cross section
      !!
      !! @param E beam energy (MeV)
      !! @param TH scattered angle (degrees) 
      !! @param W the photon energy, \f$ \nu \f$
      !! @param A 
      !!
      !!  REVISION HISTORY:
      !!  - 10/26/2011 Whitney Armstrong
      !!    Added documentation
      !! 
      !!  - 9/5/2013 David Flay 
      !!    Arbitrary scaling of the cross section components 
      !!    to fit the data, valid for JLab E01-012 data (Nitrogen and some 3He options)
      !!    were commented out.  If we want to use these adjustments, we should think of 
      !!    a more intuitive way to add them in, aside from this hard-coded approach. 
      !!  
      !! @ingroup QFS
      !!
      REAL*8 FUNCTION SIGX(E,TH,W,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      ALPH=1./137.03604
      PI=ACOS(-1.)

      SIG0=100.D-4
      SIG1=54.*1.D-1
      GAM0=650.
C     SIG0=111.*1.E-30
C     SIG1=60.*1.E-27
C     GAM0=550.
      R=0.10

      PIMASS=140.
      PM=939.
      AQ=250.    ! NEVER USED  

      THR=TH*PI/180.
C      IF(W.LT.1.E-5)GO TO 4
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      ARG0=W-QMS/2./PM-PIMASS-PIMASS**2/2./PM
      ARG1=ARG0/GAM0
      ARG=ARG1**2/2.
      IF(ARG1.GT.8.)THEN
      SHAPE=1.+SIG1/SIG0/ARG0
      ELSEIF(ARG1.LT.1.E-5)THEN
      SHAPE=0.
      ELSEIF(ARG1.LT.0.1)THEN
      SHAPE=SIG1*ARG0/2./GAM0**2/SIG0
      ELSE
      SHAPE=(1.-EXP(-ARG))*(1.+SIG1/SIG0/ARG0)
      ENDIF
      EKAPPA=W-QMS/2./PM
      SIGGAM=SIG0*SHAPE
      QS=QMS+W**2
      EPS=1./(1.+2.*QS*TAN(THR/2.)**2/QMS)
      FLUX=ALPH*EKAPPA*(E-W)/2./PI**2/QMS/E/(1.-EPS)
      IF(FLUX.LT.1.E-20)FLUX=0.
      SIGEE=FLUX*SIGGAM*FPHENOM(QMS)**2
C     SIGEE=FLUX*SIGGAM
      R=0.56*1.E6/(QMS+PM**2)
      FACTOR1=1.+EPS*R
      SIGEE=SIGEE*FACTOR1

C4     if (A.eq.14.0) then
C        SIGX=A*SIGEE*(2.7*(E-W)*sin(THR/2.0)/(QMS*1E-3))
C        SIGX=A*SIGEE*(1.0/sin(THR/2.0))*(E/(E-W))*(1.0/16.5)*(1.25)
C        SIGX=A*SIGEE*(1.0/sin(THR/2.0))*(E/(E-W))*1.0/16.5
C      endif 
c      else if (A.eq.3.0) then
C        SIGX=A*SIGEE*(4.0*(E-W)*sin(THR/2.0)/(QMS*1E-3))*0.8*(E/4730.0)
C        SIGX=A*SIGEE*(4.0/(SIN(THR/2.0)))*0.86*(1E+2/(E-W))*
C     &       (5009/E)*(SIN(32.0*PI/180.0)/SIN(THR))
C        SIGX=A*SIGEE*(4.0/(QMS*1E-6))*(1E-3*E)*(1E+2/(E-W))
C        IF (TH .EQ. 15.5) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C        ELSE IF (TH . EQ. 25 ) THEN 
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(W*1E-3)*0.65
C        ELSE IF (TH . EQ. 32 ) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(W*1E-3)*0.35
C        ELSE IF (TH . EQ. 45 ) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(1E+1/(E-W))
C        ENDIF
C        SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR)) 
C      else 
C        SIGX=A*SIGEE
C      endif
        SIGX=A*SIGEE

      RETURN
      END
*SIGR1
      REAL*8 FUNCTION SIGR1(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1500.
      EPSR=0.
      AR0=1000.
      AR1=1000.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=1.20D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1240.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1240.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
c
      WTHRESH=0.0
c
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR1=A*THRESH*SIGR
C      if (A.eq.14.0) then
C        SIGR1=SIGR1
C        SIGR1=SIGR1
C      else if (A.eq.3.0) then
C        SIGR1=SIGR1*(4.0/(E*1E-3))
C        SIGR1=SIGR1*(1.0/(SIN(THR)))
C        SIGR1=SIGR1*(3.0/(QMS*1E-6))
C        SIGR1=SIGR1*(SIN((32.0)*(PI/180.0))/SIN(THR))*0.8
C      else
C        SIGR1=SIGR1
C      endif

      RETURN
      END

*SIGR2
      REAL*8 FUNCTION SIGR2(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1700.
      EPSR=0.
      AR0=1200.
      AR1=1200.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=0.68D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1520.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1520.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR2=A*THRESH*SIGR
C      if (A.eq.14.0) then
C        SIGR2=SIGR2*(0.01*QMS*1E-6)*(E/W)*5
C        SIGR2=SIGR2
C      else if (A.eq.3.0) then
C        SIGR2=SIGR2*(0.01*QMS*1E-6)
C        SIGR2=SIGR2*(0.5*QMS*1E-6)
C        SIGR2=SIGR2*1.1*(SIN(THR)/SIN(15.5*(PI/180.0))) 
C      else
C        SIGR2=SIGR2
C      endif
      RETURN
      END
*SIG2N
      REAL*8 FUNCTION SIG2N(E,TH,W,Z,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=ACOS(-1.)
      THR=TH*PI/180.
      DM=1232.
      PIMASS=140.
      PM=940.
      A2=550.
      PFR=60.
      GAM2N=20.
      GAMQFR=40.
      GAMREF=300.
      GAMR=GAMREF
      SIGREF=0.20D-7
      QMSR=4.*596.8*(596.8-380.)*SIN(60.*PI/180./2.)**2
      QVSR=QMSR+380.**2
      SIGKIN=0.5*SIGMOT(596.8D0,60.*PI/180.)
      SIGKIN=SIGKIN*(QMSR/2./QVSR+TAN(60.*PI/180./2.)**2)
      SIGKIN=SIGKIN*QVSR*FD(QMSR,A2)**2
      SIGKIN=SIGKIN*GAMR/GAMREF
      SIGCON=SIGREF/SIGKIN
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      GAMQF=GAMQFR*(PF/PFR)*(SQRT(QVS)/SQRT(QVSR))
      EFFMASS=(PM+DM)/2.
      SIG=(Z*(A-Z)/A)*SIGMOT(E,THR)
      SIG=SIG*(QMS/2./QVS+TAN(THR/2.)**2)
      SIG=SIG*QVS*FD(QMS,A2)**2
      EKAPPA=W-QMS/2./PM
      CMTOT2=PM**2+2.*PM*EKAPPA
C     GAM=SQRT(GAMR**2+GAMQF**2)
      GAM=GAMR
      SIG=SIG*CMTOT2*GAM**2
      SIG=SIG/((CMTOT2-EFFMASS**2)**2+CMTOT2*GAM**2)
      SIG=SIG*(GAMR/GAM)*SIGCON
      SIG2N=SIG
      WTHRESH=QMS/4./PM
c      
      WTHRESH=0.0
c
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAM2N)
      ELSE
      THRESH=0.
      ENDIF
C      if (A.eq.14.0) then
c        SIG2N=SIG2N*THRESH*0.65
C        SIG2N=SIG2N*THRESH*12.0*(1.0+E/W)
C        SIG2N=SIG2N*THRESH*(68.0*E/1000)*(1.0/(E-W))*1E+2*(1.0+E/W)
C        SIG2N=SIG2N*THRESH*2.0*(E-W)*1E-2*(1.0+E/W)
C      else if (A.eq.3.0) then
C        SIG2N=SIG2N*THRESH*0.9*(E-W)*1E-2*(1.0+E/W)
C        SIG2N=SIG2N*THRESH*(E/5009)*(W*1E-3)
C     &        *(SIN((32.0/2.0)*(PI/180.0)))/SIN(THR/2.0)
C        SIG2N=SIG2N*THRESH*QMS*1E-6*(4018/E)
C        SIG2N=SIG2N*THRESH
C      else
C        SIG2N=SIG2N*THRESH
C      endif
      SIG2N=SIG2N*THRESH

      RETURN
      END
*SIGDEL
      REAL*8 FUNCTION SIGDEL(E,TH,W,A,EPSD,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PM     = 939.
      PIMASS = 140.
c      DM     = 1219.
      DM     = 1232.
      AD1    = 685.
      AD0    = 774.
      PI     = ACOS(-1.)
      ALPH   = 1./137.03604
      HBARC  = 197.32858
      GAMDP  = 110.
      GAMSPRD= 140.
      GAMR   = 120.
      GAMPI  = 5.
      QFDP   = 1.02D-7
      PFR    = 230.
      QMSR   = 4.*730.*(730.-390.)*SIN(37.1*PI/180./2.)**2
      QVSR   = QMSR+390.**2
      QMSRQ  = 4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ  = QMSRQ+115.**2

      NA=INT(A)
      IF(NA.EQ.1)THEN
        QFD=QFDP
        GSPRDA=0.
        AD=AD0
      ELSEIF(NA.LT.4)THEN
        QFD=QFDP
        GSPRDA=(A-1.)*GAMSPRD/3.
        AD=AD0+(A-1.)*(AD1-AD0)/3.
      ELSE
        AD=AD1
        GSPRDA=GAMSPRD
        QFD=QFDP
      ENDIF
      THR = TH*PI/180.
      QMS = 4.*E*(E-W)*SIN(THR/2.)**2
      QVS = QMS+W**2
      EKAPPA = W-QMS/2./PM
      CMTOT2 = PM**2+2.*PM*EKAPPA
C  BEGIN DELTA CALCULATION
      IF(NA.GT.1)THEN
        GAMQ=GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
      ELSE
        GAMQ=0.
      ENDIF

      EPD = E-(DM-PM)*(DM+PM)/2./PM
      EPD = EPD/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPD = EPD-EPSD
      WD  = E-EPD
      QMSPK = 4.*E*EPD*SIN(THR/2.)**2
      QVSPK = QMSPK+WD**2
C
C NOTE WIDTH INCLUDES E-DEPENDENCE,FERMI BROADENING,& SPREADING
C
      WTHRESH = 4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH = WTHRESH/2./PM
      THRESHD = 1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH = WTHRESH/THRESHD
c
      WTHRESH = 0.0
c
      IF(W.GT.WTHRESH)THEN
        THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
        THRESH=0.
      ENDIF
      GAMD = GAMDP
      GAM  = SQRT(GAMD**2+GAMQ**2+GSPRDA**2)
      SIGD = QFDP*(GAMDP/GAM)
      SIGD = SIGD*CMTOT2*GAM**2
      SIGD = SIGD/((CMTOT2-(DM+EPSD)**2)**2+CMTOT2*GAM**2)
      SIGD = SIGD*FD(QMS,AD)**2/FD(QMSR,AD)**2
      TEST = QVS/QVSR
      SIGD = SIGD*TEST
      SIGD = SIGD*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGD = SIGD/(QMSR/2./QVSR+TAN(37.1*PI/180./2.)**2)
      SIGD = SIGD*SIGMOT(E,THR)/SIGMOT(730.D0,37.1*PI/180.)
      SIGD = SIGD*A
      SIGD = SIGD*THRESH
C      if (A.eq.14.0) then
C        SIGDEL = SIGD*( 2.5*E/(4000.0)+ (E-W) )/W
C        SIGDEL = SIGD*(0.1-2.65*(E-W)*1E-3)
C        SIGDEL = SIGD*(2.0+0.75*(E-W)*1E-3)
C      else if (A.eq.3.0) then
C        SIGDEL = SIGD*( 4.*E/(4000.0)+ (E-W) )/W
C        SIGDEL = SIGD*0.4
C      else
C        SIGDEL = SIGD
C      endif
      SIGDEL = SIGD
      RETURN
      END
*SIGQFS
      REAL*8 FUNCTION SIGQFS(E,TH,W,Z,A,EPS,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL Qdep
      COMMON/QDEPENDENCE/XPAR0,XPAR1,Qdep
      PM   = 939.
      UP   = 2.7928456
      UN   = -1.91304184
      AP0  = 840.
      AP1  = 750.
      ALPH = 1./137.03604
      HBARC= 197.32858
      PI   = ACOS(-1.)
      GAMR = 120.
      PFR  = 230.
      QMSRQ= 4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ= QMSRQ+115.**2
      NA   = INT(A)

      IF(NA.EQ.1)THEN
        AP=AP0
      ELSEIF(NA.LT.4)THEN
        AP=AP0+(A-1.)*(AP1-AP0)/3.
      ELSE
        AP=AP1
      ENDIF
      THR = TH*PI/180.
      QMS = 4.*E*(E-W)*SIN(THR/2.)**2
      QVS = QMS+W**2

C  START QFS SECTION
      SIGNS  = SIGMOT(E,THR)*RECOIL(E,THR,PM)

      SIGEP = GEP(QMS,AP)**2 + TAU(QMS) * GMP(QMS,AP)**2
      SIGEP = SIGEP/(1.0+TAU(QMS) )
      SIGEP = SIGEP+2.0*TAU(QMS)*GMP(QMS,AP)**2 * (TAN(THR/2.))**2
      SIGEP = SIGNS*SIGEP

      SIGEN = GEN(QMS,AP)**2 + TAU(QMS) * GMN(QMS,AP)**2
      SIGEN = SIGEN/(1.0+TAU(QMS) )
      SIGEN = SIGEN+2.0*TAU(QMS)*GMN(QMS,AP)**2 * (TAN(THR/2.))**2
      SIGEN = SIGNS*SIGEN

      EPQ    = 4.*E**2*SIN(THR/2.)**2/2./PM
      EPQ    = EPQ/(1.+2.*E*SIN(THR/2.)**2/PM)+EPS
      EPQ    = E-EPQ

CDEBUG   QSEPQ=4.*E*EPQ*SIN(THR/2.)**2  ! Q^2 at Q.E. Peak

      IF(INT(A).EQ.1)THEN
        ARG = (E-W-EPQ)/SQRT(2.)/1.
        DEN = 2.51
      ELSE
        GAMQ = GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
        ARG  = (E-W-EPQ)/1.20/(GAMQ/2.)
        DEN  = 2.13*(GAMQ/2.)
      ENDIF
      NQ = INT(ARG)
      IF(ABS(NQ).GT.10)THEN
        SIGQ = 0.
      ELSE
        SIGQ = (Z*SIGEP+(A-Z)*SIGEN)*EXP(-ARG**2)/DEN
      ENDIF
      SIGQFS = SIGQ 
C      if (A.eq.14.0) then
C        SIGQFS=SIGQ*E*(E/4000.0)*(1/1178.0)
C      else if (A.eq.3.0) then
CC        SIGQFS=SIGQ*(E/W)*(E/4000.0)
CC        SIGQFS=SIGQ*(SIN((32.0/2.0)*(PI/180.0)))/SIN(THR/2.0)
C         SIGQFS=SIGQ*1.2
C      else
C        SIGQFS=SIGQ
C      endif
CDEBUG Q2 DEPENDENCE 
      if (Qdep) then
        QQ    = QMS/1.E6
        XCOR  = XPAR1/QQ + XPAR0   
        SIGQFS= SIGQFS/XCOR
      endif

      RETURN
      END

cAdd acc funcs
      REAL*8 FUNCTION TAU(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PM   = 939.
      TAU  = QMS/4.0/PM**2
      RETURN
      END

      REAL*8 FUNCTION GEP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      GEP=1./(1.+QMS/AP**2)**2
      RETURN
      END

      REAL*8 FUNCTION GEN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PM   = 939.
      UN   = -1.91304184

      GEN = -UN
      GEN = GEN * TAU(QMS)/( 1.0+5.6*TAU(QMS) )
      GEN = GEN * GEP(QMS,AP)
      RETURN
      END

      REAL*8 FUNCTION GMP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      UP   =  2.7928456
      GMP  =  UP * GEP(QMS,AP)
      RETURN
      END

      REAL*8 FUNCTION GMN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      UN   = -1.91304184
      GMN  = UN * GEP(QMS,AP)
      RETURN
      END
