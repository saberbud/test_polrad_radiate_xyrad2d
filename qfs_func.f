C#####################################################################
*FD
      REAL*8 FUNCTION FD(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FD=1./(1.+QMS/A**2)**2
      RETURN
      END
*FM
      REAL*8 FUNCTION FM(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FM=1./(1.+QMS/A**2)
      RETURN
      END
*FPHENOM
      REAL*8 FUNCTION FPHENOM(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      A1=.55
      A2=20./1.E6
      B1=.45
      B2=.45/1.E6
      C1=0.03
      C2=0.2/1.E12
      FPHENOM=A1*EXP(-A2*QMS)+B1*EXP(-B2*QMS)
      FPHENOM=FPHENOM+C1*EXP(-C2*(QMS-4.5E6)**2)
      FPHENOM=SQRT(FPHENOM)
      RETURN
      END
*FYUKAWA
      REAL*8 FUNCTION FYUKAWA(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(QMS.LT.1.E-5.OR.A.LT.1.E-5)THEN
      FYUKAWA=0.
      ELSE
      ARG=SQRT(QMS/2.)/A
      FYUKAWA=ATAN(ARG)/ARG
      ENDIF
      RETURN
      END
*SIGMOT
      REAL*8 FUNCTION SIGMOT(E,THR)
      IMPLICIT REAL*8 (A-H,O-Z)
      ALPH=1./137.03604
      HBARC=197.3286
      SIGMOT=(ALPH*HBARC*COS(THR/2.)/2./E/SIN(THR/2.)**2)**2
C  FM**2/SR
      RETURN
      END
*RECOIL
      REAL*8 FUNCTION RECOIL(E,THR,TM)
      IMPLICIT REAL*8 (A-H,O-Z)
      RECOIL=1./(1.+2.*E*SIN(THR/2.)**2/TM)
      RETURN
      END
