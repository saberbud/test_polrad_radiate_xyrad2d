*     RADIATE
      SUBROUTINE RADIATE(E,TH,W,SIGNR,SIGRAD,xF,xstype)
C     DOES NOT INCLUDE CONTRIBUTION FROM ELASTIC SCATTERING
C
C-----K. Slifer. 09/16/02
C
C     Rewrote subroutine to include external bremsstrahlung using
C     formalism from S. Stein et al. Phys. Rev. D 12 7. Equation (A82)
C     Where possible the equation number is given with each expression.
C----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'const.inp'
      INCLUDE 'data_common_const.inp'

      REAL*8 SIGNR,SIGRAD,DEL,PREC,W,W0,Z,A,SPENCE
      REAL*8 type_cxsn
      REAL*8 SCALE,UNSCALE,Tb,Ta
      REAL*8 PI,THR,XMT,ARG!,ALPH,EMASS,xb,XM
      REAL*8 QMS,D1,tr,D2,xF,Tpb,Tpa
      REAL*8 eta,xi,R,TERM1,TERM2,TERM3,TERM4
      REAL*8 X1LO,X1HI,X2LO,X2HI,ANS_Es,ANS_Ep
      REAL E,E0,TH,TH0
      REAL xyE,xyTH
      REAL*8 xyW
      REAL*8 E0gev,EPgev
      integer NSP,KF,trial,xstype
      LOGICAL extrad
      COMMON/PAR/E0,TH0,W0,Z,A,SPENCE
!      COMMON/PAR/E,TH,W,Z,A,SPENCE
      COMMON/xy/xyE,xyTH,xyW
      COMMON/ADDONS/SCALE,UNSCALE,Tb,Ta,extrad,trial

c      Z=2
c      A=3

c      Z=1
c      A=1

c      A=Aset
c      Z=Zset

      A=Asir
      Z=Zsir

      Ta=0.05949
      Tb=0.00263

      SCALE=1.0

c      extrad=.true.
      extrad=.false.

c$$$     write(*,*)'SCALE,UNSCALE,Tb,Ta,extrad,trial',SCALE,UNSCALE,Tb,Ta,extrad,trial
c      write(*,*)'Z,A',Z,A

      DEL   = 10
      PREC  = .001  ! 0.001 gets rid of glitches
      !WRITE(6,*) "IN RADIATE" 
      if (.NOT.extrad) THEN  ! don't apply external rad. cor.
        Tb=0.0
        Ta=0.0
      endif
c      write(*,*)'Z,A,Ta,Tb',Z,A,Ta,Tb

      E0=E
      TH0=TH
      W0=W

c      write(*,*)'E,TH,W,E0,TH0,W0',E,TH,W,E0,TH0,W0

      xyE=E
      xyTH=TH
      xyW=W

      E0gev=E/1000.
      EPgev=(E-W)/1000.

      SIGNR=type_cxsn(E0gev,EPgev,TH,xstype)

c      write(*,*)'xyE,xyTH,xyW',xyE,xyTH,xyW

      XMT   = A*XM   ! Mass of the target
      PI    = ACOS(-1.)
      THR   = TH    ! angle in radians
      ARG   = COS(THR/2.)**2

      SPENCE= PI**2/6.-LOG(ARG)*LOG(1.-ARG)
      DO 10 NSP=1,50
 10     SPENCE = SPENCE-ARG**NSP/FLOAT(NSP)**2

      QMS   = 4.*E*(E-W)*SIN(THR/2.)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      tr=D1/2./xb

c      write(*,*)'tr=',tr

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-17./36. ! this term dominates D2
      D2 = D2 +0.5*(PI**2/6.-SPENCE)
      D2 = D2 -1./4.*( LOG( E/(E-W) ) )**2        ! Correct. to peak. appr.
      D2 = D2*(2.*ALPH/PI)                 
      D2 = D2+0.5772*xb*(Tb+Ta)                   ! Here D2= F-1
      xF = (1.+D2)                                ! (A44)

      Tpb = tr + Tb
      Tpa = tr + Ta  
   
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) ) ! (A83)
      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )    ! (A46)
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )              ! (A52)

      SIGRAD = SIGNR * xF
      SIGRAD = SIGRAD*( (R*DEL/E  )**(xb*Tpb) )
      SIGRAD = SIGRAD*( (DEL/(E-W))**(xb*Tpa) )
      SIGRAD = SIGRAD*(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )

c      write(*,*)'SIGNR=',SIGNR
c      write(*,*)'QMS',QMS
c      write(*,*)'xF=',xF
c      write(*,*)'E,W,QMS',E,W,QMS
c      write(*,*)'W=',W
c      write(*,*)'R=',R
c      write(*,*)'DEL=',DEL
c      write(*,*)"(R*DEL/E  )**(xb*Tpb)=",(R*DEL/E  )**(xb*Tpb)
c      write(*,*)"(DEL/(E-W))**(xb*Tpa)=",(DEL/(E-W))**(xb*Tpa)
c      write(*,*)"(1.-xi/DEL/( 1.-xb*(Tpb+Tpa)))=",
c     %     (1.-xi/DEL/( 1.-xb*(Tpb+Tpa)))

c$$$      write(*,*)'SIGRAD before integral=',SIGRAD

      TERM1=(R*DEL/E  )**(xb*Tpb)
      TERM2=(DEL/(E-W))**(xb*Tpa)
      TERM3=(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )
      TERM4=xF
C
C-----Stein's 1st integral wrt dEs' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X1LO   = (E-W)*( XMT/( XMT-2.*(E-W)*(SIN(THR/2.))**2) -1.0 )
      X1HI   = W-R*DEL

c      if (mod(trial,1000).lt.10) THEN
c         write(6,*) "dEs integral: ",E,THR,W,x1lo,x1hi,R,XMT,trial
c      endif

c      write(*,*)'XMT',XMT
c      write(*,*)'X1LO,X1HI',X1LO,X1HI

      ANS_Es = 0.
      IF (X1HI.GT.X1LO) THEN
        CALL ROM(X1LO,X1HI,PREC,ANS_Es,KF,1,trial,xstype)
      ELSE
C        write(6,*) "Integral dEs. SKIPPING:nu,lower,higher ",W,X1LO,X1HI
      ENDIF

c      write(*,*)'SCALE=',SCALE

      ANS_Es = ANS_Es * SCALE
C
C-----Stein's 2nd integral wrt dEp' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X2LO   = E*( 1.0-XMT/( XMT+2.*E*(SIN(THR/2.))**2) )
      X2HI   = W-DEL 
c      if (trial.EQ.23279) THEN
c         write(6,*) "dEp integral: ",E,THR,W,x2lo,x2hi,XMT,trial
c      endif 

c      write(*,*)'X2LO,X2HI',X2LO,X2HI

      ANS_Ep = 0.
      IF (X2HI.GT.X2LO) THEN
         CALL ROM(X2LO,X2HI,PREC,ANS_Ep,KF,2,trial,xstype)
      ELSE
C         write(6,*) "Integral dEp. SKIPPING:nu,lower,higher ",W,X2LO,X2HI 
      ENDIF
      ANS_Ep = ANS_Ep * SCALE

CCDEBUG
C      WRITE(6,'(F7.1,3F8.2)') W,(SIGRAD-SIGNR)/SIGNR*100., 
C     &                          ANS_Es/SIGNR*100.,
C     &                          ANS_Ep/SIGNR*100.
C

c      write(*,*)'ANS_Es=',ANS_Es
c      write(*,*)'ANS_Ep=',ANS_Ep
c      write(*,*)'SIGRAD before sum=',SIGRAD

      SIGRAD = SIGRAD + ANS_Es + ANS_Ep

c      write(*,*)'SIGRAD final=',SIGRAD

c$$$      write(*,*)'SIGRAD after integral=',SIGRAD

CDEBUG
C      WRITE(6,'(F7.1,6F8.2)') 
C     &         W,TERM1,TERM2,TERM3,TERM4,
C     &         ANS_Es/SIGRAD*100.0,ANS_Ep/SIGRAD*100.0

   
      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


*VALY
      SUBROUTINE VALY(X,F,IFUNC,xstype)
      implicit none

      INCLUDE 'const.inp'
      
      REAL*8 X,F,W,Z,A,SPENCE
      REAL*8 SCALE,UNSCALE,Tb,Ta
      REAL*8 PI,THR,XMT
      REAL*8 eta,xi,R
      REAL*8 QMS2,tr2,Tpb,Tpa,D2,SIG2,F2
      REAL*8 QMS1,tr1,SIG1,F1
      real*8 inelas_cxsn,qfs_phi,type_cxsn
      REAL E,TH,Eps,dEs
      REAL*8 xyW
      REAL xyE,xyTH
      Integer IFUNC,trial,xstype
      LOGICAL extrad
C     REAL*8 FUNCTION F(X)
      COMMON/PAR/E,TH,W,Z,A,SPENCE
      COMMON/xy/xyE,xyTH,xyW
      COMMON/ADDONS/SCALE,UNSCALE,Tb,Ta,extrad,trial

c      write(*,*)'Z,A=',Z,A
c      write(*,*)'Ta,Tb=',Ta,Tb
c$$$
c$$$      Z=2
c$$$      A=3

c      write(*,*)'VALY:xyE,xyTH,xyW',xyE,xyTH,xyW
c      write(*,*)'VALY:E,TH,W,Z,A,X,EMASS=',E,TH,W,Z,A,X,EMASS
c      stop

      if (.NOT.extrad) THEN  ! apply external rad. cor.
        Tb=0.0
        Ta=0.0
      endif
      PI    = ACOS(-1.)
      THR   = TH
      XMT   = A*XM   ! Mass of the target

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) )
C
C-----Stein's 2nd integral dEp'
C
      QMS2  = 4.*E*(E-X)*SIN(THR/2.)**2  ! 1/15/03
      tr2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)
      Tpb   = tr2 + Tb
      Tpa   = tr2 + Ta

      D2    = 13.*(LOG(QMS2/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !KS. Correction to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta)
      
      Eps = (E-X)/1000.
c      SIG2 = inelas_cxsn(E/1000.,Eps,TH) ! in nb/MeV/Sr
      SIG2 = type_cxsn(E/1000.,Eps,TH,xstype)
      IF (SIG2.lt.0) THEN
c         print*, 'Error!  dEp-Cross section is negative or NAN',
c     +        E/1000.,Eps,TH,SIG2
c         STOP
      ENDIF
c      SIG2 = SIG2*UNSCALE
c      SIG2  = SIGQFS(E,TH,X,Z,A,EPS,PF)
      F2    = ( xb*Tpa/(W-X) ) *qfs_phi((W-X)/(E-X))
      F2    = F2 + xi/(2.*(W-X)**2)
      F2    = F2 * SIG2*(1.+D2)
      F2    = F2 * ( (W-X)/(E-X) )**(xb*Tpa)
      F2    = F2 * ( (W-X)*R/(E) )**(xb*Tpb)
C
C-----Stein's 1st integral dEs'
C
      QMS1  = 4.*(E-W+X)*(E-W)*SIN(THR/2.)**2    !    1/15/03
      tr1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)
      Tpb   = tr1 + Tb
      Tpa   = tr1 + Ta

      D2    = 13.*(LOG(QMS1/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !Corr. to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta) ! 1/14/02

      dEs = (E-W+X)/1000.
      Eps = (E-W)/1000.

c      SIG1 = inelas_cxsn(dEs,Eps,TH) ! in nb/MeV/Sr
      SIG1 = type_cxsn(dEs,Eps,TH,xstype)
c$$$      write(*,*)'SIG1=',SIG1
      IF (SIG1.lt.0) THEN
c         print*, 'Error!  dEs-Cross section is negative or NAN',
c     +        dEs,Eps,TH,SIG2
c         STOP
      ENDIF
c      if ((dEs-Eps)*1000.lt.-1.0) Then
c         write(6,*) E,Eps,X,W,dEs,(dEs-Eps)*1000.
c      endif
c      SIG1  = SIGQFS(E-W+X,TH,X,Z,A,EPS,PF)
      F1    = ( xb*Tpb/(W-X) ) *qfs_phi((W-X)/(E))   ! 
c$$$      write(*,*)'F1=',F1
      F1    = F1 + xi/(2.*(W-X)**2)
c$$$      write(*,*)'F1=',F1
      F1    = F1 * SIG1*(1.+D2)
c$$$      write(*,*)'F1=',F1
c$$$      write(*,*)'( (W-X)/((E-W)*R) )=',( (W-X)/((E-W)*R) )
c$$$      write(*,*)'( (W-X)/ (E)      )=',( (W-X)/ (E)      )
c$$$      write(*,*)'(xb*Tpa),(xb*Tpb)=',(xb*Tpa),(xb*Tpb)
c$$$      write(*,*)'fac1=',( (W-X)/((E-W)*R) )**(xb*Tpa)
c$$$      write(*,*)'fac2=',( (W-X)/ (E)      )**(xb*Tpb)
      F1    = F1 * ( (W-X)/((E-W)*R) )**(xb*Tpa)
      F1    = F1 * ( (W-X)/ (E)      )**(xb*Tpb) ! 
c$$$      write(*,*)'F1=',F1

c$$$      write(*,*)'IFUNC=',IFUNC
      IF(IFUNC.EQ.2) THEN      ! dEp'
        F=F2
c        write(6,*) "dEp = ",E/1000.,X,Eps,inelas_cxsn(E/1000.,Eps,TH)
c        write(6,*) E,TH,X,sig2
      ELSEIF (IFUNC.EQ.1) THEN ! dEs'
        F=F1
c        if (trial.EQ.23279) THEN
c           write(6,*) "cxsn dEs = ",dEs,Eps,E,W,X,sig1
c        endif
c        write(6,*) E,TH,X,sig1 
      ELSE
        write(6,*) "PROBLEM. "
        STOP
      ENDIF
      RETURN
      END

      SUBROUTINE ROM(A,B,EPS,ANS,K,IFUNC,trial,xstype)
      implicit none
      REAL*8 A,B,EPS,ANS,W,H,FA,FB,SIG
      REAL*8 X,F,E
      integer IFUNC,IU,IV,J,J1,K,L,M,trial,xstype
C  ROMBERG METHOD OF INTEGRATION
      DIMENSION W(50,50)
      H=B-A
      K=0
      CALL VALY(A,FA,IFUNC,xstype)
      CALL VALY(B,FB,IFUNC,xstype)
      W(1,1)=(FA+FB)*H/2.

c      write(*,*)'ROM:H,FA,FB=',H,FA,FB

    4 K=K+1
c$$$      write(*,*)'K=',K
      IF(K.GE.20)GO TO 5
      IF(K.GE.17)GO TO 3  !Force quit with result
      H=H/2.
      SIG=0.
      M=2**(K-1)
C      print*,K,M,H,W(1,1)
c      if (K.GT.9) THEN
c         print*,K,M,H,W(1,1),A,B,trial,IFUNC
c      ENDIF
      DO 1 J=1,M
      J1=2*J-1
      X=A+FLOAT(J1)*H
c      print*, J,J1,X
      CALL VALY(X,F,IFUNC,xstype)
c      write(6,*) "DEBUG: ",k,IFUNC,f
    1 SIG=SIG+F
c$$$      write(*,*)'F=',F
      W(K+1,1)=W(K,1)/2.+H*SIG
      DO 2 L=1,K
         IU=K+1-L
         IV=L+1
    2 W(IU,IV)=(4.**(IV-1)*W(IU+1,IV-1)-W(IU,IV-1))/(4.**(IV-1)-1.)
      If(W(IU,IV).GT.0.) Then 
         E=(W(IU,IV)-W(IU,IV-1))/W(IU,IV)
      else 
         E =1.0
      endif 
!      if (trial.EQ.23279) THEN
c      if (K.GT.9) THEN
c         write(6,*) W(IU,IV),W(IU,IV-1),(ABS(E)-EPS),IU,IV,K
c      ENDIF
      IF(ABS(E)-EPS) 3,3,4
    3 ANS=W(1,IV)
c$$$      write(*,*)'ANS=',ANS
      RETURN
    5 PRINT 100
  100 FORMAT(' K OVERFLOW')
      CALL EXIT(0)
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL*8 function qfs_phi(x)
      implicit none
      REAL*8 x
      qfs_phi=1.0-x+3./4.*x**2
      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE TARGETCELL(TAG,Tb,Ta)
      implicit none
      REAL*8 Tb,Ta
      CHARACTER TAG*6

C-----DETERMINE TARGET CELL CHARACTERISTICS
C     This info is specific to E97110

      if (TAG.EQ.'PEN_6D') then  ! Target is "Penelope" 6 degs.
         Tb     = 2.116E-3      ! Thickness before scatt., including 1/2 target
         Ta     = 6.492E-2      ! Thickness after scatt., including targ contibution
      elseif(TAG.EQ.'REF16D') then ! N2 Refrence Cell #1 6 degs.
         Tb     = 7.760E-3
         Ta     = 6.883E-2     
      elseif(TAG.EQ.'PRP_6D') then ! Priapus 6 degs.
         Tb     = 2.071E-3
         Ta     = 6.293E-2
      elseif(TAG.EQ.'REF26D') then ! N2 Refrence Cell #2 6 degs.
         Tb     = 7.811E-3
         Ta     = 6.632E-2     
      elseif(TAG.EQ.'PRP_9D') then ! Priapus 9 degs.
         Tb     = 2.073E-3
         Ta     = 4.492E-2    
      elseif(TAG.EQ.'REF29D') then ! N2 Refrence Cell #2 9 degs.
         Tb     = 7.733E-3
         Ta     = 4.716E-2
      elseif(TAG.EQ.'C12_9D') then ! Carbon 9 degs.
         Tb     = 1.119E-3
         Ta     = 9.287E-3
      else
        WRITE(6,*) TAG,':  Unknown Energy'
        STOP
      endif
C      density=density*2.6868E19 ! Convert from Amagat to (cm)^-3
C      write(8,'(A,F6.1,x,2(F9.5,x),E10.5,x,F5.1)') 
C     &         '#',E,    Tb,Ta,    density,zlength 
      RETURN 
      END
 
      SUBROUTINE GETSCALE(units,SCALE,UNSCALE)
      implicit none
      REAL*8 SCALE,UNSCALE
      CHARACTER units*2

C-----DETERMINE DESIRED UNITS OF CROSS SECTION
      if(units.EQ.'cm') then     ! cm^2/(MeV sr)
        SCALE=1.D-33
        UNSCALE=1.0
      elseif(units.EQ.'mb')then  ! mb/(MeV sr)
        SCALE=1.D-06
        UNSCALE=1.D-27
      elseif(units.EQ.'ub')then  ! ub/(MeV sr)
        SCALE=1.D-03
        UNSCALE=1.D-30
      elseif(units.EQ.'nb')then  ! nb/(MeV sr)
        SCALE=1.0
        UNSCALE=1.D-33
      elseif(units.EQ.'pb') then ! pb/(MeV sr)
        SCALE=1.D+03
        UNSCALE=1.D-36
      else
        WRITE(6,*) 'Unknown scale'
        STOP
      endif
C      write(8,'(A,A)') '#',units
      RETURN
      END

      SUBROUTINE GETZA(nucleus,Z,A)
      implicit none
      integer Z,A
      CHARACTER nucleus*3 
C-----DETERMINE nucleus charge, mass
      if (nucleus.EQ.'He3') THEN ! He-3
        Z=2
        A=3
      elseif (nucleus.EQ.'C12') THEN ! Carbon
        Z=6
        A=12
      elseif (nucleus.EQ.'Nit') THEN ! Nitrogen
        Z=7
        A=14
      elseif (nucleus.EQ.'Sil') THEN ! Silicon
        Z=14
        A=28
      elseif (nucleus.EQ.'Irn') THEN ! Iron
        Z=26
        A=56
      ELSE
        write(6,*) 'Problem'
        stop
      endif
C------
      RETURN
      END
