! POLSIG Start -----------------------------------------------------------------
!     Calculation of the cross section
!
!     INPUT
!       E0 = Incident electron energy in MeV
!       EP = Scattered electron energy in MeV
!       TH_RAD = Scattering angle in radians
!
!
!     Kinematic region
!
!       IELAS_IN = 0 : quasi-elastic and inelastic
!       IELAS_IN = 1 : elastic
!
!     Polarization specification
!
!       IPOL = 0 : parallel, unpolarized
!       IPOL = 1 : parallel, polarized
!       IPOL = 2 : perpendicular, unpolarized
!       IPOL = 3 : perpendicular, polarized
!
!     The output SIGMA_BORN and SIGMA_RAD are in nb/MeV.Sr
!
!---- For XJACOB:
!---- Factor 1.0D-3 converts nb/GeV.Sr to nb/MeV.Sr for inelastic and quasi-
!---- elastic cross section and elastic tail. For elastic cross section, it
!---- simply converts nb to ub. And also remember that elastic cross section
!---- will lack recoil factor (E'/E)
!
!----------------------------------------------------------------------
      SUBROUTINE POLSIG(E0,EP,TH_RAD,IELAS_IN,IPOL,SIGMA_BORN,SIGMA_RAD,xfp,xstype)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL E0,EP,TH_RAD,SIGMA_BORN,SIGMA_RAD,xfp
      INTEGER IFLAG,ixypol
      integer xstype,xyxstype
      common/xytype/xyxstype
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     &     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     &     sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich

      COMMON/CHOI/W2_ELAS
      COMMON/SOFTPHOTON/EXTAI2
      COMMON/ELASTIC/IELAS
      COMMON/XY/ixypol
!
!----------------------------------------------------------------------
!
      IELAS = IELAS_IN  ! CMASS
      ixypol=IPOL

      Q2    = 4.0*E0*EP*SIN(0.5*TH_RAD)**2*1.0E-6
      XNU   = (E0-EP)*1.0E-3

      xyxstype=xstype

c      write(*,*)'Q2,IPOL=',Q2,IPOL

      IF (IELAS.EQ.0) THEN
         W2_N    = AMH*(AMH + 2.0*XNU) - Q2   ! just for reference
         W_N     = SQRT(W2_N)*1000.0
         W2      = AMT*(AMT + 2.0*XNU) - Q2
         XS      = Q2/(2.0*AMT*XNU)
         YS      = XNU/(E0*1.0E-3)

         XJACOB  = 2.0D0*PI*AMH*YS/(1.0-YS)   ! XJacob = 1E-3 * (1-YS)/(2PI*Mh*YS)
         XJACOB  = 1.0D-3/XJACOB              !        = (1-YS)/2PI*Mh*YS) [1/MeV]

         Y_ELAS  = 1.0/(1.0 + AMT/(2.0*E0*1.0E-3*SIN(0.5*TH_RAD)**2))

         W2_ELAS = AMT**2                     ! CDEBUG Add breakup

      ELSE IF (IELAS.EQ.1) THEN

         W2      = AMT*(AMT + 2.0*XNU) - Q2
         XS      = Q2/(2.0*AMT*XNU)
         YS      = XNU/(E0*1.0E-3)
         XJACOB  = 2.0D0*PI*AMH*YS/(1.0-YS)
         XJACOB  = 1.0D-3/XJACOB 
         Y_ELAS  = 1.0/(1.0 + AMT/(2.0*E0*1.0E-3*SIN(0.5*TH_RAD)**2))
         W2_ELAS = AMT**2       ! use He3 mass for W2
      ENDIF

      IF ( (W2.LE.0.0).OR.(YS.LT.Y_ELAS-1.0E-02) ) THEN
         WRITE(6,*) "PROB IN POLSIG.  W<0 or Y < Y_ELAS"
         write(*,*)'W2,YS,(Y_ELAS-1.0E-02)',W2,YS,(Y_ELAS-1.0E-02)
         SIGMA_BORN = 0.0
         SIGMA_RAD  = 0.0
         RETURN
      ENDIF

      SNUC = 2.0*AMH*SQRT((E0*1.0E-3)**2+AML2)       ! 2*Mh*E0

c      write(*,*)'xy1'

      IF (IELAS.EQ.0) THEN
        Y   = SNUC*XS*YS*(AMT/AMH)
        call conkin(snuc,amT,IPOL)
      ELSE IF (IELAS.EQ.1) THEN
        Y   = SNUC*XS*YS*(AMT/AMH)
        CALL CONKIN(SNUC,AMT,IPOL)
      ENDIF

c      write(*,*)'xy2'   

!--------------------------------------------------------------------------
!-----Delta is factorizing part of virtual and real leptonic bremsstrahlung
!--------------------------------------------------------------------------
      call deltas(tr,factor1,del_sub)

c      write(*,*)'xy3'

      IF (IPOL.EQ.0.OR.IPOL.EQ.2) THEN  ! Cross section for unpolarized hadron target   
         UN   = 1.0                     ! un = 1 calculates F1 and F2 struct funct    
         PL   = 1.0                     ! Lepton polarization                       
         PN   = 0.0                     ! Unpolarized target, g1 = g2 = 0        
         QN   = 0.0                     ! Qn = 0,the hadronic tensor for spin 1/2 part.
         ISF1 = 1
         ISF2 = 2
         ISF3 = 1
      ELSEIF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN ! Diff. between 2 hadron polar. directions 
         UN   = 0.0                        ! un = 0 means F1 = F2 = 0         
         PL   = 1.0                        ! Lepton polarization                
         PN   = 1.0                        ! Pn defines hadron pol. g1&g2 non-zero
         QN   = 0.0
         ISF1 = 3
         ISF2 = 4
         ISF3 = 1
      ENDIF
           
      EXTAI2     = 1.0D0
      
      CALL BORNIN(SIB)                                  ! SIB is dsig/(dxdy) in nb

      write(*,*)'Polsig:SIGMA_BORN',(SIB*XJACOB)
      stop

c      write(*,*)'Born xs finish'

      EXTAI2     = ((SX-Y/TARA)**2/S/(S-Y/TARA))**TR    !used only in elastic calc.
      
!DEBUG
!      IF ( (W_N.GT.1232.0 .AND. W_N.LT.1233.0 ).AND.(E0.EQ.862.0) )THEN
!       IF(IPOL.EQ.1) THEN
!        WRITE(68,'(A,2F10.2)') "#",E0,W_N
!        WRITE(68,'(A)')        "#, TAU  R"
!        WRITE(69,'(A,2F10.2)') "#",E0,W_N
!        WRITE(69,'(A)')        "#, Taupass, RCUT(2)"
!        TAIL       = AN*ALFA/PI*TAIL_INTEG(TAMIN,TAMAX)  
!       ENDIF
!      ENDIF
!DEBUG

c      write(*,*)'xy4'
      TAIL       = AN*ALFA/PI*TAIL_INTEG(TAMIN,TAMAX)   
c      write(*,*)'xy5'
      SIGMA_BORN = SIB*XJACOB                           ! nb/MeV-sr

      write(*,*)'Polsig:SIGMA_BORN',SIGMA_BORN
      stop

      SIGMA_RAD  = ( SIB*FACTOR1*(1.+ALFA/PI*DEL_SUB) + TAIL )*XJACOB

      SIGMA_RADnt=( SIB*FACTOR1*(1.+ALFA/PI*DEL_SUB))*XJACOB

      IF (IPOL.EQ.1.OR.IPOL.EQ.3) THEN
        SIGMA_BORN = -SIGMA_BORN
        SIGMA_RAD  = -SIGMA_RAD
      ENDIF

      xfp=SIGMA_RADnt/SIGMA_BORN

c      write(*,*)'SIB,FACTOR1,DEL_SUB,TAIL,XJACOB',SIB,FACTOR1,DEL_SUB,TAIL,XJACOB
c      write(*,*)'SIGMA_BORN,SIGMA_RAD,SIGMA_RADnt',SIGMA_BORN,SIGMA_RAD,SIGMA_RADnt
c      write(*,*)'SIGMA_RADnt/SIGMA_BORN=',SIGMA_RADnt/SIGMA_BORN

      RETURN
      end

!$$$  EXTAI1     = EXP(ALFA/PI*DELINF)
!$$$  EXTAI2     = 1.0                                ! to check proton elastic tail
!$$$  SIGMA_RAD  = (SIB*EXTAI1*(1.+ALFA/PI*(DELTA-DELINF))+TAIL)*XJACOB
!SL   EXTAI3 = 1.0D0
!SL   EXTAI3 = ((SX-Y)**2/S/(S-Y))**TR                ! not used
! POLSIG End ###################################################################

! CONKIN Start -----------------------------------------------------------------
!     set of kinematical constants
      subroutine conkin(snuc,amtar,IPOL)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
!
!     Polarization specification
!
!       IPOL = 0 : parallel, unpolarized
!       IPOL = 1 : parallel, polarized
!       IPOL = 2 : perpendicular, unpolarized
!       IPOL = 3 : perpendicular, polarized
!

      amp  = amtar        !        M_target
      ap   = 2.*amp       !      2*M_target
      amp2 = amp**2       !       (M_target)**2
      ap2  = 2.*amp**2    !  2 .0*(M_target)**2
      s    = snuc*amp/amh ! s = S*(M_target/M_hadron) = 2*M_target*E0
      x    = s*(1.-ys)    ! X = (1-y)S                = 2*M_target*Ep
      sx   = s-x          ! Sx=                       = 2*M_target*nu
      sxp  = s+x          ! Sp=                       = 2*M_target*(E0+Ep)
      ym   = y+al2        ! Q^2 + 2m^2 = Q_m^2        ~ Q^2

      tpl  = s**2+x**2    ! 4*M_target**2(E0**2 + Ep**2)
      tmi  = s**2-x**2    ! 4*M_target**2(E0**2 - Ep**2)

                                              ! w2 is already defined in main program
!SL   w2   = amp2+s-y-x                       ! use proton mass for W2
!$$$  W2   = AMT*(AMT + (S-X)/AMH) - Y        ! use He3 mass for W2

      als  = s*s-al2*ap2                      ! Lambda_s = S^2 - 2m^2*M^2 ~ S^2 
      alx  = x*x-al2*ap2                      ! Lambda_x = X^2 - 2m^2*M^2 ~ X^2
      alm  = y*y+4.*aml2*y                    ! Lambda_y = Y^2 + 4m^2*Y   ~ Q2^2
      aly  = sx**2+4.*amp2*y                  ! Lambda_q = Sx**2 + 4M_target^2 * Q2
                                              !          = 4*M_target^2 * |qvec|^2

      sqls = dsqrt(als)                       ! sqrt(Lambda_s) ~ S
      sqlx = dsqrt(alx)                       ! sqrt(Lambda_x) ~ X
      sqly = dsqrt(aly)                       ! sqrt(Lambda_q) = 2*M_target * |qvec| 
      sqlm = dsqrt(alm)                       ! sqrt(Lambda_y) ~ Q2
 
      allm = dlog((sqlm+y)/(sqlm-y))/sqlm     ! ~ LOG( Q2/m^2 + 1) / Q2
      axy  = pi*(s-x)                         ! 2*M_target*nu
      an   = 2.*alfa**2/sqls*axy*barn*amh/amp ! 2*alph^2*Pi*(nu/E0)(M_p/M_t)*barn

c      write(*,*)'an:alfa,barn',alfa,barn

!     tamin= (sx-sqly)/ap2                    ! equivalent to expression below.
      tamax= (sx+sqly)/ap2                    ! tau_max
      tamin= -y/amp2/tamax                    ! tau_min

      as   = s/2./aml/sqls                    ! ~1/m_e
      bs   = 0.
      cs   = -aml/sqls                        ! -m_e/(2*M_t*E0)

      IF (IPOL/2.EQ.0) THEN                   ! Parallel configuration
         ae  = amp/sqls                        ! M_t/(2*M_t*E0) = 1/(2*E0)
         be  = 0.                              !
         ce  = -s/ap/sqls                      ! ~ 1/(2*M_t)
      ELSE                                    ! Perpendicular configuration
         sqn = dsqrt(s*x*y-aly*aml2-amp2*y*y)  ! ~sqrt(4M_t^2*E0*Ep*Q2 - M_t^2*Q2*Q2)
         ae  = (-s*x+ap2*ym)/sqls/sqn/2.       ! [1/GeV]
         be  = sqls/sqn/2.                     ! [1/GeV]
         ce  = -(s*y+al2*sx)/sqls/sqn/2.       ! [1/GeV]
      ENDIF

      apq   = -y*(ae-be)+ce*sx                                 ! q*eta        [GeV]
      apn   = (y+4.*aml2)*(ae+be)+ce*sxp                       ! (k1+k2)*eta  [GeV]
      dk2ks = as*ym+al2*bs+cs*x                                ! k_2*ksi      [GeV]
      dksp1 = as*s+bs*x+cs*ap2                                 ! ksi*p        [GeV]
      dapks = 2.*(al2*(as*ae+bs*be)+ap2*cs*ce+ym*(as*be+bs*ae)
     .        +s*(as*ce+cs*ae)+x*(bs*ce+cs*be))                ! ksi*eta      [ 1 ]
      return
      end
! CONKIN End ###################################################################

! BORNIN Start -----------------------------------------------------------------
!
!     sibor is born cross section with polarized initial
!     lepton and polarized target
!     siamm is contribution of anomalous magnetic moment.
!
      subroutine bornin(sibor)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/print/ipri1
      dimension sfm0(8),tm(8),SFM(8)
!
!-----First determine reaction region (elastic, quasi-elastic or inelastic)
!
      ipri1 = 1
      call strf(0d0,0d0,sfm0,SFM)
      ipri1 = 0

!-----Only tm(3) and tm(4) used for inelastic.

      tm(1) = -(2.*aml2-y)
      tm(2) = (-(amp2*y-s*x))/(2.*amp2)
      tm(3) = (2.*(apq*dk2ks-dapks*y)*aml)/amp           ! [GeV^2]
      tm(4) = apq/amp*(-(dk2ks*sx-2.*dksp1*y)*aml)/amp2  ! [GeV^2]
      tm(7) = (-(4.*aml2+3.*apn**2-3.*apq**2+y))/2.
      tm(8) = apq/amp*(-3.*(apn*sxp-apq*sx))/(2.*ap)
      ek    = (3.*apq**2-y)/amp2
      tm(5) = -ek*tm(1)
      tm(6) = -ek*tm(2)
      ssum  = 0.

      do 1 isf=isf1,isf2,isf3
c         write(*,*)'BORNIN:isf,sfm0(isf),tm(isf)',isf,sfm0(isf),tm(isf)

         ppol = 1.

         if(isf.eq.3.or.isf.eq.4) ppol = -pn
         if(isf.ge.5)             ppol = qn/6

         ssum = ssum+tm(isf)*sfm0(isf)*ppol               ! [GeV^2]

         write(*,*)'Born:tm(isf),sfm0(isf)',tm(isf),sfm0(isf)

    1 continue

c      write(*,*)'BORNIN:ssum,an,y',ssum,an,y

      sibor = ssum*2.0*an/y**2.                           ! [GeV^2]*[GeV^2 nb]/[GeV^4]      

      write(*,*)'Born:ssum,sibor,ppol',ssum,sibor,ppol

      return                                              ! = [nb]
      end
! BORNIN End ###################################################################

! DELTAS Start -----------------------------------------------------------------
!
!-----delta is factorizing part of virtual and real leptonic bremsstrahlung
!
!     subroutine deltas(delta,delinf,tr,factor1,del_sub)
      subroutine deltas(tr,factor1,del_sub)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      COMMON/CHOI/W2_ELAS
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)

c      write(*,*)'xy2_1'

      del1   = -ym*(alm*allm**2/2.+2.*fspen(2d0*sqlm/(y+sqlm))-pi2/2.)/sqlm
      del2   = (3.*y/2.+4.*aml2)*allm-2.

c      write(*,*)'xy2_2'

      sum    = vacpol(y)

c      write(*,*)'xy2_3'

      aj0    = 2.*(ym*allm-1.)
      deltai = aj0*dlog(DABS(w2-W2_ELAS)/aml/dsqrt(w2))

      ss     = x+y
      xx     = s-y
      alss   = ss**2-2.*w2*al2
      alxx   = xx**2-2.*w2*al2
      sqlss  = dsqrt(alss)
      sqlxx  = dsqrt(alxx)
      allss  = dlog((sqlss+ss)/(-sqlss+ss))/sqlss
      allxx  = dlog((sqlxx+xx)/(-sqlxx+xx))/sqlxx
      dlm    = dlog(y/aml2)
      sfpr   = dlm**2/2.-dlm*dlog(ss*xx/(aml2*w2))
     &          -(dlog(ss/xx))**2/2.+fspen((s*x-y*amp2)/(ss*xx))-pi2/3.
      delta0 = (ss*allss+xx*allxx)/2.+sfpr

      delta   = deltai+delta0+del1+del2+sum
!$$$  delinf  = (dlm-1.)*dlog((w2-amc2)**2/(ss*xx))
      delinf  = (dlm-1.)*dlog((w2-W2_ELAS)**2/(ss*xx))
      tr      = alfa/pi*(dlm-1.)

!$$$  factor1 = ((W2-W2_ELAS)**2/(SS*XX))**(alfa/pi*0.5*aj0)
      factor1 = ((W2-W2_ELAS)**2/(SS*XX))**tr

      DEL_SUB = aj0*dlog(DSQRT(SS*XX)/aml/dsqrt(w2))
     &          + DELTA0 + DEL1 + DEL2 + SUM
      return
      end
! DELTAS End ###################################################################

! VACPOL Start -----------------------------------------------------------------
      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/p/pi,pi2,alfa,i1(8),i2(8)
      dimension am2(3)
!
!     am2 : squared masses of charge leptons
!
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
         a2=2.*am2(i)
         sqlmi=dsqrt(t*t+2.*a2*t)
         allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
        aaa = -1.345d-9
        bbb = -2.302d-3
        ccc = 4.091
      elseif(t.lt.64d0)then
        aaa = -1.512d-3
        bbb =  -2.822d-3
        ccc = 1.218
      else
        aaa = -1.1344d-3
        bbb = -3.0680d-3
        ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alfa

      vacpol=suml+sumh

      end
! VACPOL End ###################################################################

! FSPENS Start -----------------------------------------------------------------
!
!    spence function
!
      double precision function fspens(x)
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end
! FSPENS End ###################################################################

! FSPEN Start ------------------------------------------------------------------
      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end
! FSPEN End ####################################################################

! TAILS Start ------------------------------------------------------------------
       subroutine tails(ta,tm)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
       dimension tm(8,6),ajm2(2),ajm3(3),ii(8)
      data ii/1,2,3,4,1,2,5,6/

      b2=(-aly*ta+sxp*sx*ta+2.*sxp*y)/2. ! B_2
      b1=(-aly*ta-sxp*sx*ta-2.*sxp*y)/2. ! B_1
      c1=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(s*ta+y)**2) ! C_1
      c2=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(ta*x-y)**2) ! C_2
      bb=1./sqly                ! F
      sc1=dsqrt(c1)
      sc2=dsqrt(c2)
      bi12=(sxp*(sx*ta+2.*y))/(sc1*sc2*(sc1+sc2)) ! F_d
      bi1pi2=1./sc2+1./sc1      ! F_1+
      bis=-b1/sc1/c1+b2/sc2/c2  ! F_2+
      bir=b2/sc2/c2+b1/sc1/c1   ! F_2-
      b1i=-b1/aly/sqly          ! F_i
      b11i=(3.*b1**2-aly*c1)/2./aly**2/sqly ! F_ii
      sps=as+bs                 ! s_ksi
      spe=ae+be                 ! s_eta
      ccpe=(ae-be)*ta+2.*ce     ! r_eta
      ccps=(as-bs)*ta+2.*cs     ! r_ksi

      sis=(2.*bi1pi2*sps+bir*sps*ta+bis*ccps)/2. ! F_{2+}^\ksi
      sir=( (2.*bi12*sps*ta+bir*ccps+bis*sps*ta))/2. ! F_{2-}^\ksi
      si12=(bi12*ccps+bi1pi2*sps)/2. ! F_d^\ksi
      eis=(2.*bi1pi2*spe+bir*spe*ta+bis*ccpe)/2. ! F_{2+}^\eta
      eir=( (2.*bi12*spe*ta+bir*ccpe+bis*spe*ta))/2. ! F_{2-}^\eta
      ei12=(bi12*ccpe+bi1pi2*spe)/2. ! F_d^\eta

      ois=((2.*bi1pi2+bir*ta)*(ccpe*sps+ccps*spe)+(ccpe*ccps+
     . spe*sps*ta**2)*bis+8.*bb*spe*sps+4.*bi12*spe*sps*ta**2)/
     .     4.                   ! F_{2+}^{\ksi\eta}
      oir=( ((2.*bi12+bis)*(ccpe*sps+ccps*spe)*ta+(ccpe*ccps+
     . spe*sps*ta**2)*bir+4.*bi1pi2*spe*sps*ta))/4. ! F_{2-}^{\ksi\eta}
      oi12=((ccpe*ccps+spe*sps*ta**2)*bi12+(ccpe*sps+ccps*spe)*
     . bi1pi2+4.*bb*spe*sps)/4. ! F_d^{\ksi\eta}
      eeis=((ccpe**2+spe**2*ta**2)*bis+8.*bb*spe**2+4.*bi12*spe
     . **2*ta**2+4.*bi1pi2*ccpe*spe+2.*bir*ccpe*spe*ta)/4. ! F_{1+}^{\eta\eta}
      eeir=( ((ccpe**2+spe**2*ta**2)*bir+4.*bi12*ccpe*spe*ta+4.
     . *bi1pi2*spe**2*ta+2.*bis*ccpe*spe*ta))/4.
      eei12=((ccpe**2+spe**2*ta**2)*bi12+4.*bb*spe**2+2.*bi1pi2
     . *ccpe*spe)/4.
      ei1pi2=(4.*bb*spe+bi12*spe*ta**2+bi1pi2*ccpe)/2.
      eei1i2=((ccpe**2+spe**2*ta**2)*bi1pi2+4.*(2.*ccpe-spe*ta)
     . *bb*spe+8.*b1i*spe**2+2.*bi12*ccpe*spe*ta**2)/4.
      eb=((ccpe-spe*ta)*bb+2.*b1i*spe)/2.
      eeb=((ccpe-spe*ta)**2*bb+4.*(ccpe-spe*ta)*b1i*spe+4.*b11i
     . *spe**2)/4.
       call ffu(1,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .,eis,eir,ei12,ei1pi2,ta)
       call ffu(2,eb,eis,eir,ei12,ei1pi2,oir,ois,oi12
     .,eeis,eeir,eei12,eei1i2,ta)
       call ffu(3,eeb,eeis,eeir,eei12,eei1i2,0d0,0d0,0d0
     .,0d0,0d0,0d0,0d0,ta)
       ajm2(1)=apq/amp
       ajm2(2)=-1./amp
       ajm3(1)=(y-3.*apq**2)/amp2
       ajm3(2)=6.*apq/amp2
       ajm3(3)=-3./amp2
       do 15 i=1,8
       do 13 l=1,6
   13  tm(i,l)=0
       do 10 k=1,i2(i)
       ajk=1.
       if(i.eq.4.or.i.eq.8)ajk=ajm2(k)
       if(i.eq.5.or.i.eq.6)ajk=ajm3(k)
       do 10 j=k,i1(i)+k-1
       tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,k)*ajk
       if((i.eq.5.or.i.eq.6).and.k.eq.2)
     . tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,1)*ta/amp2
  10   continue
  15   continue
       return
       end
! TAILS End ####################################################################

! FFU Start --------------------------------------------------------------------
       subroutine ffu(n,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .        ,eis,eir,ei12,ei1pi2,ta)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
      hi2=aml2*bis-ym*bi12
      shi2=aml2*sis-ym*si12
      ehi2=aml2*eis-ym*ei12
      ohi2=aml2*ois-ym*oi12
       goto(10,20,30)n
  10   continue
      tm3(3,1,n)=(8.*(apq*dk2ks-dapks*y)*aml*hi2)/amp
      tm3(3,2,n)=(-2.*((2.*(bi12*dk2ks*ta-2.*shi2)*apq+(2.*shi2-
     . sir*y+sis*ym)*apn+4.*dapks*hi2*ta)-4.*((2.*ei12-eis)*
     . dk2ks-(si12-sis)*apn)*aml2)*aml)/amp
      tm3(3,3,n)=(2.*(((2.*si12+sir-sis)*apn*ta-2.*dk2ks*ei12*ta
     . -6.*ohi2-oir*y+ois*ym)-4.*aml2*oi12)*aml)/amp
      tm3(3,4,n)=(2.*(2.*oi12-oir+ois)*aml*ta)/amp
      tm3(5,1,n)=-2.*(4.*aml2+3.*apn**2-3.*apq**2+y)*hi2
      tm3(5,2,n)=-2.*(6.*aml2*apn*eir-3.*apn**2*bi12*ta+3.*apn*
     . apq*bi1pi2+6.*apq*ehi2+hi2*ta)
      tm3(5,3,n)=-(24.*aml2*eei12-6.*apn*ei1pi2-6.*apq*ei12*ta-
     . 2.*bb-bi12*ta**2)
  20   continue
      tm3(4,1,n)=(-4.*(dk2ks*sx-2.*dksp1*y)*aml*hi2)/amp2
      tm3(4,2,n)=(((2.*(sxp-2.*sx)*shi2+2.*bi12*dk2ks*sx*ta+8.*
     . dksp1*hi2*ta-sir*sxp*y+sis*sxp*ym)-4.*(2.*bi12*dk2ks-bis*
     . dk2ks-si12*sxp+sis*sxp)*aml2)*aml)/amp2
      tm3(4,3,n)=((((sxp*ta-ym)*sis-(sxp*ta-y)*sir+2.*bi12*dk2ks
     . *ta+6.*shi2-2.*si12*sxp*ta)+4.*aml2*si12)*aml)/amp2
      tm3(4,4,n)=(-(2.*si12-sir+sis)*aml*ta)/amp2
      tm3(6,1,n)=(-3.*(apn*sxp-apq*sx)*hi2)/amp
      tm3(6,2,n)=(-3.*(2.*(apn*bir+eir*sxp)*aml2-(2.*bi12*sxp*ta
     . -bi1pi2*sx)*apn+(bi1pi2*sxp+2.*hi2)*apq+2.*ehi2*sx))/(2.*
     . amp)
      tm3(6,3,n)=(-3.*(8.*aml2*ei12-apn*bi1pi2-apq*bi12*ta-ei12*
     . sx*ta-ei1pi2*sxp))/(2.*amp)
  30   continue
      tm3(1,1,n)=-4.*(2.*aml2-y)*hi2
      tm3(1,2,n)=4.*hi2*ta
      tm3(1,3,n)=-2.*(2.*bb+bi12*ta**2)
      tm3(2,1,n)=(((sxp**2-sx**2)-4.*amp2*y)*hi2)/(2.*amp2)
      tm3(2,2,n)=(2.*aml2*bir*sxp-4.*amp2*hi2*ta-bi12*sxp**2*ta+
     . bi1pi2*sxp*sx+2.*hi2*sx)/(2.*amp2)
      tm3(2,3,n)=(2.*(2.*bb+bi12*ta**2)*amp2+4.*aml2*bi12-bi12*
     . sx*ta-bi1pi2*sxp)/(2.*amp2)
       return
       end
! FFU End ######################################################################

! STRF Start -------------------------------------------------------------------
!
!     the programm calculates deep inelastic (ita=1),
!     elastic (ita=2), quasielastic (ita=3) structure functions
!     in kinematical point (ta,rr).
!          rr=sx-tt,
!          ta=(t-y)/rr,
!     where tt=t+amf2-amp2, amf2 is invarint mass of final hadrons
!
      subroutine strf(ta,rr,sfm,SFM0)
      implicit real*8(a-h,o-z)
      INCLUDE 'data_common_const.inp'
      integer xstype,xyxstype
      common/xytype/xyxstype
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      COMMON/SOFTPHOTON/EXTAI2
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/print/ipri1
      dimension sfm(8),SFM0(8),FIGI(4,2)
      COMMON/ELASTIC/IELAS
      COMMON/XY/ixypol
      REAL F1_FROM_DATA,F2_FROM_DATA,G1_FROM_DATA,G2_FROM_DATA
      COMMON/DEBUG/IDEBUG1
!
!     Initialize form factors
!
      DO I = 1,8
         SFM(I)  = 0.0D0
         SFM0(I) = 0.0D0
      ENDDO

      xstype=xyxstype
c      write(*,*)'strf:xstype=',xstype
c      stop

c      write(*,*)'ixypol=',ixypol
c      stop

      t    = y+rr*ta              ! same as Q**2 when there is no real photon emitted
      tt   = sx-rr                ! same as 2*M*nu when there is no real photon
      amf2 = tt-t+amp2            ! same as W**2 when there is no real photon
      aks  = t/tt                 ! same as x_bjorken when there is no real photon
      anu  = tt/ap                ! same as nu when there is no real photon

c      write(*,*)'aks,t,tt,y,rr,ta',aks,t,tt,y,rr,ta
c      write(*,*)'amp=',sqrt(amp2),amp,amtar,amt
      
      b1   = 0.d0
      b2   = 0.d0
      b3   = 0.d0
      b4   = 0.d0
!
!-----If t<0, then inelastic (or quasi-elastic) channel is not possible
!
      IF (T.LE.0.0)    RETURN
      IF (AMF2.LT.0.0) RETURN
      IF (IELAS.GT.0)  GOTO 20

      EPSI = AP2/SX                  ! =2M^2/(2M Nu) = M/Nu
      AKS0 = Y/SX                    ! =Q2/2Mnu      = Xbjorken 
      xynu0=SX/AMT/2.0d0
      xymp=0.93827
      xyw20=xymp**2+2*xymp*xynu0-Y

c      write(*,*)'SX,AMT',SX,AMT
c      write(*,*)'AKS0,Y,W2,xynu0,xyw20',AKS0,Y,W2,xynu0,xyw20

c      F1 = F1_FROM_DATA(AKS0,Y)      ! SAME AS F1SFUN(X,Q^2)
c      F2 = F2_FROM_DATA(AKS0,Y)      ! SAME AS F2SFUN(X,Q^2)
c      G1 = G1_FROM_DATA(AKS0,Y)
c      G2 = G2_FROM_DATA(AKS0,Y)

c      F1=1.0
c      F2=1.0
c      G1=1.0
c      G2=1.0

c      Atar=3.0d0   !3He
c      Ztar=2.0d0   !3He
c      Atar=1.0d0   !H
c      Ztar=1.0d0   !H
      Atar=Aset
      Ztar=Zset
c      write(*,*)'Atar,Ztar',Atar,Ztar

      F1_in=0
      F2_in=0
      F1_qe=0
      F2_qe=0

      if(ixypol.eq.0)then
        call F1F2IN09(Ztar,Atar,Y,xyw20,F1_in,F2_in,rtemp)
        call F1F2QE09(Ztar,Atar,Y,xyw20,F1_qe,F2_qe)
      endif

      F1=F1_in+F1_qe
      F2=F2_in+F2_qe

      if(xstype.eq.2)then
        F1=F1_in
        F2=F2_in
      elseif(xstype.eq.1)then
        F1=F1_qe
        F2=F2_qe
      endif

c      write(*,*)'strf:x,Q2',AKS0,Y
c      stop

      xynp=Ztar
      xynn=Atar-Ztar
      if(ixypol.eq.1.or.ixypol.eq.3)then
        call g1g2types(AKS0,Y,xynp,xynn,xstype,G1,G2)
      endif

c      write(*,*)'F1_in,F2_in',F1_in,F2_in
c      write(*,*)'F1,F2,G1,G2',F1,F2,G1,G2

      SFM0(1) = UN*F1+QN/6.*B1
      SFM0(2) = EPSI    * (UN*F2+QN/6.*B2)

      SFM0(3) = EPSI    * (G1+G2)
      SFM0(4) = EPSI**2 * G2

      SFM0(5) = EPSI**2 * B1
      SFM0(6) = EPSI**3 * (B2/3.+B3+B4)
      SFM0(7) = EPSI    * (B2/3.-B3)
      SFM0(8) = EPSI**2 * (B2/3.-B4)

      write(*,*)'strf:G1,G2,EPSI,SFM0(3),SFM0(4)',G1,G2,EPSI,SFM0(3),SFM0(4)

      IF (TA.EQ.0.0D0.AND.RR.EQ.0.0D0) THEN
         DO I = 1,8
            SFM(I) = SFM0(I)
         ENDDO
      ELSE ! Contribution from the inelastic channel
         epsi = ap2/tt            ! same as M/nu when there is no real photon
         xynu=tt/2.0d0/amt
         xyw2=xymp**2+2*xymp*xynu-t

c         write(*,*)'strf:aks,t(Q2),tt,amf2,xyw2,ta,rr',aks,t,tt,amf2,xyw2,ta,rr

c         f1 = F1_FROM_DATA(aks,t) ! same as f1sfun(x,Q^2)
c         f2 = F2_FROM_DATA(aks,t) ! same as f2sfun(x,Q^2)
c         g1 = G1_FROM_DATA(aks,t)
c         g2 = G2_FROM_DATA(aks,t)

         if(ixypol.eq.0)then
           call F1F2IN09(Ztar,Atar,t,xyw2,F1_in,F2_in,rtemp)
           call F1F2QE09(Ztar,Atar,t,xyw2,F1_qe,F2_qe)
         endif

         F1=F1_in+F1_qe
         F2=F2_in+F2_qe

         if(xstype.eq.2)then
           F1=F1_in
           F2=F2_in
         elseif(xstype.eq.1)then
           F1=F1_qe
           F2=F2_qe
         endif

         xynp=Ztar
         xynn=Atar-Ztar
         if(ixypol.eq.1.or.ixypol.eq.3)then
           call g1g2types(AKS,t,xynp,xynn,xstype,G1,G2)
         endif

         IF(g1.ne.g1 .or. g2.ne.g2)then
           write(*,*)'f1,f2,g1,g2',f1,f2,g1,g2
           stop
         endif

c         if(g1.eq.g1 .and. g2.eq.g2)write(*,*)'f1,f2,g1,g2,IPOL',f1,f2,g1,g2,IPOL

c         IF(g1.ne.g1)g1=0.0D0
c         IF(g2.ne.g2)g2=0.0D0

c         f1=1.0
c         f2=1.0
c         g1=1.0
c         g2=1.0

c         write(*,*)'f1,f2,g1,g2',f1,f2,g1,g2

         sfm(1) = un*f1+qn/6.*b1
         sfm(2) = epsi    * (un*f2+qn/6.*b2)

         sfm(3) = epsi    * (g1+g2)
         sfm(4) = epsi**2 * g2

         sfm(5) = epsi**2 * b1
         sfm(6) = epsi**3 * (b2/3.+b3+b4)
         sfm(7) = epsi    * (b2/3.-b3)
         sfm(8) = epsi**2 * (b2/3.-b4)
      ENDIF

c      write(*,*)'F1,F2,G1,G2',F1,F2,G1,G2

 10   CONTINUE
      IF (IELAS.LE.0) RETURN
 20   CONTINUE

!
!-----Contribution from the elastic channel
!
      WRITE(6,*) "ELASTIC CONTRIBUTION"
      EPSI = AP2/T*(AMT/AMP)**2
      tau  = t/4./amt**2
      tau1 = 1.+tau

      call ffhe3(t,ge,gm)

!$$$  GE = (1.0d0 + t/0.71d0)**-2 ! Proton form factor used by Mo and Tsai
!$$$  GM = 2.793*GE

      XNU_EL = 0.5*T/AMT

      IF (ABS( ANU-XNU_EL ).GT.1.0D-3) THEN
         SPREAD = 0.0
      ELSE
         SPREAD = 1.0
      ENDIF
      
      f1 =     amt * tau    *            gm**2       * SPREAD
      f2 = 2.0*amt * tau    * (ge**2+tau*gm**2)/tau1 * SPREAD
      g1 =     amt * tau    *    gm*(ge+tau*gm)/tau1 * SPREAD
      g2 =     amt * tau**2 *    gm*(ge-gm)    /tau1 * SPREAD

      FACTOR  = AMP/AMT
      FACTOR2 = FACTOR**2
      FACTOR3 = FACTOR2*FACTOR
!-----
      sfm(1) = sfm(1) + EXTAI2 * FACTOR            * (un*f1+qn/6.*b1)
      sfm(2) = sfm(2) + EXTAI2 * FACTOR  * epsi    * (un*f2+qn/6.*b2)

      sfm(3) = sfm(3) + EXTAI2 * FACTOR2 * epsi    * (g1+g2)
      sfm(4) = sfm(4) + EXTAI2 * FACTOR3 * epsi**2 * g2

      sfm(5) = sfm(5) + EXTAI2 * FACTOR3 * epsi**2 * b1
      sfm(6) = sfm(6) + EXTAI2 * FACTOR3 * epsi**3 * (b2/3.+b3+b4)
      sfm(7) = sfm(7) + EXTAI2 * FACTOR  * epsi    * (b2/3.-b3)
      sfm(8) = sfm(8) + EXTAI2 * FACTOR2 * epsi**2 * (b2/3.-b4)
!-----
 30   CONTINUE
      return
      end
! STRF End #####################################################################

! FFHE3 Start ------------------------------------------------------------------
      subroutine ffhe3(t,ge,gm)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      tf=t/chbar**2
      IF (TF.GE.0.0) THEN
         qf=sqrt(tf)
      ELSE
         QF=-SQRT(-TF)
      ENDIF
      a=.675
      b=.366
      c=.836
      am=.654
      bm=.456
      cm=.821
      d=-6.78d-3
      p=.9
      q0=3.98
      f0=ddexp(-a**2*tf) - b**2*qf**2*ddexp(-c**2*tf)
      fm=ddexp(-am**2*tf) - bm**2*qf**2*ddexp(-cm**2*tf)
      df=d*ddexp(-((qf-q0)/p)**2)
      ge=(f0+df)*tarz
      gm=fm*tara      * (-2.13)
      end
! FFHE3 End ####################################################################

! DDEXP Start ------------------------------------------------------------------
      double precision function ddexp(x)
      implicit real*8(a-h,o-z)
        ddexp=0.
        if(x.gt.-50.)ddexp=exp(x)
      return
      end
! DDEXP End ####################################################################

! BLOCK DATA Start -------------------------------------------------------------
!          CHECK THIS
      block data
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     &     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/p/pi,pi2,alfa,i1(8),i2(8)
      data
     &amm/2.7928456d0/,amn/-1.913148d0/,chbar/.197328d0/,barn/.389379d6/
     &aml/.511000d-3/,aml2/.261112d-6/,al2/.522240d-6/,
!     &amt/2.8094d0/      ! He-3
     &amt/0.93827231d0/, ! to check proton elastic tail
     &tara/3d0/,         ! CMASS (only used in elastic and q.e. calculation)
     &tarz/2d0/,         ! CMASS (only used in elastic and q.e. calculation)
     &fermom/.164d0/,    !  - ? Fermi momentum in 3He
     &isf20/4/,
     &pi/3.1415926d0/,pi2/9.869604d0/,alfa/.729735d-2/,
     &amc2/1.151857d0/,
     &amp/.938272d0/,amh/.938272d0/,
     &i2/1,1,1,2,3,3,1,2/,
     &i1/3,3,4,4,3,3,3,3/
      end
! BLOCK DATA End ###############################################################

! DELTA_FTN Start --------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DELTA_FTN(X,XM,SIG)
      DOUBLE PRECISION X,XM,SIG,PI,SQ2PI
      DATA PI/3.1415926D0/
      DATA SQ2PI/2.506628275D0/

      DELTA_FTN = DEXP(-0.5D0*((X-XM)/SIG)**2)/(SIG*SQ2PI)
      RETURN
      END
! DELTA_FTN End ################################################################

! TAIL_INTEGRAND Start ---------------------------------------------------------
      DOUBLE PRECISION FUNCTION TAIL_INTEGRAND(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION SFM0(8),SFM(8),TM(8,6)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire,ich
      common/p/pi,pi2,alfa,i1(8),i2(8)
      COMMON/TAIL_INTEGRAL/TAU,TM

c      write(*,*)'TAIL_INTEGRAND:TAU,R',TAU,R

      CALL STRF(TAU,R,SFM,SFM0)
      SUM = 0.0D0
      DO ISF = ISF1,ISF2,ISF3
         IF (ISF.EQ.3.OR.ISF.EQ.4) THEN
            PPOL = -PN
         ELSE IF (ISF.EQ.5) THEN
            PPOL = QN/6.0D0
         ELSE
            PPOL = 1.0D0
         ENDIF
         DO IRR = 1,I1(ISF)+I2(ISF)-1
            IF (IRR.EQ.1) THEN
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*(SFM(ISF)/(Y+
     &              R*TAU)**2-SFM0(ISF)/Y**2)
            ELSE
               TEMP = -0.5*TM(ISF,IRR)*R**(IRR-2)*PPOL*SFM(ISF)/(Y+
     &              R*TAU)**2
            ENDIF
            SUM = SUM + TEMP
         ENDDO
      ENDDO
      TAIL_INTEGRAND = SUM

c      write(*,*)'TAIL_INTEGRAND',TAIL_INTEGRAND
!SL   WRITE(68,'(5F15.3)') TAU,R,TAIL_INTEGRAND,TAMIN,TAMAX
      RETURN
      END
! TAIL_INTEGRAND End ###########################################################

! TAIL_INT_OVER_R Start --------------------------------------------------------
      DOUBLE PRECISION FUNCTION TAIL_INT_OVER_R(TALN)
      DOUBLE PRECISION TALN,TAU,TAIL_INTEGRAND,TM(8,6)
      EXTERNAL TAIL_INTEGRAND
      DOUBLE PRECISION R_EL,R_QE,DR_E_1SIG,DR_Q_1SIG,RCUT(100),SUM,TEMP
      DOUBLE PRECISION TAU_PASS,DGQUAD
      COMMON/TAIL_INTEGRAL/TAU_PASS,TM
      DOUBLE PRECISION amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,amt,tara,
     &     tarz,fermom,amm,amn,chbar,barn,s,x,sx,sxp,y,ym,w2,als,
     &     alx,alm,aly,sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,
     &     tmi
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     &     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     &     sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      COMMON/ELASTIC/IELAS

      ! C030405
      COMMON/XGOOF/ESMIN,ESMAX,EPMIN,EPMAX,E_COPY,EP_COPY,W2_START 

      DATA NENTRY/0/
      NENTRY = NENTRY+1

      TAU = DEXP(TALN)-XS
      CALL TAILS(TAU,TM)
      TAU_PASS = TAU       ! pass Tau argument to TAIL_INTEGRAND

c      write(*,*)'TAIL_INT_OVER_R:TAU,TAU_PASS',TAU,TAU_PASS

      IF (IELAS.EQ.0) THEN
         NINT = 1
         RCUT(1) = 1.0D-10

!$$$     w2      = amp2+s-y-x                        ! use proton mass for W2
!$$$     RCUT(2) = (W2-AMP2-1D-10)/(1.0D0+TAU)

!030405  RCUT(2) = ( (S-X)-Y -1.0D-10) /(1.0D0+TAU)  ! starting from Elastic
         W2_START=(amt+0.137)**2
c         write(*,*)'TAIL_INT_OVER_R:W2_START',W2_START

         w2      = amp2+s-y-x                        ! W in terms of target mass.
         RCUT(2) = (W2-W2_START)/(1.0D0+TAU)         ! starts integral at lowest input W.

c         write(*,*)'TAIL_INT_OVER_R:W2,TAU',W2,TAU

         IF (RCUT(2).LT.RCUT(1)) THEN 
            RCUT(2) = RCUT(1)
         ENDIF
         IF (TAU.LT.0.0) THEN
            R_MAX = -Y/TAU
            IF (RCUT(2).GT.R_MAX) RCUT(2) = R_MAX
         ENDIF
         IF (RCUT(2).GT.SX) THEN
            RCUT(2) = SX
         ENDIF
!
!------- Inner Integral wrt R.
!------- for DGQUAD N can be 2-16,20,24,32,40,48,64,80,or 96
!
c         write(*,*)'Num integral 1.'

         SUM = 0.0D0

         DO I = 1,NINT
c            write(*,*)'I,RCUT(I),RCUT(I+1)',I,RCUT(I),RCUT(I+1)
            IF (RCUT(I+1).GT.RCUT(I)) THEN
! slifer recommends 040307: 
!               TEMP = DGQUAD(TAIL_INTEGRAND,RCUT(I),RCUT(I+1),96)     
! original version: 
              !TEMP = DGAUSS(TAIL_INTEGRAND,RCUT(I),RCUT(I+1),1.0D-6)
! for testing:
!               TEMP = DGQUAD(TAIL_INTEGRAND,RCUT(I),RCUT(I+1),96,1)
!               TEMP = DGQUADT2(TAIL_INTEGRAND,RCUT(I),RCUT(I+1),96,1)
               CALL simpsx(RCUT(I),RCUT(I+1),20,0.1,TAIL_INTEGRAND,TEMP)
c               write(*,*)'TAIL_INT_OVER_R:TEMP=',TEMP
               SUM  = SUM + TEMP
            ENDIF
         ENDDO
      ELSE IF (IELAS.EQ.1) THEN ! Add contrib. from the QE peak assuming delta function
!$$$     TEMP = TAIL_INTEGRAND(R_QE)*2.0*AMP/(1.0+TAU)

!        Add contribution from the elastic peak
!        First identify elastic and quasielastic R
!        Elastic R
!
         R_EL = (SX-Y)/(1.0 + TAU)
         SUM  = TAIL_INTEGRAND(R_EL)*2.0*AMT/(1.0+TAU)
      ENDIF

      TAIL_INT_OVER_R = SUM*(XS+TAU)
c      write(*,*)'TAIL_INT_OVER_R,XS,TAU',TAIL_INT_OVER_R,XS,TAU
      RETURN
      END
! TAIL_INT_OVER_R End ##########################################################

! TAIL_INTEG Start -------------------------------------------------------------
      DOUBLE PRECISION FUNCTION TAIL_INTEG(TAU_MIN,TAU_MAX)
      DOUBLE PRECISION TAU_MIN,TAU_MAX
      DOUBLE PRECISION TAIL_INT_OVER_R,TEMP,DGQUAD
      EXTERNAL TAIL_INT_OVER_R
      DOUBLE PRECISION TCUT(100)
      DIMENSION NGAUSSPT(7)
      DATA NGAUSSPT/8,16,16,8,16,16,8/
      COMMON/DEBUG/IDEBUG1
      COMMON/XGOOF2/IMARK,JMARK
      DOUBLE PRECISION amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,amt,tara,
     &     tarz,fermom,amm,amn,chbar,barn,s,x,sx,sxp,y,ym,w2,als,
     &     alx,alm,aly,sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,
     &     tmi
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     &     amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     &     sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi

      DOUBLE PRECISION TAU_PASS,TM(8,6),R,TAIL_INTEGRAND
      COMMON/TAIL_INTEGRAL/TAU_PASS,TM
      COMMON/ELASTIC/IELAS
      REAL TEMP1,TEMP2,TEMP3,WEIGHT
      IF (IELAS.EQ.0) THEN                  ! quasi-elastic and inelastic
         TCUT(1) = DLOG(TAU_MIN+XS)
!$$$     TCUT(2) = DLOG(XS)
!$$$     TCUT(3) = DLOG(TAU_MAX+XS)
         TCUT(3) = DLOG(XS-Y/S)              
         TCUT(6) = DLOG(XS+Y/X)              
         DGUESS  = (TCUT(6)-TCUT(3))*1.0D-1
         TCUT(2) = TCUT(3) - DGUESS
         TCUT(4) = TCUT(3) + DGUESS
         TCUT(5) = TCUT(6) - DGUESS
         TCUT(7) = TCUT(6) + DGUESS
         TCUT(8) = DLOG(XS+TAU_MAX)
         ISTART  = 1                    
         IEND    = 7                    

         ! C030405
         ! This should cover the canonical RC triangle fine.
         ! Not sure what the previous choice of limits was for.

corig         TCUT(10) = DLOG(XS-Y/S)
corig         TCUT(11) = DLOG(XS+Y/X)

         !in FUNCTION TAIL_INT_OVER_R, TAU = DEXP(TALN)-XS

         TCUT(10) = DLOG(TAU_MIN+XS)
         TCUT(11) = DLOG(XS+TAU_MAX)

c         DGUESS=sqrt((s-x)**2+4.0*amt**2*y)
c         TCUT(10)=(s-x-DGUESS)/2.0/amt**2
c         TCUT(11)=(s-x+DGUESS)/2.0/amt**2
c         write(*,*)'sx,y(Q2),amt(M),lamda',(s-x),y,amt,DGUESS
c         write(*,*)'tamin,tamax,TAU_MIN,TAU_MAX',tamin,tamax,TAU_MIN,TAU_MAX
c         write(*,*)'Y/S,Y/X',Y/S,Y/X

         ISTART   = 10                   
         IEND     = 10

      ELSE IF (IELAS.EQ.1) THEN
         ISTART  = 1             
         TCUT(1) = DLOG(TAU_MIN+XS)
         TCUT(2) = DLOG(XS)
         TCUT(3) = DLOG(TAU_MAX+XS)
         IEND = 2
      ENDIF

!      
!---- Outer Integral wrt tau.
!---- for DGQUAD N can be 2-16,20,24,32,40,48,64,80,or 96
!----
!---- DGAUSS and DGQUAD give almost identical results (for E94010 kinematics)
!---- with 1.0D-3 precision.  At 1.0D-6 precision it is just way too slow.

c      write(*,*)'Num integral 2'

      TAIL_INTEG = 0.0D0

      DO I = ISTART,IEND
c        write(*,*)'I,TCUT(I),TCUT(I+1)',I,TCUT(I),TCUT(I+1)
! slifer recommends 040307:
        !TEMP       = DGAUSS(TAIL_INT_OVER_R,TCUT(I),TCUT(I+1),1.0D-1)
! what was orirginally in the version that i got:
cxy        TEMP       = DGQUAD(TAIL_INT_OVER_R,TCUT(I),TCUT(I+1),96,2)
!        TEMP       = DGQUADT2(TAIL_INT_OVER_R,TCUT(I),TCUT(I+1),96,2)
c        TEMP=TAIL_INT_OVER_R(TCUT(I))
!        TEMP=0.0D0
! js for testing 
        !TEMP       = DGAUSS(TAIL_INT_OVER_R,TCUT(I),TCUT(I+1),1.0D-1)
        CALL simptx(TCUT(I),TCUT(I+1),100,0.005,TAIL_INT_OVER_R,TEMP)
         TAIL_INTEG = TAIL_INTEG + TEMP
      ENDDO

      RETURN
      END
! TAIL_INTEG End ###############################################################
