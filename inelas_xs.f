        real*8 function type_cxsn(e,ep,theta,xstype)

        implicit none
        INCLUDE 'const.inp'
        INCLUDE 'data_common_const.inp'

        real*8 A,Z,W2,Q2,M,x,nu,xs
        real E,Ep,theta
        real*8 EPS,PF
        real*8 inelas_cxsn,qf_cxsn,in_cxsn,SIGQFS
        real*8 F1,F2,r,temp,F1qe,F2qe
        real*8 E0,TH_deg,W
        real*8 xs_pol_perp,xs_pol_para
        integer xstype,g12_type
        common/xsset_g1g2/g12_type

        PF=130.0
        EPS=10.0
        E0=e*1000.
        TH_deg=theta*180.0/3.14159265354
        W=(e-ep)*1000.0
        Z=Zset
        A=Aset

c        write(*,*)'typexs:e,ep,theta,xstype',e,ep,theta,xstype

        if(xstype.eq.0)then
           type_cxsn=inelas_cxsn(e,ep,theta)
        elseif(xstype.eq.1)then
           type_cxsn=qf_cxsn(e,ep,theta)
           if(type_cxsn.lt.0)type_cxsn=0.0
        elseif(xstype.eq.2)then
           type_cxsn=in_cxsn(e,ep,theta)
        elseif(xstype.eq.3)then
           type_cxsn=SIGQFS(E0,TH_deg,W,Z,A,EPS,PF)
           type_cxsn=type_cxsn*1.D+07    !nb/MeV/Sr
        elseif(xstype.eq.4)then
           call cross_tot_mod(Z,A,e,ep,theta,type_cxsn,f1,f2,q2,w2)
           type_cxsn=type_cxsn*1.D-03    !nb/MeV/Sr
        elseif(xstype.eq.5)then
           call cross_qe_mod(Z,A,e,ep,theta,type_cxsn,f1,f2,q2,w2)
           type_cxsn=type_cxsn*1.D-03    !nb/MeV/Sr
           if(type_cxsn.lt.0)type_cxsn=0.0
        elseif(xstype.eq.6)then
           type_cxsn=xs_pol_perp(e,ep,theta,g12_type)
        elseif(xstype.eq.7)then
           type_cxsn=xs_pol_para(e,ep,theta,g12_type)
        else
           write(*,*)'wrong xs type. stop'
           stop
        endif

        if(type_cxsn.ne.type_cxsn)type_cxsn=0.0

        return
        END

        subroutine xsg1g2_init(g12_type_in)
        implicit none
        integer g12_type_in,g12_type
        common/xsset_g1g2/g12_type

        g12_type=g12_type_in

        end


!!!!!!!!function to compute inelastic cross section
        real*8 function inelas_cxsn(e,ep,theta)

        implicit none
        INCLUDE 'const.inp'
        INCLUDE 'data_common_const.inp'

!       INELASTIC CROSS SECTION FOR nucleus
!       P. Bosted from A. Deur

!
!       E      - INCIDENT ENERGY (GeV)
!       THETA  - SCATTERING ANGLE (RADIANS)
!       M      - TARGET MASS (GeV)
!       A      - NUCLEUS NUMBER OF NUCLEONS
!       Z      - NUMBER OF PROTONS IN NUCLEUS
!       RECOIL - RECOIL FACTOR
!       EP     - SCATTERED ENERGY (GEV)
!       NU     - ENERGY TRANSFER (GeV)
!       Q2     - MOMENTUM TRANSFER SQUARED (GeV**2)
!       MOTT   - MOTT CROSS SECTION
!       F1,F2  - are structure functions per nucleus(not per nucleon)
!       R      - ratio sigl/sigt
!

        real*8 A,Z,W2,Q2,M,x,nu,xs
        real E,Ep,theta,MOTT
        real*8 F1,F2,r,temp,F1qe,F2qe
        integer debug

        parameter (debug=0)
        PARAMETER (M=0.93827)
c        PARAMETER (A=1.)
c        PARAMETER (Z=1.)
c        PARAMETER (A=3.)
c        PARAMETER (Z=2.)
        
!
        A=Aset
        Z=Zset
        
        nu=E-Ep
        Q2=4.*E*Ep*(sin(theta/2.))**2.
        w2=M**2.+2.*M*nu-Q2
        x=Q2/2./M/nu
!       print*, e,ep,theta,Q2,w2,x,M

!       call F1F2IN06(Z,A,q2,w2,F1,F2,r)
        call F1F2IN09(Z,A,q2,w2,F1,F2,r)    ! DELTA &DIS TOGETHER
        call F1F2QE09(Z,A,q2,w2,F1qe,F2qe)  ! QUASI-ELASTIC ONLY

c        write(*,*)'inelas:q2,w2,F1,F2',q2,w2,(F1+F1qe),(F2+F2qe)

c$$$        call F1F2IN07(Z,A,q2,w2,F1,F2,r)    ! DELTA &DIS TOGETHER
c$$$        call F1F2QE07(Z,A,q2,w2,F1qe,F2qe)  ! QUASI-ELASTIC ONLY

c$$$        F1=1.3*F1
c$$$        F2=1.2*F2

        xs=MOTT(E,THETA) ! in ub/Sr
        temp=(2/M*(F1+F1qe)*tan(theta/2)**2+(F2+F2qe)/nu)
        inelas_cxsn=xs*temp ! in nb/MeV/Sr

c        write(*,*)'inelas:inelas_cxsn,()',inelas_cxsn,temp

        if (debug.eq.1) then
           print*,'For nucleus of A=',A,' and Z=',Z
           print*,'W=',sqrt(W2),' GeV ','x=',x
           print*,'Mott XS:',xs,' ub/Sr'
           print*,'inelastic proton XS:',inelas_cxsn,' nb/MeV/Sr'
           print*,'F1=',F1,' F2=',F2
        endif
!PRINT*,e,theta*180./3.1416,F1,F2,xs,temp
!print*,'Xsec= ',inelas_cxsn,' microbarn/sr/GeV'
        return
        END


        real*8 function qf_cxsn(e,ep,theta)

        implicit none
        INCLUDE 'const.inp'
        INCLUDE 'data_common_const.inp'

        real*8 A,Z,W2,Q2,M,x,nu,xs
        real E,Ep,theta,MOTT
        real*8 F1,F2,r,temp,F1qe,F2qe
        integer debug

        parameter (debug=0)
        PARAMETER (M=0.93827)

        A=Aset
        Z=Zset

        nu=E-Ep
        Q2=4.*E*Ep*(sin(theta/2.))**2.
        w2=M**2.+2.*M*nu-Q2
        x=Q2/2./M/nu

        call F1F2QE09(Z,A,q2,w2,F1qe,F2qe)  ! QUASI-ELASTIC ONLY

        xs=MOTT(E,THETA) ! in ub/Sr
        temp=(2/M*(F1qe)*tan(theta/2)**2+(F2qe)/nu)
        qf_cxsn=xs*temp ! in nb/MeV/Sr

        return
        END

        real*8 function in_cxsn(e,ep,theta)

        implicit none
        INCLUDE 'const.inp'
        INCLUDE 'data_common_const.inp'

        real*8 A,Z,W2,Q2,M,x,nu,xs
        real E,Ep,theta,MOTT
        real*8 F1,F2,r,temp,F1qe,F2qe
        integer debug

        parameter (debug=0)
        PARAMETER (M=0.93827)

        A=Aset
        Z=Zset

        nu=E-Ep
        Q2=4.*E*Ep*(sin(theta/2.))**2.
        w2=M**2.+2.*M*nu-Q2
        x=Q2/2./M/nu

        call F1F2IN09(Z,A,q2,w2,F1,F2,r)    ! DELTA &DIS TOGETHER

        xs=MOTT(E,THETA) ! in ub/Sr
        temp=(2/M*(F1)*tan(theta/2)**2+(F2)/nu)
        in_cxsn=xs*temp ! in nb/MeV/Sr

        return
        END


!
!     MOTT CROSS SECTION
!
        REAL FUNCTION MOTT(ES, THETA)
        IMPLICIT NONE

        INCLUDE 'const.inp'

!       
!       ES    - INCIDENT ENERGY (GEV)
!       THETA - SCATTERING ANGLE (RADIANS)
!       Z     - TARGET ATOMIC NUMBER
!       ALPHA - FINE STRUCTURE CONSTANT
!       
        REAL ES, THETA, ALPHA, C2, S,Z,ST
        PARAMETER (ALPHA = 1./137.035989561)
        PARAMETER (Z = 6)
!
        C2 = COS(0.5*THETA)
        ST = SIN(0.5*THETA)
        S = ALPHA*C2/(2.*ES*ST*ST) 
        MOTT = S*S*hbc2*1000.     ! microbarn
c        write(*,*)'MOTT:MOTT,hbc2',MOTT,hbc2

        RETURN
        END

