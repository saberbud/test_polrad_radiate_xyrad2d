      real*8 function xs_pol_perp(e,ep,theta,xstype)

      implicit none
      INCLUDE 'const.inp'
      INCLUDE 'data_common_const.inp'

      real*8 A,Z,W2,Q2,M,x,nu,xs
      real*8 E,Ep,theta
      real*8 xynp,xynn
      real*8 g1,g2
      integer xstype

      Z=Zset
      A=Aset
      xynp=Z
      xynn=A-Z

      M=0.93827
      nu=e-ep
      Q2=4.*E*Ep*(sin(theta/2.))**2.
      w2=M**2.+2.*M*nu-Q2
      x=Q2/2./M/nu

c      write(*,*)'xs_pol_perp:Z,A,xynp,xynn,ALPH',Z,A,xynp,xynn,ALPH
c      write(*,*)'xs_pol_perp:e,ep,theta,x,Q2',e,ep,theta,x,Q2

      call g1g2types(x,Q2,xynp,xynn,xstype,g1,g2)
c      write(*,*)'xs_pol_perp:g1,g2',g1,g2
c      write(*,*)g1,g2
c      write(22,*)g1,g2

      xs=4.0*ALPH*ALPH*Ep*Ep*sin(theta)*(g1+g2*2.0*E/nu)/(M*Q2*E*nu)

      xs_pol_perp=xs*389.379
      return
      END

      real*8 function xs_pol_para(e,ep,theta,xstype)

      implicit none
      INCLUDE 'const.inp'
      INCLUDE 'data_common_const.inp'

      real*8 A,Z,W2,Q2,M,x,nu,xs
      real*8 E,Ep,theta
      real*8 xynp,xynn
      real*8 g1,g2
      integer xstype

      Z=Zset
      A=Aset
      xynp=Z
      xynn=A-Z

      M=0.93827
      nu=e-ep
      Q2=4.*E*Ep*(sin(theta/2.))**2.
      w2=M**2.+2.*M*nu-Q2
      x=Q2/2./M/nu

      call g1g2types(x,Q2,xynp,xynn,xstype,g1,g2)

c      write(*,*)g1,g2
c      write(22,*)g1,g2

      xs=4.0*ALPH*ALPH*Ep*((E+Ep*cos(theta))*g1-g2*2.0*M*x)/(M*Q2*E*nu)

      xs_pol_para=xs*389.379
      return
      END






