      SUBROUTINE xyrad2d(Exy,THxy,nuxy,SIGNR,SIGRAD,xstype,din1,din2)
      IMPLICIT NONE

      INCLUDE 'const.inp'
      INCLUDE 'data_common_const.inp'

      REAL*8 SIGNR,SIGRAD,DEL,nu,Z,A,SPENCExy,TH
      REAL*8 type_cxsn
      REAL*8 Esmin,Mxy,Exy,THxy,nuxy,Ep
      real*8 e0gev,epgev
      integer xstype,cxstype
      real*8 din1,din2
      real*8 del1,del2
      real*8 Es_int
      real*8 Iprob,tr,Es_xy,Ep_xy,TH_xy
      real*8 eta
      real*8 Ta,Tb,xia,xib,ba,bb
      COMMON/TAIL_INT/Es_xy,Ep_xy,cxstype,tr,TH_xy,SPENCExy
     &                ,del1,del2
      COMMON/TAIL_INT2/Es_int,Ta,Tb,xia,xib,ba,bb
      real*8 int_es,int_ep,int_esdp,int_epds,temp,temp1,temp2,low,high
      real*8 D1,D2,PI,QMS,ARG,xF,bt,temp3
      integer NSP
      real*8 coef1,coef2

      Mxy=938.27
      PI=acos(-1.0)
      Ep=Exy-nuxy
      Esmin=Ep/(1.0-Ep*(1.0-cos(THxy))/Mxy)

      cxstype=xstype
      Es_xy=Exy
      Ep_xy=Ep
      TH_xy=THxy

      A=Asir
      Z=Zsir

      Ta=Taset
      Tb=Tbset

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xia=(3.141592653*EMASS/2./ALPH)*Ta/((Z+eta)*LOG(183.*Z**(-1./3.)))
      xib=(3.141592653*EMASS/2./ALPH)*Tb/((Z+eta)*LOG(183.*Z**(-1./3.)))

      ba=(4.0/3.0)*(1.0+(Z+1)/9.0/(Z+eta)/LOG(183.0*Z**(-1./3.)))
      bb=ba

c      write(*,*)'ba=',ba

      del1=din1
      del2=din2

      QMS=4.0*Exy*Ep*sin(THxy/2.0)**2

      ARG=cos(THxy/2.0)**2
      SPENCExy= PI**2/6.-LOG(ARG)*LOG(1.-ARG)
      DO 10 NSP=1,50
 10     SPENCExy = SPENCExy-ARG**NSP/FLOAT(NSP)**2

c      write(*,*)'Ep,Esmin',Ep,Esmin

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
c      tr=D1/2./xb
      tr=D1/2.
c      write(*,*)'tr,xb,SPENCExy=',tr,xb,SPENCExy

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-14./9. ! vac pol, vertex corr
      D2 = D2 +0.5*(PI**2/6.-SPENCExy)    !Schwinger correction
      D2 = D2 -1./4.*( LOG( Exy/(Ep) ) )**2        ! Correct. to angle peak. appr.
      D2 = D2*(2.*ALPH/PI)
c      D2 = D2+0.5772*xb*(Tb+Ta+2.0*tr)                   ! ~1/Gamma
      xF = (1.+D2)

c      write(*,*)'xF',xF

      SIGNR=type_cxsn(Exy/1000.,Ep/1000.,TH_xy,cxstype)   !Born XS
c      write(*,*)'xyrad2d:SIGNR,cxstype',SIGNR,cxstype
c      write(*,*)'xyrad2d:Exy,Ep,TH_xy',Exy,Ep,TH_xy
c      stop


      !Singlar-integral: Es singular, Ep regular
      low=Ep_xy+del2
      high=Exy/(1.0+Exy*(1.0-cos(TH_xy))/938.27)
      temp=int_epds(low) !intialize before putting in Simpson int.
      CALL simpsx(low,high,150,0.01,int_epds,temp) !proper

      !Singlar-integral: Ep singular, Es regular
      low=Exy-100.
      temp1=int_esdp(low)   !intialize before putting in Simpson int.
      low=Esmin
      high=Exy-del1
      CALL simpsx(low,high,150,0.01,int_esdp,temp1)


      !2D integral part
      low=Esmin+10.0
      temp2=int_es(temp1)  !initialization before putting in Simpson int.
      low=Esmin
      high=Exy-del1
      CALL simptx(low,high,150,0.01,int_es,temp2)
c      write(*,*)'intg int_es',temp2
c      cxstype=xstype


      !double-singular part
      bt=tr+Tb*bb !Ts
      coef1=(del1/Exy)**bt/(gamma(1.0+bt))*(1.0-xib/del1/(1-bt))

      bt=tr+Ta*ba !Tp
      coef2=(del2/Ep )**bt/(gamma(1.0+bt))*(1.0-xia/del2/(1-bt))

      temp3=SIGNR*xF*coef1*coef2

      SIGRAD=temp+temp1+temp2+temp3
c      write(*,*)'xyrad2d:SIGRAD',SIGRAD

c      write(*,*)'xyrad2d:both singular',temp3
c      write(*,*)'xyrad2d: Es singular',temp
c      write(*,*)'xyrad2d: Ep singular',temp1
c      write(*,*)'xyrad2d: regular 2d',temp2
c      write(*,*)'xyrad2d:SIGNR',SIGNR
c      write(*,*)'xyrad2d:SIGRAD',SIGRAD

      RETURN
      END


      real*8 function Iprob(E1,E2,b,t,xi)
      implicit none
      real*8 E1,E2,t,b,gamma,v,bt,xi

      if(E1.le.E2)then
        write(*,*)'Iprob:E1<E2;stop',E1,E2
        stop

        Iprob=0.0
        return
      endif

c      b=4.0/3.0
      bt=b*t
      v=(E1-E2)/E1
      Iprob=(1.0/(gamma(1.0+bt)))*(v**bt)
     %  *(bt/(E1-E2)*(1.0-v+3.0*v*v/4.0)+xi/(E1-E2)/(E1-E2))

c      write(*,*)'bt,E1,E2,v,v**bt',bt,E1,E2,v,v**bt

      return
      END

      real*8 function int_es(Es)
      implicit none
      real*8 Iprob,Es,Es_int
      integer xstype,cxstype
      real*8 del1,del2
      real*8 tr,Es_xy,Ep_xy,TH_xy,SPENCExy
      real*8 int_ep,epmax,epmin
      real*8 Ta,Tb,xia,xib,ba,bb
      COMMON/TAIL_INT/Es_xy,Ep_xy,cxstype,tr,TH_xy,SPENCExy
     &                ,del1,del2
      COMMON/TAIL_INT2/Es_int,Ta,Tb,xia,xib,ba,bb
      real*8 xs,del,low


      Es_int=Es
c      write(*,*)'Es_xy,cxstype,tr,TH_xy=',Es_xy,cxstype,tr,TH_xy

      epmax=Es/(1.0+Es*(1.0-cos(TH_xy))/938.27)

      xs=int_ep(Ep_xy+del2) !initialize before Simp int

      del=del2
      low=Ep_xy+del
c      write(*,*)'int_ep',xs
      if(low.lt.epmax)then
        CALL simpsx(Ep_xy+del,epmax,100,0.01,int_ep,xs)
      else
        xs=0.0
c        write(*,*)'Ep_xy,epmax',Ep_xy,epmax
c        stop
      endif
      int_es=xs

      return
      END

      real*8 function int_ep(Ep)
      implicit none
      INCLUDE 'const.inp'

      real*8 Iprob,Ep,Es_int
      integer xstype,cxstype
      real*8 del1,del2
      real*8 tr,Es_xy,Ep_xy,TH_xy,SPENCExy
      real*8 Ta,Tb,xia,xib,ba,bb
      COMMON/TAIL_INT/Es_xy,Ep_xy,cxstype,tr,TH_xy,SPENCExy
     &                ,del1,del2
      COMMON/TAIL_INT2/Es_int,Ta,Tb,xia,xib,ba,bb
      real*8 type_cxsn,xs
      real*8 D1,D2,PI,QMS,xF,trl
      real*8 ip1,ip2,trla,trlb

      PI=acos(-1.0)

      QMS=4.0*Es_int*Ep*sin(TH_xy/2.0)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      trl=D1/2.
      trla=trl/ba+Ta
      trlb=trl/bb+Tb

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-14./9. ! vac pol, vertex corr
      D2 = D2 +0.5*(PI**2/6.-SPENCExy)
      D2 = D2 -1./4.*( LOG( Es_int/(Ep) ) )**2        ! Correct. to angle peak. appr.
      D2 = D2*(2.*ALPH/PI)
c      D2 = D2+0.5772*xb*(2.0*trl+Ta+Tb)                   ! ~ 1/Gamma
      xF = (1.+D2)

      ip1=Iprob(Es_xy,Es_int,bb,trlb,xib)
      ip2=Iprob(Ep,Ep_xy,ba,trla,xia)
      xs=type_cxsn(Es_int/1000.,Ep/1000.,TH_xy,cxstype)

      int_ep=xF*xs*ip1*ip2

      if(int_ep.ne.int_ep)then
        write(*,*)'Es_int,Ep,Ep_xy,int_ep=',Es_int,Ep,Ep_xy,int_ep
        stop
      endif

      return
      END

      real*8 function int_esdp(Es)
      implicit none
      INCLUDE 'const.inp'

      real*8 Iprob,Es,Es_int
      integer xstype,cxstype
      real*8 del1,del2
      real*8 tr,Es_xy,Ep_xy,TH_xy,SPENCExy
      real*8 int_ep,epmax,epmin
      real*8 Ta,Tb,xia,xib,ba,bb
      COMMON/TAIL_INT/Es_xy,Ep_xy,cxstype,tr,TH_xy,SPENCExy
     &                ,del1,del2
      COMMON/TAIL_INT2/Es_int,Ta,Tb,xia,xib,ba,bb
      real*8 xs,del,low,type_cxsn
      real*8 D1,D2,PI,QMS,xF,b,bt,trl

      PI=acos(-1.0)

      QMS=4.0*Es*Ep_xy*sin(TH_xy/2.0)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      trl=D1/2./bb+Tb

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-14./9. ! vac pol, vertex corr
      D2 = D2 +0.5*(PI**2/6.-SPENCExy)
      D2 = D2 -1./4.*( LOG( Es/(Ep_xy) ) )**2        ! Correct. to angle peak. appr.
      D2 = D2*(2.*ALPH/PI)
c      D2 = D2+0.5772*xb*(2.0*trl)                   ! ~1/Gamma
      xF = (1.+D2)

      b=4.0/3.0
c      bt=b*(trl-Tb+Ta)
      bt=(trl-Tb)*bb+Ta*ba

      int_esdp=Iprob(Es_xy,Es,bb,trl,xib)
      xs=type_cxsn(Es/1000.,Ep_xy/1000.,TH_xy,cxstype)
      int_esdp=int_esdp*xs
      int_esdp=int_esdp*xF
      int_esdp=int_esdp*(del2/Ep_xy)**bt/(gamma(1.0+bt))
     %     *(1.0-xia/del2/(1-bt)) !proper

c      write(*,*)'int_esdp:Q2,xF,int_esdp,xs',QMS,xF,int_esdp,xs

      return
      END

      real*8 function int_epds(Ep)
      implicit none
      INCLUDE 'const.inp'

      real*8 Iprob,Es,Es_int,Ep
      integer xstype,cxstype
      real*8 del1,del2
      real*8 tr,Es_xy,Ep_xy,TH_xy,SPENCExy
      real*8 int_ep,epmax,epmin
      real*8 Ta,Tb,xia,xib,ba,bb
      COMMON/TAIL_INT/Es_xy,Ep_xy,cxstype,tr,TH_xy,SPENCExy
     &                ,del1,del2
      COMMON/TAIL_INT2/Es_int,Ta,Tb,xia,xib,ba,bb
      real*8 xs,del,low,type_cxsn
      real*8 D1,D2,PI,QMS,xF,b,bt,trl

      PI=acos(-1.0)

      QMS=4.0*Es_xy*Ep*sin(TH_xy/2.0)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      trl=D1/2./ba+Ta

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-14./9. ! vac pol, vertex corr
      D2 = D2 +0.5*(PI**2/6.-SPENCExy)
      D2 = D2 -1./4.*( LOG( Es_xy/(Ep) ) )**2        ! Correct. to peak. appr.
      D2 = D2*(2.*ALPH/PI)
c      D2 = D2+0.5772*xb*(2.0*trl)                   ! ~1/Gamma
      xF = (1.+D2)

      b=4.0/3.0
c      bt=b*(trl-Ta+Tb)
      bt=ba*(trl-Ta)+bb*Tb

      int_epds=Iprob(Ep,Ep_xy,ba,trl,xia)
      xs=type_cxsn(Es_xy/1000.,Ep/1000.,TH_xy,cxstype)
      int_epds=int_epds*xs
      int_epds=int_epds*xF
      int_epds=int_epds*(del1/Es_xy)**bt/(gamma(1.0+bt))
     %     *(1.0-xib/del1/(1-bt))

      return
      END



