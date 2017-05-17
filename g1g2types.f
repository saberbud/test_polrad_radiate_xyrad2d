      subroutine g1g2types(x,Q2,np,nn,xstype,g1,g2)
      implicit none
      integer xstype,xstype0,i,findg1g2f,he3comb
      integer findg1g2dssvf,findg1g2grsvf
      real*8 xyn1,xyn2,typein
      real*8 x,Q2,np,nn,g1,g2
      real*8 g1p,g2p,g1n,g2n
      real*8 g1cc,g2cf,reg1pf,reg2pf,reg1nf,reg2nf,g2_nlo_ww
      real*8 reg1pdssvf,reg2pdssvf,reg1ndssvf,reg2ndssvf
      real*8 reg1pgrsvf,reg2pgrsvf,reg1ngrsvf,reg2ngrsvf
      real*8 poln,polp
      integer xy_maid_on,xy_he3g_on
      common/g1g2set/xy_maid_on,xy_he3g_on

      xyn1=1.0
      xyn2=0.0
      poln=0.87
      polp=-0.027

      g1p=0
      g2p=0
      g1n=0
      g2n=0

      xstype0=xstype

      he3comb=0
      if(xstype.ge.10)then
        he3comb=1
        xstype=xstype-10
      endif

      if(Q2.le.1.0)then
c        xstype0=xstype
        xstype=1

        g1=0
        g2=0
        if(xy_maid_on.eq.0)then
         xstype=xstype0
         return
        endif
      endif

      if(x.gt.0.9999.or.x.le.0.0001)then
        g1=0
        g2=0

        xstype=xstype0
        return
      endif

      if(xstype.eq.0)then
        i=findg1g2f(Q2,x)
        g1p=reg1pf()
        g2p=reg2pf()
        g1n=reg1nf()
        g2n=reg2nf()
      elseif(xstype.eq.1)then
        call GET_G(x,Q2,g1p,g2p)
        call GET_GN(x,Q2,g1n,g2n)
      elseif(xstype.eq.2)then
        typein=0 !dssv
        g1p=g1cc(x,Q2,xyn1,xyn2,typein)
        g2p=g2cf(x,Q2,xyn1,xyn2,typein)
        g1n=g1cc(x,Q2,xyn2,xyn1,typein)
        g2n=g2cf(x,Q2,xyn2,xyn1,typein)
      elseif(xstype.eq.3)then
        typein=1 !grsv
        g1p=g1cc(x,Q2,xyn1,xyn2,typein)
        g2p=g2cf(x,Q2,xyn1,xyn2,typein)
        g1n=g1cc(x,Q2,xyn2,xyn1,typein)
        g2n=g2cf(x,Q2,xyn2,xyn1,typein)
      elseif(xstype.eq.4)then
        i=findg1g2f(Q2,x)
        g1p=reg1pf()
        g1n=reg1nf()
        i=0
        g2p=g2_nlo_ww(x,Q2,i)
        i=1
        g2n=g2_nlo_ww(x,Q2,i)
      elseif(xstype.eq.5)then
        i=findg1g2dssvf(Q2,x)
        g1p=reg1pdssvf()
        g2p=reg2pdssvf()
        g1n=reg1ndssvf()
        g2n=reg2ndssvf()
      elseif(xstype.eq.6)then
        i=findg1g2grsvf(Q2,x)
        g1p=reg1pgrsvf()
        g2p=reg2pgrsvf()
        g1n=reg1ngrsvf()
        g2n=reg2ngrsvf()
      else
         write(*,*)'wrong xs g1g2 type. stop',xstype
         stop
      endif

      if(Q2.le.1.0)then
        xstype=xstype0
      endif

      if(he3comb.lt.1)then
        g1=g1p*np+g1n*nn
        g2=g2p*np+g2n*nn
      else
        xstype=xstype0
c        write(*,*)'he3comb'
        g1=g1p*(np*polp-0.014)+g1n*(nn*poln+0.056)
        g2=g2p*(np*polp-0.014)+g2n*(nn*poln+0.056)
      endif

      return

      end

      subroutine g1g2types_init(maid_on,he3g_on)
      implicit none
      integer maid_on,he3g_on
      integer xy_maid_on,xy_he3g_on
      common/g1g2set/xy_maid_on,xy_he3g_on

      xy_maid_on=maid_on
      xy_he3g_on=he3g_on

      end

