        real*8 function g2_nlo_ww(x,Q2,pn)

        implicit none
        real*8 x,Q2,g1intg_nlo,temp,high,reg1pf,reg1nf
        real*8 xyintQ2
        integer i,pn,xyintpn,findg1g2f
        COMMON/g2nloint/xyintQ2,xyintpn

        xyintQ2=Q2
        xyintpn=pn

        i=findg1g2f(Q2,x)
        if(pn.eq.0)then
          g2_nlo_ww=-reg1pf()
        elseif(pn.eq.1)then
          g2_nlo_ww=-reg1nf()
        endif
c        write(*,*)'g2_nlo_ww before intg',g2_nlo_ww

        temp=g1intg_nlo(x)
c        write(*,*)'g2_nlo_ww:g1intg_nlo',temp

        temp=0.0

        high=1.0
        CALL simpux(x, high, 150, 0.001, g1intg_nlo, temp)
c        write(*,*)'intg=',temp

        g2_nlo_ww=g2_nlo_ww+temp

c        write(*,*)'x,g1intg(x)=',x,g1intg(x)

        return
        END

        real*8 function g1intg_nlo(x)

        implicit none
        real*8 x,reg1pf,reg1nf,xyintQ2
        integer i,pn,xyintpn,findg1g2f
        COMMON/g2nloint/xyintQ2,xyintpn

        pn=xyintpn

        i=findg1g2f(xyintQ2,x)
        if(pn.eq.0)then
          g1intg_nlo=reg1pf()
        elseif(pn.eq.1)then
          g1intg_nlo=reg1nf()
        endif
        g1intg_nlo=g1intg_nlo/x

        return
        END
