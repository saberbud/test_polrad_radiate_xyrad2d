        real*8 function g2cf(x,Q2,np,nn,typein)

        implicit none
        real*8 x,Q2,np,nn,typein,g1cc,g1intg,temp,high
        real*8 xyQ2,xynp,xynn,xytypein
        COMMON/g2trans/xyQ2,xynp,xynn,xytypein

        xyQ2=Q2
        xynp=np
        xynn=nn
        xytypein=typein

        g2cf=-g1cc(x,Q2,np,nn,typein)
c        write(*,*)'g2cf before intg',g2cf

        temp=g1intg(x)
        temp=0.0

        high=1.0
        CALL simpux(x, high, 150, 0.001, g1intg, temp)
c        write(*,*)'intg=',temp

        g2cf=g2cf+temp

c        write(*,*)'x,g1intg(x)=',x,g1intg(x)

        return
        END

        real*8 function g1intg(x)

        implicit none
        real*8 x,g1cc
        real*8 xyQ2,xynp,xynn,xytypein
        COMMON/g2trans/xyQ2,xynp,xynn,xytypein

        g1intg=g1cc(x,xyQ2,xynp,xynn,xytypein)
        g1intg=g1intg/x

        return
        END
