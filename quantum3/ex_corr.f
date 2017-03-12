       subroutine vexcorr(cdd,cdu, r, vxcd,vxcu, ecxt)
       implicit double precision (a-h,o-z)
c  ************************************************************************
c
c      vexcorr generates the exchange-correlation potential from 
c      the electron charge density at radius r.
c
c      All quantities in atomic units.
c
c  INPUT:
c    cdd  = spin-down charge at radius r
c    cdu  = spin-up charge at radius r
c    r    = radius
c
c  OUTPUT:
c    vxcd = spin-down exchange-correlation potential
c    vxcu = spin-up exchange-correlation potential
c    exct = total exchange-correlation energy per elec
c
c  ************************************************************************
c
c
      real*8 cdd, cdu, r, vxcd,vxcu, exct
c
c
       pi = 4*atan(1.D0)
c
       trd = 1.D0/3.D0
       ftrd = 4*trd
       tftm = 2**ftrd-2
       a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
       alp = 2 * trd
c
       cdsum = cdd + cdu             ! spin-up + spin-down charge density
       if (cdsum .le. 0.D0) then
         vxcd = 0.0d0
         vxcu = 0.0d0
         exct = 0.0d0
         return
       endif
c
       rs = (3*r**2/cdsum)**trd
c
       z = (cdd-cdu) / cdsum
       fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
       fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
c
c      exchange potential and energy density
c
       vxp = -3*alp/(pi*a0*rs)
       exp = 3*vxp/4
c
       vxf = 2**trd*vxp
       exf = 2**trd*exp
c
c   correlation potential and energy density
c
c      The Perdew-Zunger parameterization is used.
c      See Phys. Rev. B 23, 5048 (1981) - Appendix C (pg. 5075).
c
       if (rs .gt. 1.D0) then
c
       te = 1.D0+(7.D0/6.D0)*1.0529D0*sqrt(rs)+(4.D0/3.D0)*0.3334D0*rs
       be = 1.D0+1.0529D0*sqrt(rs)+0.3334D0*rs
       ecp = -0.2846D0/be
       vcp = -0.2846D0*te/be**2
       te = 1.D0+(7.D0/6.D0)*1.3981D0*sqrt(rs)+(4.D0/3.D0)*0.2611D0*rs
       be = 1.D0+1.3981D0*sqrt(rs)+0.2611D0*rs
       ecf = -0.1686D0/be
       vcf = -0.1686D0*te/be**2
c
       else
c
       ecp = +2*((0.0311D0+0.0020D0*rs)*log(rs)
     1      - 0.048D0-0.0116D0*rs)
       vcp = +2*((0.0311D0+2.D0/3.D0*0.0020D0*rs)*log(rs)
     1      - (0.048D0+0.0311D0/3.D0)
     2      - (2.D0/3.D0*0.0116D0+0.0020D0/3.D0)*rs)
       ecf = +2*((0.01555D0+0.0007D0*rs)*log(rs)
     1      - 0.0269D0-0.0048D0*rs)
       vcf = +2*((0.01555D0+2.D0/3.D0*0.0007D0*rs)*log(rs)
     1      - (0.0269D0+0.01555D0/3.D0)
     2      - (2.D0/3.D0*0.0048D0+0.0007D0/3.D0)*rs)
c
       endif
c
       vxcp = vxp + vcp
       vxcf = vxf + vcf
       vxcd = vxcp
       vxcu = vxcp
       excp = exp + ecp
       excf = exf + ecf
       vcd = vcp
       vcu = vcp
       exct = excp
       ect = ecp
       if (z .ne. 0.D0) then
         vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
         vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
         vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
         vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
         exct = exct + fz*(excf-excp)
         ect = ect + fz*(ecf-ecp)
       endif
c
       return
       end
