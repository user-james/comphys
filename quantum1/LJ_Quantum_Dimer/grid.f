c *******************************************************************
       subroutine grid(rmax, aa, bb, nrmax, r, rab, nr, a, b)
c *******************************************************************
c
c  This subroutine sets up a radial grid of points which are roughly
c  logarithmically distributed in r. This means that the separation
c  of grid of points is very fine near r=0 and coarse at large radii.
c
c  The grid points are:
c
c        r(i) = a * [ exp(b*(i-1)) -1 ]     for i=1,nr
c
c  This strategy is used because atomic potentials are typically much
c  more rapidly varying for small radii than for large radii. It is
c  a crude example of an "adaptive grid" approach.
c
c  Input Parameters:
c
c    rmax  =  largest value of r(i)  (roughly - in fact, largest r(i) .le. rmax)
c             Note: input rmax=0 gives default rmax = 80
c    aa, bb = parameters determining grid spacing (see below)
c             Note: input aa=0 and bb=0 gives default values of aa=6 and bb=40
c    nrmax  = array limits on radial array
c
c  Output Parameters:
c
c    nr  =   number of points on grid
c    r(i)  =  ith point on grid
c    rab(i) =  dr/di  at ith grid point
c              (dr/di is calculated from above expression for r(i),
c               treating i as a real number in definition of r(i) )
c               Note: rab(i) is roughly r(i+1)-r(i)  
c
c    a,b  =  grid spacing parameters used elsewhere
c 
c *******************************************************************
       implicit none
c
       integer i, nr, nrmax
c 
c
       real*8 r(*), rab(*)
       real*8 rmax, aa, bb, a, b
c
c   set default values
c
       if (rmax .eq. 0.D0) rmax=80.D0   
       if (aa .eq. 0.D0) aa=6.D0
       if (bb .eq. 0.D0) bb=40.D0
c
c   output grid parameters
c
       a = exp(-aa)     
       b = 1/bb
c
c   set up grid
c
       i = 1
       r(1) = 0.0d0
       do while (i .lt. nrmax .and. r(i) .le. rmax)
         i = i + 1
         r(i) = a*(exp(b*(i-1))-1)
         rab(i) = (r(i)+a)*b
       enddo
c
       if ( i .eq. nrmax ) then
         write(6,50)
 50      format(/,38h arraylimits for radial array exceeded,/)
         stop
       endif
c
       nr = i-1
c
       return
       end
