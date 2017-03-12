c
      double precision function dran1(iseed)
c
c Pseudo-random number generator by D. R. Hamann, AT&T Bell
c Laboratories, 1991, based in part on
c ran3 and ran2 from W. H. Press, et al.,"Numerical Recipes" 
c (Cambridge U. Press, Cambridge, 1986), p. 199,
c and on Table 1 and Algorithm B from D. E. Knuth,
c "Seminumerical Algorithms", 2nd ed. Vol 2 
c (Addison-Wesley, Reading, Mass., 1981), pp28 and 32.
c
c Don't even think of using this routine without testing it
c thoroughly yourself.  See Knuth or M. H. Kalos and P. A. Whitlock,
c "Monte Carlo Methods," Vol. 1, ( John Wiley, New York, 198?),
c Appendix A.
c
c Don't make single precision.
c
c USEAGE:
c       CALL WITH ISEED > 0 TO INITIALIZE
c       CALL WITH ISEED = 0 FOR SUCCESSIVE RANDOM NUMBERS
c
      implicit none
      integer m,ia,ic,ishuf,kc,lc
c these linear congruential "start-up" generator parameters
c shouldn't have much to do with the performance of the real generator
      parameter (m = 714025, ia = 1366, ic = 150889)
      parameter (kc = 55, lc = 31)
c any ishuf value should be ok here, but a prime several times larger
c than kc seems like it should be best
      parameter (ishuf = 499)
      save iff,da,inext,inextp,dshuf,y,v,kc1
      integer iseed,i,j,ii,inext,inextp,iff,ix,kc1
      double precision dmi
      double precision dj,dk,da(kc)
      double precision dshuf,y,v(ishuf)
      data iff /0/
      if(iseed .gt. 0 .or. iff .eq. 0)then
        iff=1
c
c fill additive generator array using linear congruential generator
c 3 times per number to get plenty of digits
c
        ix = mod(iabs(iseed), m)
        dmi = 1.0d0 / dble(m)
        kc1 = kc + 1
        do i = 1,kc
          ix = mod(ia * ix + ic, m)
          y = dmi * dble(ix)
          ix = mod(ia * ix + ic, m)
          y = dmi * (dble(ix) + y)
          ix = mod(ia * ix + ic, m)
          y = dmi * (dble(ix) + y)
          da(i) = y
        end do
c
c use all given digits of iseed to modify first entry of additive
c generator
c
        y = dble(iabs(iseed))
10      y = 0.1d0 * y
        if(y .ge. 1.0d0) go to 10
        dj=da(1) - y
        if(dj .lt. 0.0d0) dj = dj + 1.0d0
        da(1) = dj
c
c warm up knuth subtractive generator and fill shuffle table
c
        inext=0
        inextp=lc
        y = da(kc)
        dshuf = dble(ishuf)
c
        do i = 1, 10 * ishuf
          v(mod(i - 1, ishuf) + 1) = y
          inext = inext + 1
          if(inext .eq. kc1) inext = 1
          inextp = inextp + 1
          if(inextp .eq. kc1) inextp = 1
          dj = da(inext) - da(inextp)
          if(dj .lt. 0.0d0) dj = dj + 1.0d0
          da(inext) = dj
          y = dj
        end do
c
      endif
c
c entry point except on first call
c main Knuth subtractive generator plus algorithm B shuffle
c
      j = idint(dshuf * y) + 1
      y = v(j)
      dran1 = y
      inext = inext + 1
      if(inext .eq. kc1) inext = 1
      inextp = inextp + 1
      if(inextp .eq. kc1) inextp = 1
      dj = da(inext) - da(inextp)
      if(dj .lt. 0.0d0) dj = dj + 1.0d0
      da(inext) = dj
      v(j) = dj
      return
      end
c
