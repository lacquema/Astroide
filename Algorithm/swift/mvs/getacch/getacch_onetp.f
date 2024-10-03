c*************************************************************************
c                        GETACCH_ONETP.F
c*************************************************************************
c This subroutine calculates the acceleration on 1 test particle
c in the HELIOCENTRIC frame. 
c             Input:
c                  nbod        ==>  number of massive bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                  j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                  xh,yh,zh    ==>  massive part position in helio coord 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position in heliocentric coord 
c                                     (real scalars)
c                  istat       ==>  status of the test paricles
c                                      (integer scalar)
c                                      istat = 0 ==> active:  = 1 not
c             Output:
c               axht,ayht,azht ==>  tp acceleration in helio coord 
c                                   (real arrays)
c
c Author:  Hervé Beust
c Date:    Apr 20, 2023
c Adapted from getacch_tp.f

      subroutine getacch_onetp(nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     xht,yht,zht,istat,axht,ayht,azht)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod,istat
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod)
      real*8 xht(1),yht(1),zht(1),j2rp2,j4rp4

c...  Outputs:
      real*8 axht,ayht,azht
                
c...  Internals:
      integer i,ntp
      real*8 ir3h(NPLMAX),irh(NPLMAX),ir3ht(1),irht(1)
      real*8 axh0,ayh0,azh0,axh3,ayh3,azh3,aoblxt,aoblyt,aoblzt
      real*8 aoblx(NPLMAX),aobly(NPLMAX),aoblz(NPLMAX)

c----
c...  Executable code 
      ntp = 1

c...  get thr r^-3's  for the planets and test paricles
      call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
      call getacch_ir3(ntp,1,xht,yht,zht,ir3ht,irht)

c...  calc the ah0's:  recall that they are the same for all particles
      call getacch_ah0(2,nbod,mass,xh,yh,zh,ir3h,axh0,ayh0,azh0) 

c...  the first terms are 0

c...  the second terms are 0

c...  now the third terms
      call getacch_ah3_tp(nbod,1,mass,xh,yh,zh,xht,yht,zht,
     &                     istat,axh3,ayh3,azh3)

c...  add them all together
      if (istat.eq.0) then
            axht = axh0 + axh3
            ayht = ayh0 + ayh3
            azht = azh0 + azh3
      endif

      return
      end      ! getacch_onetp

c---------------------------------------------------------------------




