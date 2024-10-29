c*************************************************************************
c                        HELIOQ.INC
c*************************************************************************
C Include file for the helioq integrator
c
c Author:  Hal Levison  
c Date:    2/2/2K
c Last revision: 

c...   Size of internal arrays
       integer N6DBS
       parameter(N6DBS=6*NTPMAX)
       integer NTEMP
       parameter(NTEMP=72*NTPMAX)

c...   maximum numbers of internations
       integer NTRYS
       parameter(NTRYS=11)

c...   maximum numbers of cnts
       integer NCNTS
       parameter(NCNTS=10000)

c...   Fraction increase in rtrans to start to do the BS
      real*8 OVERSIZE
      parameter(OVERSIZE=1.2)

c...   BS accuracy
      real*8 EPS
c      PARAMETER(EPS=4.0d-15)
      PARAMETER(EPS=1.0d-14)
c      PARAMETER(EPS=1.0d-15)

c...  The number of terms to use in a fit
      integer NUSE
      PARAMETER (NUSE=7)

c...  How to grow or shring the timestep
      real*8 SHRINK,GROW
      PARAMETER (SHRINK=.95D0,GROW=1.2D0)

c.... Number of Columns used by the extrapolation routines
      integer NCOL
c      PARAMETER (NCOL=7)
      PARAMETER (NCOL=NUSE)

c------------------------------------------------------------------------


