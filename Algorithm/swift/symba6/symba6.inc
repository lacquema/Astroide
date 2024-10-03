c*************************************************************************
c                        SYMBA6.INC
c*************************************************************************
C Include file for the SyMBA subroutines
c
c Remarks:  copied from symba.inc
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 7/20/99, Search stuff added by Craig Angor

c...  Maximum number of encounters
	integer NENMAX
c	parameter(NENMAX=NTPMAX*100)   ! Must be less than 32767
	parameter(NENMAX=32766)   ! Must be less than 32767

c...  Maximum number of encounters per particle
        integer NENMAX1
        parameter(NENMAX1=100)  

c...	scale factor for hill's sphere to take shorter time step
       real*8 RHSCALE
       parameter (RHSCALE=6.5)

c...   Ratio of shell radii squared
	real*8 RSHELL
	parameter (RSHELL = 0.48075d0)   ! RSHELL ~ NTENC^(-2/3) 

c..    ratio of the number of time steps in the adjoining shells 
       integer NTENC
       parameter (NTENC=3)

c..    Maximum number of radial and azimuthal search zones
       integer NZNMAX
       parameter (NZNMAX=100)

c...   Maximum number of search sectors
       integer  NSECTMAX
       parameter (NSECTMAX = NZNMAX*NZNMAX)   ! max number of search zones 
c-------------------------------------------------------------------------

