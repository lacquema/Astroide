c*************************************************************************
c                        SWIFT.INC
c*************************************************************************
C Include file for SWIFT
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision:  3/11/2021 H Beust ! chlen

       implicit NONE    ! you got it baby

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Length of character lines
       integer*2, parameter :: chlen=150

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun 
       parameter  (NPLMAX = 11)   ! max number of planets, including the Sun 
       parameter  (NTPMAX = 200000) ! max number of test particles

c...   Size of the test particle integer status flag
	integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)  
	integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet  

c...   Size of the test particle integer status flag
	integer NSTATR    
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
	real*8 DANBYAC , DANBYB
	parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
	integer NLAG1, NLAG2
	parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
	real*8 TINY
        PARAMETER(TINY=4.D-15)
       
        integer INDICEMAX
        PARAMETER(INDICEMAX=5000)


c...    trig stuff
	real*8 PI,TWOPI,PIBY2,DEGRAD
	parameter (PI = 3.14159265358979D0)
	parameter (TWOPI = 2.0D0 * PI)
	parameter (PIBY2 = PI/2.0D0)
	parameter (DEGRAD = 180.0D0 / PI)

        integer NDT
        parameter (NDT=200)
        
        integer NB_DT
        parameter (NB_DT=200) 

     
c-------------------------------------------------------------------------
