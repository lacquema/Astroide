c*************************************************************************
c                        KICKVH_TID.F
c*************************************************************************
c To kick the velocity components rot(*) by drot(*)*dt 
c
c             Input:
c                 nbod          ==>  number of bodies (int scalar)
c                 rotx,roty,rotz   ==>  initial spin vectors
c                                    (real arrays)
c                 drotx,droty,drotz   ==>  acceleration in helio coord
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 rotx,roty,rotz   ==>  final spin vectors 
c                                    (real arrays)
c
c     ALGORITHM: Obvious  
c       
c     AUTHOR:  H. Beust
c     DATE WRITTEN:  Jun 12, 2008
c     REmark : diff with kickvh = we stars at n=1 (not 2)

      subroutine kickvh_tid(nbod,rotx,roty,rotz,drotx,droty,drotz,dt) 


      include '../swift.inc'

c...  Inputs Only: 
	integer nbod
	real*8 drotx(nbod),droty(nbod),drotz(nbod)
	real*8 dt

c...   Inputs and Output:
	real*8 rotx(nbod),roty(nbod),rotz(nbod)

c...  Internals:
	integer n

c----
c...  Executable code 

	do n= 1, nbod
	   rotx(n) = rotx(n) + drotx(n)*dt
	   roty(n) = roty(n) + droty(n)*dt
	   rotz(n) = rotz(n) + drotz(n)*dt
	enddo

        return
        end    ! kickvh_tid
c-----------------------------------------------------------------------------