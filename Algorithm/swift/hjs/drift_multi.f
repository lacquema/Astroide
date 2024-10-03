c*************************************************************************
c                        DRIFT_MULTI.F
c*************************************************************************
c This subroutine does the danby-type drift for all orbits 
c in a hierarchical Genealized Jacobi case , with Poincare time regularization
c         Source : Mikkola, 1997, CeMDA 67, 147
c appropriate vbles and redoing a drift if the accuracy is too poor 
c (as flagged by the integer iflg).
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mu,eta        ==>  masses (sats & cens) of orbits
c                 coef          ==> The coefficients of time tarnsformation
c                 x,y,z         ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 vx,vy,vz      ==>  initial position in jacobi coord 
c                                    (real scalar)
c                 h             ==>  variable step ds
c             Output:
c                 x,y,z         ==>  final position in jacobi coord 
c                                       (real scalars)
c                 vx,vy,vz      ==>  final position in jacobi coord 
c                                       (real scalars)
c                 iflg          ==>  integer (zero for successful step)
c                 dt            ==>  time step
c
c Author: Hervé Beust
c Date:    Sep 22, 20006
c Last revision: 2/10/93
c Based on drift_hjs.f & drift_one.f

      subroutine drift_multi(nbod,coef,eta,mu,eps,
     &                           xj,yj,zj,vxj,vyj,vzj,h,dt)

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod 
      real*8 h,eta(nbod),mu(nbod),eps,coef(nbod)

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Output
      real*8 dt

c...  Internals:
      integer i,j,ialpha
      real*8 htmp,gm(NPLMAX),dttmp
      real*8 a,e,inc,capom,omega,capm,dm1(nplmax),dm2(nplmax)
      integer iflg

c----
c...  Executable code 

c...  Compute the modified central masses
c...          total mass = M = eta+mu
c...          reduced mass = m' = eta*mu/M
c...          GMeff = G*M + coef*eps/m'  (with G=1)   

        do i=2,nbod
          gm(i) = (eta(i)+mu(i))*(1.d0+coef(i)*eps/eta(i)/mu(i))
c          print*,gm(i)/(eta(i)+mu(i))
        end do

        iflg = 0
        do i=2,nbod
         call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
     &     gm(i),ialpha,a,e,inc,capom,omega,capm)
         if (e.gt.0.998) iflg=10
c           dm2(i) = capm
c           print*,'drift_multi orb ',i,sngl(a),sngl(e),
c     &      sngl(inc*180d0/PI),sngl(capom*180d0/PI),
c     &      sngl(omega*180d0/PI),sngl(capm*180d0/PI)
        end do


c...    Avance the orbits simultaneously and compute the time step dt
        call kepu_multi_new(nbod-1,h,gm(2:nbod),coef(2:nbod),
     &         xj(2:nbod),yj(2:nbod),zj(2:nbod),
     &         vxj(2:nbod),vyj(2:nbod),vzj(2:nbod),dt,iflg)
c        print*,'dt',dt
c        do i=2,nbod
c          call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
c     &     gm(i),ialpha,a,e,inc,capom,omega,capm)
c          dm2(i) = capm-dm2(i)
c          dm1(i) = dt*sqrt(gm(i)/a**3)
c          print*,'drift_multi orb ',i,'m=',capm*180./PI

c           print*,'drift_multi orb ap',i,sngl(a),sngl(e),
c     &      sngl(inc*180d0/PI),sngl(capom*180d0/PI),
c     &      sngl(omega*180d0/PI),sngl(capm*180d0/PI)
c          call drift_one(gm(i),xj(i),yj(i),zj(i),
c     &      vxj(i),vyj(i),vzj(i),dt,iflg)
c          call orbel_xv2el(xj(i),yj(i),zj(i),vxj(i),vyj(i),vzj(i),
c     &     eta(i)+mu(i),ialpha,a,e,inc,capom,omega,capm)
c           print*,'drift_multi orb enc',i,sngl(a),sngl(e),
c     &      sngl(inc*180d0/PI),sngl(capom*180d0/PI),
c     &      sngl(omega*180d0/PI),sngl(capm*180d0/PI)
c
c       end do

c        stop
c...   Redo it if it fails
	   if(iflg .ne. 0) then
	     do i = 1,10
c               print*,'gloub !!!!!!!!!!!!!!!!!!!!!!!!!!'
	       htmp = h/10.d0
               dt = 0.d0
               call kepu_multi_new(nbod-1,htmp,gm(2:nbod),
     &           coef(2:nbod),xj(2:nbod),yj(2:nbod),zj(2:nbod),
     &           vxj(2:nbod),vyj(2:nbod),vzj(2:nbod),dttmp,iflg)
               if(iflg .ne. 0) then
                 write(*,*) ' No convergence in Danby !!!!!!!!!'
                 write(*,*) dt
                 do j=2, nbod
                   write(*,*) xj(j),yj(j),zj(j)
                   write(*,*) vxj(j),vyj(j),vzj(j)
                 end do
                 write(*,*) ' STOPPING '
                 call util_exit(1)
               end if
               dt = dt + dttmp
             end do
	   endif

        return
        end    ! drift_multi.f
c-------------------------------------------------------------------
