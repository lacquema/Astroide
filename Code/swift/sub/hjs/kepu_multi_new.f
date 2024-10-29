c*************************************************************************
c                        KEPU_MULTI_NEW.F
c*************************************************************************
c This subroutine solves kepler's equation in universal variables for
c all the orbits simultaneously in HJS, with Poincare time regularization
c         Source : Mikkola, 1997, CeMDA 67, 147
c It also gives the timestep dt
c using NEWTON'S METHOD
c
c             Input:
c                 n             ==>  size of the array (numbre of orbits)
c                 h             ==> Universal variable step (ds)    
c                 gm            ==>  dynamical masses of the orbits
c                 coef          ==> Coefs of the time transformation
c                 x0,y0,z0,     ==> Initial positions and velocities
c                 vx0,vy0,vz0        of the orbits.
c             Output:
c                 dt            ==>  time step (real scalor)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c Author:  Hervé Beust
c Date:    Sep 22, 2006
c Based on drift_kepu_new.f

      subroutine kepu_multi_new(n,h,gm,coef,x0,y0,z0,vx0,vy0,vz0,dt,
     &                                                         iflgn)

      include '../swift.inc'

c...  Inputs: 
      integer n
      real*8 h,gm(n),coef(n)

c...  Inputs and Outputs 
      real*8 x0(n),y0(n),z0(n),vx0(n),vy0(n),vz0(n)

c...  Outputs:
      real*8 dt
      integer iflgn

c...  Internals:
      integer nc,i,j
      logical ok,vb
      real*8 xx,ss,v0s,f,g,fdot,gdot
      real*8 r0(NPLMAX),u(NPLMAX),alpha(NPLMAX)
      real*8 c0(NPLMAX),c1(NPLMAX),c2(NPLMAX),c3(NPLMAX)
      real*8 x(NPLMAX),y(NPLMAX),z(NPLMAX)
      real*8 vx(NPLMAX),vy(NPLMAX),vz(NPLMAX)
      real*8 ufp(NPLMAX),cufp(NPLMAX)
      real*8 s(0:NPLMAX),ff(0:NPLMAX),ds(0:NPLMAX)
      real*8 ujac(0:NPLMAX,0:NPLMAX)
c      real*8 prod(0:NPLMAX,0:NPLMAX),jac(0:NPLMAX,0:NPLMAX)
c----
c...  Executable code 

c...  Compute distances, energies, and angular momenta
c      vb = .false.
c      if (iflgn.eq.10) vb= .true.
      do i=1,n
        r0(i) = sqrt(x0(i)*x0(i) + y0(i)*y0(i) + z0(i)*z0(i))
        v0s = vx0(i)*vx0(i) + vy0(i)*vy0(i) + vz0(i)*vz0(i)
        u(i) = x0(i)*vx0(i) + y0(i)*vy0(i) + z0(i)*vz0(i)
        alpha(i) = 2.0d0*gm(i)/r0(i) - v0s
      end do

      dt = 0.d0

c...  Compute initial guess
      ss = 0.d0
      do i=1,n
        ss = ss+coef(i)/r0(i)
      end do
      s(0) = h/ss
      do i=1,n
        call drift_kepu_guess(s(0),r0(i),gm(i),alpha(i),u(i),s(i))
      end do
         
c... Here is the big loop
      ok = .false.
      iflgn = 1
      nc = 0
c      if (vb) print*,'*********************'
      do while (.not.ok)
        ff(0) = -h
        do i=1,n
          ff(0) = ff(0) + coef(i)*s(i)
        end do
        do i=1,n
          xx = s(i)*s(i)*alpha(i)
          call drift_kepu_stumpff(xx,c0(i),c1(i),c2(i),c3(i))
          c1(i) = c1(i)*s(i) 
          c2(i) = c2(i)*s(i)*s(i) 
          c3(i) = c3(i)*s(i)*s(i)*s(i)
          ff(i) = r0(i)*c1(i) + u(i)*c2(i) + gm(i)*c3(i) - s(0)
          ufp(i) = 1.d0/(r0(i)*c0(i) + u(i)*c1(i) + gm(i)*c2(i))
          cufp(i) = coef(i)*ufp(i)
        end do
        ss = 0.d0
        do i=1,n
          ss = ss+cufp(i)
        end do
        ss = 1.d0/ss

c...   Build inverse Jacobian matrix
        ujac(0,0) = ss
        do i=1,n
          ujac(i,0) = ufp(i)*ss
        end do
        do j=1,n
          ujac(0,j) = -ss*cufp(j)
        end do
        do i=1,n
          do j=1,n
            ujac(i,j) = -ufp(i)*ss*cufp(j)
          end do
          ujac(i,i) = ujac(i,i) + ufp(i)
        end do 

c        jac = 0.d0
c        do i=1,n
c          jac(i,0) = -1.d0
c          jac(i,i) = 1.d0/ufp(i)
c          jac(0,i) = coef(i)
c        end do
c
c        prod(0:n,0:n) = matmul(jac(0:n,0:n),ujac(0:n,0:n))
c        print*,'****************'
c        do i=0,n
c          print*,sngl(prod(i,0:n))
c        end do
c        print*,'****************'

c...    Compute displacement
        do i=0,n
          ds(i) = 0.d0
          do j=0,n
            ds(i) = ds(i)-ujac(i,j)*ff(j)
          end do
        end do
        do i=0,n
          s(i) = s(i) + ds(i)
        end do
        ok = .true.
c       if (vb) print*,sngl(ff(0:n)*h)
        do i=0,n
          ok = (ok.and.(abs(ff(i))*h.lt.DANBYB))
        end do         
        nc = nc+1
        if (ok) iflgn = 0
        ok = (ok.or.(nc.gt.100))
      end do
c      if (vb) print*,'-------------------------------'

c...  Move the orbits
      if (iflgn.eq.0) then
        dt = s(0)
        do i=1,n
          f = 1.0d0 - (gm(i)/r0(i))*c2(i)
          g = s(0) - gm(i)*c3(i)
          fdot = -(gm(i)*ufp(i)/r0(i))*c1(i)
          gdot = 1.d0 - (gm(i)*ufp(i))*c2(i)

           x(i) = x0(i)*f + vx0(i)*g
           y(i) = y0(i)*f + vy0(i)*g
           z(i) = z0(i)*f + vz0(i)*g
           vx(i) = x0(i)*fdot + vx0(i)*gdot
           vy(i) = y0(i)*fdot + vy0(i)*gdot
           vz(i) = z0(i)*fdot + vz0(i)*gdot

           x0(i) = x(i)
           y0(i)  = y(i)
           z0(i) = z(i)
           vx0(i) = vx(i)
           vy0(i) = vy(i)
           vz0(i) = vz(i)
         end do
       end if

       return

       end  ! kepu_multi_new.f
c----------------------------------------------------------------------
