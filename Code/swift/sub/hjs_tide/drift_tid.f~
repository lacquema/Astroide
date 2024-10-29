c*************************************************************************
c                        DRIFT_TID.F
c*************************************************************************
c This subroutine moves the ROTATION AXES of the solid bodies (averaged)
c
c                Source for formulas : Mardling & LIn, 2002, ApJ 573, 829 
c                         + rotation fomulas : rotation.mws
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 tid         ==> 0 (no tide), 1 (tide) or 2 (averaged tide)  
c                 oloc        ==> Link between bodies and orbits
c                                  (2D integer array)
c                 inert       ==> Moments of inertia (real array)
c                 krpl5       ==> k(apsidal consants) * r(radius)^5
c                 Qlov        ==> Q constant's (array)
c                 eta         ==> Masses of centers for orbits (real array)
c                 mu          ==> Masses of satellites for orbits (real arr.)
c                 mass        ==>  mass of bodies (real array)
c                 xj,yj,zj    ==>  position in jacobi coord (real arrays)
c                 vxj,vyj,vzj ==>  Velocity in Jacobi coord (real arrays)
c                 dt          ==> Timestep
c                 rotx,roty,rotz ==> Spin vectors of the bodies (real arrays)
c             Output:
c                 rotx,roty,rotz => Rotation vectors of bodies
c
c Author:  H. Beust  
c Date:    Jun 23, 2008
c Modified : Oct 20, 2008

      subroutine drift_tid(nbod,tid,oloc,mass,inert,krpl5,Qlov,
     &     eta,mu,xj,yj,zj,vxj,vyj,vzj,dt,rotx,roty,rotz)

      include '../swift.inc'

c...  Inputs: 
      integer nbod
      integer oloc(NPLMAX,NPLMAX)
      integer tid(nbod)
      real*8 inert(nbod),krpl5(nbod),Qlov(nbod)
      real*8 mass(nbod),eta(nbod),mu(nbod)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 dt

c...  Inputs & outputs
      real*8 rotx(nbod),roty(nbod),rotz(nbod)

c...  Internals:
      integer i,j,k,ialpha
      real*8 rr2(NPLMAX),irj(NPLMAX)
      real*8 accijz(NPLMAX,NPLMAX),rij2,irij3,a,e,q
      real*8 f1,f2,ux,uy,uz,nn,ebt,cw,sw,inc,capom,omega,capm
      real*8 na,exq,exq9,exq12,ex2,e4,e6,a2,a3,a4
      real*8 exq2,exq3,c,b,alpha,delta,r3m,w,w2,xx,ff
      real*8 ep,ep2,ep3,ep6,q2,q3,q4,q7
      real*8 e1(3),e2(3),e3(3),r1(2),r2(2),r3(2),rot2
      real*8 hh5,hh6,hhb,hhc,r1f,r2f,r3f
      real*8 vdotr,rx,ry,vx,vy,vz,cu,su

c----
c...  Executable code 

c...  get thr r^-3's

      do i = 2,nbod
         rr2(i) = xj(i)*xj(i) + yj(i)*yj(i) + zj(i)*zj(i)
         irj(i) = 1.d0/sqrt(rr2(i))
      end do
      
      do k=2,nbod

c...   Now compute tidal terms for orbit k
c...   i and j will be satellite & center for this orbit. We assume 
c...   there is only i and j
c...   Source Mardling & Lin, 2002, ApJ 573, 829

c        tid(k)=0
        if (tid(k).gt.0) then
          i = 1
          do while (oloc(k,i).ne.1)
            i = i + 1
          end do
          j = 1
          do while (oloc(k,j).ne.-1)
            j = j + 1
          end do
c 	  call orbel_xv2el(xj(k),yj(k),zj(k),vxj(k),vyj(k),vzj(k),
c     &                 mu(k)+eta(k),ialpha,a,e,inc,capom,omega,capm)
 
          call orbel_xv2aeq(xj(k),yj(k),zj(k),vxj(k),vyj(k),vzj(k),
     &                 mu(k)+eta(k),ialpha,a,e,q) 
          na = sqrt((eta(k)+mu(k))/a)
          nn = na/a
          ux = xj(k)*irj(k)
          uy = yj(k)*irj(k)
          uz = zj(k)*irj(k)
          exq2 = 1.d0-e*e
          exq = sqrt(exq2)
          exq3 = exq2*exq
          exq9 = exq3*exq3*exq3
          exq12 = exq9*exq3
          ex2 = e*e
          e4 = ex2*ex2
          e6 = e4*ex2
          ep = 1.d0+e
          ep2 = ep*ep
          ep3 = ep2*ep
          ep6 = ep3*ep3

          q2 = q*q
          q3 = q2*q
          q4 = q2*q2
          q7 = q4*q3
            
          e3(1) = (yj(k)*vzj(k)-zj(k)*vyj(k))*exq3/(nn*q2*ep2)
          e3(2) = (zj(k)*vxj(k)-xj(k)*vzj(k))*exq3/(nn*q2*ep2)
          e3(3) = (xj(k)*vyj(k)-yj(k)*vxj(k))*exq3/(nn*q2*ep2)

          if (ex2.ge.TINY) then
            vdotr = xj(k)*vxj(k)+yj(k)*vyj(k)+zj(k)*vzj(k)
            cu = (1.d0-1.d0/a/irj(k))/e
            if (cu.gt.1.d0) cu = 1.d0
            if (cu.lt.-1.d0) cu =-1.d0
            su = sqrt(1.d0-cu*cu)
            if (vdotr.lt.0.d0) su = -su
            vx = e3(2)*uz-e3(3)*uy
            vy = e3(3)*ux-e3(1)*uz
            vz = e3(1)*uy-e3(2)*ux
            ry = su*exq*a*irj(k)
            rx = (cu-e)*a*irj(k)
            e1(1) = -ry*vx+rx*ux
            e1(2) = -ry*vy+rx*uy
            e1(3) = -ry*vz+rx*uz

c            e1(1) = ((vyj(k)*e3(3)-vzj(k)*e3(2))*exq/na-ux)/e
c            e1(2) = ((vzj(k)*e3(1)-vxj(k)*e3(3))*exq/na-uy)/e
c            e1(3) = ((vxj(k)*e3(2)-vyj(k)*e3(1))*exq/na-uz)/e
c            e1(1) = cos(omega)*cos(capom)-sin(omega)*cos(inc)*sin(capom)
c            e1(2) = cos(omega)*sin(capom)+sin(omega)*cos(inc)*cos(capom)
c            e1(3) = sin(omega)*sin(inc)
c            print*,'drift',e,e1(1)**2+e1(2)**2+e1(3)**2,cu
          else
            e1(1) = ux
            e1(2) = uy
            e1(3) = uz
          end if
       
          e2(1) = e3(2)*e1(3)-e3(3)*e1(2)
          e2(2) = e3(3)*e1(1)-e3(1)*e1(3)
          e2(3) = e3(1)*e1(2)-e3(2)*e1(1)

          hh5 = exq3/ep3

          hh6 = (16.d0+120.d0*ex2+90.d0*e4+5.d0*e6)/16.d0/ep6
          hhb = -(8.d0+24.d0*ex2+3.d0*e4)/16.d0*exq3/ep6
          hhc =  ex2*(6.d0+ex2)/8.d0*exq3/ep6

c........
          r1(1) = rotx(j)*e1(1)+roty(j)*e1(2)+rotz(j)*e1(3)
          r2(1) = rotx(j)*e2(1)+roty(j)*e2(2)+rotz(j)*e2(3)
          r3(1) = rotx(j)*e3(1)+roty(j)*e3(2)+rotz(j)*e3(3)
          r1(2) = rotx(i)*e1(1)+roty(i)*e1(2)+rotz(i)*e1(3)
          r2(2) = rotx(i)*e2(1)+roty(i)*e2(2)+rotz(i)*e2(3)
          r3(2) = rotx(i)*e3(1)+roty(i)*e3(2)+rotz(i)*e3(3)

          f1 = -(eta(k)*mu(k))/(eta(k)+mu(k))
c.......  Rotation center
          f2 = (1.d0+mu(k)/eta(k))*krpl5(j)/q3

          alpha = f1*hh5*f2/inert(j)

          f2 = -6.d0*(eta(k)+mu(k))/nn
     &                   *krpl5(j)/Qlov(j)*mu(k)/eta(k)/q7
          b = f1*f2*hhb*q/inert(j)
          c = f1*f2*hhc*q/inert(j)
          delta = 0.5d0*nn*hh6/hhb
          ebt = exp(b*dt)
          xx = 2.d0*b*dt
          if (abs(xx).lt.1d-2) then      ! Pade approximant
            ff =  -xx*(60.d0+xx*(-10.d0+xx))
     &                      /(-120.d0+xx*(60.d0+xx*(-12.d0+xx)))
          else
            ff = (exp(xx)-1.d0)/xx-1.d0
          end if
          r3m = delta*ff+r3(1)*(ff+1.d0)
c          r3m = -(delta+r3(1))*(1.d0-ebt*ebt)/(2.d0*b*dt)-delta
          w2 = -alpha*alpha*r3m*r3m+c*c
          if (w2.lt.-tiny) then
            w = sqrt(-w2)
            cw = cos(w*dt)
            sw = sin(w*dt)
            r1f = (cw*r1(1)+(c*r1(1)-alpha*r3m*r2(1))*sw/w)*ebt
            r2f = (cw*r2(1)-(c*r2(1)-alpha*r3m*r1(1))*sw/w)*ebt
c            r3f = ebt*ebt*(delta+r3(1))-delta
            r3f = delta*(ff+1.d0)*xx+r3(1)*((ff+1.d0)*xx+1.d0)
          else if (w2.gt.tiny) then
            w = sqrt(w2)
            cw = cosh(w*dt)
            sw = sinh(w*dt)
            r1f = (cw*r1(1)+(c*r1(1)-alpha*r3m*r2(1))*sw/w)*ebt
            r2f = (cw*r2(1)-(c*r2(1)-alpha*r3m*r1(1))*sw/w)*ebt
            r3f = delta*(ff+1.d0)*xx+r3(1)*((ff+1.d0)*xx+1.d0)
c            r3f = ebt*ebt*(delta+r3(1))-delta
          else
            r1f = (r1(1)+dt*r1(1)*c-dt*r2(1)*alpha*r3m)*ebt
            r2f = (r2(1)+dt*r1(1)*alpha*r3m-dt*r2(1)*c)*ebt
            r3f = delta*(ff+1.d0)*xx+r3(1)*((ff+1.d0)*xx+1.d0)
c            r3f =  ebt*ebt*(delta+r3(1))-delta
          end if

          rotx(j) = r1f*e1(1)+r2f*e2(1)+r3f*e3(1)
          roty(j) = r1f*e1(2)+r2f*e2(2)+r3f*e3(2)
          rotz(j) = r1f*e1(3)+r2f*e2(3)+r3f*e3(3)

c.......  Rotation satellite
          f2 = (1.d0+eta(k)/mu(k))*krpl5(i)/q3

c

c          fx = -f1*r3(1)*r2(1)*hh5*f2/inert(j) 
c          fy = +f1*r3(1)*r1(1)*hh5*f2/inert(j)
c          fz = 0.d0
c            f2 = -6.d0*nn*krpl5(j)/Qlov(j)*mu(k)/eta(k)/a4
c            fz = fz + f1*f2*hh6*na/inert(j)
c            fx = fx + f1*f2*r1(1)*hh7*a/inert(j)
c            fy = fy + f1*f2*r2(1)*hh8*a/inert(j)
c            fz = fz + f1*f2*r3(1)*hh9*a/inert(j)

c         fx = dr1/dt = -alpha*r2*r3+(b+c)*r1
c         fy = dr2/dt = alpha*r1*r3+(b-c)*r2
c         fz = dr3/dt = 2*b*(r3+delta) => r3 = -delta+(delta+r30)*e^(2*b*t)

          alpha = f1*hh5*f2/inert(i)
          f2 = -6.d0*(eta(k)+mu(k))/nn
     &                   *krpl5(i)/Qlov(i)*eta(k)/mu(k)/q7
          b = f1*f2*hhb*q/inert(i)
          c = f1*f2*hhc*q/inert(i)
          delta = 0.5d0*nn*hh6/hhb
          ebt = exp(b*dt)
          xx = 2.d0*b*dt
          if (abs(xx).lt.1d-2) then      ! Pade approximant
            ff =  -xx*(60.d0+xx*(-10.d0+xx))
     &                      /(-120.d0+xx*(60.d0+xx*(-12.d0+xx)))
          else
            ff = (exp(xx)-1.d0)/xx-1.d0
          end if
          r3m = delta*ff+r3(2)*(ff+1.d0)
c          r3m = -(delta+r3(2))*(1.d0-ebt*ebt)/(2.d0*b*dt)-delta
          w2 = -alpha*alpha*r3m*r3m+c*c
          if (w2.lt.-tiny) then
c            print*,'cos'
            w = sqrt(-w2)
            cw = cos(w*dt)
            sw = sin(w*dt)
            r1f = (cw*r1(2)+(c*r1(2)-alpha*r3m*r2(2))*sw/w)*ebt
            r2f = (cw*r2(2)-(c*r2(2)-alpha*r3m*r1(2))*sw/w)*ebt
            r3f = delta*(ff+1.d0)*xx+r3(2)*((ff+1.d0)*xx+1.d0)
c            r3f = ebt*ebt*(delta+r3(2))-delta
          else if (w2.gt.tiny) then
c            print*,'cosh'
            w = sqrt(w2)
            cw = cosh(w*dt)
            sw = sinh(w*dt)
            r1f = (cw*r1(2)+(c*r1(2)-alpha*r3m*r2(2))*sw/w)*ebt
            r2f = (cw*r2(2)-(c*r2(2)-alpha*r3m*r1(2))*sw/w)*ebt
            r3f = delta*(ff+1.d0)*xx+r3(2)*((ff+1.d0)*xx+1.d0)
c            r3f = ebt*ebt*(delta+r3(2))-delta
          else
c            print*,'zero'
            r1f = (r1(2)+dt*r1(2)*c-dt*r2(2)*alpha*r3m)*ebt
            r2f = (r2(2)+dt*r1(2)*alpha*r3m-dt*r2(2)*c)*ebt
            r3f = delta*(ff+1.d0)*xx+r3(2)*((ff+1.d0)*xx+1.d0) 
c           r3f = ebt*ebt*(delta+r3(2))-delta
          end if

          rotx(i) = r1f*e1(1)+r2f*e2(1)+r3f*e3(1)
          roty(i) = r1f*e1(2)+r2f*e2(2)+r3f*e3(2)
          rotz(i) = r1f*e1(3)+r2f*e2(3)+r3f*e3(3)
        end if
      end do

      return
      end      ! drift_tid

c---------------------------------------------------------------------




