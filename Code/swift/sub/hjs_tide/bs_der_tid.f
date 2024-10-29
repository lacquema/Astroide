C**************************************************************************
C	    		        BS_DER_TID
C**************************************************************************
c This is the subroutine that calculates the tidal derivatives 
c
c             Input:
c              ff(10)  ==> various constants
c              ybs(13)   ==> Variables 
c
c             Output:
c 	         dy  ==> derivatives of the independant var (real array)
c
c Authors:  Herve Beust
c Date:    2/17/2009 (Based on bs_der)
c

      subroutine bs_der_tid(ff,ybs,dy)
       
      include '../swift.inc'

c...  Inputs Only: 
      real*8 ff(10)
      real*8 ybs(13)

c...  Output
      real*8 dy(13)

c...  Internals:
      integer i,j,k,kk,ialpha,iflg
      real*8 a,e,q,gm,f2m2,f2m1
      real*8 f1,f2,f3,alpha,ux,uy,uz,nq,nn
      real*8 f2f3,f4,f5,f6,f7,f8,f9,f10,f11
      real*8 exq,exq2,exq3,ex4,ex6,ex2,exq4
      real*8 ep,ep2,ep3,ep5,ep6,ep12,ep52,ep72,ep92,ep132
      real*8 q2,q3,q5,q8
      real*8 gqdx(2),gqdy(2),gqdz(2),gtfx(2),gtfy(2),gtfz(2)
      real*8 dqdx(2),dqdy(2),dqdz(2),dtfx(2),dtfy(2),dtfz(2)
      real*8 grelx,grely,grelz
      real*8 e1(3),e2(3),e3(3),rx(2),ry(2),rz(2)
      real*8 vecth(3),vecte(3),h
      real*8 hh1,hh2,hh3,hh4,hh5,hh6,hh7
      real*8 vdotr,posx,posy,vx,vy,vz,cu,su,capm
      real*8 gtotx,gtoty,gtotz,dtotx,dtoty,dtotz,mrel,dcapm
      real*8 rotx(2),roty(2),rotz(2),mqd(2) 
      real*8 dvecte(3),dvecth(3),drotx(2),droty(2),drotz(2)
      real*8, parameter :: 
     &           cl = 2.9979248d8*86400.d0*365.2422d0/1.4959787061d11
      real*8, parameter :: cl2 = cl*cl

c----
c...  Executable code 

c...  move things so that I can deal with it
      vecte(1:3) = ybs(1:3)
      vecth(1:3) = ybs(4:6)
      capm = ybs(7)
      rotx(1) = ybs(8)
      roty(1) = ybs(9)
      rotz(1) = ybs(10)
      rotx(2) = ybs(11)
      roty(2) = ybs(12)
      rotz(2) = ybs(13)

      gm = ff(7)
      e = sqrt(vecte(1)*vecte(1) + vecte(2)*vecte(2)
     &                                  + vecte(3)*vecte(3))
      h = sqrt(vecth(1)*vecth(1) + vecth(2)*vecth(2)
     &                                  + vecth(3)*vecth(3))
      q = h*h/gm/(1.d0+e)
      q2 = q*q
      q3 = q2*q
      q5 = q3*q2
      q8 = q5*q3
c      print*,'yb',ybs

      a = q/(1.d0-e)
      nn = sqrt(gm/(a*a*a))
      nq = sqrt(gm/q3)

      ex2 = e*e
      exq2 = 1.d0-ex2
      exq = sqrt(exq2)
      exq3 = exq2*exq
      exq4 = exq2*exq2
      ex4 = ex2*ex2
      ex6 = ex4*ex2
      ep = 1.d0+e
      ep12 = sqrt(ep)
      ep2 = ep*ep
      ep3 = ep2*ep
      ep5 = ep2*ep3
      ep6 = ep3*ep3 
      ep52 = ep2*ep12
      ep72 = ep3*ep12
      ep92 = ep2*ep52
      ep132 = ep6*ep12

      e3(1) = vecth(1)/h
      e3(2) = vecth(2)/h
      e3(3) = vecth(3)/h

      if (e*e.gt.TINY) then
        e1(1) = vecte(1)/e
        e1(2) = vecte(2)/e
        e1(3) = vecte(3)/e
      else
        i = 1
        if (abs(vecth(2)).gt.abs(vecth(i))) i = 2
        if (abs(vecth(3)).gt.abs(vecth(i))) i = 3
        j = mod(i,3)+1
        k = mod(j,3)+1
        f2 = 1.d0/sqrt(1.d0-e3(k)*e3(k))
        e1(i) = -e3(j)*f2
        e1(j) = e3(i)*f2
        e1(k) = 0.d0
      end if 

      e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
      e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
      e2(3) = e3(1)*e1(2) - e3(2)*e1(1)
     
      hh1 = 30.d0*(1.d0+1.5d0*ex2+0.125d0*ex4)
      hh2 = 54.d0*(1.d0+3.75d0*ex2+1.875d0*ex4+0.078125d0*ex6)
      hh3 = 3.d0*(1.d0+1.5d0*ex2+0.125d0*ex4)
      hh4 = 6.d0*(1.d0+7.5d0*ex2+5.625d0*ex4+0.3125d0*ex6)
      hh5 = 3.d0*(1.d0+4.5d0*ex2+0.625d0*ex4)
      hh6 = 6.d0*(1.d0+3.d0*ex2+0.375d0*ex4)
      hh7 = 6.d0*(1.d0-4.5d0*ex2-0.875d0*ex4)

c........

      rx(1) = rotx(1)*e1(1) + roty(1)*e1(2) + rotz(1)*e1(3)
      ry(1) = rotx(1)*e2(1) + roty(1)*e2(2) + rotz(1)*e2(3)
      rz(1) = rotx(1)*e3(1) + roty(1)*e3(2) + rotz(1)*e3(3)
      rx(2) = rotx(2)*e1(1) + roty(2)*e1(2) + rotz(2)*e1(3)
      ry(2) = rotx(2)*e2(1) + roty(2)*e2(2) + rotz(2)*e2(3)
      rz(2) = rotx(2)*e3(1) + roty(2)*e3(2) + rotz(2)*e3(3)

      f4 = 1.d0/(nq*ep72*q5)
      f5 = 1.d0/(nq*ep132*q8)
      f7 = nq/(ep92*q3)
      f8 = e*hh3/(ep5*q5)
      f9 = nq*e*hh2/(ep132*q5)
      f10 = nq*nq*hh4/(ep6*q3)
      f11 = e*exq3*f4

      f2 = ff(1)
      f2m2 = ff(9)
      f3 = (rx(1)*rx(1)+ry(1)*ry(1)-2.d0*rz(1)*rz(1))
      f2f3 = f2*f3
      f6 = f2*exq3/(q3*ep3)

      gqdx(1) = 0.d0
      gqdy(1) = -0.5d0*f2f3*f11
      gqdz(1) = f2*ry(1)*rz(1)*f11
      gqdy(1) = gqdy(1) + f2m2*hh1*exq3*e*f5

      dqdx(1) = -f6*rz(1)*ry(1)
      dqdy(1) = f6*rz(1)*rx(1)
      dqdz(1) = 0.d0

      mqd(1) = -(0.5d0*f2f3*f4 + f2m2*hh7*f5)*exq4

      f2 = ff(2)
      f2m1 = ff(10)
      f3 = (rx(2)*rx(2)+ry(2)*ry(2)-2.d0*rz(2)*rz(2))
      f2f3 = f2*f3
      f6 = f2*exq3/(q3*ep3)

      gqdx(2) = 0.d0
      gqdy(2) = -0.5d0*f2f3*f11
      gqdz(2) = f2*ry(2)*rz(2)*f11
      gqdy(2) = gqdy(2) + f2m1*hh1*exq3*e*f5
          
      dqdx(2) = -f6*rz(2)*ry(2)
      dqdy(2) = f6*rz(2)*rx(2)
      dqdz(2) = 0.d0
       
      mqd(2) = -(0.5d0*f2f3*f4 + f2m1*hh7*f5)*exq4

      f2 = ff(3)

      gtfx(1) = -f2*f9
      gtfx(1) = gtfx(1) + 11.d0*rz(1)*f2*f8
      gtfz(1) = -rx(1)*f2*f8
      gtfy(1) = 0.d0

      dtfz(1) = -f10*f2
      dtfx(1) = f2*rx(1)*hh3*f7
      dtfy(1) = f2*ry(1)*hh5*f7
      dtfz(1) = dtfz(1) + f2*rz(1)*hh6*f7

      f2 = ff(4)

      gtfx(2) = -f2*f9    
      gtfx(2) = gtfx(2) + 11.d0*rz(2)*f2*f8
      gtfz(2) = -rx(2)*f2*f8
      gtfy(2) = 0.d0

      dtfz(2) = -f10*f2  
      dtfx(2) = f2*rx(2)*hh3*f7
      dtfy(2) = f2*ry(2)*hh5*f7
      dtfz(2) = dtfz(2) + f2*rz(2)*hh6*f7

      f4 = nq*nq*nq*q2/(cl2*ep52)

      grelx = 0.d0
      grely = 3.d0*f4*exq3*e
      grelz = 0.d0
        
      f2 = ff(5)

      mrel = -(15.d0-6.d0*exq-9.d0*f2+7.d0*exq*f2)*f4*exq4

c          print*,sngl(gqdx),sngl(dqdx)
c          print*,sngl(gqdy),sngl(dqdy)
c          print*,sngl(gqdz),sngl(dqdz)
c          print*,sngl(gtfx),sngl(dtfx)
c          print*,sngl(gtfy),sngl(dtfy)
c          print*,sngl(gtfz),sngl(dtfy)
c          print*,sngl(grelx),sngl(grely),sngl(grelz)
c
      gtotx = gqdx(1) + gqdx(2) + gtfx(1) + gtfx(2) + grelx
      gtoty = gqdy(1) + gqdy(2) + gtfy(1) + gtfy(2) + grely
      gtotz = gqdz(1) + gqdz(2) + gtfz(1) + gtfz(2) + grelz

      dtotx = dqdx(1) + dqdx(2) + dtfx(1) + dtfx(2)
      dtoty = dqdy(1) + dqdy(2) + dtfy(1) + dtfy(2)
      dtotz = dqdz(1) + dqdz(2) + dtfz(1) + dtfz(2)

      dcapm = mqd(1) + mqd(2) + mrel
c          print*,gtotx,e1(1),gtoty,e2(1),gtotz,e3(1)
c          print*,gtotx,e1(2),gtoty,e2(2),gtotz,e3(2)
c          print*,gtotx,e1(3),gtoty,e2(3),gtotz,e3(3)
          
      dvecte(1) = gtotx*e1(1) + gtoty*e2(1) + gtotz*e3(1)
      dvecte(2) = gtotx*e1(2) + gtoty*e2(2) + gtotz*e3(2)
      dvecte(3) = gtotx*e1(3) + gtoty*e2(3) + gtotz*e3(3)

      dvecth(1) = dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1)
      dvecth(2) = dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2)
      dvecth(3) = dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3)

      dtotx = dqdx(1) + dtfx(1)
      dtoty = dqdy(1) + dtfy(1)
      dtotz = dqdz(1) + dtfz(1)

      alpha = ff(6)
      drotx(1) = (dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1))*alpha
      droty(1) = (dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2))*alpha
      drotz(1) = (dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3))*alpha

      dtotx = dqdx(2) + dtfx(2)
      dtoty = dqdy(2) + dtfy(2)
      dtotz = dqdz(2) + dtfz(2)

      alpha = ff(8)
      drotx(2) = (dtotx*e1(1) + dtoty*e2(1) + dtotz*e3(1))*alpha
      droty(2) = (dtotx*e1(2) + dtoty*e2(2) + dtotz*e3(2))*alpha
      drotz(2) = (dtotx*e1(3) + dtoty*e2(3) + dtotz*e3(3))*alpha

      dy(1:3) = dvecte(1:3)
      dy(4:6) = dvecth(1:3)
      dy(7) = dcapm
      dy(8) = drotx(1)
      dy(9) = droty(1)
      dy(10) = drotz(1)
      dy(11) = drotx(2)
      dy(12) = droty(2)
      dy(13) = drotz(2)
c      print*,'dy',dy

      return
      end     ! bs_der_tid

c--------------------------------------------------------------------
