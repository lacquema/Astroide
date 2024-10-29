c*************************************************************************
c                            HELIOQ_MIDD.F
c*************************************************************************
c This is the Midpoing integrator that BS uses
c
c             Input:
c 	           y            ==> initial value dependent variables  
c                                     (real array)
c 	           dy           ==> initial value dy/dx
c                                     (real array)
c                 nbod          ==>  number of massive bodies (int scalar)
c                  xs           ==> initial value independent variable 
c                                        (real scalar)
c 	           htot         ==> total time to integrate  (real scalar)
c 	           nstep        ==> number of steps to take  (int scalar)
c                 iwhich        ==> determines which H to integrate
c                                       (int scalar)
c                 rtrans        ==>  solar encounter radii (real array)
c                 mass          ==>  mass of bodies (real array)
c             Output:
c 	          yout          ==> final value dependent variables  
c                                       (real array)
c
c
c Remarks: Based on numerical recipe routine midd.f
c Authors:  Hal Levison 
c Date:    2/8/2K
c Last revision: 

      subroutine helioq_mmid(y,dy,nbod,xs,htot,nstep,yout,iwhich,
     &     rtrans,mass)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer iwhich,nstep,nbod
      real*8 rtrans(2)
      real*8 y(N6DBS),xs,htot,dy(N6DBS)
      real*8 mass(NTPMAX)

c...  Inputs Only: 
      real*8 yout(N6DBS)

c...  Internals:
      integer nvar,i,n
      real*8 ym(N6DBS),yn(N6DBS)
      real*8 h,x,h2,swap

c----
c...  Executable code 

      nvar = 6*(nbod-1)

      h=htot/nstep

      do i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dy(i)
      enddo

      x=xs+h

      if(iwhich.eq.1) then
         call helioq_bs_der1(rtrans,nbod,mass,yn,yout)
      else
         call helioq_bs_der2(rtrans,nbod,mass,yn,yout)
      endif

      h2=2.0d0*h

      do n=2,nstep
         do i=1,nvar
            swap=ym(i)+h2*yout(i)
            ym(i)=yn(i)
            yn(i)=swap
         enddo
         x=x+h

         if(iwhich.eq.1) then
            call helioq_bs_der1(rtrans,nbod,mass,yn,yout)
         else
            call helioq_bs_der2(rtrans,nbod,mass,yn,yout)
         endif

      enddo

      do i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
      enddo

      return
      end   ! helioq_midd.f
c-------------------------------------------------------------------

