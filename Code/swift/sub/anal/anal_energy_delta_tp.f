c*************************************************************************
c                          ANAL_ENERGY_TP.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.

      subroutine anal_energy_delta_tp(nbod,mass,istat,j2rp2,j4rp4,xh,yh,
     &           zh,
     &           vxh,vyh,vzh,ntp,masst,xht,yht,zht,vxht,vyht,vzht,HandE,
     &            energy,eltot)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 masst,ialpha,a,e,inc,capom,omega,capm
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 HandE(NTPMAX,4)
      integer istat(NTPMAX,NSTAT)

c...  Output
      real*8 energy,eltot(3),ke,pot

c...  Internals
      real*8 elx,ely,elz
      real*8 xx,yy,zz,rr2,oblpot,msys,irh(NTPMAX),ir3h(NTPMAX)
      real*8 xb(NTPMAX),yb(NTPMAX),zb(NTPMAX)       ! Used NTPMAX for symba
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
      real*8 xbt(NTPMAX),ybt(NTPMAX),zbt(NTPMAX)       ! Used NTPMAX for symba
      real*8 vxbt(NTPMAX),vybt(NTPMAX),vzbt(NTPMAX)

      integer i,j

c----
c...  Executable code 

      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &           xb,yb,zb,vxb,vyb,vzb,msys)

      call coord_h2b_tp(ntp,xht,yht,zht,vxht,vyht,vzht,
     &                  xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &                  xbt,ybt,zbt,vxbt,vybt,vzbt)
      eltot(1)=0.d0
      eltot(2)=0.d0
      eltot(3)=0.d0
      energy=0.0
      ke = 0.d0
      pot= 0.d0
c      Etheo=0.0d0
c      Htheo=0.0d0

      do i = 1,ntp
      if(istat(i,1).eq.0) then 
c          call orbel_xv2el(xht(2),yht(2),zht(2),vxht(2),vyht(2),vzht(2),
c     &    mass(1),ialpha,a,e,inc,capom,omega,capm)

c         Etheo=Etheo - mass(1)*masst/(2*a)
c         Htheo=Htheo + sqrt(mass(1)*a*(1-e**2))*masst

         elx=(ybt(i)*vzbt(i)-zbt(i)*vybt(i))*masst
         ely=(zbt(i)*vxbt(i)-xbt(i)*vzbt(i))*masst
         elz=(xbt(i)*vybt(i)-ybt(i)*vxbt(i))*masst
         eltot(1) = eltot(1) + elx - HandE(i,1)
         eltot(2) = eltot(2) + ely - HandE(i,2)
         eltot(3) = eltot(3) + elz - HandE(i,3)
         
         HandE(i,1)=elx
         HandE(i,2)=ely
         HandE(i,3)=elz


         ke =  0.5*masst*(vxbt(i)**2 + vybt(i)**2 + vzbt(i)**2)
         pot = 0.d0

         do j = 1,nbod
	    xx = xbt(i) - xb(j)
	    yy = ybt(i) - yb(j)
	    zz = zbt(i) - zb(j)
	    rr2 = xx**2 + yy**2 + zz**2 
            if((masst.ne.0.0d0).and.(mass(j).ne.0.0d0)) then
               pot = pot - masst*mass(j)/(sqrt(rr2))
            endif
         enddo
          
         energy = energy + ke + pot - HandE(i,4)
         HandE(i,4) = ke + pot  

      endif
      enddo

      return	
      end      ! anal_energy
c-----------------------------------------------------------------------

