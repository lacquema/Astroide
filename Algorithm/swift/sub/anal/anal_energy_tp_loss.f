c*************************************************************************
c                          ANAL_ENERGY_TP.F
c*************************************************************************
c Calculates the energy of the total system (massive bodies) wrt time.
c returns the total energy of n objects by direct pairwise summation
c G = 1., and we allow diff. masses.  Also returns square of total ang. mom.

      subroutine anal_energy_tp_loss(nbod,mass,istat,
     &           j2rp2,j4rp4,xh,yh,
     &            zh,
     &           vxh,vyh,vzh,ntp,masst,xht,yht,zht,vxht,vyht,vzht,
     &           energy,eltotal)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,ntp
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 masst
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer istat(NTPMAX,NSTAT)

c...  Output
      real*8 energy,eltot(3),eltotal

c...  Internals
      real*8 elx,ely,elz,ke,pot
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


      ke = 0.d0
      pot= 0.d0
      energy=0.d0
 
      eltot(1)=0.d0
      eltot(2)=0.d0
      eltot(3)=0.d0
 
      do i = 1,ntp
      if((istat(i,1).eq.1)) then
        write(*,*) istat(i,4)
      endif
      enddo 
        write (*,*) '------------' 
      do i = 1,ntp
      if((istat(i,1).eq.1).and.(istat(i,4).eq.0)) then
         write(*,*) istat(i,4) 
         elx=(ybt(i)*vzbt(i)-zbt(i)*vybt(i))*masst
         ely=(zbt(i)*vxbt(i)-xbt(i)*vzbt(i))*masst
         elz=(xbt(i)*vybt(i)-ybt(i)*vxbt(i))*masst
         eltot(1) = eltot(1) + elx
         eltot(2) = eltot(2) + ely
         eltot(3) = eltot(3) + elz
         
         ke = ke + 0.5*masst*(vxbt(i)**2 + vybt(i)**2 + vzbt(i)**2)
         do j = 1,nbod
	    xx = xbt(i) - xb(j)
	    yy = ybt(i) - yb(j)
	    zz = zbt(i) - zb(j)
	    rr2 = xx**2 + yy**2 + zz**2 
            if((rr2.ne.0.0d0).and.(mass(j).ne.0.0d0)) then
               pot = pot - masst*mass(j)/(sqrt(rr2))
            endif
         enddo
        istat(i,4)=1
      endif
      enddo
       write (*,*) '------------' 
      do i = 1,ntp
      if((istat(i,1).eq.1)) then
        write(*,*) istat(i,4)
      endif
      enddo 
        write (*,*) '------------' 
       energy = ke + pot
       eltotal =sqrt(eltot(1)**2+eltot(2)**2+eltot(3)**2)



      return	
      end      ! anal_energy
c-----------------------------------------------------------------------
