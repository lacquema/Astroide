c*************************************************************************
c                        SYMBATR_HELIO_GETACCH.F
c*************************************************************************
c This subroutine calculates the acceleration on the massive particles
c in the HELIOCENTRIC frame. 
c             Input:
c                 nbod        ==>  number of massive bodies (int scalor)
c                 nbodm       ==>  The last massive particle
c                                  (int scalor)
c                 mass        ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh    ==>  position in heliocentric coord (real arrays)
c             Output:
c                 axh,ayh,azh ==>  acceleration in helio coord (real arrays)
c
c Remarks Based on helio_getacch.f
c Author:  Hal Levison  
c Date:    3/20/97
c Last revision: 

      subroutine symbatr_helio_getacch(nbod,nbodm,mass,j2rp2,
     &     j4rp4,xh,yh,zh,axh,ayh,azh)

      include '../swift.inc'

c...  Inputs: 
      integer nbod,nbodm
      real*8 mass(nbod),j2rp2,j4rp4
      real*8 xh(nbod),yh(nbod),zh(nbod)

c...  Outputs:
      real*8 axh(nbod),ayh(nbod),azh(nbod)
                
c...  Internals:
      integer i,j
      real*8 aoblx(NTPMAX),aobly(NTPMAX),aoblz(NTPMAX) 
      real*8 ir3h(NTPMAX),irh(NTPMAX)
      real*8 dx,dy,dz,rji2,faci,facj,irij3

c----
c...  Executable code 

      do i=1,nbod
         axh(i) = 0.0
         ayh(i) = 0.0
         azh(i) = 0.0
      enddo

c...  now the third terms
      do i=2,nbodm
         do j=i+1,nbod

             dx = xh(j) - xh(i)
             dy = yh(j) - yh(i)
             dz = zh(j) - zh(i)
             rji2 = dx*dx + dy*dy + dz*dz

             irij3 = 1.0d0/(rji2*sqrt(rji2))
             faci = mass(i)*irij3
             facj = mass(j)*irij3

             axh(j) = axh(j) - faci*dx
             ayh(j) = ayh(j) - faci*dy
             azh(j) = azh(j) - faci*dz

             axh(i) = axh(i) + facj*dx
             ayh(i) = ayh(i) + facj*dy
             azh(i) = azh(i) + facj*dz

         enddo
      enddo

c     New: Do the interactions among the bodies 
c     with mass < mtiny


      call symbatr_interface(xh,yh,zh,axh,ayh,azh,
     &     mass,nbod,nbodm)


c...  Now do j2 and j4 stuff
      if(j2rp2.ne.0.0d0) then
         call getacch_ir3(nbod,2,xh,yh,zh,ir3h,irh)
         call obl_acc(nbod,mass,j2rp2,j4rp4,xh,yh,zh,irh,
     &        aoblx,aobly,aoblz)
         do i = 2,nbod
            axh(i) = axh(i) + aoblx(i) - aoblx(1)
            ayh(i) = ayh(i) + aobly(i) - aobly(1)
            azh(i) = azh(i) + aoblz(i) - aoblz(1)
         enddo
      endif

      return
      end      ! symbatr_helio_getacch

c---------------------------------------------------------------------




