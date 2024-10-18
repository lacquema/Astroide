C**************************************************************************
C	    		        HELIOQ_BS_DER1
C**************************************************************************
c This is the subroutine that calculates the derivatives of the independant 
c variable for the Sun part of the Hamiltonian
c
c             Input:
c              nbod  ==> number of planets  (int scalar)
c              mass  ==>  mass of bodies (real array)
c 	        y    ==> values dependent variables  (real array)
c             rtrans ==>  solar encounter radii (real array)
c
c             Output:
c 	         dy  ==> derivatives of the independant var (real array)
c
c Remarks:  
c Authors:  Hal Levison
c Date:    2/2/2K
c Last revision: 

      subroutine helioq_bs_der1(rtrans,nbod,mass,y,dy)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod)
      real*8 rtrans(2),y(6,NTPMAX)

c...  Outputs:
      real*8 dy(6,NTPMAX)

c...  Internals:
      integer i
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
      real*8 ptx,pty,ptz,psun2
      real*8 ffunc,dffunc(NTPMAX)

c----
c...  Executable code 

c...  move things so that I can deal with it
      do i=1,nbod-1
         xh(i+1) = y(1,i)
         yh(i+1) = y(2,i)
         zh(i+1) = y(3,i)
         vxb(i+1) = y(4,i)
         vyb(i+1) = y(5,i)
         vzb(i+1) = y(6,i)
      enddo

c...  Calculate the momentum of the Sun
      ptx = mass(2)*vxb(2)
      pty = mass(2)*vyb(2)
      ptz = mass(2)*vzb(2)
      do i=3,nbod
         ptx = ptx + mass(i)*vxb(i)
         pty = pty + mass(i)*vyb(i)
         ptz = ptz + mass(i)*vzb(i)
      enddo

      ptx = ptx/mass(1)
      pty = pty/mass(1)
      ptz = ptz/mass(1)

      psun2 = ptx**2+pty**2+ptz**2

c....  get the F
      call helioq_getf(nbod,xh,yh,zh,rtrans,ffunc,dffunc)

      ffunc = 1.0d0 - ffunc 

c.... Calculate the derivatives
      do i=2,nbod
         dy(1,i-1) = ptx*ffunc
         dy(2,i-1) = pty*ffunc
         dy(3,i-1) = ptz*ffunc

         dy(4,i-1) = psun2*dffunc(i)*xh(i)
         dy(5,i-1) = psun2*dffunc(i)*yh(i)
         dy(6,i-1) = psun2*dffunc(i)*zh(i)
      enddo

      return
      end       ! helioq_bs_der1
c---------------------------------------------------------------------
