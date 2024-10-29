C**************************************************************************
C	    		        HELIOQ_BS_DER2
C**************************************************************************
c This is the subroutine that calculates the derivatives of the independant 
c variable for the modified Kep part of the Hamiltonian
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

      subroutine helioq_bs_der2(rtrans,nbod,mass,y,dy)

      include '../swift.inc'
      include 'helioq.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod)
      real*8 rtrans(2),y(6,NTPMAX)

c...  Outputs:
      real*8 dy(6,NTPMAX)

c...  Internals:
      integer i,j
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)
      real*8 ptx,pty,ptz,psun2,rr2,fac
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

c...  calculate the momentum of the Sun
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

c...  get the F
      call helioq_getf(nbod,xh,yh,zh,rtrans,ffunc,dffunc)

c.... Calculate the derivatives
      do i=2,nbod
         dy(1,i-1) = vxb(i) + ptx*ffunc
         dy(2,i-1) = vyb(i) + pty*ffunc
         dy(3,i-1) = vzb(i) + ptz*ffunc

         rr2 = xh(i)**2 + yh(i)**2 + zh(i)**2
         fac = -1.0*mass(1)/(rr2*sqrt(rr2))

         dy(4,i-1) = fac*xh(i) - psun2*dffunc(i)*xh(i)
         dy(5,i-1) = fac*yh(i) - psun2*dffunc(i)*yh(i)
         dy(6,i-1) = fac*zh(i) - psun2*dffunc(i)*zh(i)

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c         do j=1,6
c            if(abs(dy(j,i-1)).gt.1.0d6) then
c               write(*,*) 'Here Der2_1a',i,j
c               write(*,*) 'Here Der2_1b',ffunc,ptx,pty,ptz
c               write(*,*) 'Here Der2_1c',psun2,dffunc(i)
c            endif
c         enddo
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo

      return
      end            ! helioq_bs_der2
c---------------------------------------------------------------------
