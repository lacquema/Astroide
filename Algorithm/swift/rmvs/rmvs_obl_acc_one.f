c*************************************************************************
c                            RMVS_OBL_ACC_ONE.F
c*************************************************************************
c This subroutine adds the differential accel due to J2 and J4 (for 1 tp)
c
c             Input:
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ipl            ==>  the planet that is currently in the 
c                                      center (integer scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xpl,ypl,zpl    ==>  massive part position at 
c                                       (real arrays)
c              aoblx,aobly,aoblz ==> acceleration of the Sun on the central pl
c                                     due to J2 anf J4
c                                      (real scalars)
c                 xpt,ypt,zpt    ==>  initial part position in planet-coord 
c                                      (real scalars)
c                 istat           ==>  status of the test particle
c                                      (integer array)
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(2) = -1 ==> Danby did not work
c              axpt,aypt,azpt    ==>  accel on tp w/o J2 and J4
c                                      (real arrays)
c
c             Output:
c                 axpt,aypt,azpt ==>  accel on tp WITH J2 and J4 added.
c                                       (real arrays)
c
c Remarks: Taken from rmvs_obl_acc.f
c Authors:  Hervé Beust
c Date:    Apr 20, 2023
c Last revision: 

      subroutine rmvs_obl_acc_one(nbod,ipl,mass,j2rp2,j4rp4,xpl,ypl,
     &     zpl,aoblx,aobly,aoblz,xpt,ypt,zpt,istat,axpt,aypt,azpt)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,ipl
      integer istat(NSTAT)
      real*8 mass(nbod),j2rp2,j4rp4  
      real*8 xpl(nbod),ypl(nbod),zpl(nbod)
      real*8 xpt,ypt,zpt
      real*8 aoblx,aobly,aoblz

c...  Inputs and Outputs:
      real*8 axpt,aypt,azpt

c...  Internals:
      integer i
      real*8 xht,yht,zht,irht,aoblxt,aoblyt,aoblzt


c----
c...  Executable code 

c...  Do we need to do this?
      if(j2rp2.eq.0.0d0) then
         return                  !!!!!!! NOTE 
      endif

c...  first get barycentric accel 
      xht = xpt-xpl(ipl)
      yht = ypt-ypl(ipl)
      zht = zpt-zpl(ipl)
      irht = 1.0d0/sqrt(xht**2 + yht**2 + zht**2)

      call obl_acc_tp(1,istat(1),mass(ipl),j2rp2,j4rp4,xht,yht,zht,
     &     irht,aoblxt,aoblyt,aoblzt)

      axpt = axpt + aoblxt - aoblx
      aypt = aypt + aoblyt - aobly
      azpt = azpt + aoblzt - aoblz

      return
      end     ! rmvs_obl_acc_one
c------------------------------------------------------------------------
