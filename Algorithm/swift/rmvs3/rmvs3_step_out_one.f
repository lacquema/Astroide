c*************************************************************************
c                            RMVS3_STEP_OUT_ONE.F
c*************************************************************************
c     This subroutine takes a full dt step in helio coord for 1 test particle
c in the outer region of an encounter.  It also sets up and call
c the inner region inetgration if necessary
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 r2hill         ==>  Hill radii (rael array)
c                 j2rp2,j4rp4        ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c             xbeg,ybeg,zbeg     ==>  initial planet position in helio coord 
c                                         (real arrays)
c             xtmp,ytmp,ztmp     ==>  position of planet wrt time
c                                       (2d real arrays)
c             vxtmp,vytmp,vztmp   ==>  velocity of planet wrt time
c                                       (2d real arrays)
c             ienc               ==> ienc = 0 if tp not involved in enc 
c                                       in outer region: = planet# if it is. 
c                                         (integer scalar)
c                                              =  0 No (integer scalar)  
c                 dt            ==>  time step (real sclar)
c             Output:
c                 xht,yht,zht    ==>  tp position in helio coord 
c                                       (real scalars)
c                 vxht,vyht,vzht ==>  tp velocity in helio coord 
c                                       (real scalars) 
c                 axht,ayht,azht ==>  tp acceleration in helio coord 
c                                       (real scalars) 
c                 istat              ==>  status of the test particles
c                                      (1d integer array)
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(2) = -1 ==> Danby did not work
c                 iencio,ienc  ==> Multiplied by -1 if particle entered inner
c                                     region (integer scalar)
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer scalar)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real scalar)
c
c Remarks: Adopted from rmvs3_step_out.f
c Authors:  Hervé Beust
c Date:    Apr 18, 2023

      subroutine rmvs3_step_out_one(i1st,nbod,mass,r2hill,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,xtmp,ytmp,ztmp,vxtmp,
     &     vytmp,vztmp,xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht,
     &     istat,dt,iencio,ienc,isperi,peri)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4,r2hill(nbod)
      real*8 xtmp(NPLMAX,NTENC),ytmp(NPLMAX,NTENC)
      real*8 ztmp(NPLMAX,NTENC)
      real*8 vxtmp(NPLMAX,NTENC),vytmp(NPLMAX,NTENC)
      real*8 vztmp(NPLMAX,NTENC)
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 vxbeg(nbod),vybeg(nbod),vzbeg(nbod)

c...  Inputs and Outputs:
      integer istat(NSTAT),ienc,iencio
      real*8 xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht

c...  Outputs:
      real*8 peri
      integer isperi

c...  Internals
      integer i1sttp,i,j,istattmp(NSTAT)
      real*8 dto
      real*8 xbegi(NPLMAX),ybegi(NPLMAX),zbegi(NPLMAX)
      real*8 xendi(NPLMAX),yendi(NPLMAX),zendi(NPLMAX)
      real*8 vxbegi(NPLMAX),vybegi(NPLMAX),vzbegi(NPLMAX)
      real*8 vxendi(NPLMAX),vyendi(NPLMAX),vzendi(NPLMAX)
      real*8 peril
      integer isperil

c----
c...  Executable code 

      i1sttp = i1st
      dto = dt/float(NTENC)

      istattmp(1) = 0 ! integrate the particle
      istattmp(2:NSTAT) = istat(2:NSTAT)

c...  do integration of outer loop
      do i=1,NTENC
c...      remember the current position of the planets
         if (i.eq.1) then
             xbegi(1:nbod) = xbeg(1:nbod)
             ybegi(1:nbod) = ybeg(1:nbod)
             zbegi(1:nbod) = zbeg(1:nbod)
             vxbegi(1:nbod) = vxbeg(1:nbod)
             vybegi(1:nbod) = vybeg(1:nbod)
             vzbegi(1:nbod) = vzbeg(1:nbod)
         else
             xbegi(1:nbod) = xtmp(1:nbod,i-1)
             ybegi(1:nbod) = ytmp(1:nbod,i-1)
             zbegi(1:nbod) = ztmp(1:nbod,i-1)
             vxbegi(1:nbod) = vxtmp(1:nbod,i-1)
             vybegi(1:nbod) = vytmp(1:nbod,i-1)
             vzbegi(1:nbod) = vztmp(1:nbod,i-1)
         endif
         xendi(1:nbod) = xtmp(1:nbod,i)
         yendi(1:nbod) = ytmp(1:nbod,i)
         zendi(1:nbod) = ztmp(1:nbod,i)
         vxendi(1:nbod) = vxtmp(1:nbod,i)
         vyendi(1:nbod) = vytmp(1:nbod,i)
         vzendi(1:nbod) = vztmp(1:nbod,i)

         call rmvs3_step_out2_one(i1sttp,nbod,mass(1:nbod),
     &        r2hill,j2rp2,j4rp4,
     &        xbegi(1:nbod),ybegi(1:nbod),zbegi(1:nbod),
     &        xendi(1:nbod),yendi(1:nbod),zendi(1:nbod),
     &        vxbegi(1:nbod),vybegi(1:nbod),vzbegi(1:nbod),
     &        vxendi(1:nbod),vyendi(1:nbod),vzendi(1:nbod),
     &        xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht,
     &        istattmp,dto,iencio,ienc,isperil,peril)	

         if (isperil.eq.0) then
                peri = peril
                isperi = isperil
         endif
         if ((i.eq.1).or.(isperi.ne.0)) isperi = isperil

      end do

c...   Have to update istat just in case Danby had problems
      if ((ienc.ne.0).and.(istattmp(1).ne.0)) then
            istat(1) = 1       ! it had problems
            istat(2:NSTAT) = istattmp(2:NSTAT)
      end if
c...   Now update the encounter stuff
      if (ienc.lt.0) istat(4:NSTAT) = istattmp(4:NSTAT)

      return
      end     ! rmvs3_step_out_one
c----------------------------------------------------------------
