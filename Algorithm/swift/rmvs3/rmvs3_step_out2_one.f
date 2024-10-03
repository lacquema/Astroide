c*************************************************************************
c                            RMVS3_STEP_OUT2_ONE.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for 1 test particle
c in the outer region of an encounter.  It also sets up and call
c the inner region inetgration if necessary
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c             j2rp2,j4rp4        ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c             xbeg,ybeg,zbeg     ==>  initial planet position in helio coord 
c                                         (real arrays)
c             vxbeg,vybeg,vzbeg  ==>  initial planet velcoity in helio coord 
c                                         (real arrays)
c             xend,yend,zend     ==>  final planet position in helio coord 
c                                         (real arrays)
c             vxend,vyend,vzend  ==>  final planet velcoity in helio coord 
c                                         (real arrays)
c             istat              ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c             ienc0               ==> ienc0(j) = 0 if tp j not involved in enc 
c                                       in outer region: = planet# if it is. 
c                                         (integer array)
c                                              =  0 No (integer scalar)  
c                 dt            ==>  time step (real sclar)
c             Input and Output:
c                 xht,yht,zht    ==>  tp position in helio coord 
c                                       (real scalars)
c                 vxht,vyht,vzht ==>  tp velocity in helio coord 
c                                       (real scalars)
c                 axht,ayht,azht ==>  tp accel's in helio coord 
c                                       (real arrays)      
c             istat              ==>  status of the test paricles
c                                      (1d integer array)
c                                      istat(1) = 0 ==> active:  = 1 not
c                                      istat(2) = -1 ==> Danby did not work
c             iencio,ienc0       ==> Multiplied by -1 if particle entered inner
c                                     region (integer scalar)
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer scalar)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real scalar)
c
c Remarks: Adapted from rmvs3_step_out2.f
c Authors:  Hervé Beust
c Date:    Apr 18, 2023
c Last revision: 

      subroutine rmvs3_step_out2_one(i1st,nbod,mass,r2hill,
     &      j2rp2,j4rp4,xbeg,ybeg,zbeg,xend,yend,zend,vxbeg,vybeg,
     &      vzbeg,vxend,vyend,vzend,xht,yht,zht,
     &      vxht,vyht,vzht,axht,ayht,azht,istat,dt,
     &      iencio,ienc0,isperi,peri)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4,r2hill(nbod)
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 vxbeg(nbod),vybeg(nbod),vzbeg(nbod)
      real*8 xend(nbod),yend(nbod),zend(nbod)
      real*8 vxend(nbod),vyend(nbod),vzend(nbod)

c...  Inputs and Outputs:
      integer istat(NSTAT),ienc0,iencio
      real*8 xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht

c...  Outputs:
      real*8 peri
      integer isperi

c...  Internals
      integer i1sttp,i,j,icflg,i1sto
      real*8 rts
      Integer nenc(NPLMAX),itpenc(NPLMAX),ienc
      real*8 xtmp(NPLMAX,NTPHENC),ytmp(NPLMAX,NTPHENC)
      real*8 ztmp(NPLMAX,NTPHENC)
      real*8 vxtmp(NPLMAX,NTPHENC),vytmp(NPLMAX,NTPHENC)
      real*8 vztmp(NPLMAX,NTPHENC)


c----
c...  Executable code 

      i1sttp = i1st

c...  are there any encounters?
      rts = RHPSCALE*RHPSCALE
      call rmvs3_chk_one(nbod,mass,r2hill,xbeg,ybeg,zbeg,
     &     vxbeg,vybeg,vzbeg,xht,yht,zht,vxht,vyht,vzht,
     &     istat(1),dt,rts,icflg,nenc,ienc)

c...  keep track of encounters in the inner region
      call rmvs3_elog_one(icflg,iencio,ienc,istat)

c.... if not just do a normal step and leave
      if(icflg.eq.0) then
          call step_kdk_onetp(i1st,nbod,mass,j2rp2,j4rp4,
     &           xbeg,ybeg,zbeg,xend,yend,zend,
     &           xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht,istat,dt)	

         return       !  NOTE AN EXIT
      endif

c...  INNER ENCOUNTER STUFF FROM HERE ON!!!!!
c...  Now do the interpolation for intermediate steps
      call rmvs3_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,dt,mass(1),NTPHENC,
     &     xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)
c...  Do the inner integration
      call rmvs3_step_in_one(i1sttp,nbod,mass,j2rp2,j4rp4,xtmp,
     &     ytmp,ztmp,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,xht,yht,zht,
     &     vxht,vyht,vzht,istat,nenc(1:nbod),isperi,peri,dt)
      
      if (istat(1).eq.0) istat(2) = 0

c...  put the enc info into istat
      if(istat(1).eq.0) then
          if(ienc.ne.0) then
             ienc0 = -1 * abs(ienc0)
          endif
      endif

c...  we MUST make sure that the saved accel arrays are ok
c...  calculate them again
      i1st = 0 

      return
      end   ! rmvs3_step_out2_one
c------------------------------------------------------------------------

