c*************************************************************************
c                            RMVS3_STEP_PAR.F
c*************************************************************************
c This subroutine takes a *dtout* step in helio coord.  
c both massive and test particles. Prallelized treatment of test particles
c INCLUDES close approuches between test particles and planets
c
c VERSION 3 of RMVS
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 tout,tstop    ==> time+dtout, stopping time
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xh,yh,zh      ==>  initial planet position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial planet velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final planet position in helio coord 
c                                       (real arrays)
c                 vxh,vyh,vzh   ==>  final planet velocity in helio coord 
c                                       (real arrays)
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c
c Remarks: Adapted from rmvs3_step.f
c Authors:  Herv� Beust 
c Date:    Apr. 18, 2023

      subroutine rmvs3_step_par(i1st,time,tout,tstop,nbod,
     &     ntp,mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,
     &     xht,yht,zht,vxht,vyht,vzht,istat,rstat,iencio,dt)	

      include '../swift.inc'
      include '../rmvs/rmvs.inc'
      
c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4,tout,tstop

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT),iencio(ntp)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals
      integer, parameter :: NSTEPMAX=150000
      integer i1sttp,icflg,i,j,imax,ii
      real*8 tpla,titp
      real*8 xhtp,yhtp,zhtp,vxhtp,vyhtp,vzhtp,axhtp,ayhtp,azhtp
      real*8 xtmp(NPLMAX,NTENC),ytmp(NPLMAX,NTENC)
      real*8 ztmp(NPLMAX,NTENC),rts,peri(NTPMAX)
      real*8 vxtmp(NPLMAX,NTENC),vytmp(NPLMAX,NTENC)
      real*8 vztmp(NPLMAX,NTENC)
c      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
c      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
c      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
c      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)
      Integer nenc(NPLMAX),ienc(NTPMAX)
      integer istattmp(NTPMAX,NSTAT),isperi(NTPMAX)
      real*8, dimension(:,:), target, allocatable :: xhs,yhs,zhs
      real*8, dimension(:,:), target, allocatable :: vxhs,vyhs,vzhs
      real*8, dimension(:), pointer :: xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg
      real*8, dimension(:), pointer :: xend,yend,zend,vxend,vyend,vzend
c     real*8 xhs(NPLMAX,NSTEPMAX),yhs(NPLMAX,NSTEPMAX)
c      real*8 zhs(NPLMAX,NSTEPMAX),vxhs(NPLMAX,NSTEPMAX)
c      real*8 vyhs(NPLMAX,NSTEPMAX),vzhs(NPLMAX,NSTEPMAX)
      real*8 r2hill(NPLMAX)
      
c-----
c...  Executable code; first compute hill radii 

      call util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill)
      i1sttp = i1st
      imax = nint((tout-time)/dt)
      allocate(xhs(nbod,imax+1))
      allocate(yhs(nbod,imax+1))
      allocate(zhs(nbod,imax+1)) 
      allocate(vxhs(nbod,imax+1))
      allocate(vyhs(nbod,imax+1))
      allocate(vzhs(nbod,imax+1))       

c...  First do the planets up to next ``tout''
      
      tpla = time
c      i = 1
      xhs(:,1) = xh
      yhs(:,1) = yh
      zhs(:,1) = zh
      vxhs(:,1) = vxh
      vyhs(:,1) = vyh
      vzhs(:,1) = vzh

c...  In this version we recalculate accel arrays at each beginning.
      i1st = 0
c      do while ( (tpla.lt.tout) .and. (tpla.lt.tstop) )
      do ii = 1,imax
         call step_kdk_pl(i1st,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,dt)	
c     if (i.gt.NSTEPMAX) then
c            write(*,*) 'NSTEPMAX = ',NSTEPMAX, ' exceeded !'
c            call util_exit(1)
c         end if
c... Store planet positions at intermediate step
         xhs(:,ii+1) = xh
         yhs(:,ii+1) = yh
         zhs(:,ii+1) = zh         
         vxhs(:,ii+1) = vxh
         vyhs(:,ii+1) = vyh
         vzhs(:,ii+1) = vzh
         tpla = tpla+dt
      end do

      rts = RHSCALE*RHSCALE
c...  Now the parallel loop over the test particles up to next ``tout''

!$omp parallel do
!$omp& default(shared)
!$omp& private(i,titp,i1sttp,icflg,nenc)
!$omp& private(xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg)
!$omp& private(xend,yend,zend,vxend,vyend,vzend)
!$omp& private(xhtp,yhtp,zhtp,vxhtp,vyhtp,vzhtp)
!$omp& private(axhtp,ayhtp,azhtp)
!$omp& private(xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)

      do j = 1,ntp        
         if (istat(j,1).eq.0) then
            i1sttp = 0
            i = 1
            titp = time
            xhtp = xht(j)
            yhtp = yht(j)
            zhtp = zht(j)
            vxhtp = vxht(j)
            vyhtp = vyht(j)
            vzhtp = vzht(j)
c...   now the loop up to the next 'tout'
            do while ( (titp.lt.tout).and.(titp.lt.tstop).and.
     &                     (istat(j,1).eq.0).and.(i.le.imax) )
            
c...  remember the begin and end position & velocities of the planets
               xbeg => xhs(:,i)
               ybeg => yhs(:,i)
               zbeg => zhs(:,i)
               vxbeg => vxhs(:,i)
               vybeg => vyhs(:,i)
               vzbeg => vzhs(:,i)
               xend => xhs(:,i+1)
               yend => yhs(:,i+1)
               zend => zhs(:,i+1)
               vxend => vxhs(:,i+1)
               vyend => vyhs(:,i+1)
               vzend => vzhs(:,i+1)

c...  are there any encounters?
               call rmvs3_chk_one(nbod,mass,r2hill,xbeg,ybeg,zbeg,
     &              vxbeg,vybeg,vzbeg,xhtp,yhtp,zhtp,vxhtp,vyhtp,
     &              vzhtp,istat(j,1),dt,rts,icflg,nenc,ienc(j))
c...     nenc and itpenc not used here!
c... if not just do a normal step
               if(icflg.eq.0) then

                 call step_kdk_onetp(i1sttp,nbod,mass,j2rp2,j4rp4,
     &                 xbeg,ybeg,zbeg,xend,yend,zend,
     &                 xhtp,yhtp,zhtp,vxhtp,vyhtp,vzhtp,
     &                 axhtp,ayhtp,azhtp,istat(j,:),dt)	
                  if (istat(j,1).eq.0) istat(j,2) = 0
               else

c...  ENCOUNTER STUFF FROM HERE ON!!!!!

c...  Do the interpolation for intermediate steps
                  call rmvs3_interp(nbod,xbeg,ybeg,zbeg,
     &                 vxbeg,vybeg,vzbeg,xend,yend,zend,
     &                 vxend,vyend,vzend,dt,mass(1),NTENC,
     &                 xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)
c...  do the integration
                  call rmvs3_step_out_one(i1sttp,nbod,mass,r2hill,
     &                 j2rp2,j4rp4,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &                 xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp,xhtp,yhtp,zhtp,
     &                 vxhtp,vyhtp,vzhtp,axhtp,ayhtp,azhtp,istat(j,:),
     &                 dt,iencio(j),ienc(j),isperi(j),peri(j))
               end if  

               xht(j) = xhtp
               yht(j) = yhtp
               zht(j) = zhtp
               vxht(j) = vxhtp
               vyht(j) = vyhtp
               vzht(j) = vzhtp

               i = i+1
               titp = titp+dt
            end do   
         end if
      end do

!$omp end parallel do
      
      time = tpla    ! � voir

      deallocate(xhs)
      deallocate(yhs)
      deallocate(zhs) 
      deallocate(vxhs)
      deallocate(vyhs)
      deallocate(vzhs)   
      return
      end   ! rmvs3_step_par
c------------------------------------------------------------------------