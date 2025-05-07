c*************************************************************************
c                            STEP_S6B_HJS_PAR.F
c*************************************************************************
c This subroutine takes a *dtout* step in Jacobi coord.  
c both massive and test particles. Prallelized treatment of test particles
c    Uses S6B 6th order symplectic algorithm (Chambers & Murison 2000)
c     
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 time          ==>  current time (real scalar)
c                 tout,tstop    ==> time+dtout, stopping time      
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 oloc          ==>  Link bodies - orbits (2D int array)
c                 oloct         ==>  Link orbits - tp's (2D int array)
c                 mat,umat      ==>  Conversion matrixes for bodies
c                                     (2D real arrays)
c                 matp,umatp    ==>  Conversion vectors for tp's
c                                     (2D real arrays)
c                 eta,mu        ==>  Masses for center & satellites for bods.
c                                     (real arrays)
c                 etatp         ==>  Masses for centers for tp's (real array)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj ==> initial pos. in Jacobi coord (real arrays)
c                 vxj,vyj,vzj ==> initial vel. in Jacobi coord (real arrays)
c                 xjt,yjt,zjt ==> initial tp pos in Jacobi coord (real arrays)
c                 vxjt,vyjt,vzjt==>initial tp vel in Jacobi coord (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 rstat           ==>  status of the test paricles
c                                      (2d real array)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==> final position in Jacobi coord (real arrays)
c                 vxj,vyj,vzj   ==> final velocity in Jacobi coord (real arrays)
c                 xjt,yjt,zjt  ==> final position in Jacobi coord (real arrays)
c                 vxjt,vyjt,vzjt ==> final vels. in Jacobi coord (real arrays)
c
c
c Remarks: Adopted from step_s6b_hjs.f & rmvs3_step_par.f
c Authors:  Hervé Beust
c Date:    A
c Revision: Feb 25, 2025

      subroutine step_s6b_hjs_par(i1st,time,tout,tstop,nbod,ntp,
     &     oloc,oloct,mat,umat,matp,umatp,mass,eta,mu,etatp,
     &     xj,yj,zj,vxj,vyj,vzj,xjt,yjt,zjt,vxjt,vyjt,vzjt,
     &     istat,rstat,dt)	

      include '../swift.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,tout,tstop
      real*8 mu(nbod),eta(nbod),etatp(ntp)
      real*8 umat(NPLMAX,NPLMAX),mat(NPLMAX,NPLMAX)
      real*8 umatp(NPLMAX,NTPMAX),matp(NPLMAX,NTPMAX)
      integer oloc(NPLMAX,NPLMAX),oloct(NPLMAX,NTPMAX)

      
c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals
      integer, parameter :: NSTEPMAX=1000000
      integer i1sttp,i,j,k,imax
      real*8 tpla,titp
      real*8 xjtp,yjtp,zjtp,vxjtp,vyjtp,vzjtp,axjtp,ayjtp,azjtp
      real*8, pointer :: xjpl(:,:),yjpl(:,:),zjpl(:,:),ir3j(:,:)
      real*8, pointer :: vxjpl(:,:),vyjpl(:,:),vzjpl(:,:)
      real*8, pointer :: xbpl(:,:),ybpl(:,:),zbpl(:,:)
      real*8, pointer :: axbpl(:,:),aybpl(:,:),azbpl(:,:)
      real*8, dimension(:,:), allocatable :: xjplc,yjplc,zjplc,ir3jc
      real*8, dimension(:,:), allocatable :: vxjplc,vyjplc,vzjplc
      real*8, dimension(:,:), allocatable :: xbplc,ybplc,zbplc
      real*8, dimension(:,:), allocatable :: axbplc,aybplc,azbplc
      real*8, dimension(:,:,:), target, allocatable :: xjs,yjs,zjs,ir3js
      real*8, dimension(:,:,:), target, allocatable :: vxjs,vyjs,vzjs
      real*8, dimension(:,:,:), target, allocatable :: xbs,ybs,zbs
      real*8, dimension(:,:,:), target, allocatable :: axbs,aybs,azbs
      
c----
c...  Executable code 

      i1sttp = i1st
      imax = nint((tout-time)/dt)
      allocate(xjplc(nbod,4))
      allocate(yjplc(nbod,4))
      allocate(zjplc(nbod,4))
      allocate(vxjplc(nbod,4))
      allocate(vyjplc(nbod,4))
      allocate(vzjplc(nbod,4))
      allocate(xbplc(nbod,4))
      allocate(ybplc(nbod,4))
      allocate(zbplc(nbod,4))
      allocate(axbplc(nbod,4))
      allocate(aybplc(nbod,4))
      allocate(azbplc(nbod,4))
      allocate(ir3jc(nbod,4))
      allocate(xjs(nbod,4,imax+1))
      allocate(yjs(nbod,4,imax+1))
      allocate(zjs(nbod,4,imax+1))      
      allocate(vxjs(nbod,4,imax+1))
      allocate(vyjs(nbod,4,imax+1))
      allocate(vzjs(nbod,4,imax+1))  
      allocate(xbs(nbod,4,imax+1))
      allocate(ybs(nbod,4,imax+1))
      allocate(zbs(nbod,4,imax+1))  
      allocate(axbs(nbod,4,imax+1))
      allocate(aybs(nbod,4,imax+1))
      allocate(azbs(nbod,4,imax+1))  
      allocate(ir3js(nbod,4,imax+1))
      
c...  First do the planets up to next ``tout''
      
      tpla = time
c      i = 1

c...  In this version we recalculate accel arrays at each beginning.
      i1st = 0
c      do while ( (tpla.lt.tout) .and. (tpla.lt.tstop) )
      do i = 1,imax
         call step_s6b_pl_hjs(i1st,nbod,oloc,umat,mat,mass,eta,mu,
     &        xjplc,yjplc,zjplc,vxjplc,vyjplc,vzjplc,ir3jc,
     &        xbplc,ybplc,zbplc,axbplc,aybplc,azbplc,
     &        xj,yj,zj,vxj,vyj,vzj,dt)	
         do k = 1,4
           xjs(:,k,i) = xjplc(:,k)
           yjs(:,k,i) = yjplc(:,k)
           zjs(:,k,i) = zjplc(:,k)         
           vxjs(:,k,i) = vxjplc(:,k)
           vyjs(:,k,i) = vyjplc(:,k)
           vzjs(:,k,i) = vzjplc(:,k)   
           xbs(:,k,i) = xbplc(:,k)
           ybs(:,k,i) = ybplc(:,k)
           zbs(:,k,i) = zbplc(:,k)         
           axbs(:,k,i) = axbplc(:,k)
           aybs(:,k,i) = aybplc(:,k)
           azbs(:,k,i) = azbplc(:,k)   
           ir3js(:,k,i) = ir3jc(:,k)
         end do
         tpla = tpla+dt
      end do
      
c...  Now the parallel loop over the test particles up to next ``tout''

!$omp parallel do
!$omp& default(shared)
!$omp& private(i,k,titp,i1sttp)
!$omp& private(xjpl,yjpl,zjpl,ir3j,vxjpl,vyjpl,vzjpl)
!$omp& private(xbpl,ybpl,zbpl,axbpl,aybpl,azbpl)
!$omp& private(xjtp,yjtp,zjtp,vxjtp,vyjtp,vzjtp,axjtp,ayjtp,azjtp)
    
      do j = 1,ntp
         if (istat(j,1).eq.0) then
            i1sttp = 0
            i = 1
            titp = time
            xjtp = xjt(j)
            yjtp = yjt(j)
            zjtp = zjt(j)
            vxjtp = vxjt(j)
            vyjtp = vyjt(j)
            vzjtp = vzjt(j)
c...   now the loop up to the next 'tout'
            do while ( (titp.lt.tout).and.(titp.lt.tstop).and.
     &                            (istat(j,1).eq.0) )
            
c...  remember intermediate positions of the planets
               xjpl => xjs(:,:,i)
               yjpl => yjs(:,:,i)
               zjpl => zjs(:,:,i)
               ir3j => ir3js(:,:,i)
               vxjpl => vxjs(:,:,i)
               vyjpl => vyjs(:,:,i)
               vzjpl => vzjs(:,:,i)
               xbpl => xbs(:,:,i)
               ybpl => ybs(:,:,i)
               zbpl => zbs(:,:,i)
               axbpl => axbs(:,:,i)
               aybpl => aybs(:,:,i)
               azbpl => azbs(:,:,i)

               call step_s6b_onetp_hjs(i1sttp,nbod,matp,umatp,oloct,
     &              mass,eta,mu,etatp(j),xjpl,yjpl,zjpl,ir3j,
     &              vxjpl,vyjpl,vzjpl,xbpl,ybpl,zbpl,axbpl,aybpl,azbpl,
     &              xjtp,yjtp,zjtp,vxjtp,vyjtp,vzjtp,
     &              axjtp,ayjtp,azjtp,istat(j,:),dt)	
               if (istat(j,1).eq.0) istat(j,2) = 0

               xjt(j) = xjtp
               yjt(j) = yjtp
               zjt(j) = zjtp
               vxjt(j) = vxjtp
               vyjt(j) = vyjtp
               vzjt(j) = vzjtp
               i = i+1
               titp = titp+dt
            end do   
         end if
      end do

!$omp end parallel do
      
      time = tout    ! à voir
c...
      deallocate(xjplc)
      deallocate(yjplc)
      deallocate(zjplc)
      deallocate(vxjplc)
      deallocate(vyjplc)
      deallocate(vzjplc)
      deallocate(xbplc)
      deallocate(ybplc)
      deallocate(zbplc)
      deallocate(axbplc)
      deallocate(aybplc)
      deallocate(azbplc)
      deallocate(ir3jc)
      deallocate(xjs)
      deallocate(yjs)
      deallocate(zjs)      
      deallocate(vxjs)
      deallocate(vyjs)
      deallocate(vzjs)  
      deallocate(xbs)
      deallocate(ybs)
      deallocate(zbs)  
      deallocate(axbs)
      deallocate(aybs)
      deallocate(azbs)  
      deallocate(ir3js)
      return
      end   ! step_s6b_par
c------------------------------------------------------------------------

