c*************************************************************************
c                            RMVS_STEP_IN_MIG_TP.F
c*************************************************************************
c This subroutine takes a step for tp in inner region.
c
      subroutine rmvs_step_in_mig_tp(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istat,dt,ipl,
     &              aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze)

      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st,ipl
      real*8 mass(nbod),dt,j2rp2,j4rp4  
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 aoblxb,aoblyb,aoblzb,aoblxe,aoblye,aoblze

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)

c...  Internals:
      real*8 dth,j2rp2i,j4rp4i
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

c...  disable J2 and J4 terms in the getacch_tp calls
      j2rp2i = 0.0d0
      j4rp4i = 0.0d0

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,ntp,mass,j2rp2i,j4rp4i,xbeg,ybeg,zbeg,
     &                  xht,yht,zht,istat,axht,ayht,azht)
          call rmvs_obl_acc(nbod,ntp,ipl,mass,j2rp2,j4rp4,xbeg,
     &         ybeg,zbeg,aoblxb,aoblyb,aoblzb,xht,yht,zht,istat,
     &         axht,ayht,azht)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_mig_tp(ntp,mass(1),xht,yht,zht,vxht,vyht,vzht,dt,istat)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,ntp,mass,j2rp2i,j4rp4i,xend,yend,zend,
     &                  xht,yht,zht,istat,axht,ayht,azht)
      call rmvs_obl_acc(nbod,ntp,ipl,mass,j2rp2,j4rp4,xend,yend,
     &     zend,aoblxe,aoblye,aoblze,xht,yht,zht,istat,
     &     axht,ayht,azht)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

      return
      end   ! rmvs_step_in_tp
c---------------------------------------------------------------------
