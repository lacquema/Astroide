c*************************************************************************
c                        DRIFT_MIG_TP.F
c*************************************************************************
c This subroutine loops thorugh the TEST particles and calls the danby routine

      subroutine drift_mig_tp(ntp,msun,xjt,yjt,zjt,vxjt,vyjt,vzjt,dt
     &                    ,istat)	

      include '../../swift.inc'

c...  Inputs Only: 
      integer ntp
      real*8 msun,dt

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp),yjt(ntp),zjt(ntp)
      real*8 vxjt(ntp),vyjt(ntp),vzjt(ntp)

c...  Internals:
	integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth
	do j = 1,ntp
           if(istat(j,1).eq.0) then
	      call drift_one(msun,xjt(j),yjt(j),zjt(j),
     &             vxjt(j),vyjt(j),vzjt(j),dt,iflg)
              if(iflg.ne.0) then
                 istat(j,1) = 1
                 istat(j,2) = -1
              endif
           endif          
	enddo
	return
	end
c--------------------------------------------------------------------------

