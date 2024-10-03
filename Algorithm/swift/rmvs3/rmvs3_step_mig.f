c*************************************************************************
c                            RMVS3_STEP_MIG.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c both massive and test particles
c INCLUDES close approuches between test particles and planets
c

      subroutine rmvs3_step_mig(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,
     &     xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,
     &     istat,rstat,dt,masst,energy,eltot,HandE,ttime
     &     , LastEnergy,LastH,Hvartot,Envartot)
     	
      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,ntp,i1st
      real*8 mass(nbod),dt,time,j2rp2,j4rp4

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 rstat(NTPMAX,NSTATR)
      real*8 LastEnergy(NB_DT),LastH(NB_DT),ttime(NB_DT)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      real*8 Envartot,Hvartot,LostEnergy,LostH,Etheo,Htheo
      real*8 HandE(NTPMAX,4),rand


c...  Internals
      integer i1sttp,i1stpl,i1sto,icflg,i,j,np
      real*8 xtmp(NPLMAX,NTENC),ytmp(NPLMAX,NTENC)
      real*8 ztmp(NPLMAX,NTENC),rts,peri(NTPMAX)
      real*8 vxtmp(NPLMAX,NTENC),vytmp(NPLMAX,NTENC)
      real*8 vztmp(NPLMAX,NTENC)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)
      Integer ienc(NTPMAX),indice
      integer istattmp(NTPMAX,NSTAT),isperi(NTPMAX)
      real*8 masst,ke,pot,energy,energyp,energylast,Hlast
      real*8 eltot(3),eltotp(3),energy_loss,eltot_loss
      real *8 mu,ialpha,a,e,q,inc,capom,omega,capm
      real*8 siga,sigb,chi2



c----
c...  Executable code 

      i1sttp = i1st
      i1sto = i1st

c      call anal_energy_delta_tp(nbod,mass,istat,j2rp2,j4rp4,xh,yh,zh,
c     &           vxh,vyh,vzh,ntp,masst,xht,yht,zht,vxht,vyht,vzht,
c     &           HandE,energy,eltot)

c      write(*,*)time,energy,Etheo
c      write(*,*)time,sqrt(eltot(1)**2+eltot(2)**2+eltot(3)**2),Htheo

c      call anal_energy_tp_loss(nbod,mass,istat,
c     &           j2rp2,j4rp4,xh,yh,zh,
c     &           vxh,vyh,vzh,ntp,masst,xht,yht,zht,vxht,vyht,vzht,
c     &           energy_loss,eltot_loss)

c        energy = -2.5d-10
c        eltot(3) = 0.0

c        if(time.gt.dt) then
c          Envartot=energy/dt
c          Hvartot= eltot(3)/dt
c        else 
c          Envartot=0.0
c          Hvartot=0.0  
c        endif
c        mu=mass(2)+mass(1)


      indice=int(mod(time/dt,dble(NB_DT)))+1
      ttime(indice)=time 
      LastEnergy(indice) = energy       
      LastH(indice) = eltot(3)

c      if(time.lt.((NB_DT+5)*dt)) then
c        Envartot=0.0
c        Hvartot=0.0
c      else  
c         Hvartot=0.0
c        Envartot=0.0
c         do i=1,NB_DT
c             Hvartot = Hvartot + LastH(i)
c             Envartot = Envartot + LastEnergy(i)
c          enddo 
c          Hvartot=Hvartot/(NB_DT*dt)
c          Envartot=Envartot/(NB_DT*dt)
c      endif

  


      if((indice.eq.1).and.(time.ne.0.0)) then 
               Hvartot=0.0
               Envartot=0.0
               do i=1,NB_DT
                 Hvartot = Hvartot + LastH(i)
                 Envartot = Envartot + LastEnergy(i)
               enddo 
               Hvartot=Hvartot/(NB_DT*dt)
               Envartot=Envartot/(NB_DT*dt)
      endif

      if(time.lt.(2*NB_DT*dt)) then
        Envartot=0.0
        Hvartot=0.0
      endif  

c         call orbel_xv2el(xh(2),yh(2),zh(2),vxh(2),vyh(2),vzh(2),
c     &     mu,ialpha,a,e,inc,capom,omega,capm)

c      if(time.lt.0.25d6) then
c        Envartot=-2.5d-10 
c        Hvartot=0.0 
c      endif 



c             if(mod(time,1.0d4).eq.0.) then
c         call orbel_xv2el(xh(2),yh(2),zh(2),vxh(2),vyh(2),vzh(2),
c     &     mu,ialpha,a,e,inc,capom,omega,capm)
c        write(*,*)time, Envartot,a
c        write(*,*)time, Hvartot,e
c          endif    


c...  are there any encounters?
      rts = RHSCALE*RHSCALE
      call rmvs3_chk_par(nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,xht,yht,
     &       zht,vxht,vyht,vzht,istat,dt,rts,icflg,ienc)

c.... if not just do a normal step and leave
      if(icflg.eq.0) then
         call step_kdk_mig(i1st,time,nbod,ntp,mass,j2rp2,j4rp4,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,dt,
     &        Envartot,Hvartot)	

         do i=1,ntp
            if(istat(i,1).eq.0) then
               istat(i,2) = 0
            endif
         enddo

         return       !  NOTE AN EXIT
      endif

c...  ENCOUNTER STUFF FROM HERE ON!!!!!

c...  save initial x and v of planets if there are planocentric enc
      do i=1,nbod
            xbeg(i) = xh(i)
            ybeg(i) = yh(i)
            zbeg(i) = zh(i)
            vxbeg(i) = vxh(i)
            vybeg(i) = vyh(i)
            vzbeg(i) = vzh(i)
      enddo

c... do a full step for the planets
      i1stpl = i1st
      call step_kdk_mig_pl(i1stpl,nbod,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,dt,-Envartot,-Hvartot)


c...  save the final position and velocity of planets
      do i=1,nbod
         xend(i) = xh(i)
         yend(i) = yh(i)
         zend(i) = zh(i)
         vxend(i) = vxh(i)
         vyend(i) = vyh(i)
         vzend(i) = vzh(i)
      enddo

c...  Now do the interpolation for intermediate steps
      call rmvs3_interp(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &     xend,yend,zend,vxend,vyend,vzend,dt,mass(1),NTENC,
     &     xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp)

c...  do the integration
      call rmvs3_step_mig_out(i1st,nbod,ntp,mass,j2rp2,j4rp4,
     &     xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,xtmp,ytmp,ztmp,vxtmp,
     &     vytmp,vztmp,xht,yht,zht,vxht,vyht,vzht,istat,ienc,dt,
     &     isperi,peri)


c...  As of this point all the test particles that are involved in an
c...  encounter have been moved.  But not the ones that have not.
c...  so move those,  BUT NOT the onces in the encounter

 
c...  make a temporary istat array so only they are active

      do i=1,ntp
         if(istat(i,1).eq.0) then
            if(ienc(i).ne.0) then
               istattmp(i,1) = 1 ! don't integrate it
            else
               istattmp(i,1) = 0 ! integrate it
            endif
            do j=2,NSTAT
               istattmp(i,j) = 0
            enddo
         else
            istattmp(i,1) = 1   ! don't integrate it
         endif
       enddo

      

c...  do a full step
      i1sto = 0      ! we need to recalculate accel arrays
      call step_kdk_mig_par_tp(i1sto,nbod,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,
     &              xht,yht,zht,vxht,vyht,vzht,istattmp,dt)

c      write (*,*) '***********' 
c      do i = 1,ntp
c      if((istat(i,1).eq.1)) then
c        write(*,*) istat(i,4)
c      endif
c      enddo 
c        write (*,*) '**********' 



c...  fix up the istat array
c      do i=1,ntp
c         if(istattmp(i,2) .ne. 0) then   ! danby screwed up
c            istat(i,1) = 1
c            do j=2,NSTAT
c               istat(i,j) = istattmp(i,j)
c            enddo
c          endif
c       enddo

c...  put the enc info into istat

      do i=1,ntp
         if(istattmp(i,2) .ne. 0) then   ! danby screwed up
            istat(i,1) = 1
            do j=2,NSTAT
               istat(i,j) = istattmp(i,j)
            enddo
          endif
         if(istat(i,1).eq.0) then
            if(ienc(i).lt.0) then
               istat(i,2) = ienc(i)
               istat(i,3) = abs(ienc(i))
               if(isperi(i).eq.0) then
                  rstat(i,3) = peri(i)
                  np = NSTATP + abs(ienc(i)) - 1
                  if(rstat(i,np).eq.0) then
                     rstat(i,np) = peri(i)
                  else
                     rstat(i,np) = min(rstat(i,np),peri(i))
                  endif
               else
                  rstat(i,3) = 0.0d0
               endif
            else if(ienc(i).gt.0) then
               istat(i,2) = ienc(i)
            else
               istat(i,2) = 0
            endif
         endif
      enddo

 

c...  we MUST make sure that the saved accel arrays are ok
c...  calculate them again
       i1st = 0 

      return
      end   ! step_enc
c------------------------------------------------------------------------






