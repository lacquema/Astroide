c*************************************************************************
c                            RMVS3_CHK_ONE.F
c*************************************************************************
c This subroutine checks to see if there are encounters for 1 tp
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 r2hill        ==> Hill radii (real array) 
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test particles (integer)
c                                      istat = 0 ==> active:  = 1 not
c                 dt            ==>  time step  (real scalar)
c                 rts           ==>  The fraction of hill sphere for encounter
c                                               (real scalar)
c             Output:
c                 icflg         ==> ecounters? = 1 Yes
c                                              =  0 No (integer scalar)  
c                 nenc          ==> nenc is the number of tp enc planet i
c                                   (integer array)
c                 itpenc        ==> itpenc(i) is a list of tp enc planet i
c                                   (integer array)
c                 ienc          ==> ienc = 0 if tp not involved in enc 
c                                           = planet# if it is. 
c                                           (integer scalar)
c
c Remarks: Based on RMVS3_CHK.F
c Authors:  Hervé Beust 
c Date:    Apr 20, 2023
c
c Last revision: 

      subroutine rmvs3_chk_one(nbod,mass,r2hill,xh,yh,zh,vxh,vyh,vzh,
     &      xht,yht,zht,vxht,vyht,vzht,istat,dt,rts,icflg,nenc,ienc)

      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs: 
      integer nbod,ntp,istat
      real*8 mass(nbod),xh(nbod),yh(nbod),zh(nbod),dt
      real*8 xht,yht,zht,vxht,vyht,vzht
      real*8 vxh(nbod),vyh(nbod),vzh(nbod),r2hill(nbod)
      real*8 rts

c...  Outputs
       integer icflg
       Integer nenc(nbod)
       integer ienc

c...  Internals
       integer i,j,iflag
       real*8 r2crit,r2critp
       real*8 xr,yr,zr,vxr,vyr,vzr

c-----
c...  Executable code 

c...    clear everything out
        icflg = 0
        nenc(1:nbod) = 0
        ienc = 0

        if (istat.eq.0) then

	    i = 2
	    iflag = 0				! precaution
 
c... Check for close encounters until we find one or run out of planets. BG

            do while ((iflag.eq.0).and.(i.le.nbod))

               r2crit = r2hill(i)*rts
               r2critp = 0.0           ! dummy so we can use rmvs routine

               xr = xht - xh(i)
               yr = yht - yh(i)
               zr = zht - zh(i)
               vxr = vxht - vxh(i)
               vyr = vyht - vyh(i)
               vzr = vzht - vzh(i)
               call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &                         r2crit,r2critp,iflag)
               if (iflag.gt.0) then
                  icflg = 1
                  ienc = i
                  nenc(i) = nenc(i) + 1
               endif
		 
 	       i = i + 1    		! next planet
            enddo
        endif

        return
        end  ! rmvs3_chk_one
c------------------------------------------------------

