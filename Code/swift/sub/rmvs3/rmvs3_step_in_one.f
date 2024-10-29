c*************************************************************************
c                            RMVS3_STEP_IN_ONE.F
c*************************************************************************
c This subroutine takes a full dt step in helio coord for 1 test particle
c in the inner region of an encounter.
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xpl,ypl,zpl    ==>  position of planet wrt time
c                                       (2d real arrays)
c                 xbeg,ybeg,zbeg ==>  position of planet @ beginning of dt
c                                       (real arrays)
c                 vxbeg,vybeg,vzbeg ==>  vel of planet @ beginning of dt
c                                       (2d real arrays)
c                 xend,yend,zend ==>  position of planet @ end of dt
c                                       (2d real arrays)
c                 vxend,vyend,vzend ==>  vel of planet @ end of dt
c                                       (2d real arrays)
c                 xht,yht,zht    ==>  initial tp position in helio coord 
c                                      (real scalars)
c                 vxht,vyht,vzht ==>  initial tp velocity in helio coord 
c                                        (real scalars)
c                 istat          ==>  status of the test paricles
c                                      (1d integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                      istat(i) = -1 ==> Danby did not work
c                 nenci          ==> nenci(i) is the number of tp enc planet i
c                                    in inner region (integer array)
c                 itpenci       ==> itpenci(*,i) is a list of tp enc planet i
c                                   in inner region (2d integer array)
c                 dt            ==>  time step (real sclar)
c             Output:
c                 xht,yht,zht    ==>  final tp position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final tp position in helio coord 
c                                       (real arrays)
c                                      NOTE: only the tp in the inner region
c                                            will have their x and v's changed 
c                 isperi         ==> = 0 if tp went through peri
c                                    =-1 if tp pre peri
c                                    = 1 if tp post peri
c                                         (integer array)
c                 peri           ==> set to pericenter dist. if isperi=0
c                                         (real array)
c
c
c Remarks: Adapted from rmvs3_step_in
c Authors:  Hervé Beust
c Date:    Apr 18, 2023

      subroutine rmvs3_step_in_one(i1st,nbod,mass,j2rp2,j4rp4,
     &           xpl,ypl,zpl,xbeg,ybeg,zbeg,vxbeg,vybeg,vzbeg,
     &           xend,yend,zend,vxend,vyend,vzend,xht,yht,zht,
     &           vxht,vyht,vzht,istat,nenci,isperi,peri,dt)


      include '../swift.inc'
      include '../rmvs/rmvs.inc'

c...  Inputs Only: 
      integer nbod,i1st
      real*8 mass(nbod),dt,j2rp2,j4rp4
      real*8 xpl(NPLMAX,NTPHENC),ypl(NPLMAX,NTPHENC)
      real*8 zpl(NPLMAX,NTPHENC)
      Integer nenci(nbod)
      real*8 xbeg(nbod),ybeg(nbod),zbeg(nbod)
      real*8 vxbeg(nbod),vybeg(nbod),vzbeg(nbod)
      real*8 xend(nbod),yend(nbod),zend(nbod)
      real*8 vxend(nbod),vyend(nbod),vzend(nbod)

c...  Inputs and Outputs:
      integer istat(NSTAT)
      real*8 xht,yht,zht,vxht,vyht,vzht,axht,ayht,azht

c...  Outputs:
      real*8 peri
      integer isperi

c...  Internals
      integer i1sttp,i,j,k,link,istattmp(NSTAT)
      real*8 xtpt,ytpt,ztpt,dti
      real*8 vxtpt,vytpt,vztpt,axtpt,aytpt,aztpt
      real*8 xpltb(NPLMAX),ypltb(NPLMAX),zpltb(NPLMAX)
      real*8 xplte(NPLMAX),yplte(NPLMAX),zplte(NPLMAX)
      real*8 masst(NPLMAX),perit
      integer isperit
      logical*2 lperit

      real*8 irh(NPLMAX),aoblx(NPLMAX,0:NTPHENC)
      real*8 aobly(NPLMAX,0:NTPHENC),aoblz(NPLMAX,0:NTPHENC)

c----
c...  Executable code 

      dti = dt/float(NTPHENC)

c...  First get the accel due to J2 and J4 in barycentric frame on planets
      if(j2rp2.ne.0.0d0) then
         do i=2,nbod
            irh(i) = 1.0d0/sqrt( xbeg(i)**2 + ybeg(i)**2 + zbeg(i)**2 )
         enddo
         call obl_acc(nbod,mass,j2rp2,j4rp4,xbeg,ybeg,zbeg,
     &        irh,aoblx(1,0),aobly(1,0),aoblz(1,0))
         do j=1,NTPHENC
            do i=2,nbod
               irh(i) = 1.0d0/sqrt( xpl(i,j)**2 + ypl(i,j)**2 + 
     &              zpl(i,j)**2 )
            enddo
            call obl_acc(nbod,mass,j2rp2,j4rp4,xpl(1,j),ypl(1,j),
     &        zpl(1,j),irh,aoblx(1,j),aobly(1,j),aoblz(1,j))
         enddo
      endif

c...  Must do each planet one at a time
      do i=2,nbod
         if(nenci(i).ne.0) then

c...          set up test particle
            xtpt = xht - xbeg(i)
            ytpt = yht - ybeg(i)
            ztpt = zht - zbeg(i)
            vxtpt = vxht - vxbeg(i)
            vytpt = vyht - vybeg(i)
            vztpt = vzht - vzbeg(i)
            istattmp(1) = 0
            lperit = .false.

C...          set up planets at t=0
            call rmvs_step_in_mvpl(i,nbod,mass,xbeg,ybeg,zbeg,
     &                masst,xpltb,ypltb,zpltb,xplte,yplte,zplte)
            i1sttp = 0

            call util_peri(0,1,xtpt,ytpt,ztpt,vxtpt,
     &             vytpt,vztpt,masst(1),isperit,perit,lperit)

            do j=1,NTPHENC
               call rmvs_step_in_mvpl(i,nbod,mass,
     &              xpl(1:nbod,j),ypl(1:nbod,j),zpl(1:nbod,j),
     &              masst,xpltb,ypltb,zpltb,xplte,yplte,zplte)
               call rmvs_step_in_onetp(i1sttp,nbod,masst,j2rp2,j4rp4,
     &              xpltb(1:nbod),ypltb(1:nbod),zpltb(1:nbod),
     &              xplte(1:nbod),yplte(1:nbod),zplte(1:nbod),
     &           xtpt,ytpt,ztpt,vxtpt,vytpt,vztpt,axtpt,aytpt,aztpt,
     &           istattmp,dti,i,aoblx(i,j-1),aobly(i,j-1),aoblz(i,j-1),	
     &           aoblx(i,j),aobly(i,j),aoblz(i,j))

               call util_peri(1,nenci(i),xtpt,ytpt,ztpt,vxtpt,
     &                vytpt,vztpt,masst(1),isperit,perit,lperit)

            end do

c...          put things back
            xht = xtpt + xend(i)
            yht = ytpt + yend(i)
            zht = ztpt + zend(i)
            vxht = vxtpt + vxend(i)
            vyht = vytpt + vyend(i)
            vzht = vztpt + vzend(i)
            if (lperit) then
                isperi = 0
            else
                isperi = isperit
            endif
            peri = perit
            if (istattmp(1).ne.0) then
               istat(1) = 1
               istat(2:NSTAT) = istattmp(2:NSTAT)
            endif
          endif
        enddo

        return
        end  ! rmvs3_step_in_one
c------------------------------------------------------------------




