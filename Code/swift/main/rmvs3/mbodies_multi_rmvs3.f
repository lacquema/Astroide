C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C
        PROGRAM MULTI_EXTRACT

        include '../../sub/swift.inc'
       
        REAL*8, PARAMETER :: BINR=1.0d0
        REAL*8, PARAMETER :: RMAXBIN=3000.0d0
        INTEGER*4, PARAMETER :: NMAXBIN = NINT(RMAXBIN/BINR)
        INTEGER*4, PARAMETER :: NSMAX = 100
        INTEGER*4 NBIN(NMAXBIN)
        REAL*4, DIMENSION(:,:), ALLOCATABLE :: DATA4
        REAL*4, DIMENSION(-NPLMAX:10*NTPMAX) :: TA,TINC,TE,TR,TEX,
     &                   TEY,TEZ,TX,TY,TZ,TCAPOM,TOMEGA,TCAPM,
     &                   TEPX,TEPY,TEPZ

        INTEGER*4 NT,IPL,K,IBIN,ISNAPTOT,NSNAP,IREAD,SIZ,ISNAP
        REAL*8 XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),MASS(NPLMAX)
        REAL*8 VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX)
        REAL*8 VXH(NPLMAX),VYH(NPLMAX),VZH(NPLMAX)
        REAL*8 XH(NPLMAX),YH(NPLMAX),ZH(NPLMAX)
        REAL*8 VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX)
        REAL*8 XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX)
        REAL*8 EX(NPLMAX),EY(NPLMAX),EZ(NPLMAX),AINIT(NTPMAX*NSMAX)
        REAL*8 EPX(NPLMAX),EPY(NPLMAX),EPZ(NPLMAX)
        REAL*8 MAT(NPLMAX,NPLMAX),UMAT(NPLMAX,NPLMAX)
        real*8 xht(10*NTPMAX),yht(10*NTPMAX),zht(10*NTPMAX)
        real*8 vxht(10*NTPMAX),vyht(10*NTPMAX),vzht(10*NTPMAX)
        real*8 matp(NPLMAX,10*NTPMAX),umatp(NPLMAX,10*NTPMAX)
        REAL*8 ETA(NPLMAX),MU(NPLMAX),WIMPMFP(10*NTPMAX)
        REAL*8 EPXT,EPYT,EPZT,EXT,EYT,EZT,MPAS(3,3)
        INTEGER ISTAT(10*NTPMAX,NSTAT),OLOC(NPLMAX,NPLMAX)
        REAL*8 RSTAT(10*NTPMAX,NSTATR),etatp(10*NTPMAX)
        REAL*8 DTSNAP,RR,TIM,XT,YT,ZT,VXT,VYT,VZT,TDEB,q
        INTEGER NBOD,NTP,IERR,IFOL,ISTEP,KBIN,OLOCT(NPLMAX,10*NTPMAX)
        INTEGER IFLGCHK,IU,NLEFT,I,ID,J,ISIZ,ISNAPST,SNAPST,REFNS
        INTEGER IO_READ_HDR,IO_READ_LINE,NPO,NBOD0,NTP0(NSMAX)
        INTEGER IO_READ_HDR_R,IO_READ_LINE_R,II,ialpha,NS,NBSIMU
        REAL*8 CAPMB0(NPLMAX),CAPMB(NPLMAX)

        REAL*8 T0,TSTOP,DT,DTOUT,DTDUMP,DR,CX,CY,CZ,MATR(3,3),GM
        REAL*8 T,TMAX,EMAX(10*NTPMAX),VARPI,KR,KT,RCOLL,RX,RY,R,U

        REAL*8 RMIN,RMAX,RMAXU,QMIN,RPLSQ(NPLMAX),VX,VY,VZ,MSYS
        LOGICAL ERROR,WIMPS,EVAP,COLL,REFOVER
        LOGICAL*2 LCLOSE,OK
        REAL*8 A,E,INC,CAPOM,OMEGA,CAPM,CXB,CYB,CZB,CVXB,CVYB,CVZB
         REAL*8 ELH,ELK,ELP,ELQ,Z0,ALPHEVAP,REVAP,RADMIN,NN,TS
        REAL*8 SUMX,SUMY,SUMZ,SUMVX,SUMVY,SUMVZ,ETAL,j2rp2,j4rp4

        CHARACTER*(CHLEN) OUTFILE(NSMAX),INPARFILE,INPLFILE,
     &           INTPFILE(NSMAX),FOPENSTAT,DIRO(NSMAX),DIRS,
     &           GNAME(NSMAX),DEV,EXTRFILE

        NPO=200
        print*,' '
        KBIN = NMAXBIN/10
        COLL=.FALSE.
        WRITE(*,*)'Enter number of individual simulations : '
        READ(*,*) NBSIMU
        DO I=1,NBSIMU
          WRITE(*,*)'Enter name of parameter data file : ',I
          READ(*,'(a)')inparfile
          call io_init_param_hb(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro(i),dirs,
     &         gname(i),outfile(i),fopenstat)
        END DO

        WRITE(*,*) 'Enter beginning time (estimated) :'
        READ(*,*)TDEB
        IF (TDEB.GT.TSTOP) THEN
          PRINT*,'Time too large'
          STOP
        END IF

        WRITE(*,*) 'Enter final time (estimated) :'
        READ(*,*)TIM
        IF (TIM.GT.TSTOP) THEN
          PRINT*,'Time too large'
          STOP
        END IF
        WRITE(*,*) 'Enter timestep between snapshots :'
        READ(*,*)DTSNAP
        NSNAP = NINT(DTSNAP/DTOUT)
        SNAPST = FLOOR(TDEB/DTSNAP)

        WRITE(*,*) ' '
        WRITE(*,*) 'Enter name of planet data file : '
        READ(*,'(A)')INPLFILE

        CALL IO_INIT_PL(INPLFILE,LCLOSE,IFLGCHK,NBOD0,
     &       MASS,XH,YH,ZH,VXH,VYH,VZH,RPLSQ,j2rp2,j4rp4)

        WRITE(*,*) ' '
        WRITE(*,*) 'Enter name of extract file: '
        READ(*,'(A)') EXTRFILE
        IU = 20
        DR = 180.0/PI
c        DEV=TRIM(DIRO(1))//'/'//TRIM(EXTRFILE)//'.gdf'
        ! DEV=TRIM(DIRS)//'/'//TRIM(EXTRFILE)//'.dat'
        DEV='./'//TRIM(EXTRFILE)//'.dat'

        IF (BTEST(IFLGCHK,0)) THEN
           WRITE(*,*) ' Reading INTEGER*2 binary files '
        ELSE IF(BTEST(IFLGCHK,1)) THEN
           WRITE(*,*) ' Reading REAL*4 binary files '
        ELSE
           WRITE(*,*) ' Error: no binary file format specified '
           WRITE(*,*) '        in param file '
           STOP
        ENDIF
        OLOC=0
        DO J=2,nbod0
            OLOC(J,1)=-1
            OLOC(J,J)=1
        END DO

        DO NS = 1,NBSIMU
          write(*,*) 'Enter name of test particle data file : '
          read(*,'(a)') intpfile(ns)
          call io_init_tp(intpfile(ns),ntp0(ns),xht,yht,zht,vxht,vyht,
     &     vzht,istat,rstat)
          do i=1,ntp0(ns)
            j=sum(ntp0(1:ns-1))
            call orbel_xv2aeq(xht(i),yht(i),zht(i),
     &        vxht(i),vyht(i),vzht(i),mass(1),ialpha,ainit(j+i),e,q)
          end do
        END DO

        SIZ = MAX(NBOD0,3)+NBOD0
        REFNS = 1
        print*,'ici'

        TMAX = 0.d0
        DO NS = 1,NBSIMU
          OPEN(UNIT=IU, FILE=trim(diro(ns))//'/'//OUTFILE(NS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
          mpas =0.0d0
          IREAD = 0
          ISNAP = 0
          ISNAPST = 0
          NTP = 0
          T = -1.0d0
          IPL = 2
c          OLOC=0
          TIMELOOP1 : DO WHILE(T.LT.TIM)
            IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
              IERR = IO_READ_HDR(IU,T,NBOD,NLEFT)
            ELSE
              IERR = IO_READ_HDR_R(IU,T,NBOD,NLEFT)
            ENDIF
            IF (IERR.NE.0) THEN
c              print*,'T=',T,'Simu',NS,' finie, SIZ=',SIZ
              EXIT TIMELOOP1
            END IF
            IF (T.GT.TMAX) THEN
              TMAX = T
              IF (NS.NE.REFNS) THEN
                REFNS = NS
                print*,'====> refns=',refns
              END IF
            END IF
            IF (MOD(IREAD,200).EQ.0) 
     &              print*,'T=',T,'Simu',NS,' NLEFT=',nleft,siz
            OK=(MOD(IREAD,NSNAP).EQ.0)
            IF (OK.AND.(ISNAPST.LT.SNAPST)) THEN
              ISNAPST = ISNAPST+1
              OK = .FALSE.
            END IF
            BODLOOP1 : DO I=2,NBOD
              IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              IF(IERR.NE.0) EXIT BODLOOP1
            END DO BODLOOP1
            NTP=0
            LEFTLOOP1 : DO I=1,NLEFT
              IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              END IF
              IF(IERR.NE.0) THEN
                PRINT*,'Pb avec la taille'
                EXIT LEFTLOOP1
              END IF
              IF (OK) THEN
                IF ((E.LT.1.0d0).and.(e.ge.0.0d0)) THEN
                  NTP=NTP+1
                END IF
              END IF
            END DO LEFTLOOP1
            IF (OK) THEN
              IF (NS.EQ.REFNS) THEN
                SIZ = SIZ+NBOD+NTP+1
              ELSE
                SIZ = SIZ+NTP
              END IF
            END IF
            IREAD=IREAD+1
          END DO TIMELOOP1
          CLOSE(IU)
        END DO
        IF (TIM.GT.TMAX) TIM = TMAX
        ISNAPTOT = FLOOR(TIM/DTSNAP)+1-SNAPST

        ALLOCATE(DATA4(SIZ,11))
        DATA4 = SNGL(0.0d0)
        ISIZ = 0

        DO I=1,NBOD0
          DATA4(I,1) = SNGL(MASS(I))
        END DO
        DO I=2,NBOD0
          DATA4(I,2) = SNGL(MASS(1))
        END DO
        DO I=2,NBOD0
          DATA4(I,3) = SNGL(MASS(I))
        END DO
        DO I=1,3
          DO J=1,3
            DATA4(J,I+3) = SNGL(MPAS(J,I))
          END DO
        END DO
        DATA4(1,7) = REAL(ISNAPTOT)
        print*,nbod0,sum(ntp0(1:nbsimu))
        DATA4(2,7) = REAL(NBOD0)
        DATA4(3,7) = REAL(SUM(NTP0(1:NBSIMU)))
        DO I=1,NBOD0
          DO J=1,NBOD0
            DATA4(MAX(NBOD0,3)+J,I) = REAL(OLOC(J,I))
          END DO
        END DO
        ISIZ = MAX(NBOD0,3)+NBOD0

        DO NS = 1,NBSIMU
          OPEN(UNIT=IU+NS, FILE=TRIM(DIRO(NS))//'/'//OUTFILE(NS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
        END DO
        IREAD = 0
        ISNAP = 0
        ISNAPST = 0
        NTP = 0
        T = -1.0d0

c        print*,oloc(1:nbod0,1:nbod0)
        DO WHILE(T.LT.TIM)
c
          OK = (MOD(IREAD,NSNAP).EQ.0)
          IF (OK.AND.(ISNAPST.LT.SNAPST)) THEN
            ISNAPST = ISNAPST+1
            OK = .FALSE.
          END IF
c
          NTP = 0
          SIMLOOP : DO NS = 1,NBSIMU
            IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
              IERR = IO_READ_HDR(IU+NS,TS,NBOD,NLEFT)
            ELSE
              IERR = IO_READ_HDR_R(IU+NS,TS,NBOD,NLEFT)
            ENDIF
            IF (NS.EQ.REFNS) T = TS
            IF (IERR.NE.0) THEN
               IF (MOD(IREAD,200).EQ.0) 
     &              print*,'T=',SNGL(T),'Simu',NS,' finie'
               CYCLE SIMLOOP
            END IF
            IF (MOD(IREAD,200).EQ.0) THEN
                   print*,'T=',SNGL(T),'Simu',NS,' NLEFT=',nleft
            END IF
            XH(1) = 0.0d0
            YH(1) = 0.0d0
            ZH(1) = 0.0d0
            VXH(1) = 0.0d0
            VYH(1) = 0.0d0
            VZH(1) = 0.0d0
            BODLOOP2 : DO I=2,NBOD
              IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
                IERR = IO_READ_LINE(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              IF(IERR.NE.0) EXIT BODLOOP2
              IF ((OK).AND.(NS.EQ.REFNS)) THEN
                TA(ID)=SNGL(A)
                TE(ID)=SNGL(E)
                TINC(ID)=SNGL(INC*DR)
                TCAPOM(ID)=SNGL(CAPOM*DR)
                TOMEGA(ID)=SNGL(OMEGA*DR)
                TCAPM(ID)=SNGL(CAPM*DR)
                CAPMB0(-ID)=CAPM
                GM=MASS(-ID)+MASS(1)
                CALL KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &    R,XH(-ID),YH(-ID),ZH(-ID),VXH(-ID),VYH(-ID),VZH(-ID),
     &       EX(-ID),EY(-ID),EZ(-ID),EPX(-ID),EPY(-ID),EPZ(-ID))
              END IF
              CAPMB(-ID) = CAPM
c             if (ok.and.(ns.eq.38)) 
c     &             print*,ns,t,ts,sngl(capmb0(2)),sngl(capmb(2))
c              IF (CAPMB(-ID).NE.CAPMB0(-ID)) print*,-id,
c     &          (CAPMB(-ID)-CAPMB0(-ID))   ! *180.d0/PI
            END DO BODLOOP2
            IF ((OK).AND.(NS.EQ.REFNS)) THEN
              call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xb,yb,zb,vxb,vyb,vzb,msys)
              DO I=1,NBOD
                TR(-I) = SNGL(SQRT(XB(I)*XB(I)+YB(I)*YB(I)+ZB(I)*ZB(I)))
                TX(-I) = SNGL(XB(I))
                TY(-I) = SNGL(YB(I))
                TZ(-I) = SNGL(ZB(I))
                TEX(-I) = SNGL(EX(I))
                TEY(-I) = SNGL(EY(I))
                TEZ(-I) = SNGL(EZ(I))
                TEPX(-I) = SNGL(EPX(I))
                TEPY(-I) = SNGL(EPY(I))
                TEPZ(-I) = SNGL(EPZ(I))
              END DO
            END IF
c
            LEFTLOOP2 : DO I=1,NLEFT
              IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
                IERR = IO_READ_LINE(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU+NS,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              ID = ID + SUM(NTP0(1:NS-1))
              IF(IERR.NE.0) EXIT LEFTLOOP2
              IF (OK) THEN
                IF ((E.LT.1.0d0).and.(e.ge.0.0d0)) THEN
                  NTP=NTP+1
                  IF (T.NE.TS) THEN
                    NN = SQRT(MASS(1)/(A*A*A))
c                    print*,'t-ts',sngl(t),sngl(t-ts),sngl(a)
c     &                   sngl(capmb0(2)),sngl(capmb(2))
                    CAPM = CAPM+(T-TS)*NN
                  END IF
                  TA(NTP)=SNGL(A)
                  TE(NTP)=SNGL(E)
                  TINC(NTP)=SNGL(INC*DR)
                  TCAPOM(NTP)=SNGL(CAPOM*DR)
                  TOMEGA(NTP)=SNGL(OMEGA*DR)
                  TCAPM(NTP)=SNGL(CAPM*DR)
                  GM = MASS(1)
                  CALL KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &               R,CX,CY,CZ,VX,VY,VZ,EXT,EYT,EZT,EPXT,EPYT,EPZT)
                  call coord_h2b_tp(1,cx,cy,cz,vx,vy,vzh,
     &               xb(1),yb(1),zb(1),vxb(1),vyb(1),vzb(1),
     &               cxb,cyb,czb,cvxb,cvyb,cvzb)
                  TR(NTP) = SNGL(AINIT(ID))
                  TX(NTP) = SNGL(CXB)
                  TY(NTP) = SNGL(CYB)
                  TZ(NTP) = SNGL(CZB)
                END IF
              END IF
            END DO LEFTLOOP2
          END DO SIMLOOP

          IF (OK) THEN
            ISNAP = ISNAP+1
            DATA4(ISIZ+1,1) = SNGL(T)
            DATA4(ISIZ+1,2) = REAL(NBOD)
            DATA4(ISIZ+1,3) = REAL(NTP)
            DATA4(ISIZ+1,4) = REAL(ISNAP)
            DATA4(ISIZ+1,5) = REAL(ISNAPTOT)
            ISIZ = ISIZ+1
            DO I=-1,-NBOD,-1
              DATA4(ISIZ-I,1)=TA(I)
              DATA4(ISIZ-I,2)=TE(I)
              DATA4(ISIZ-I,3)=TEX(I)
              DATA4(ISIZ-I,4)=TEY(I)
              DATA4(ISIZ-I,5)=TEZ(I)
              DATA4(ISIZ-I,6)=TEPX(I)
              DATA4(ISIZ-I,7)=TEPY(I)
              DATA4(ISIZ-I,8)=TEPZ(I)
              DATA4(ISIZ-I,9)=TX(I)
              DATA4(ISIZ-I,10)=TY(I)
              DATA4(ISIZ-I,11)=TZ(I)
            END DO
            DO I=1,NTP
              DATA4(ISIZ+I+NBOD,1)=TA(I)
              DATA4(ISIZ+I+NBOD,2)=TE(I)
              DATA4(ISIZ+I+NBOD,3)=TINC(I)
              DATA4(ISIZ+I+NBOD,4)=TCAPOM(I)
              DATA4(ISIZ+I+NBOD,5)=TOMEGA(I)
              DATA4(ISIZ+I+NBOD,6)=TCAPM(I)
              DATA4(ISIZ+I+NBOD,7)=TR(I)
              DATA4(ISIZ+I+NBOD,9)=TX(I)
              DATA4(ISIZ+I+NBOD,10)=TY(I)
              DATA4(ISIZ+I+NBOD,11)=TZ(I)
            END DO
            ISIZ = ISIZ+NBOD+NTP
            print*,'Creating snapshot at T=',SNGL(T),' npts=',NBOD+NTP
          END IF
          IREAD=IREAD+1
        END DO
        DO NS = 1,NBSIMU
          CLOSE(IU+NS)
        END DO
        print*,siz,isiz,isnaptot
c        CALL IMWRITE(DEV,
c     &             DATA4,SIZ,BINR,RMAXBIN,DBLE(NBOD),11,0.,0.,0.,
c     &             1,0.,0.,0.,SNGL(-10000.),SNGL(0.001),
c     &             0.,0.,0.,ERROR)
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)SIZ,11,1
        WRITE(18,*)DATA4
        CLOSE(18)
        DEALLOCATE(DATA4)

        print*,'Entering followbodies tim=',tim
        CALL FOLLOWBODIES(TIM,DTOUT,LCLOSE,IFLGCHK,NBSIMU,REFNS,
     &   DIRS,TRIM(INPLFILE),DIRO(1:NBSIMU),OUTFILE(1:NBSIMU))
        END

C
C--------------------------------------------------------------------
C               Keplerian orbit stuff
C--------------------------------------------------------------------
C
        SUBROUTINE KEPLER_STUFF(GM,A,E,INC,CAPOM,OMEGA,CAPM,
     &                   R,CX,CY,CZ,VX,VY,VZ,EX,EY,EZ,EPX,EPY,EPZ)

        include '../../sub/swift.inc'

        REAL*8 GM,A,E,INC,CAPOM,OMEGA,CAPM
        REAL*8 R,RX,RY,CX,CY,CZ,VRX,VRY,VX,VY,VZ
        REAL*8 EX,EY,EZ,EPX,EPY,EPZ,U,F,F1,F2,F3,DELA
        LOGICAL OK
        INTEGER I

        EX = COS(OMEGA)*COS(CAPOM)-SIN(OMEGA)*COS(INC)*SIN(CAPOM)
        EY = COS(OMEGA)*SIN(CAPOM)+SIN(OMEGA)*COS(INC)*COS(CAPOM)
        EZ = SIN(OMEGA)*SIN(INC)
        EPX = -SIN(OMEGA)*COS(CAPOM)-COS(OMEGA)*COS(INC)*SIN(CAPOM)
        EPY = -SIN(OMEGA)*SIN(CAPOM)+COS(OMEGA)*COS(INC)*COS(CAPOM)
        EPZ = COS(OMEGA)*SIN(INC)
        U = CAPM
        OK = .FALSE.
        I=0
        DO WHILE (.NOT.OK)
          F = U-E*SIN(U)-CAPM
          F1 = 1.-E*COS(U)
          F2 = U-CAPM-F
          F3 = 1.-F1
          DELA = -F/F1
          DELA = -F/(F1+0.5*DELA*F2)
          DELA = -F/(F1+0.5*DELA*F2+DELA*DELA*F3/6.)
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.1e-12).OR.(I.GE.500))
        END DO
        R = A*(1.-E*COS(U))
        RX = A*(COS(U)-E)
        RY = A*SQRT(1.-E*E)*SIN(U)
        VRX = -SQRT(GM/A)*A/R*SIN(U)
        VRY = SQRT(GM/A)*A/R*COS(U)
        CX = EX*RX+EPX*RY
        CY = EY*RX+EPY*RY
        CZ = EZ*RX+EPZ*RY
        VX = EX*VRX+EPX*VRY
        VY = EY*VRX+EPY*VRY
        VZ = EZ*VRX+EPZ*VRY
        END

C--------------------------------------
        SUBROUTINE KEPLER(CAPM,E,U)

        IMPLICIT NONE

        REAL*8 U,DELA,CAPM,F1,F2,F3,F,E
        LOGICAL OK
        INTEGER*4 I

        U = CAPM
        OK = .FALSE.
        I=0
        DO WHILE (.NOT.OK)
          F = U-E*SIN(U)-CAPM
          F1 = 1.-E*COS(U)
          F2 = U-CAPM-F
          F3 = 1.-F1
          DELA = -F/F1
          DELA = -F/(F1+0.5*DELA*F2)
          DELA = -F/(F1+0.5*DELA*F2+DELA*DELA*F3/6.)
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.1e-12).OR.(I.GE.500))
        END DO
        END

C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ,A)

        include '../../sub/swift.inc'

        INTEGER I,NBOD

        REAL*8 MASS(NBOD),X(NBOD),Y(NBOD),Z(NBOD),
     &                    VX(NBOD),VY(NBOD),VZ(NBOD),
c     &                    XT(NPLMAX),YT(NPLMAX),ZT(NPLMAX),
c     &                    VXT(NPLMAX),VYT(NPLMAX),VZT(NPLMAX),
     &                    C1,C2,C3,C,BOT,A(3,3)

        C1 = SUM(MASS(1:NBOD)*
     &           (Y(1:NBOD)*VZ(1:NBOD)-Z(1:NBOD)*VY(1:NBOD)))
        C2 = SUM(MASS(1:NBOD)*
     &           (Z(1:NBOD)*VX(1:NBOD)-X(1:NBOD)*VZ(1:NBOD)))
        C3 = SUM(MASS(1:NBOD)*
     &           (X(1:NBOD)*VY(1:NBOD)-Y(1:NBOD)*VX(1:NBOD)))

        C = SQRT( C1*C1 + C2*C2 + C3*C3 )
        C1 = C1/C
        C2 = C2/C
        C3 = C3/C
        BOT = 1.0D0/(1.0D0 + C3)

        A(1,1) = 1.0D0 - C1*C1*BOT
        A(1,2) = -1.0D0*C1*C2*BOT
        A(1,3) = -1.0D0*C1

        A(2,1) = -1.0D0*C1*C2*BOT
        A(2,2) = 1.0D0 - C2*C2*BOT
        A(2,3) = -1.0D0*C2

        A(3,1) = C1
        A(3,2) = C2
        A(3,3) = C3


c       DO I=1,NBOD
c          XT(I) = A(1,1)*X(I) + A(1,2)*Y(I) + A(1,3)*Z(I)
c          YT(I) = A(2,1)*X(I) + A(2,2)*Y(I) + A(2,3)*Z(I)
c          ZT(I) = A(3,1)*X(I) + A(3,2)*Y(I) + A(3,3)*Z(I)
c          VXT(I) = A(1,1)*VX(I) + A(1,2)*VY(I) + A(1,3)*VZ(I)
c          VYT(I) = A(2,1)*VX(I) + A(2,2)*VY(I) + A(2,3)*VZ(I)
c          VZT(I) = A(3,1)*VX(I) + A(3,2)*VY(I) + A(3,3)*VZ(I)
c        END DO
c
C..    Check
c
c        C1 = SUM(MASS(1:NBOD)*
c     &           (YT(1:NBOD)*VZT(1:NBOD)-ZT(1:NBOD)*VYT(1:NBOD)))
c        C2 = SUM(MASS(1:NBOD)*
c     &           (ZT(1:NBOD)*VXT(1:NBOD)-XT(1:NBOD)*VZT(1:NBOD)))
c        C3 = SUM(MASS(1:NBOD)*
c     &           (XT(1:NBOD)*VYT(1:NBOD)-YT(1:NBOD)*VXT(1:NBOD)))
c        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3


        END
c-------------------------------------------------------------------------
C
C 
      SUBROUTINE FOLLOWBODIES(TSTOP,DTOUT,LCLOSE,IFLGCHK,
     &                NBSIMU,REFNS,DIRS,INPLFILE,DIRO,OUTFILE)

        include '../../sub/swift.inc'

        INTEGER*4 :: NBSIMU,REFNS,NS
        CHARACTER*(*) DIRS,INPLFILE
        CHARACTER*(*),DIMENSION(1:NBSIMU) :: DIRO,OUTFILE
        REAL*8 :: TSTOP,DTOUT
        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4, PARAMETER :: NFOLMAX = 20
        REAL*4, DIMENSION(500000,-NPLMAX:NFOLMAX) :: TT,TA,TINC,TE,
     &                                         TCAPOM,TOMEGA,TCAPM
	
        INTEGER*4 NB(-NPLMAX:NFOLMAX),NBMAX,NT
        INTEGER*4 NBC(-NPLMAX:NFOLMAX),NCG,J,JFOL
	REAL*8 XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),MASS(NPLMAX)
	REAL*8 VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX),TTINC(NTPMAX)
        REAL*8 TTX(NTPMAX),TTY(NTPMAX),TTZ(NTPMAX),TTH(NTPMAX)
        REAL*8 TTR(NTPMAX) 
	INTEGER ISTAT(NTPMAX,NSTAT),NPFOL,IDFOL(NFOLMAX),NSFOL(NFOLMAX)
        REAL*8 RSTAT(NTPMAX,NSTATR),RR,TTA(NTPMAX),TTE(NTPMAX)
	REAL*8 MOYM,XX,YY,ZZ,XP,DTT,TEV,EH,TEVII,CEPS,SEPS,EPS
	INTEGER NBOD,NTP,IERR,IFOL,ISTEP
	INTEGER IFLGCHK,IU,NLEFT,I,ID,OLOC(NPLMAX,NPLMAX)
        INTEGER IO_READ_HDR,IO_READ_LINE,IREAD
        INTEGER IO_READ_HDR_R,IO_READ_LINE_R,II

	REAL*8 T0,DT,DTDUMP,DR,TTCAPOM(NTPMAX),IMOY
	REAL*8 T,TMAX,EMAX(NTPMAX),VARPI,KR,KT,RCOLL,AMOY,EMOY
        REAL*8 MAT(NPLMAX,NPLMAX),UMAT(NPLMAX,NPLMAX),OR2,IR2
        REAL*8 ETA(NPLMAX),MU(NPLMAX),IRMOY,OMOY,KXMOY,KYMOY,KZMOY
        REAL*8 PPP(3),SIGP(3),CHI2,DQA,DQAA,RHA,SIGA

	REAL*8 RMIN,RMAX,RMAXU,QMIN,RPLSQ(NPLMAX),KX,KY,KZ,KXR,KYR
        LOGICAL ERROR,WIMPS,EVAP,COLL
        LOGICAL*2 LCLOSE,OK
        REAL*8 A,E,INC,CAPOM,OMEGA,CAPM,J2RP2,J4RP4,KZR,SIGI,SIGO
        REAL*8 ELH,ELK,ELP,ELQ,Z0,ALPHEVAP,REVAP,RADMIN,ORMOY,RHO

        CHARACTER*(CHLEN) :: DEV
        COLL=.FALSE.
c	WRITE(*,*) 'Enter name of parameter data file : '
c	READ(*,'(a)') inparfile
c        call io_init_param_hb(inparfile,t0,tstop,dt,dtout,dtdump,
c     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,dirs,gname,
c     &         outfile,fopenstat)

c        WRITE(*,*) 'Enter final time (estimated) :'
c        READ(*,*)TSTOP
        NT=NINT(TSTOP/DTOUT)
        print*,nt
        WRITE(*,*) 'One data over :'
        READ(*,*)NCG
	WRITE(*,*) ' '
c	WRITE(*,*) 'Enter name of planet data file : '
c	READ(*,'(A)') INPLFILE
        call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xj,yj,zj,
     &     vxj,vyj,vzj,rplsq,j2rp2,j4rp4)
        IU = 20
        NB = 0
        NBC = 0
        DR = 180.0D0/PI
        T=-1.0d0

        NPFOL = 0
	WRITE(*,*) ' '
	WRITE(*,*) 'Enter number of particles to follow : '
        READ(*,*)NPFOL
        IF (NPFOL.GT.NFOLMAX) THEN
          WRITE(*,*)'Maximum number of particles allowed = ',NFOLMAX
          STOP
        END IF
        IF (NPFOL.GT.0) THEN
  	  WRITE(*,*) 'Enter SIMU/ID''s of particles to follow : '
          DO I=1,NPFOL
            READ(*,*)NSFOL(I),IDFOL(I)
          END DO
         END IF


        IF (BTEST(IFLGCHK,0)) THEN
           WRITE(*,*) ' Reading an INTEGER*2 binary file '
        ELSE IF(BTEST(IFLGCHK,1)) THEN
           WRITE(*,*) ' Reading a REAL*4 binary file '
        ELSE
           WRITE(*,*) ' Error: no binary file format specified '
           WRITE(*,*) '        in param file '
           STOP
        ENDIF
        NTP=0

        OPEN(UNIT=IU, FILE=TRIM(DIRO(REFNS))//'/'//OUTFILE(REFNS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
        IREAD = 0
        MAIN_LOOP: DO WHILE (T.LT.TSTOP)
          IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
            IERR = IO_READ_HDR(IU,T,NBOD,NLEFT) 
          ELSE
            IERR = IO_READ_HDR_R(IU,T,NBOD,NLEFT) 
          ENDIF
          IF (IERR.NE.0) EXIT MAIN_LOOP
          if (mod(iread,500).eq.0) print*,'T = ',t,' NBOD=',nbod
          iread = iread+1
c          write(*,*)sngl(T),nleft
          DO I=2,NBOD
            IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
              IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM) 
            ELSE
              IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM) 
            ENDIF
            IF(IERR.NE.0) EXIT MAIN_LOOP
            IF (NBC(ID).EQ.0) THEN
              NB(ID)=NB(ID)+1
              TT(NB(ID),ID)=SNGL(T)
              TA(NB(ID),ID)=SNGL(A)
              TE(NB(ID),ID)=SNGL(E)
              TINC(NB(ID),ID)=SNGL(INC*DR)
              TCAPOM(NB(ID),ID)=SNGL(CAPOM*DR)
              TOMEGA(NB(ID),ID)=SNGL(OMEGA*DR)
              TCAPM(NB(ID),ID)=SNGL(CAPM*DR)
            END IF
            NBC(ID)=NBC(ID)+1
            IF (NBC(ID).EQ.NCG) NBC(ID)=0
          END DO
          DO I=1,NLEFT
            IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
              IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM) 
            ELSE
              IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM) 
            ENDIF
            IF(IERR.NE.0) EXIT MAIN_LOOP
          END DO
        END DO MAIN_LOOP
        CLOSE(IU)

c        open(unit=15,file='particles.dat',status='unknown')
        SIMLOOP1 : DO NS = 1,NBSIMU
          OPEN(UNIT=IU, FILE=trim(diro(ns))//'/'//OUTFILE(NS),
     &                    STATUS='OLD',FORM='UNFORMATTED')
          IREAD = 0
          T = -1.0d0
          TIMELOOP1 : DO WHILE(T.LT.TSTOP)
            IF (BTEST(IFLGCHK,0))  THEN        ! Bit 0 i set = Integer*2
              IERR = IO_READ_HDR(IU,T,NBOD,NLEFT)
            ELSE
              IERR = IO_READ_HDR_R(IU,T,NBOD,NLEFT)
            ENDIF
            IF (IERR.NE.0) EXIT TIMELOOP1
            if (mod(iread,500).eq.0) print*,'Simu',NS,'T = ',t
            IREAD = IREAD+1
            BODLOOP1 : DO I=2,NBOD
              IF (BTEST(IFLGCHK,0))  THEN    ! Bit 0 i set = Integer*2
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ENDIF
              IF(IERR.NE.0) EXIT BODLOOP1
            END DO BODLOOP1
            NTP=0
            LEFTLOOP1 : DO I=1,NLEFT
              IF(BTEST(IFLGCHK,0))  THEN ! BIT 0 IS SET
                IERR = IO_READ_LINE(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              ELSE
                IERR = IO_READ_LINE_R(IU,ID,A,E,INC,CAPOM,OMEGA,CAPM)
              END IF
              IF(IERR.NE.0) EXIT LEFTLOOP1
c              if ((t.gt.0.995*tstop).and.(e.gt.0.8)) then
c                 write(15,*)ns,id,sngl(t),sngl(a),sngl(e)
c              end if
              DO J=1,NPFOL
                IF ((ID.EQ.IDFOL(J)).AND.(NS.EQ.NSFOL(J))) THEN 
                  JFOL = J
                  IF (NBC(JFOL).EQ.0) THEN
                    NB(JFOL)=NB(JFOL)+1
                    TT(NB(JFOL),JFOL)=SNGL(T)
                    TA(NB(JFOL),JFOL)=SNGL(A)
                    TE(NB(JFOL),JFOL)=SNGL(E)
                    TINC(NB(JFOL),JFOL)=SNGL(INC*DR)
                    TCAPOM(NB(JFOL),JFOL)=SNGL(CAPOM*DR)
                    TOMEGA(NB(JFOL),JFOL)=SNGL(OMEGA*DR)
                    TCAPM(NB(JFOL),JFOL)=SNGL(CAPM*DR)
                  END IF
                  NBC(JFOL)=NBC(JFOL)+1
                  IF (NBC(JFOL).EQ.NCG) NBC(JFOL)=0
                END IF
              END DO
            END DO LEFTLOOP1
          END DO TIMELOOP1
          CLOSE(IU)
        END DO SIMLOOP1
        close(15)
        NBMAX=MAXVAL(NB)
        print*,'nbmax=',nbmax,nb
        ALLOCATE(DATA4(NBMAX,7,NBOD+NPFOL))
        DATA4=0.
        DO I=-2,-NBOD,-1
          DATA4(1:NB(I),1,-I-1)=TT(1:NB(I),I)
          DATA4(1:NB(I),2,-I-1)=TA(1:NB(I),I)
          DATA4(1:NB(I),3,-I-1)=TE(1:NB(I),I)
          DATA4(1:NB(I),4,-I-1)=TINC(1:NB(I),I)
          DATA4(1:NB(I),5,-I-1)=TCAPOM(1:NB(I),I)
          DATA4(1:NB(I),6,-I-1)=TOMEGA(1:NB(I),I)
          DATA4(1:NB(I),7,-I-1)=TCAPM(1:NB(I),I)
        END DO
        DATA4(1:NB(0),1,NBOD)=TT(1:NB(0),0)
        DATA4(1:NB(0),2,NBOD)=TA(1:NB(0),0)
        DATA4(1:NB(0),3,NBOD)=TE(1:NB(0),0)
        DATA4(1:NB(0),4,NBOD)=TINC(1:NB(0),0)
        DATA4(1:NB(0),5,NBOD)=TCAPOM(1:NB(0),0)
        DATA4(1:NB(0),6,NBOD)=TOMEGA(1:NB(0),0)
        DATA4(1:NB(0),7,NBOD)=TCAPM(1:NB(0),0)
        DO I=1,NPFOL
          print*,I,NB(I)
          DATA4(1:NB(I),1,NBOD+I)=TT(1:NB(I),I)
          DATA4(1:NB(I),2,NBOD+I)=TA(1:NB(I),I)
          DATA4(1:NB(I),3,NBOD+I)=TE(1:NB(I),I)
          DATA4(1:NB(I),4,NBOD+I)=TINC(1:NB(I),I)
          DATA4(1:NB(I),5,NBOD+I)=TCAPOM(1:NB(I),I)
          DATA4(1:NB(I),6,NBOD+I)=TOMEGA(1:NB(I),I)
          DATA4(1:NB(I),7,NBOD+I)=REAL(NB(I))
        END DO

c        print*,DATA4(1:NB(0),7,NBOD)
        ! DEV = TRIM(DIRS)//'/followbodies.dat'
        DEV = './followbodies.dat'

c        CALL IMWRITE(DEV,DATA4,NBMAX,DBLE(NBOD),0.,0.,7,
c     &                    0.,0.,0.,NBOD,0.,0.,0.,
c     &                    SNGL(-10000.),SNGL(0.001),0.,0.,0.,ERROR)
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)NBMAX,7,NBOD+NPFOL
        WRITE(18,*)DATA4
        CLOSE(18)
        DEALLOCATE(DATA4)

        END
C--------------------------------------------------------------------
C               Keplerian orbit stuff
C--------------------------------------------------------------------
C
        SUBROUTINE KEPLER_SMALL(A,E,INC,CAPOM,OMEGA,CAPM,CX,CY,CZ)

        include '../../sub/swift.inc'
     
        REAL*8 GM,A,E,INC,CAPOM,OMEGA,CAPM
        REAL*8 R,RX,RY,CX,CY,CZ,VRX,VRY,VX,VY,VZ
        REAL*8 EX,EY,EZ,EPX,EPY,EPZ,U,F,F1,F2,F3,DELA
        LOGICAL OK
        INTEGER I

        EX = COS(OMEGA)*COS(CAPOM)-SIN(OMEGA)*COS(INC)*SIN(CAPOM)
        EY = COS(OMEGA)*SIN(CAPOM)+SIN(OMEGA)*COS(INC)*COS(CAPOM)
        EZ = SIN(OMEGA)*SIN(INC)
        EPX = -SIN(OMEGA)*COS(CAPOM)-COS(OMEGA)*COS(INC)*SIN(CAPOM)
        EPY = -SIN(OMEGA)*SIN(CAPOM)+COS(OMEGA)*COS(INC)*COS(CAPOM)
        EPZ = COS(OMEGA)*SIN(INC)
        U = CAPM
        OK = .FALSE.
        I=0
        DO WHILE (.NOT.OK)
          F = U-E*SIN(U)-CAPM
          F1 = 1.-E*COS(U)
          F2 = U-CAPM-F
          F3 = 1.-F1
          DELA = -F/F1
          DELA = -F/(F1+0.5*DELA*F2)
          DELA = -F/(F1+0.5*DELA*F2+DELA*DELA*F3/6.)
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.1e-12).OR.(I.GE.500))
        END DO
        R = A*(1.-E*COS(U))
        RX = A*(COS(U)-E)
        RY = A*SQRT(1.-E*E)*SIN(U)
c        VRX = -SQRT(GM/A)*A/R*SIN(U)
c        VRY = SQRT(GM/A)*A/R*COS(U)
        CX = EX*RX+EPX*RY  
        CY = EY*RX+EPY*RY 
        CZ = EZ*RX+EPZ*RY
c        VX = EX*VRX+EPX*VRY  
c        VY = EY*VRX+EPY*VRY 
c        VZ = EZ*VRX+EPZ*VRY
        END 

C
