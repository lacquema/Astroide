C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C 
        PROGRAM GEN_MULTI_RMVS3

        include '../../sub/swift.inc'

        REAL*8, PARAMETER :: SMASSYR = TWOPI*TWOPI,
     &          DR = 1.7453292519943294d-2,      ! 180/pi 
     &          EPOCH = 2449101.0d0,
     &          YEAR = 365.2422d0, 
     &          AU = 1.495978707d8               ! km

        REAL*8 MASS(NPLMAX),PERIOD,JDPERI,Q,RPLSQ(NPLMAX),
     &         XH(NPLMAX),YH(NPLMAX),ZH(NPLMAX),
     &         VXH(NPLMAX),VYH(NPLMAX),VZH(NPLMAX),
     &         XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX),
     &         VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX),EBOD(NPLMAX),
     &         XCM,YCM,ZCM,DT,DTTRY,DTMIN,DT0,ABOD(NPLMAX),
     &         VXCM,VYCM,VZCM,FDT,X,CPUTIME,CPUREF,
     &         MASSO(NPLMAX),RPLSQO(NPLMAX),J2RP2,J4RP4,
     &         XHO(NPLMAX),YHO(NPLMAX),ZHO(NPLMAX),
     &         VXHO(NPLMAX),VYHO(NPLMAX),VZHO(NPLMAX),
     &         XBO(NPLMAX),YBO(NPLMAX),ZBO(NPLMAX),
     &         VXBO(NPLMAX),VYBO(NPLMAX),VZBO(NPLMAX),
     &         MSYS,FACHILL,R2HILL(NPLMAX),DTBOD(NPLMAX)

	INTEGER*4	NDATA          ! Taille du tableau de donnees.
        REAL*8 XHT(NTPMAX),YHT(NTPMAX),ZHT(NTPMAX),
     &         VXHT(NTPMAX),VYHT(NTPMAX),VZHT(NTPMAX),
     &         BETAHT(NTPMAX),BETAMIN,BETAMAX,GME,
     &         GM,A(10*NTPMAX),E,INC,CAPOM,OMEGA,CAPM,AMIN,AMAX,DA,
     &         RSTAT(NTPMAX,NSTATR),MSUN,RMIN,RMAX,RMIN25,RMAX25,
     &         CIMAX,IMAX,DE,EMAX,EMIN,WIMPMFP(NTPMAX),rand
        REAL*4 SEED

        INTEGER ISTAT(NTPMAX,NSTAT),IPOINT,IPAR
        INTEGER IT(10*NTPMAX),MT(10*NTPMAX),IFILE
        INTEGER NTP,NTPFIRST,I,J,K,ISEED,LLM,IND,NTPT,CNTFILE,NTPFILE
        LOGICAL OK,WIMPS,OKC
        CHARACTER*(CHLEN) INTPFILE,TPFILE,INPARFILE,PARFILE, OUTFILE
        CHARACTER*(CHLEN) DIRO,DIRS,GNAME,gname_old,FOPENSTAT,GOFILE
        CHARACTER*(CHLEN) MVSFILE
        CHARACTER*(CHLEN) BINDIR,MEXTRFILE,CONTFILE,STARTFILE, str
        CHARACTER*(CHLEN) GOCFILE,GOCMD,INTEG, CLEARFILE, CORRFILE
        CHARACTER*(CHLEN) OPTIONFILE, STATESFILE
        CHARACTER*1 PR
        real*8 t0,tstop
        real*8 dtout,dtdump
        real*8 rmaxu,qmin

C
C...  Inverse mass of planets in 1/solar masses and physical radii in km
C
        INTEGER, PARAMETER :: NPLR = 9           ! Nombre de planetes
        REAL*8, DIMENSION(NPLR), PARAMETER ::
     &                      RMASS = (/ 6023600.0d0,     ! Mercury
     &                                 408524.0d0,      ! Venus
     &                                 328900.0d0,      ! Earth
     &                                 3098710.0d0,     ! Mars
     &                                 1047.355d0,      ! Jupiter
     &                                 3498.5d0,        ! Saturn
     &                                 22869.0d0,       ! Uranus
     &                                 19424.0d0,       ! Neptune
     &                                 135300000.0d0 /),! Pluto
     &                      RPHY =  (/ 2339.0d0,        ! Mercury
     &                                 6052.0d0,        ! Venus
     &                                 6378.0d0,        ! Earth
     &                                 3397.0d0,        ! Mars
     &                                 71390.0d0,       ! Jupiter
     &                                 60000.0d0,       ! Saturn
     &                                 25400.0d0,       ! Uranus
     &                                 25225.0d0,       ! Neptune
     &                                 600.0d0 /)       ! Pluto

        LOGICAL*2 LCLOSE
        INTEGER IFLGCHK,NCOR
        CHARACTER*1 :: REP
        CHARACTER*32 NAME
        CHARACTER*10 INPLFILE
        CHARACTER*65 FPLANETS
        INTEGER NBOD,IALPHA,IP1,IP2,IUFLG,ICFLG,IRFLG
        INTEGER NNN1,NNN

        CHARACTER*256 LINE, TRIMMED_LINE
        INTEGER POS_COM, pos_slash
        CHARACTER*100 dirs_temp, diro_temp, gname_temp
        CHARACTER*1 BOOLS(6)


0009    FORMAT('    rm -rf ./run_',I2.2) ! added by antoine

0001    FORMAT((a),'/mbodies_multi_rmvs3 <<!') ! added by antoine
2204    FORMAT((a),'/swift_',(a),' < ./run_',I2.2,'/files.in') ! modified by antoine
0002    FORMAT((a),'/swift_',(a),' < ',(a),'/',(a),'/', ! added by antoine
     &                'run_',I2.2,'/files_',I2.2,'.in')
0003    FORMAT((a),'/swift_',(a),' < ./files.in') ! added by antoine
0004    FORMAT((a),'/swift_',(a),' < ./files_',I2.2,'.in') ! added by antoine
0005    FORMAT((a),'/swift_',(a),' < ',(a),'/',(a), ! added by antoine
     &                '/files_',I2.2,'.in')
0006    FORMAT('start_',I2.2,'.sh')
0007    FORMAT('continue_',I2.2,'.sh')
0008    FORMAT('cd ',(a)) ! added by antoine

1000    FORMAT(I6,6(1X,F7.3))   !étiquettes de format
2000    FORMAT((a),'_',I2.2,'.in')
3000    format(100(a1,1x))
4000    FORMAT((a),'_',I2.2)
5000    FORMAT('files_',I2.2,'.in')
6000    FORMAT('go_',(a),'_',I2.2,'.sh')
6001    FORMAT('cont_',(a),'_',I2.2,'.sh')
2200    FORMAT((a),'/run_',I2.2,'/',(a),'_',I2.2,'.in') ! modified by antoine
2201    FORMAT(I0) ! modified by antoine
2202    FORMAT('tfin="',f10.1,'"')  
22021    FORMAT('t_0=${1:-"',f0.0,'"}')   ! added by antoine
22022    FORMAT('t_f=${2:-"',f0.0,'"}')   ! added by antoine
22023    FORMAT('dt=${3:-"',f0.0,'"}')   ! added by antoine
2203    FORMAT(f10.1)

2205    FORMAT('order=${1:-"oarsub -l /nodes=1/cpu=1/core=',I1,
     &       ',walltime=48 --project dynapla"}') ! modified by antoine
2206    FORMAT('export OMP_NUM_THREADS=',I1)
c     ........entree des parametres

        READ(*,'(a)')BINDIR  ! added by antoine
        

        MEXTRFILE='mextract_multi.sh'
        CONTFILE='continues.sh'
        STARTFILE='starts.sh'
        CLEARFILE='clear.sh'
        CORRFILE='corr_multi.dat'
        OPTIONFILE='options.in'
        STATESFILE='states.sh'

        WRITE(*,*)'Remove old files' ! added by antoine
        CALL SYSTEM('[ -f ./'//CLEARFILE//' ] && ./'//CLEARFILE) ! added by antoine

        WRITE(*,*)'Generate sub-simulations:' ! added by antoine

        OPEN(51,FILE=MEXTRFILE,STATUS='UNKNOWN')
        OPEN(52,FILE=CONTFILE,STATUS='UNKNOWN')
        OPEN(53,FILE=STARTFILE,STATUS='UNKNOWN')  ! added by antoine
        OPEN(54,FILE=CLEARFILE,STATUS='UNKNOWN')  ! added by antoine
        OPEN(55,FILE=OPTIONFILE,STATUS='UNKNOWN')  ! added by antoine

        WRITE(53,'(a)')'#! /bin/bash' ! added by antoine
        WRITE(54,'(a)')'#! /bin/bash' ! added by antoine

        WRITE(51,'(a)')'#! /bin/bash' ! modified by antoine
        ! WRITE(51,'(a)')'unlimit'
        WRITE(51,'(a)')'export OMP_NUM_THREADS=1'
        WRITE(51,'(a)')'export STACKSIZE=1000000'
        
        WRITE(*,*)' Name of integrator'
        WRITE(*,*)' (rmvs3,whm,whm_s6b,rmvs3_parp,whm_s6b_parp) '
        READ(*,'(A)')INTEG


        ! Ne garder que ce qui est avant le '#' dans INTEG:     ! added by antoine
        pos_com = INDEX(INTEG, '#')

        IF (pos_com > 0) THEN 
          INTEG = ADJUSTL(TRIM(INTEG(1:pos_com-1)))
        ELSE
          INTEG = ADJUSTL(TRIM(INTEG))
        ENDIF
        ! Supprimer les tabulations et espaces restants
        DO I = LEN(INTEG), 1, -1
          IF (INTEG(I:I) .EQ. CHAR(9) .OR. INTEG(I:I) .EQ. '\t') THEN
            INTEG(I:I) = ''
          ENDIF
        END DO

        ! Ecrire le fichier options.in  ! added by antoine
        DO I = 1, 3 ! 3 premieres lignes d'options
          READ(*, '(A)') LINE
          POS_COM = INDEX(LINE, '#')
          IF (POS_COM > 0) THEN
            TRIMMED_LINE = TRIM(LINE(1:POS_COM-1))
          ELSE
            TRIMMED_LINE = TRIM(LINE)
          END IF
          WRITE(55, '(A)') TRIM(TRIMMED_LINE)
          WRITE(*, '(A)') TRIM(TRIMMED_LINE)
          IF (I .EQ. 3) THEN ! si dans la 3e
            READ(TRIMMED_LINE, *) (BOOLS(J), J=1,6)
            IF (BOOLS(2).EQ.'T') THEN ! si le bool associe au removal des limites, alors ecrire la ligne suivante
              WRITE(*, '(A)') 'Removing limits activated'
              READ(*, '(A)') LINE
              POS_COM = INDEX(LINE, '#')
              IF (POS_COM > 0) THEN
                TRIMMED_LINE = TRIM(LINE(1:POS_COM-1))
              ELSE
                TRIMMED_LINE = TRIM(LINE)
              END IF
              WRITE(55, '(A)') TRIM(TRIMMED_LINE)
              WRITE(*, '(A)') TRIM(TRIMMED_LINE)
            END IF
          END IF
        END DO
        READ(*, '(A)') LINE ! directory and sub-directory
        POS_COM = INDEX(LINE, '#')
        IF (POS_COM > 0) THEN
          TRIMMED_LINE = TRIM(LINE(1:POS_COM-1))
        ELSE
          TRIMMED_LINE = TRIM(LINE)
        END IF
        diro_temp = TRIM(TRIMMED_LINE)
        WRITE(*,*) 'diro_temp : ', diro_temp
        pos_slash = INDEX(diro_temp, '/', .TRUE.)  ! .TRUE. pour chercher depuis la fin
        if (pos_slash > 0) then
          dirs_temp = TRIM(diro_temp(1:pos_slash-1))  ! Extraire la partie avant le dernier '/'
          gname_temp = TRIM(diro_temp(pos_slash+1:))  ! Extraire la partie après le dernier '/'
        else
          write(*,*) 'Erreur : "/" non trouvé dans diro'
          dirs_temp = ''
          gname_temp = ''
        endif
        WRITE(55, '(A)') dirs_temp
        WRITE(*, '(A)') dirs_temp
        WRITE(55, '(A)') gname_temp
        WRITE(*, '(A)') gname_temp
        outfile = 'simulation'
        WRITE(55, '(A)') outfile
        WRITE(55, '(A)') 'new'
        CLOSE(55)
        CALL SYSTEM('chmod ogu+x '//TRIM(OPTIONFILE))  ! added by antoine

        

        WRITE(*,*)' Number of cores'
        READ(*,*)NCOR

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Units menu (G=1) : '
          WRITE(*,*) '       0 ==> Solar masses and AU '
          WRITE(*,*) '       1 ==> AU and years '
          READ(*,*) IUFLG
          OK = ((IUFLG.EQ.0).OR.(IUFLG.EQ.1))
        END DO

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Coordinate Menu: '
          WRITE(*,*) '       0 ==> Ecliptic '
          WRITE(*,*) '       1 ==> Invariable plane '
          READ(*,*) ICFLG
          OK = ((ICFLG.EQ.0).OR.(ICFLG.EQ.1))
        END DO

        OK=.FALSE.
        DO WHILE(.NOT.OK)
          WRITE(*,*) ' Planetary radius menu: '
          WRITE(*,*) '       0 ==> Do not include a radius '
          WRITE(*,*)
     &'       1 ==> Include the physical radius of planet '
          WRITE(*,*)
     &'       2 ==> Include a multiple of the Hill radius ' 
          READ(*,*) IRFLG
          OK=((IRFLG.EQ.0).OR.(IRFLG.EQ.1).OR.(IRFLG.NE.2))
        END DO

        IF(IRFLG.EQ.2) THEN
          WRITE(*,*) ' What is that multiple ? '
          READ(*,*) FACHILL
        END IF

        WRITE(*,*)' Give mass of the star in solar masses '
        READ(*,*)MASS(1)
        WRITE(*,*)'Give number of massive planets '
        READ(*,*)NBOD
        NBOD=NBOD+1

        IF(IUFLG.EQ.1) MASS(1) = MASS(1)*SMASSYR

        WRITE(*,*) ' Mass of the Sun is : ',MASS(1)

        
  
        XH(1) = 0.0
        YH(1) = 0.0
        ZH(1) = 0.0
        VXH(1) = 0.0
        VYH(1) = 0.0
        VZH(1) = 0.0

        DO I=2,NBOD
          WRITE(*,*)' Information relative to planet #',I-1,':'
          WRITE(*,*)' Give mass in Jupiter masses '
          READ(*,*)MASS(I)
          MASS(I)=MASS(I)/RMASS(5) ! donne la masse de la planète en Msol
          WRITE(*,*) 'Give now the orbital parameters '
          READ(*,*)ABOD(I),EBOD(I),INC,CAPOM,OMEGA,CAPM  
          print*,ABOD(I),EBOD(I),INC,CAPOM,OMEGA,CAPM
          IF(IUFLG.EQ.1) MASS(I) = MASS(I)*SMASSYR
          CAPM = DR*CAPM    ! conversion des angles en degrés
          OMEGA = DR*OMEGA
          CAPOM = DR*CAPOM
          INC = DR*INC
          GM = MASS(1)+MASS(I)  
          IALPHA = -1 
          CALL ORBEL_EL2XV(GM,IALPHA,ABOD(I),EBOD(I),INC,CAPOM,
     &          OMEGA,CAPM,XH(I),YH(I),ZH(I),VXH(I),VYH(I),VZH(I))
        END DO

C...  Calc the radii
        RPLSQ(1) = 0.0D0   
        CALL UTIL_HILLS(NBOD,MASS,XH,YH,ZH,VXH,VYH,VZH,R2HILL)
        IF(IRFLG.EQ.0) THEN !IRFLG sélection type rayon 
          LCLOSE = .FALSE.
          RPLSQ(1:NBOD) = 0.d0 ! car pas de rayon
        ELSE IF(IRFLG.EQ.1) THEN
          LCLOSE = .TRUE.
          DO I=2,NBOD
            RPLSQ(I) = (RPHY(I-1)/AU)**2 
c ! ici on utilise la table de rayons donnée au début, ms que fait-on? **2 ?
          END DO
        ELSE 
          LCLOSE = .TRUE.
c           subroutine ds le cas où R est un multiple du rayon de Hill
          RPLSQ(2:NBOD) = (FACHILL*SQRT(R2HILL(2:NBOD)))**2
        ENDIF

        CALL COORD_H2B(NBOD,MASS,XH,YH,ZH,VXH,VYH,VZH,
     &     XB,YB,ZB,VXB,VYB,VZB,MSYS)

        IF(ICFLG.EQ.1) CALL INVAR(NBOD,MASS,XB,YB,ZB,VXB,VYB,VZB)
c         en lien avec le choix du réf de coord (écliptique ou invariant)

        CALL COORD_B2H(NBOD,MASS,XB,YB,ZB,VXB,VYB,VZB, 
     &      XH,YH,ZH,VXH,VYH,VZH)

	      WRITE(*,*) ' '
	      WRITE(*,*) 'Enter name of planet data file : '
	      ! READ(*,'(a)')INPLFILE  ! commented by antoine
        INPLFILE = 'bodies.in'  ! added by antoine
        WRITE(*,*)  INPLFILE
        IFLGCHK = 0
        J2RP2 = 0.0D0
        J4RP4 = 0.0D0
	      CALL IO_DUMP_PL(INPLFILE,NBOD,MASS,XH,YH,ZH,  
     &       VXH,VYH,VZH,LCLOSE,IFLGCHK,RPLSQ,J2RP2,J4RP4)

        

        WIMPMFP=0.0d0


        OKC = .TRUE.
        WRITE(*,*) ' Input ISEED: large and odd '
        READ(*,*) ISEED
        SEED=FLOAT(ISEED)
        CALL RANDOM_SEED()
        print*,seed
        WRITE(*,*) 'Enter name of test particle data file : '
        ! READ(*,'(a)')INTPFILE   ! commented by antoine   
        INTPFILE = 'debris.in'  ! added by antoine   
        IPOINT = INDEX(INTPFILE,'.')
        WRITE(*,*) 'Enter name of parameter file : '
        ! READ(*,'(a)')INPARFILE
        INPARFILE = OPTIONFILE   ! added by antoine  
         
        call io_init_param_hb(inparfile,t0,tstop,dtmin,dtout,
     &         dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,dirs,
     &         gname,outfile,fopenstat)

        ! WRITE(*,*)TRIM(DIRO)
        ! WRITE(*,*)TRIM(DIRS)
        ! WRITE(*,*)TRIM(gname)

        

        DO I=2,NBOD
          DTBOD(I) = DTMIN*2.d0*PI*SQRT(ABOD(I)**3/(MASS(I)+MASS(1)))
          DTBOD(I) = DTBOD(I)*(1.d0-EBOD(I))**1.5d0/SQRT(1.d0+EBOD(I))
          print*,'dtmbod',dtbod(i),abod(i),mass(1)
        END DO

        
        WRITE(52,2205)NCOR ! added by antoine
        WRITE(53,2205)NCOR ! added by antoine

        ! WRITE(53,'(a)')'success_count=0' ! added by antoine

        ! WRITE(51,'(a)')'simname="'//TRIM(GNAME)//'"'
        WRITE(PARFILE,22021)SNGL(T0) ! added by antoine
        WRITE(51, '(a)')TRIM(PARFILE) ! added by antoine
        WRITE(PARFILE,22022)SNGL(TSTOP)
        WRITE(51,'(a)')TRIM(PARFILE)
        dt = 0.1d0*TSTOP ! added by antoine
        WRITE(PARFILE,22023)SNGL(dt) ! added by antoine
        WRITE(51, '(a)')TRIM(PARFILE) ! added by antoine

        ! WRITE(51,0008)TRIM(DIRO)  ! added by antoine
        
        WRITE(51,0001)TRIM(BINDIR)

        IPAR = INDEX(INPARFILE,'.')
        CNTFILE = 0
        CPUTIME = 0.d0
        NTPT = 0


        DO WHILE(OKC)
          WRITE(*,*) ' Input size of grid for particles '
          READ(*,*) NTP
          OKC = (NTP.NE.0)
          IF (OKC) THEN
            WRITE(*,*) ' Input number of tp of first file '
            READ(*,*) NTPFIRST
            WRITE(*,*) ' Input limits for eccentricities'
            READ(*,*)RMIN25,RMAX25
            WRITE(*,*) ' Input max inc of particles (in deg) '
            READ(*,*) IMAX
            WRITE(*,*) ' Input min and max values of A'
            READ(*,*) AMIN,AMAX
            DT0 = DTMIN*2.d0*PI*SQRT(AMIN**3/MASS(1))
            
            E=EMIN
            CIMAX = COS(IMAX/DEGRAD)
            GM = MASS(1)
            WRITE(*,*) ' Test particles: '
            WRITE(*,*) '   A      E      I    CAPOM  OMEGA    M '
            DO I=1,NTP
              call random_number(rand)
              A(I) = AMIN + (AMAX - AMIN)*rand
            END DO
c            CALL ORDREB(A,IT,MT,NTP)
            IFILE = 1
            OPEN(12,FILE=CORRFILE,status='UNKNOWN')
            DO I=1,NTP
              IALPHA = -1
              call random_number(rand)
              INC = ACOS( 1.0D0 - (1.0D0-CIMAX)*rand)
              call random_number(rand)
              CAPM = 0.*TWOPI*rand
              call random_number(rand)
              OMEGA = TWOPI*rand
              call random_number(rand)
              CAPOM = TWOPI*rand
              call random_number(rand)
              J = NTPT+1
              call random_number(rand) 
              E = RMIN25 + (RMAX25-RMIN25)*RAND
c              WRITE(*,1000) J,A,E,INC,CAPOM,OMEGA,CAPM
c              WRITE(12,*) A(IT(I)),E,INC,CAPOM,OMEGA,CAPM
              GME=GM
              CALL ORBEL_EL2XV(GME,IALPHA,A(I),E,INC,CAPOM,
     &          OMEGA,CAPM,XHT(J),YHT(J),ZHT(J),VXHT(J),
     &          VYHT(J),VZHT(J))
c              WRITE(*,1000) I,A,E,INC,XHT(J),YHT(J),ZHT(J)
              DO K=1,NSTAT
                ISTAT(J,K) = 0
              END DO
              DO K=1,NSTATR
                RSTAT(J,K) = 0.0D0
              END DO
              NTPT = NTPT+1 
c Estimation du pas de temps a une distance donnee
              DTTRY = DT0*(A(I)/AMIN)**1.5d0
              DO K=2,NBOD
                DTTRY = DTTRY*(1.d0
     &         -0.9d0*EXP(-0.5d0*(A(I)-ABOD(K))**2/R2HILL(K)))
              END DO
              DO K=2,NBOD
                DTTRY = MIN(DTTRY,DTBOD(K))
              END DO
              DTTRY = DTOUT/NINT(DTOUT/DTTRY)
c
              IF ((NTPT.EQ.1).OR.(DTTRY.LT.DT)) THEN
                DT = DTTRY
              END IF
              CPUTIME = CPUTIME + 1.d0/DT
c
              IF (CNTFILE.EQ.0) THEN
                CPUREF = CPUTIME
                OK = (NTPT.EQ.NTPFIRST)
              ELSE
                OK = ((CPUTIME.GT.CPUREF).OR.(I.EQ.NTP))
              END IF
              IF (OK) THEN
                CNTFILE = CNTFILE+1
                IF (CNTFILE.GT.99) STOP
                IFILE = IFILE+1
                WRITE(TPFILE,2000)INTPFILE(1:IPOINT-1),CNTFILE
                CALL IO_DUMP_TP(TPFILE,NTPT,XHT,YHT,ZHT,
     &              VXHT,VYHT,VZHT,ISTAT,RSTAT)
                WRITE(PARFILE,2000)INPARFILE(1:IPAR-1),CNTFILE
                WRITE(str,4000)'/run',CNTFILE ! modified by antoine
                gname_old = gname
                gname = TRIM(gname) // TRIM(str) ! modified by antoine
                CALL io_dump_param_hb(parfile,t0,tstop,dt,
     &           dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &           lclose,dirs,gname,outfile,fopenstat)
                gname = gname_old ! added by antoine
                WRITE(MVSFILE,5000)CNTFILE
                OPEN(31,FILE=MVSFILE,STATUS='UNKNOWN') ! mvs file
                WRITE(31,'(a)')'gen_multi_rmvs3.sh'
                WRITE(31,'(a)')TRIM(PARFILE)
                WRITE(31,'(a)')TRIM(INPLFILE)
                WRITE(31,'(a)')TRIM(TPFILE)
                WRITE(31,'(a)')'1d-8'
                CLOSE(31)
                ! WRITE(GOFILE,6000)TRIM(INTEG),CNTFILE  ! commented by antoine
                ! WRITE(GOCFILE,6001)TRIM(INTEG),CNTFILE  ! commented by antoine
                WRITE(GOFILE,0006)CNTFILE  ! added by antoine
                WRITE(GOCFILE,0007)CNTFILE  ! added by antoine
                OPEN(31,FILE=GOFILE,STATUS='UNKNOWN') ! go file
                WRITE(31,'(a)')'#! /bin/bash' ! modified by antoine
                ! WRITE(31,'(a)')'unlimit' ! commented by antoine
                WRITE(GOCMD,2206)NCOR
                WRITE(31,'(a)')TRIM(GOCMD)
                WRITE(31,'(a)')'export STACKSIZE=1000000'
                ! WRITE(31, 0008)TRIM(DIRO)      ! added by antoine
                WRITE(31,'(a)')TRIM(BINDIR)//'/swift_'//TRIM(INTEG)//
     &                ' < ./'//TRIM(MVSFILE) ! modified by antoine
                ! WRITE(31,0004)TRIM(BINDIR),TRIM(INTEG),CNTFILE  ! added by antoine
                CLOSE(31)           
                OPEN(32,FILE=GOCFILE,STATUS='UNKNOWN')  ! continuation file
                WRITE(32,'(a)')'#! /bin/bash' ! modified by antoine
                ! WRITE(32,'(a)')'unlimit' ! commented by antoine
                WRITE(GOCMD,2206)NCOR
                WRITE(32,'(a)')TRIM(GOCMD)
                WRITE(32,'(a)')'export STACKSIZE=1000000'
                ! WRITE(32,'(a)')' '
                ! WRITE(32, 0008)TRIM(DIRO)      ! added by antoine
                
                WRITE(GOCMD,2204)TRIM(BINDIR),TRIM(INTEG),CNTFILE
                WRITE(32,'(a)')TRIM(GOCMD)
                ! WRITE(GOCMD,0003)TRIM(BINDIR),TRIM(INTEG) ! added by antoine
                CLOSE(32)           
                CALL SYSTEM('chmod ogu+x '//TRIM(GOFILE))
                CALL SYSTEM('chmod ogu+x '//TRIM(GOCFILE))   

                WRITE(53,'(a)')'${order} ./'//TRIM(GOFILE)   ! modied by antoine
                WRITE(52,'(a)')'${order} ./'//TRIM(GOCFILE)   ! modied by antoine

                print*,tpfile,i,ntpt,
     &               sngl(a(i)),sngl(dt),sngl(cputime)
                WRITE(12,*)'Fichier ',TRIM(tpfile),' nb parts ',ntpt,
     &               'derniere part :',i,'a=',sngl(a(i)),
     &               'step :',sngl(dt),'cpu : ',sngl(cputime)
                NTPT = 0
                CPUTIME = 0.d0
                
                ! Clear file
                WRITE(54,'(a)')'rm -f ./'//TRIM(GOFILE)
                WRITE(54,'(a)')'rm -f ./'//TRIM(GOCFILE)
                WRITE(54,'(a)')'rm -f ./'//TRIM(MVSFILE)
                WRITE(54,'(a)')'rm -f ./'//TRIM(PARFILE)
                WRITE(54,'(a)')'rm -f ./'//TRIM(TPFILE)


              END IF      
            END DO
          endif
        enddo
        CLOSE(12)
        CLOSE(52)

    !     WRITE(53,'(a)')'[ "$(ls start_* 2>/dev/null | wc -l)" -eq 0 ]'//
    !  &      ' && rm -f ./start.sh' ! added by antoine
        close(53) ! added by antoine

        WRITE(51,2201)CNTFILE
        DO I=1,CNTFILE
          WRITE(PARFILE,2200)TRIM(DIRO),I,INPARFILE(1:IPAR-1),I
          WRITE(51,'(a)')TRIM(PARFILE)
        END DO
        WRITE(51,'(a)')'${t_0}'
        WRITE(51,'(a)')'${t_f}'
        ! WRITE(PARFILE,2203)SNGL(T0) ! commented by antoine
        ! WRITE(51,'(a)')TRIM(PARFILE)
        ! WRITE(PARFILE,2203)SNGL(0.1d0*TSTOP)
        WRITE(51,'(a)')'${dt}' ! modified by antoine
        ! WRITE(51,'(a)')TRIM(PARFILE) ! commented by antoine
        WRITE(51,'(a)')TRIM(DIRO)//'/run_01/'//TRIM(INPLFILE)
        WRITE(51,'(a)')'mextract'
        DO I=1,CNTFILE
          WRITE(TPFILE,2200)TRIM(DIRO),I,INTPFILE(1:IPOINT-1),I
          WRITE(51,'(a)')TRIM(TPFILE)
        END DO
        WRITE(51,'(a)')'1'
        WRITE(51,'(a)')'0'
        WRITE(51,'(a)')'!'
        ! WRITE(51,'(a)')'exit' ! commented by antoine
        ! WRITE(51,'(a)')'!'
        CLOSE(51)

        WRITE(54,'(a)')'rm -f ./'//TRIM(OPTIONFILE)
        WRITE(54,'(a)')'rm -f ./'//TRIM(INPLFILE)
        WRITE(54,'(a)')'rm -f ./'//TRIM(CONTFILE)
        WRITE(54,'(a)')'rm -f ./'//TRIM(STARTFILE)
        WRITE(54,'(a)')'rm -f ./fort.7'
        WRITE(54,'(a)')'rm -f ./OAR*'

        WRITE(54,'(a)')'if [ "$1" == "all" ]; then'
        WRITE(54,'(a)')'  if [ "$(pwd)" != "'//TRIM(DIRO)//'" ]; then'
        WRITE(54,'(a)')'    rm -rf '//TRIM(DIRO)
        WRITE(54,'(a)')'  else'
        DO I=1,CNTFILE
          WRITE(54,0009)I
        END DO
        WRITE(54,'(a)')'  fi'
        WRITE(54,'(a)')'  rm -f ./'//TRIM(MEXTRFILE)
        WRITE(54,'(a)')'  rm -f ./'//TRIM(CORRFILE)
        WRITE(54,'(a)')'  rm -f ./'//TRIM(CLEARFILE)
        WRITE(54,'(a)')'  rm -f ./'//TRIM(STATESFILE)
        WRITE(54,'(a)')'  rm -f ./jacobi.out'
        WRITE(54,'(a)')'fi'
        CLOSE(54)  ! added by antoine

        CALL SYSTEM('chmod ogu+x '//TRIM(MEXTRFILE))
        CALL SYSTEM('chmod ogu+x '//TRIM(CONTFILE))  
        CALL SYSTEM('chmod ogu+x '//TRIM(STARTFILE))  ! added by antoine
        CALL SYSTEM('chmod ogu+x '//TRIM(CLEARFILE))  ! added by antoine
        
        
        ! WRITE(*,2204)BINDIR,'1','2','3',01
        ! WRITE(*,'(a)')TRIM(DIRS)//'/${simname}'
        ! WRITE(*,'(a)')DIRO
        ! WRITE(*,0001)TRIM(BINDIR)
        ! WRITE(*,'(a)')TRIM(GOCMD)

        ! write(*,*)NCOR

        ! Création du fichier 'states.sh' pour afficher l'état des sous-simulations
        OPEN(56, FILE=STATESFILE, STATUS='UNKNOWN')
        WRITE(56,'(a)') '#! /bin/bash'
        WRITE(56,'(a)') ''
        WRITE(56,'(a)') 'WORKPATH='//TRIM(DIRO)
        WRITE(56,'(a,i0)') 'NB_SIMU=',CNTFILE
        WRITE(56,'(a)') ''
        WRITE(56,'(a)') 'for ((i=1; i<=$NB_SIMU; i++)); do'
        WRITE(56,'(a)') 'dir="$WORKPATH/run_$(printf "%02d" $i)"'
        WRITE(56,'(a)') 'file="$dir/dump_param.dat"'
        WRITE(56,'(a)') 'if [ ! -f "$file" ]; then'
        WRITE(56,'(a)') 'echo "sub-simulation $(printf "%02d" $i): 0%"'
        WRITE(56,'(a)') 'continue'
        WRITE(56,'(a)') 'fi'
        WRITE(56,'(a)') ''
        WRITE(56,'(a)') 'read a b c < "$file"'
        WRITE(56,'(a)') ''
        WRITE(56,'(a)') 'if (( $(echo "$a < $b" | bc -l) )); then'
        WRITE(56,'(a)') 'p=$(echo "($a/$b)*100" | bc -l)'
        WRITE(56,'(a)') 'p=$(echo "$p/1" | bc)'
        WRITE(56,'(a)') 'else'
        WRITE(56,'(a)') 'p=100'
        WRITE(56,'(a)') 'fi'
        WRITE(56,'(a)') 'echo "sub-simulation $(printf "%02d" $i): $p%"'
        WRITE(56,'(a)') 'done'
        CLOSE(56)
        CALL SYSTEM('chmod ogu+x states.sh')


        END PROGRAM GEN_MULTI_RMVS3
        ! end subroutine GEN_MULTI_RMVS3

C
C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ)

        include '../../sub/swift.inc'

        INTEGER I,NBOD

        REAL*8 MASS(NBOD),X(NBOD),Y(NBOD),Z(NBOD),
     &                    VX(NBOD),VY(NBOD),VZ(NBOD),
     &                    XT,YT,ZT,VXT,VYT,VZT,
     &                    C1,C2,C3,C,BOT,A(3,3)

        C1 = SUM(MASS(1:NBOD)*
     &           (Y(1:NBOD)*VZ(1:NBOD)-Z(1:NBOD)*VY(1:NBOD)))
        C2 = SUM(MASS(1:NBOD)*
     &           (Z(1:NBOD)*VX(1:NBOD)-X(1:NBOD)*VZ(1:NBOD)))
        C3 = SUM(MASS(1:NBOD)*
     &           (X(1:NBOD)*VY(1:NBOD)-Y(1:NBOD)*VX(1:NBOD)))

        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3

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

        DO I=1,NBOD
          XT = A(1,1)*X(I) + A(1,2)*Y(I) + A(1,3)*Z(I) 
          YT = A(2,1)*X(I) + A(2,2)*Y(I) + A(2,3)*Z(I) 
          ZT = A(3,1)*X(I) + A(3,2)*Y(I) + A(3,3)*Z(I) 
          X(I) = XT
          Y(I) = YT
          Z(I) = ZT
          VXT = A(1,1)*VX(I) + A(1,2)*VY(I) + A(1,3)*VZ(I) 
          VYT = A(2,1)*VX(I) + A(2,2)*VY(I) + A(2,3)*VZ(I) 
          VZT = A(3,1)*VX(I) + A(3,2)*VY(I) + A(3,3)*VZ(I) 
          VX(I) = VXT
          VY(I) = VYT
          VZ(I) = VZT
        END DO

C..    Check

        C1 = SUM(MASS(1:NBOD)*
     &           (Y(1:NBOD)*VZ(1:NBOD)-Z(1:NBOD)*VY(1:NBOD)))
        C2 = SUM(MASS(1:NBOD)*
     &           (Z(1:NBOD)*VX(1:NBOD)-X(1:NBOD)*VZ(1:NBOD)))
        C3 = SUM(MASS(1:NBOD)*
     &           (X(1:NBOD)*VY(1:NBOD)-Y(1:NBOD)*VX(1:NBOD)))
        WRITE(*,*) ' C1,C2,C3 ',C1,C2,C3

        END

c**********************************************************************
c		      ORDREB.F
c**********************************************************************
c
c  Sorting routine.  Given an array a(n) of n numbers, it sorts 
c  it by ascending order. a itself is not modified, but the output
c  is two integer arrays it(n) and mt(n)
c      it(m) = # in array a() of number classified at rank m
c      mt(i) = rank of number a(i)
c
c Author:  Michel Henon,
c Note : Modified by Herve Beust
c Date:    10/23/2000

      subroutine ordreb(a,it,mt,n)

      implicit none

      integer it(n),mt(n)
      integer ngi,k,j,m,i2,i1,ng,i,n,jmax,imax,lg
      real*8 a(n)
      logical ok

      if (n.le.0) return
      if (n.eq.1) then
        it(1) = 1
        mt(1) = 1
        return
      end if        

c classement initial par paires.
      i=1
      it(n)=n
      do while(i.lt.n)
        if (a(i).gt.a(i+1)) then
          it(i)=i+1
          it(i+1)=i
        else
          it(i)=i
          it(i+1)=i+1
        end if
        i=i+2
      end do

      if(n.eq.2) then
        do m=1,n
          i=it(m)
          mt(i)=m
        end do
        return
      end if

c... interclassement.la table it contient ng groupes de longueur lg(sauf le
c... dernier qui peut etre plus court). ces groupes sont interclasses deux
c... a deux et les groupes resultants sont transferes dans la table mt.
c... imax et jmax sont les adresses finales des groupes en cours d'inter-
c... classement. i et j sont les adresses des deux nombres presentement
c... compares dans la table it. k est l'adresse dans la nouvelle table mt.

      ng=(n+1)/2
      lg=2
      do
        imax=lg
        jmax=min0(2*lg,n)
        i=1
        j=lg+1
        k=0
        ngi=0

        do while (ng-ngi.ge.2)
          i1=it(i)
          i2=it(j)
          ok=.false.
          do while (.not.ok)
            k=k+1
            if (a(i1).gt.a(i2)) then
              mt(k)=i2
              if (j.eq.jmax) then
                do while (i.le.imax)
                  k=k+1
                  mt(k)=it(i)
                  i=i+1
                end do
                ok=.true.
              else
                j=j+1
                i2=it(j)
              end if
            else
              mt(k)=i1
              if(i.eq.imax) then
                do while (j.le.jmax)
                  k=k+1
                  mt(k)=it(j)
                  j=j+1
                end do
                ok=.true.
              else
                i=i+1
                i1=it(i)
              end if
            end if
          end do
          ngi=ngi+2

c... on passe aux deux groupes suivants
          if (ng-ngi.ge.2) then
            i=jmax+1
            j=i+lg
            imax=jmax+lg
            jmax=min0(imax+lg,n)
          end if
        end do

c... ng est impair. le dernier groupe est transfere tel quel.
        if (ng.ne.ngi) then
          i=jmax
          do while (i.lt.n)
            i=i+1
            mt(i)=it(i)
          end do
        end if

        if (ng.eq.2) then
          do m=1,n
            it(m)=mt(m)
          end do
          do m=1,n
            i=it(m)
            mt(i)=m
          end do
          exit
        end if

c... on passe a l'interclassement suivant.transfert de mt a it.
        ng=(ng+1)/2
        lg=2*lg
        imax=lg
        jmax=min0(2*lg,n)
        i=1
        j=lg+1
        k=0
        ngi=0
        do while (ng-ngi.ge.2) 
          i1=mt(i)
          i2=mt(j)
          ok=.false.
          do while (.not.ok)
            k=k+1
            if (a(i1).gt.a(i2)) then
              it(k)=i2
              if (j.eq.jmax) then
                do while (i.le.imax)
                  k=k+1
                  it(k)=mt(i)
                  i=i+1
                end do
                ok=.true.
              else
                j=j+1
                i2=mt(j)
              end if
            else
              it(k)=i1
              if (i.eq.imax) then
                do while (j.le.jmax)
                  k=k+1
                  it(k)=mt(j)
                  j=j+1
                end do
                ok=.true.
              else
                i=i+1
                i1=mt(i)
              end if
            end if
          end do
          ngi=ngi+2

          if (ng-ngi.ge.2) then
            i=jmax+1
            j=i+lg
            imax=jmax+lg
            jmax=min0(imax+lg,n)
          end if
        end do
        
        if (ng.ne.ngi) then
          i=jmax
          do while (i.lt.n)
            i=i+1
            it(i)=mt(i)
          end do
        end if

        if (ng.eq.2) then

c... formation de la table inverse mt(i).
          do m=1,n
            i=it(m)
            mt(i)=m
          end do
          exit
        end if

c... on passe a l'interclassement suivant.
        ng=(ng+1)/2
        lg=2*lg
      end do
      return
      end    ! ordreb.f
c---------------------------------------------------------------------
c************************************************************************
c                          IO_DUMP_PARAM_HB.F
c************************************************************************
c IO_DUMP_PARAM dumps out the parameters for the integration. 
c
c      Input:
c       dparfile      ==>  Name of file to write to (character*(chlen))
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            iflgchk  ==>  =0 don't run diagnostic routines
c                         !=0 run them
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c       rmaxu         ==>  maximum distance from Sun in not bound
c                                 if <0  then don't check
c                                      (real scalar)
c       qmin          ==> Smallest perihelion distance
c                                 if <0  then don't check
c                                      (real scalar)
c       lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c       outfile       ==>  Name of binary output file (character*(chlen))

c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  5/10/94 HFL

	subroutine io_dump_param_hb(dparfile,t,tstop,dt,dtout,
     &          dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,dirs,
     &           gname,outfile,fopenstat)	

	include '../../sub/swift.inc'
	include '../../sub/io/io.inc'

 4000    FORMAT((a),'_',I2.2)

c...   Inputs
	real*8 t,tstop,dt
	integer iflgchk
	real*8 dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*(*) outfile,dparfile,fopenstat,dirs,gname

c...  Internals
        character*1 lflg(0:IO_NBITS-1),cclose
        integer i,ierr

c-----
c...  Executable code 

      ! WRITE(*,*)'CNTFILE=',CNTFILE

c Open parameter data file for the dump
        call io_open(7,dparfile,'unknown','formatted',ierr)

	write(7,*) sngl(t),tstop,dt
	write(7,*) dtout,dtdump

        do i=0,IO_NBITS-1
           if(btest(iflgchk,i))  then 
              lflg(i) = 'T'
           else
              lflg(i) = 'F'
           endif
        enddo

        write(7,3000) (lflg(i),i=IO_NBITS-1,0,-1)
 3000   format(100(a1,1x))

        if(btest(iflgchk,4))  then ! bit 4 is set
           if(lclose) then
              cclose = 'T'
           else
              cclose = 'F'
           endif
           write(7,*) sngl(rmin),sngl(rmax),sngl(rmaxu),
     &                sngl(qmin),' ',cclose
        endif

        if(btest(iflgchk,0).or.btest(iflgchk,1))  then 
          ! write(7,'(a)') trim(dirs)//'/'//trim(gname) ! added by antoine
           write(7,'(a)') trim(dirs) ! commented by antoine
           write(7,'(a)') trim(gname) ! commented by antoine
           write(7,'(a)') trim(outfile) ! commented by antoine
        endif


        write(7,'(a)') trim(fopenstat)

	close(unit = 7)

	return
	end     ! io_dump_param_hb
c____________________________________________________________________________
c
c




