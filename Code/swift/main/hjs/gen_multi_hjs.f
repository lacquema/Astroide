C
C-------------------------------------------------------------------------
C                   Main program
C-------------------------------------------------------------------------
C 
        PROGRAM GEN_MULTI_HJS

        include '../../sub/swift.inc'

        REAL*8, PARAMETER :: SMASSYR = TWOPI*TWOPI,
     &          DR = 1.7453292519943294d-2,      ! 180/pi 
     &          EPOCH = 2449101.0d0,
     &          YEAR = 365.2422d0, 
     &          AU = 1.495978707d8               ! km

        REAL*8 MASS(NPLMAX),PERIOD,JDPERI,Q,RPLSQ(NPLMAX),
     &         ETA(NPLMAX),MU(NPLMAX),MAT(NPLMAX,NPLMAX),
     &         UMAT(NPLMAX,NPLMAX),MTOT,VSAT,VCEN,
     &         XJ(NPLMAX),YJ(NPLMAX),ZJ(NPLMAX),
     &         VXJ(NPLMAX),VYJ(NPLMAX),VZJ(NPLMAX),
     &         XB(NPLMAX),YB(NPLMAX),ZB(NPLMAX),
     &         VXB(NPLMAX),VYB(NPLMAX),VZB(NPLMAX),EBOD(NPLMAX),
     &         XCM,YCM,ZCM,DT,DTTRY,DTMIN,DT0,ABOD(NPLMAX),
     &         VXCM,VYCM,VZCM,FDT,X,CPUTIME,CPUREF,
     &         MASSO(NPLMAX),RPLSQO(NPLMAX),J2RP2,J4RP4,
     &         XJO(NPLMAX),YJO(NPLMAX),ZJO(NPLMAX),
     &         VXJO(NPLMAX),VYJO(NPLMAX),VZJO(NPLMAX),
     &         XBO(NPLMAX),YBO(NPLMAX),ZBO(NPLMAX),
     &         VXBO(NPLMAX),VYBO(NPLMAX),VZBO(NPLMAX),
     &         MSYS,FACHILL,R2HILL(NPLMAX),DTJ(NPLMAX)
        REAL*8 ETATPT,MATPT(NPLMAX),UMATPT(NPLMAX),XJTEMP(NPLMAX),
     &         YJTEMP(NPLMAX),ZJTEMP(NPLMAX),VXJTEMP(NPLMAX),
     &         VYJTEMP(NPLMAX),VZJTEMP(NPLMAX),MASSTEMP(NPLMAX),
     &         XBTEMP(NPLMAX),YBTEMP(NPLMAX),ZBTEMP(NPLMAX),
     &         VXBTEMP(NPLMAX),VYBTEMP(NPLMAX),VZBTEMP(NPLMAX),
     &         UMPART(NPLMAX,NPLMAX),MATP(NPLMAX,NTPMAX),
     &         UMATP(NPLMAX,NTPMAX)

	INTEGER*4	NDATA          ! Taille du tableau de donnees.
        REAL*8 XJT(NTPMAX),YJT(NTPMAX),ZJT(NTPMAX),
     &         VXJT(NTPMAX),VYJT(NTPMAX),VZJT(NTPMAX),
     &         MATR(3,3),MATI(3,3),
     &         GME,AJ(NPLMAX),EJ(NPLMAX),IDISK,ODISK,
     &         GM,KA(10*NTPMAX),A,E,INC,CAPOM,OMEGA,CAPM,AMIN,AMAX,DA,
     &         RSTAT(NTPMAX,NSTATR),MSUN,RMIN,RMAX,RMIN25,RMAX25,
     &         CIMAX,IMAX,DE,EMAX,EMIN,WIMPMFP(NTPMAX),rand,
     &         XT,YT,ZT,VXT,VYT,VZT
        INTEGER ISTAT(NTPMAX,NSTAT),OLOC(NPLMAX,NPLMAX),IFLGCHK
        INTEGER ORBCT(NPLMAX),NC,KK,JJ,OLOCTT(NPLMAX),
     &          OLOCT(NPLMAX,NTPMAX),NTP,NTPT
        INTEGER II,K,IUFLG,IND,DKLOC,NCOR

        REAL*4 SEED

        INTEGER IPOINT,IPAR
        INTEGER IT(10*NTPMAX),MT(10*NTPMAX),IFILE
        INTEGER NTPFIRST,I,J,ISEED,LLM,CNTFILE,NTPFILE
        LOGICAL OK,SSYS,OKC,SAT,CEN,STARTDISK
        LOGICAL*2 LCLOSE
        CHARACTER*(CHLEN) INTPFILE,TPFILE,INPARFILE,PARFILE,OUTFILE
        CHARACTER*(CHLEN) DIRO,DIRS,GNAME,FOPENSTAT,GOFILE,MVSFILE
        CHARACTER*(CHLEN) BINDIR,GOCFILE,GOCMD,INTEG, CLEARFILE
        CHARACTER*(CHLEN) MEXTRFILE,CONTFILE,STARTFILE, CORRFILE
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

        CHARACTER*1 :: REP
        CHARACTER*32 NAME
        CHARACTER*14 INPLFILE
        CHARACTER*65 FPLANETS
        INTEGER NBOD,IALPHA,IP1,IP2,ICFLG,IRFLG
        INTEGER NNN1,NNN2

0009    FORMAT('  rm -rf ./run_',I2.2) ! added by antoine

0008    FORMAT('cd ',(a),'/',(a), '/') ! added by antoine
0001    FORMAT('(a)',f10.1) ! added by antoine
1000    FORMAT(I6,6(1X,F7.3))   !étiquettes de format
2000    FORMAT((a),'_',I2.2,'.in')
3000    format(100(a1,1x))
4000    FORMAT((a),'_',I2.2)
5000    FORMAT('files_',I2.2,'.in') ! modified by antoine
6000    FORMAT('start_',I2.2,'.sh') ! modified by antoine
6001    FORMAT('continue_',I2.2,'.sh') ! modified by antoine
2200    FORMAT('./run_',I2.2,'/',(a),'_',I2.2,'.in') ! modified by antoine
2201    FORMAT(I0) ! modified by antoine
2202    FORMAT('tfin=',f10.1)
22021    FORMAT('t_0=',f0.0)   ! added by antoine
22022    FORMAT('t_f=',f0.0)   ! added by antoine
22023    FORMAT('dt=',f0.0)   ! added by antoine
2203    FORMAT(f10.1)
2204    FORMAT((a),'swift_',(a),' < ./run_',I2.2,'/files.in') ! modified by antoine
2205    FORMAT('order=${1:-"oarsub -l /nodes=1/cpu=1/core=',I1,
     &       ',walltime=48 --project dynapla"}') ! modified by antoine
2206    FORMAT('export OMP_NUM_THREADS=',I1)

c ........entree des parametres

        READ(*,'(a)')BINDIR  ! added by antoine

        MEXTRFILE='mextract_multi.sh'
        CONTFILE='continue.sh'
        STARTFILE='start.sh'
        CLEARFILE='clear.sh'
        CORRFILE='corr_multi.dat'

        WRITE(*,*)'Remove old files' ! added by antoine
        CALL SYSTEM('[ -f ./'//CLEARFILE//' ] && ./'//CLEARFILE) ! added by antoine

        WRITE(*,*)'Generate sub-simulations:' ! added by antoine

        OPEN(51,FILE=MEXTRFILE,STATUS='UNKNOWN')
        OPEN(52,FILE=CONTFILE,STATUS='UNKNOWN')
        OPEN(53,FILE=STARTFILE,STATUS='UNKNOWN') ! added by antoine
        OPEN(54,FILE=CLEARFILE,STATUS='UNKNOWN')  ! added by antoine


        WRITE(53,'(a)')'#! /bin/bash' ! added by antoine
        WRITE(54,'(a)')'#! /bin/bash' ! added by antoine



        WRITE(51,'(a)')'#! /bin/bash' ! modified by antoine
        ! WRITE(51,'(a)')'unlimit' ! commented par antoine
        WRITE(51,'(a)')'export OMP_NUM_THREADS=1' ! modified by antoine
        WRITE(51,'(a)')'export STACKSIZE=1000000' ! modified by antoine
c
        WRITE(*,*)' Name of integrator'
        WRITE(*,*)' (hjs,hjs_parp) '
        READ(*,'(a)')INTEG
        WRITE(*,*)' Number of cores'
        READ(*,*)NCOR

        WRITE(52,2205)NCOR ! added by antoine
        WRITE(53,2205)NCOR ! added by antoine

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

        WRITE(*,*)' Give number of massive bodies in the system '
        READ(*,*)NBOD
        DO I=1,NBOD
          WRITE(*,*)' Give mass of body #',i,' in solar masses '
          READ(*,*)MASS(I)
        END DO

        IF(IUFLG.EQ.1) MASS(1:NBOD) = MASS(1:NBOD)*SMASSYR

        WRITE(*,*) ' Mass of the Sun is : ',MASS(1)
  
        XJ(1) = 0.0
        YJ(1) = 0.0
        ZJ(1) = 0.0
        VXJ(1) = 0.0
        VYJ(1) = 0.0
        VZJ(1) = 0.0

        DO I=2,NBOD
          WRITE(*,*)' Information relative to orbit #',I-1,':'
          WRITE(*,*)' Give situation of bodies relative to that orbit.'
          WRITE(*,*)' 0 = foreign,   1 = satellite,  -1 = center'
          READ(*,*)(OLOC(I,J),J=1,NBOD)

          WRITE(*,*) 'Give now the orbital parameters '
          READ(*,*)AJ(I),EJ(I),INC,CAPOM,OMEGA,CAPM

          CAPM = DR*CAPM
          OMEGA = DR*OMEGA
          CAPOM = DR*CAPOM
          INC = DR*INC
         
          GM = 0.
          DO J = 1,NBOD
            IF (OLOC(I,J).NE.0) GM = GM+MASS(J)
          END DO
          IALPHA = -1
          CALL ORBEL_EL2XV(GM,IALPHA,AJ(I),EJ(I),INC,CAPOM,OMEGA,CAPM,
     &          XJ(I),YJ(I),ZJ(I),VXJ(I),VYJ(I),VZJ(I))
        END DO

C...  Calc the radii
        RPLSQ(1) = 0.0D0   
        IF(IRFLG.EQ.0) THEN !IRFLG sélection type rayon 
          LCLOSE = .FALSE.
          RPLSQ(1:NBOD) = 0.d0 ! car pas de rayon
        ELSE IF(IRFLG.EQ.1) THEN
          LCLOSE = .TRUE.
          DO I=2,NBOD
            RPLSQ(I) = (RPHY(I-1)/AU)**2 
c ! ici on utilise la table de rayons donnée au début, ms que fait-on? **2 ?
          END DO
c        ELSE 
c          LCLOSE = .TRUE.
c           subroutine ds le cas où R est un multiple du rayon de Hill
c          RPLSQ(2:NBOD) = (FACHILL*SQRT(R2HILL(2:NBOD)))**2
        ENDIF

c... Compute eta's and mu's: center and satellite masses for orbits
        eta = 0
        mu = 0
        do j = 2,nbod
          do i = 1,nbod
            if (oloc(j,i).eq.1) mu(j) = mu(j)+mass(i)
            if (oloc(j,i).eq.-1) eta(j) = eta(j)+mass(i)
          end do
        end do

c... Build transform matrix Barycentric --> Jacobi
        mat = 0.0d0
        mtot = sum(mass(1:nbod))
        do i = 1,nbod
          mat(1,i) = mass(i)/mtot
        end do
        do j = 2,nbod
          do i = 1,nbod
            if (oloc(j,i).eq.1) mat(j,i) = mass(i)/mu(j)
            if (oloc(j,i).eq.-1) mat(j,i) = -mass(i)/eta(j)
          end do
        end do    

c...    Build inverse transform matrix Jacobi --> Barycentric
        umat = 0.0d0
        do i = 1,nbod
          umat(i,1) = 1.0d0
        end do
        do j = 2,nbod
          vsat = eta(j)/(mu(j)+eta(j))
          vcen = -mu(j)/(mu(j)+eta(j))
          do i = 1,nbod
            if (oloc(j,i).eq.1) umat(i,j) = vsat
            if (oloc(j,i).eq.-1) umat(i,j) = vcen
          end do
        end do         
c        do i=1,nbod
c          print*,sngl(umat(i,1:nbod))
c        end do

        CALL COORD_G2B(NBOD,UMAT,MASS(1:NBOD),
     &     XJ(1:NBOD),YJ(1:NBOD),ZJ(1:NBOD),VXJ(1:NBOD),VYJ(1:NBOD),
     &     VZJ(1:NBOD),XB(1:NBOD),YB(1:NBOD),ZB(1:NBOD),VXB(1:NBOD),
     &     VYB(1:NBOD),VZB(1:NBOD))

        IF(ICFLG.EQ.1) CALL INVAR(NBOD,MASS,XB,YB,ZB,VXB,VYB,VZB,MATI)

        CALL COORD_B2G(NBOD,MAT,MASS(1:NBOD),
     &     XB(1:NBOD),YB(1:NBOD),ZB(1:NBOD),VXB(1:NBOD),VYB(1:NBOD),
     &     VZB(1:NBOD),XJO(1:NBOD),YJO(1:NBOD),ZJO(1:NBOD),
     &     VXJO(1:NBOD),VYJO(1:NBOD),VZJO(1:NBOD))


	WRITE(*,*) ' '
	WRITE(*,*) 'Enter name of planet data file : '
	READ(*,'(a)')INPLFILE
        IFLGCHK = 0
	CALL IO_DUMP_PL_HJS(INPLFILE,NBOD,OLOC,
     &       MASS(1:NBOD),UMAT,XJO(1:NBOD),
     &       YJO(1:NBOD),ZJO(1:NBOD),VXJO(1:NBOD),VYJO(1:NBOD),
     &       VZJO(1:NBOD),LCLOSE,IFLGCHK,RPLSQ(1:NBOD))

        OKC = .TRUE.
        WRITE(*,*) ' Input ISEED: large and odd '
        READ(*,*) ISEED
        SEED=FLOAT(ISEED)
        CALL RANDOM_SEED()
        print*,seed
        WRITE(*,*) 'Enter name of test particle data file : '
        READ(*,'(a)')INTPFILE         
        IPOINT = INDEX(INTPFILE,'.')
        WRITE(*,*) 'Enter name of parameter file : '
        READ(*,'(a)')INPARFILE         
        call io_init_param_hb(inparfile,t0,tstop,dtmin,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,diro,dirs,
     &         gname,outfile,fopenstat)
        DO I=2,NBOD
          DTJ(I) = DTMIN*2.d0*PI*SQRT(AJ(I)**3/(ETA(I)+MU(1)))
          DTJ(I) = DTJ(I)*(1.d0-EJ(I))**1.5d0/SQRT(1.d0+EJ(I))
          print*,'dtj',dtj(i),aj(i),mass(1)
        END DO

        ! WRITE(51,'(a)')'simname='//TRIM(GNAME)
        WRITE(PARFILE,22021)SNGL(T0) ! added by antoine
        WRITE(51, '(a)')TRIM(PARFILE) ! added by antoine
        WRITE(PARFILE,22022)SNGL(TSTOP) ! modified by antoine
        WRITE(51, '(a)')TRIM(PARFILE)
        dt = 0.1d0*TSTOP
        WRITE(PARFILE,22023)SNGL(dt) ! added by antoine
        WRITE(51, '(a)')TRIM(PARFILE)
        
        WRITE(51, 0008)TRIM(DIRS),TRIM(GNAME)      ! added by antoine
        WRITE(51,'(a)')TRIM(BINDIR)//'mbodies_multi_hjs <<!' ! modified by antoine

        IPAR = INDEX(INPARFILE,'.')
        CNTFILE = 0
        CPUTIME = 0.d0
        NTPT = 0

        DO WHILE(OKC)
          WRITE(*,*) ' Input size of grid for particles '
          READ(*,*) NTP
          OKC = (NTP.NE.0)
          IF (OKC) THEN
            STARTDISK = .TRUE.
            WRITE(*,*) ' Input number of tp of first file '
            READ(*,*) NTPFIRST
            WRITE(*,*)' Give orbct : centers for tp '
            READ(*,*) (ORBCT(I),I=1,NBOD)
            write(*,*)'Centers : ',(orbct(j),j=1,nbod)
            WRITE(*,*)
     &'Disk in ecliptic (0),centers (1),invariable (2),or orbit (<0) ?'
            READ(*,*)DKLOC 
c... Check the validity
            call hierarchktp(nbod,oloc,orbct,ok)
            write(*,*)' '
            if (.not.ok) then
              write(*,*)'stopping...'
              stop
            end if
            write(*,*)'Hierarchical structure ok !'
            write(*,*)' '

c... Computation of the oloct array
            do j=1,nbod
              oloctt(j) = 0
            end do
            do j = 2,nbod
              sat = .true.
              cen = .true.
              do i = 1,nbod
                if (orbct(i).eq.-1) then                  
                  sat = sat.and.(oloc(j,i).eq.1)
                  cen = cen.and.(oloc(j,i).eq.-1)
                end if
              end do
              if (sat) oloctt(j) = 1
              if (cen) oloctt(j) = -1
            end do
            print*,'oloct',(oloctt(j),j=1,nbod)

            etatpt = 0.0d0   
            matpt(:) = 0.0d0
            umatpt(:) = 0.0d0
            do j = 1,nbod
              if (orbct(j).eq.-1) etatpt = etatpt+mass(j)
            end do
            do j = 1,nbod
              if (orbct(j).eq.-1) matpt(j) = -mass(j)/etatpt
            end do
            do k = 1,nbod
              if (orbct(k).eq.-1) then
                do j=1,nbod     
                  umatpt(j) = umatpt(j) - matpt(k)*umat(k,j)
                end do
              end if
            end do

            WRITE(*,*) ' Input limits for eccentricities'
            READ(*,*)RMIN25,RMAX25
            WRITE(*,*) ' Input max inc of particles (in deg) '
            READ(*,*) IMAX
            WRITE(*,*) ' Input min and max values of A'
            READ(*,*) AMIN,AMAX
            DT0 = DTMIN*2.d0*PI*SQRT(AMIN**3/ETATPT)
            
            E=EMIN
            CIMAX = COS(IMAX/DEGRAD)
            GM = ETATPT
            WRITE(*,*) ' Test particles: '
            WRITE(*,*) '   A      E      I    CAPOM  OMEGA    M '
            DO I=1,NTP
              call random_number(rand)
              KA(I) = AMIN + (AMAX - AMIN)*rand
            END DO
            CALL ORDREB(KA,IT,MT,NTP)
            IFILE = 1
            OPEN(12,FILE=CORRFILE,status='UNKNOWN')

c... Here we need to calculate the plane perp. to the local angular momentum 
c...    of the centers of the tp's. The orbital elements of the tp's are given 
c...    relative to this base plane.
            nc = 0
            do j = 1,nbod
              if (orbct(j).eq.-1) then
                nc = nc+1
              end if
            end do
    
c... Of course we calculate the midplane if there is more than one center...
c... and if demanded
            if ((nc.gt.1).and.(dkloc.eq.1)) then 

c... Now we extract a submatrix from the umat matrix 
              umpart = 0.0d0
              do j = 1,nc
                umpart(j,1) = 1.0d0
              end do
              jj = 1
              do j = 2,nbod
                ok = .true.
                do k = 1,nbod
                  if ((matpt(k).eq.0.0d0).and.(oloc(j,k).ne.0))
     &                                                ok=.false.
                end do
                if (ok) then
                  jj = jj+1
                  kk = 0
                  xjtemp(jj) = xjo(j)
                  yjtemp(jj) = yjo(j)
                  zjtemp(jj) = zjo(j)
                  vxjtemp(jj) = vxjo(j)
                  vyjtemp(jj) = vyjo(j)
                  vzjtemp(jj) = vzjo(j)
                  do k = 1,nbod
                    if (matpt(k).ne.0.0d0) then
                      kk = kk+1
                      masstemp(kk) = mass(k)
                      umpart(kk,jj) = umat(k,j)
                    end if
                  end do
                end if
              end do
                   
              call coord_g2b(nc,umpart,masstemp(1:nc),
     &          xjtemp(1:nc),yjtemp(1:nc),zjtemp(1:nc),vxjtemp(1:nc),
     &          vyjtemp(1:nc),vzjtemp(1:nc),xbtemp(1:nc),ybtemp(1:nc),
     &          zbtemp(1:nc),vxbtemp(1:nc),vybtemp(1:nc),vzbtemp(1:nc))
                print*,xbtemp(1:nc),xjtemp(1:nc)
              CALL INVAR(NC,MASSTEMP(1:NC),XBTEMP(1:NC),YBTEMP(1:NC),
     &          ZBTEMP(1:NC),VXBTEMP(1:NC),VYBTEMP(1:NC),VZBTEMP(1:NC),
     &          MATR)
            END IF
c
c.. Disque dans le plan d'une orbite donnée
c 
            IF (DKLOC.LT.0) THEN
               MASSTEMP(1) = 1.d0
               CALL INVAR(1,MASSTEMP(1:1),
     &                     XJ(-DKLOC),YJ(-DKLOC),ZJ(-DKLOC),
     &                     VXJ(-DKLOC),VYJ(-DKLOC),VZJ(-DKLOC),MATR)
            END IF
            IF (DKLOC.EQ.3) THEN
              WRITE(*,*) ' Input inclination in degrees '
              READ(*,*) IDISK
              WRITE(*,*) ' Input Omega (in deg) '
              READ(*,*) ODISK
              IDISK = IDISK*DR
              ODISK = ODISK*DR
    
              MATR(1,1) = COS(ODISK)
              MATR(1,2) = SIN(ODISK)
              MATR(1,3) = 0.0d0

              MATR(2,1) = -COS(IDISK)*SIN(ODISK)
              MATR(2,2) = COS(IDISK)*COS(ODISK)
              MATR(2,3) = SIN(IDISK)

              MATR(3,1) = SIN(IDISK)*SIN(ODISK)
              MATR(3,2) = -SIN(IDISK)*COS(ODISK)
              MATR(3,3) = COS(IDISK)
            END IF

c...... Now draw tp's

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
              J = NTPT+1
              call random_number(rand) 
              E = RMIN25 + (RMAX25-RMIN25)*RAND
c              WRITE(*,1000) J,A,E,INC,CAPOM,OMEGA,CAPM
c              WRITE(12,*) A(IT(I)),E,INC,CAPOM,OMEGA,CAPM

              GME=GM
              CALL ORBEL_EL2XV(GME,IALPHA,KA(IT(I)),E,INC,CAPOM,
     &          OMEGA,CAPM,XJT(J),YJT(J),ZJT(J),VXJT(J),
     &          VYJT(J),VZJT(J))


              IF (((NC.GT.1).AND.(DKLOC.EQ.1))
     &           .OR.(DKLOC.LT.0).OR.(DKLOC.EQ.3)) THEN
                XT = MATR(1,1)*XJT(J)+MATR(2,1)*YJT(J)
     &                                      +MATR(3,1)*ZJT(J) 
                YT = MATR(1,2)*XJT(J)+MATR(2,2)*YJT(J)
     &                                      +MATR(3,2)*ZJT(J) 
                ZT = MATR(1,3)*XJT(J)+MATR(2,3)*YJT(J)
     &                                      +MATR(3,3)*ZJT(J) 
                XJT(J) = XT
                YJT(J) = YT
                ZJT(J) = ZT
                VXT = MATR(1,1)*VXJT(J)+MATR(2,1)*VYJT(J)
     &                                      +MATR(3,1)*VZJT(J) 
                VYT = MATR(1,2)*VXJT(J)+MATR(2,2)*VYJT(J)
     &                                      +MATR(3,2)*VZJT(J) 
                VZT = MATR(1,3)*VXJT(J)+MATR(2,3)*VYJT(J)
     &                                      +MATR(3,3)*VZJT(J) 
                VXJT(J) = VXT
                VYJT(J) = VYT
                VZJT(J) = VZT
              ELSE IF ((DKLOC.EQ.0).AND.(ICFLG.EQ.1)) THEN
                XT = MATI(1,1)*XJT(J)+MATI(1,2)*YJT(J)
     &                                      +MATI(1,3)*ZJT(J) 
                YT = MATI(2,1)*XJT(J)+MATI(2,2)*YJT(J)
     &                                      +MATI(2,3)*ZJT(J) 
                ZT = MATI(3,1)*XJT(J)+MATI(3,2)*YJT(J)
     &                                      +MATI(3,3)*ZJT(J) 
                XJT(J) = XT
                YJT(J) = YT
                ZJT(J) = ZT
                VXT = MATI(1,1)*VXJT(J)+MATI(1,2)*VYJT(J)
     &                                      +MATI(1,3)*VZJT(J) 
                VYT = MATI(2,1)*VXJT(J)+MATI(2,2)*VYJT(J)
     &                                      +MATI(2,3)*VZJT(J) 
                VZT = MATI(3,1)*VXJT(J)+MATI(3,2)*VYJT(J)
     &                                      +MATI(3,3)*VZJT(J) 
                VXJT(J) = VXT
                VYJT(J) = VYT
                VZJT(J) = VZT
              END IF

              oloct(1:nbod,J) = OLOCTT(1:NBOD)
              matp(1:nbod,J) = MATPT(1:NBOD)
              umatp(1:nbod,J) = UMATPT(1:NBOD)


              DO K=1,NSTAT
                ISTAT(J,K) = 0
              END DO
              DO K=1,NSTATR
                RSTAT(J,K) = 0.0D0
              END DO
              NTPT = NTPT+1 
c Estimation du pas de temps a une distance donnee
              DTTRY = DT0*(KA(IT(I))/AMIN)**1.5d0
c              DO K=2,NBOD
c                DTTRY = DTTRY*(1.d0
c     &         -0.9d0*EXP(-0.5d0*(A(IT(I))-ABOD(K))**2/R2HILL(K)))
c              END DO
              DO K=2,NBOD
                DTTRY = MIN(DTTRY,DTJ(K))
              END DO
              DTTRY = DTOUT/NINT(DTOUT/DTTRY)
c
              IF ((NTPT.EQ.1).OR.(DTTRY.LT.DT)) THEN
                DT = DTTRY
              END IF
              CPUTIME = CPUTIME + 1.d0/DT
c
              IF (STARTDISK) THEN
                CPUREF = CPUTIME
                OK = (NTPT.EQ.NTPFIRST)
              ELSE
                OK = ((CPUTIME.GT.CPUREF).OR.(I.EQ.NTP))
              END IF
              IF (OK) THEN
                CNTFILE = CNTFILE+1
                IF (CNTFILE.GT.99) STOP
                IFILE = IFILE+1
                STARTDISK = .FALSE.
                WRITE(TPFILE,2000)INTPFILE(1:IPOINT-1),CNTFILE
                CALL IO_DUMP_TP_HJS(TPFILE,NBOD,NTPT,MATP,
     &              XJT,YJT,ZJT,VXJT,VYJT,VZJT,ISTAT,RSTAT)
                WRITE(PARFILE,2000)INPARFILE(1:IPAR-1),CNTFILE
                WRITE(DIRO,4000)'run',CNTFILE
                CALL io_dump_param_hb(parfile,t0,tstop,dt,
     &               dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &               lclose,trim(dirs)//'/'//trim(gname),
     &               diro,outfile,fopenstat)
                WRITE(MVSFILE,5000)CNTFILE
                OPEN(31,FILE=MVSFILE,STATUS='UNKNOWN')
                WRITE(31,'(a)')'gen_multi_hjs.sh'
                WRITE(31,'(a)')TRIM(PARFILE)
                WRITE(31,'(a)')TRIM(INPLFILE)
                WRITE(31,'(a)')TRIM(TPFILE)
                WRITE(31,'(a)')'1d-8'
                CLOSE(31)
                WRITE(GOFILE,6000)CNTFILE
                WRITE(GOCFILE,6001)CNTFILE
                OPEN(31,FILE=GOFILE,STATUS='UNKNOWN')
                WRITE(31,'(a)')'#! /bin/bash' ! modified by antoine
                ! WRITE(31,'(a)')'unlimit' ! commented by antoine
                WRITE(GOCMD,2206)NCOR
                WRITE(31,'(a)')TRIM(GOCMD)
                WRITE(31,'(a)')'export STACKSIZE=1000000' ! modified by antoine
                WRITE(31, 0008)TRIM(DIRS),TRIM(GNAME)      ! added by antoine
                WRITE(31,'(a)')TRIM(BINDIR)//'swift_'//TRIM(INTEG)//
     &                ' < ./'//TRIM(MVSFILE) ! modified by antoine
                CLOSE(31)       
                OPEN(32,FILE=GOCFILE,STATUS='UNKNOWN')
                WRITE(32,'(a)')'#! /bin/bash' ! modified by antoine
                ! WRITE(32,'(a)')'unlimit'
                WRITE(GOCMD,2206)NCOR
                WRITE(32,'(a)')TRIM(GOCMD)
                WRITE(32,'(a)')'export STACKSIZE=1000000'
                ! WRITE(32,'(a)')' ' ! commented by antoine
                WRITE(32, 0008)TRIM(DIRS),TRIM(GNAME)      ! added by antoine
                WRITE(GOCMD,2204)TRIM(BINDIR),TRIM(INTEG),CNTFILE
                WRITE(32,'(a)')TRIM(GOCMD)
                CLOSE(32)           
                CALL SYSTEM('chmod ogu+x '//TRIM(GOFILE))
                CALL SYSTEM('chmod ogu+x '//TRIM(GOCFILE)) 

                WRITE(53,'(a)')'${order} ./'//TRIM(GOFILE)   ! modied by antoine
                WRITE(52,'(a)')'${order} ./'//TRIM(GOCFILE)   ! modied by antoine

                print*,tpfile,i,ntpt,
     &               sngl(ka(it(i))),sngl(dt),sngl(cputime)
                WRITE(12,*)'Fichier ',TRIM(tpfile),' nb parts ',ntpt,
     &               'derniere part :',i,'a=',sngl(ka(it(i))),
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
          end if
        enddo
        CLOSE(12)
        CLOSE(52)

    !     WRITE(53,'(a)')'[ "$(ls start_* 2>/dev/null | wc -l)" -eq 0 ]'//
    !  &      ' && rm -f ./start.sh' ! added by antoine
        close(53) ! added by antoine
        
        WRITE(51,2201)CNTFILE
        DO I=1,CNTFILE
          WRITE(PARFILE,2200)I,INPARFILE(1:IPAR-1),I
          WRITE(51,'(a)')TRIM(PARFILE)
        END DO
        WRITE(51,'(a)')'${t_0}'
        WRITE(51,'(a)')'${t_f}'
        ! WRITE(PARFILE,2203)SNGL(T0) ! commented by antoine
        ! WRITE(51,'(a)')TRIM(PARFILE)
        ! WRITE(PARFILE,2203)SNGL(0.1d0*TSTOP)
        WRITE(51,'(a)')'${dt}' ! modified by antoine
        ! WRITE(51,'(a)')TRIM(PARFILE) ! commented by antoine
        WRITE(51,'(a)')'./run_01/'//TRIM(INPLFILE)
        WRITE(51,'(a)')'mextract'
        DO I=1,CNTFILE
          WRITE(TPFILE,2200)I,INTPFILE(1:IPOINT-1),I
          WRITE(51,'(a)')TRIM(TPFILE)
        END DO
        WRITE(51,'(a)')'1'
        WRITE(51,'(a)')'0'
        WRITE(51,'(a)')'!'
        ! WRITE(51,'(a)')'exit' ! commented by antoine
        ! WRITE(51,'(a)')'!'
        CLOSE(51)

        WRITE(54,'(a)')'rm -f ./'//TRIM(INPLFILE)
        WRITE(54,'(a)')'rm -f ./'//TRIM(CONTFILE)
        WRITE(54,'(a)')'rm -f ./'//TRIM(STARTFILE)
        WRITE(54,'(a)')'rm -f ./fort.7'

        WRITE(54,'(a)')'if [ "$1" == "all" ]; then'
        DO I=1,CNTFILE
          WRITE(54,0009)I
        END DO
        WRITE(54,'(a)')'  rm -f ./'//TRIM(MEXTRFILE)
        WRITE(54,'(a)')'  rm -f ./'//TRIM(CORRFILE)
        WRITE(54,'(a)')'  rm -f ./'//TRIM(CLEARFILE)
        WRITE(54,'(a)')'fi'
        CLOSE(54)  ! added by antoine

        CALL SYSTEM('chmod ogu+x '//TRIM(MEXTRFILE))
        CALL SYSTEM('chmod ogu+x '//TRIM(CONTFILE))
        CALL SYSTEM('chmod ogu+x '//TRIM(STARTFILE))  ! added by antoine 
        CALL SYSTEM('chmod ogu+x '//TRIM(CLEARFILE))  ! added by antoine    
        
        END PROGRAM GEN_MULTI_HJS
        ! END SUBROUTINE GEN_MULTI_HJS
C
C--------------------------------------------------------------------
C               For invariable plane
C--------------------------------------------------------------------
C
        SUBROUTINE INVAR(NBOD,MASS,X,Y,Z,VX,VY,VZ,A)

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
c       dparfile      ==>  Name of file to write to (character*80)
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
c       outfile       ==>  Name of binary output file (character*80)

c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  5/10/94 HFL

	subroutine io_dump_param_hb(dparfile,t,tstop,dt,dtout,dtdump,
     &           iflgchk,rmin,rmax,rmaxu,qmin,lclose,dirs,
     &           gname,outfile,fopenstat)	


        include '../../sub/swift.inc'
        include '../../sub/io/io.inc'

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
           write(7,'(a)') trim(dirs)
           write(7,'(a)') trim(gname)
           write(7,'(a)') trim(outfile)
        endif


        write(7,'(a)') fopenstat

	close(unit = 7)

	return
	end     ! io_dump_param_hb
c____________________________________________________________________________
c
c


c-------------------------------------------------------------------------

C
C-----------------------------------------------------------------------
C        Subroutine dedicated to verify the hierarchical structure
C        of the tp's we add. Returns logical ok=true if the hierarchical
C        structure is valid
C-----------------------------------------------------------------------
C
        subroutine hierarchktp(nbod,oloc,orbct,ok)

        include '../../sub/swift.inc'


        integer i,k,l,nnn,nbod,oloc(NPLMAX,NPLMAX),orbct(NPLMAX)
        logical test,ok

        write(*,*)' '
        OK = .true.
c... First check that the orbct array is made only of 0's and -1's
        do i=1,nbod
          ok = ok.and.((orbct(i).eq.0).or.(orbct(i).eq.-1))
        end do
        if (.not.ok) then
          write(*,*)'The orbit given is not a valid tp orbit'
          return
        end if
c... Adding a tp means adding an orbit. We check the validity of the
c... the new orbit in the hierarchical structure of the system, i.e.
c... with all orbits i between massive bodies
        DO I = 2,nbod
          test = .true.
c... We first test whether tp orbit and orbit i are foreign
          do k = 1,nbod
            test = test.and.(oloc(i,k)*orbct(k).eq.0)
          end do
          if (test) then
            write(*,*)'tp orbit and orbit',i-1,' are foreign'
          else
c... From here on they are not foreign. Now check whether tp orbit
c... is inner to orbit i, i.e. all centers of tp orbit
c... must fall in the centers or the satellites of orbit i.
            l = 1
            do while(orbct(l).eq.0)
              l = l+1
            end do
            nnn = oloc(i,l)
            if (nnn.ne.0) then
              test = .true.
              do k = l+1,nbod
                if (orbct(k).eq.-1) test =
     &                               test.and.(oloc(i,k).eq.nnn)
              end do
              if (test) then
                write(*,*)'tp orbit is inner to orbit',i-1
              end if
            end if
            if ((nnn.eq.0).or.(.not.test)) then
c... Now tp orbit is not inner to orbit i. We check whether orbit i
c... is inner to tp orbit
              l = 1
              do while(oloc(i,l).eq.0)
                l = l+1
              end do
              nnn = orbct(l)
              if (nnn.ne.0) then
                test = .true.
                do k = l+1,nbod
                  if (oloc(i,k).ne.0) test = 
     &                               test.and.(orbct(k).eq.nnn)
                end do
                if (test) then
                  write(*,*)'Orbit',i-1,' is inner to tp orbit'
                else
c... If all fails then we are sure that orbits i and j do not
c... fulfill the hierarchy rules
                  write(*,*)'Structure problem relative to orbit',i-1
                  ok = .false.
                end if
              else
                write(*,*)'Structure problem relative to orbit',i-1
                ok=.false. 
              end if
            end if
          end if
        end do
        
        end



