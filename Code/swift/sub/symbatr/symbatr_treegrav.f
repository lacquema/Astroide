C =============================================================
C TREEGRAV: routines to perform hierarchical force calculation.
C =============================================================



C --------------------------------------------------------
C SYMBATR_TREE_ACCEL: routine to compute the acceleration and potential
C for all bodies in the system.
C --------------------------------------------------------

      SUBROUTINE SYMBATR_TREE_ACCEL
 
	INCLUDE './inc/treedefs.f'
        DOUBLE PRECISION CPUT1, CPUT2
	INTEGER P

C       ----------------------------------------
C       Initialize timing for force calculation.
C       ----------------------------------------
C        CALL SECOND_(CPUT1)

C       --------------------------------------------
C       Initialize the force evaluation diagnostics.
C       --------------------------------------------
        NBTOT = 0
	NCTOT = 0
        NTMAX = 0
C	---------------------------------------------------
C	Loop over bodies, computing force for each in turn.
C	---------------------------------------------------
	DO 10 P = 1, NBODY
	  CALL symbatr_tree_TRWALK(P)

 10	CONTINUE
C       -----------------------------------------------
C       Compute average number of force terms per body.
C       -----------------------------------------------
        NTAVG = (NBTOT + NCTOT) / NBODY 
C       -------------------------------------
C       Compute timing for force calculation.
C       -------------------------------------
C        CALL SECOND_(CPUT2)
        CPUACC = CPUT2 - CPUT1
      END





C --------------------------------------------------------------
C TRWALK: recursive routine to walk the tree computing forces on
C particle P.
C --------------------------------------------------------------
 
      SUBROUTINE symbatr_tree_TRWALK(P)
      INTEGER P

	INCLUDE './inc/treedefs.f'
	INTEGER MXSPTR
	PARAMETER(MXSPTR = 256)
	REAL*8 EPS2, PHI0, ACC0(NDIM), POS0(NDIM)
        REAL*8 DX, DY, DZ, DR2, DR2INV, DRINV, PHIM, DR5INV, PHIQ
	INTEGER Q, NBTERM, NCTERM, SPTR, STACK(MXSPTR), K
	LOGICAL SKPSLF

	EPS2 = EPS*EPS
C	----------------------------------------------------------
C	Zero potential and acceleration for subsequent summations.
C	----------------------------------------------------------
        PHI0 = 0.0
        ACC0(1) = 0.0
        ACC0(2) = 0.0
        ACC0(3) = 0.0
C       -----------------------------------------
C       Copy position of P for quicker reference.
C       -----------------------------------------
        POS0(1) = POS(P,1)
        POS0(2) = POS(P,2)
        POS0(3) = POS(P,3)
C	----------------------------------------------
C	Initialize counters for number of force terms.
C	----------------------------------------------
        NBTERM = 0
        NCTERM = 0
	SKPSLF = .FALSE.
C	----------------------------------
C	Push the root cell onto the stack.
C	----------------------------------
        SPTR = 1
        STACK(SPTR) = ROOT
C	-------------------------------------
C	Loop while nodes on stack to process.
C	-------------------------------------
  10    IF (SPTR .GT. 0) THEN
C	  --------------------------
C	  Pop node off top of stack.
C	  --------------------------
          Q = STACK(SPTR)
          SPTR = SPTR - 1
C         ---------------------------------------------
C         Compute distance to center-of-mass of node Q.
C         ---------------------------------------------
          DX = POS0(1) - POS(Q,1)
          DY = POS0(2) - POS(Q,2)
          DZ = POS0(3) - POS(Q,3)
          DR2 = DX*DX + DY*DY + DZ*DZ
C	  -------------------------------
C	  Classify Q as a body or a cell.
C	  -------------------------------
          IF (Q .LT. INCELL) THEN
C	    -----------------------------------
C	    A body: check for self-interaction.
C	    -----------------------------------
	    IF (Q .NE. P) THEN
C             ------------------------------
C             Compute body-body interaction.
C             ------------------------------
              DR2INV = 1.0 / (DR2 + EPS2)
              PHIM = MASS_T(Q) * SQRT(DR2INV)
              PHI0 = PHI0 - PHIM
              PHIM = PHIM * DR2INV
              ACC0(1) = ACC0(1) - PHIM * DX
              ACC0(2) = ACC0(2) - PHIM * DY
              ACC0(3) = ACC0(3) - PHIM * DZ
	      NBTERM = NBTERM + 1
	    ELSE
C	      -------------------------------------------
C	      Remember that self-interaction was skipped.
C	      -------------------------------------------
	      SKPSLF = .TRUE.
            ENDIF
          ELSE
C           --------------------------------------------
C           A cell: test if interaction can be accepted.
C           --------------------------------------------
	    IF (DR2 .GE. RCRIT2(Q)) THEN
C             ----------------------------------------
C             Accepted: compute body-cell interaction.
C             ----------------------------------------
              DR2INV = 1.0 / (DR2 + EPS2)
              DRINV = SQRT(DR2INV)
              PHIM = MASS_T(Q) * DRINV
              PHI0 = PHI0 - PHIM
              PHIM = PHIM * DR2INV
              ACC0(1) = ACC0(1) - PHIM * DX
              ACC0(2) = ACC0(2) - PHIM * DY
              ACC0(3) = ACC0(3) - PHIM * DZ
C	      ------------------------------------
C	      Optionally include quadrupole terms.
C	      ------------------------------------
	      IF (USQUAD) THEN
                DR5INV = DR2INV * DR2INV * DRINV
                PHIQ = DR5INV *
     &              (0.5 * ((DX*DX - DZ*DZ) * QUAD(Q,1) +
     &                      (DY*DY - DZ*DZ) * QUAD(Q,4)) +
     &               DX*DY * QUAD(Q,2) + DX*DZ * QUAD(Q,3) +
     &               DY*DZ * QUAD(Q,5))
                PHI0 = PHI0 - PHIQ
                PHIQ = 5.0 * PHIQ * DR2INV
                ACC0(1) = ACC0(1) - PHIQ*DX + DR5INV *
     &              (DX*QUAD(Q,1) + DY*QUAD(Q,2) + DZ*QUAD(Q,3))
                ACC0(2) = ACC0(2) - PHIQ*DY + DR5INV *
     &              (DX*QUAD(Q,2) + DY*QUAD(Q,4) + DZ*QUAD(Q,5))
                ACC0(3) = ACC0(3) - PHIQ*DZ + DR5INV *
     &              (DX*QUAD(Q,3) + DY*QUAD(Q,5) -
     &               DZ*(QUAD(Q,1) + QUAD(Q,4)))
	      ENDIF
	      NCTERM = NCTERM + 1
            ELSE
C	      -----------------------------------
C	      Rejected: examine children of cell.
C	      -----------------------------------
              DO 20 K = 1, NSUBC
C		--------------------------------------
C		Push existing children onto the stack.
C		--------------------------------------
                IF (SUBP(Q,K) .NE. NULL) THEN
		  IF (SPTR .GE. MXSPTR)
     &		    CALL symbatr_io_TERROR(' TRWALK: STACK OVERFLOW')
                  SPTR = SPTR + 1
                  STACK(SPTR) = SUBP(Q,K)
                ENDIF
 20           CONTINUE
	    ENDIF
          ENDIF
          GOTO 10
        ENDIF
C       ---------------------------------------------------
C       Check that self-interaction was explicitly skipped.
C       ---------------------------------------------------
        IF (.NOT. SKPSLF)
     &    CALL symbatr_io_TERROR(' TRWALK: MISSED SELF-INTERACTION')
C       ----------------------------------------------
C       Copy total potential and acceleration to body.
C       ----------------------------------------------
        PHI(P) = PHI0
        ACC(P,1) = ACC0(1)
        ACC(P,2) = ACC0(2)
        ACC(P,3) = ACC0(3)
C       -------------------------------------
C       Update force calculation diagnostics.
C       -------------------------------------
        NBTOT = NBTOT + NBTERM
        NCTOT = NCTOT + NCTERM
        NTMAX = MAX(NTMAX, NCTERM + NBTERM)
      END












