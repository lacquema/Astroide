C ================================================================
C TREELOAD: routines to handle tree construction and verification.
C ================================================================
 
C ----------------------------------------------------------------
C MKTREE: initialize the tree structure for the force calculation.
C ----------------------------------------------------------------
 
      SUBROUTINE SYMBATR_TREE_MKTREE
 
        INCLUDE './inc/treedefs.f'
        DOUBLE PRECISION CPUT1, CPUT2
 
C       ----------------------------------
C       Start timing of tree construction.
C       ----------------------------------
C        CALL SECOND_(CPUT1)
C       --------------------------------------
C       Expand root volume to hold all bodies.
C       --------------------------------------
	CALL symbatr_tree_EXPBOX
C	--------------------------
C       Load bodies into the tree.
C	--------------------------
	CALL symbatr_tree_LDTREE
C       ------------------------------------------------
C       Compute masses, center of mass coordinates, etc.
C       ------------------------------------------------
        CALL symbatr_tree_HACKCM
C       -------------------------------------
C       Compute timing for tree construction.
C       -------------------------------------
C        CALL SECOND_(CPUT2)
        CPUMKT = CPUT2 - CPUT1
      END
 
C -----------------------------------------------------
C EXPBOX: enlarge system cube to include all particles.
C -----------------------------------------------------
 
      SUBROUTINE symbatr_tree_EXPBOX
 
        INCLUDE './inc/treedefs.f'
	REAL*8 XYZMAX
        INTEGER J, I
 
C       ------------------------------
C       Find maximum coordinate value.
C       ------------------------------
	XYZMAX = 0
        DO 10 J = 1, NDIM
	  DO 20 I = 1, NBODY
	    XYZMAX = MAX(XYZMAX, ABS(POS(I,J)))
  20	  CONTINUE
  10	CONTINUE
C	-----------------------------------------
C	Expand box by factors of 2 until it fits.
C	-----------------------------------------
  30	IF (XYZMAX .GE. RSIZE/2.0) THEN
	  RSIZE = 2.0 * RSIZE
	  GOTO 30
	ENDIF
      END



C ---------------------------------------------------------------------
C LDTREE: construct tree body-by-body.  This phase initializes the SUBP
C array and loads the geometric midpoint into the MID of each cell.
C ---------------------------------------------------------------------

      SUBROUTINE symbatr_tree_LDTREE

	INCLUDE './inc/treedefs.f'
        INTEGER symbatr_tree_mkcell, K, P

C       ---------------------------------------
C       Deallocate current tree, begin new one.
C       ---------------------------------------
        NCELL = 0
        ROOT = symbatr_tree_mkcell()
C	------------------------------------------
C	Initialize midpoint and size of root cell.
C	------------------------------------------
	DO 10 K = 1, NDIM
	  MID(ROOT,K) = 0.0
  10	CONTINUE
	CLSIZE(ROOT) = RSIZE
C       ---------------------------------------------
C       Load bodies into the new tree, one at a time.
C       ---------------------------------------------
        DO 20 P = 1, NBODY
          CALL symbatr_tree_LDBODY(P)
  20    CONTINUE
      END
 
C ----------------------------------------
C LDBODY: load body P into tree structure.
C ----------------------------------------
 
      SUBROUTINE symbatr_tree_LDBODY(P)
      INTEGER P
 
        INCLUDE './inc/treedefs.f'
        INTEGER Q, QIND, symbatr_tree_SBINDX, symbatr_tree_mkcell, 
     &       C, K, P0

C       ---------------------------------------------
C       Start Q,QIND pair in correct subcell of root.
C       ---------------------------------------------
        Q = ROOT
        QIND = symbatr_tree_SBINDX(P, Q)
C       -----------------------------------------------------
C       Loop descending tree until an empty subcell is found.
C       -----------------------------------------------------
  20    IF (SUBP(Q, QIND) .NE. NULL) THEN
C         --------------------------------------
C         On reaching another body, extend tree.
C         --------------------------------------
          IF (SUBP(Q, QIND) .LT. INCELL) THEN
C           -------------------------------------------
C           Allocate an empty cell to hold both bodies.
C           -------------------------------------------
            C = symbatr_tree_mkcell()
C           ------------------------------------------------------
C           Locate midpoint of new cell wrt. parent, and set size.
C           ------------------------------------------------------
            DO 30 K = 1, NDIM
              IF (POS(P,K) .GE. MID(Q,K)) THEN
                MID(C,K) = MID(Q,K) + CLSIZE(Q)/4.0
              ELSE
                MID(C,K) = MID(Q,K) - CLSIZE(Q)/4.0
              ENDIF
  30        CONTINUE
	    CLSIZE(C) = CLSIZE(Q) / 2.0
C           ------------------------------------------------------
C           Store old body in appropriate subcell within new cell.
C           ------------------------------------------------------
            P0 = SUBP(Q, QIND)
            SUBP(C, symbatr_tree_SBINDX(P0, C)) = P0
C           ---------------------------------------------
C           Link new cell into tree in place of old body.
C           ---------------------------------------------
            SUBP(Q, QIND) = C
          ENDIF
C         --------------------------------------------------------
C         At this point, the node indexed by Q,QIND is known to be
C	  a cell, so advance to the next level of tree, and loop.
C         --------------------------------------------------------
          Q = SUBP(Q, QIND)
          QIND = symbatr_tree_SBINDX(P, Q)
          GOTO 20
        ENDIF
C       ---------------------------------------------
C       Found place in tree for P, so store it there.
C       ---------------------------------------------
        SUBP(Q, QIND) = P
      END
 
C -------------------------------------------------------
C symbatr_tree_SBINDX: compute subcell index for node P within cell Q.
C -------------------------------------------------------
 
      INTEGER FUNCTION symbatr_tree_SBINDX(P, Q)
      INTEGER P, Q
 
        INCLUDE './inc/treedefs.f'
        INTEGER K
 
C       ---------------------------------------------------
C       Initialize subindex to point to lower left subcell.
C       ---------------------------------------------------
        symbatr_tree_SBINDX = 1
C       ---------------------------------
C       Loop over all spatial dimensions.
C       ---------------------------------
        DO 10 K = 1, NDIM
          IF (POS(P,K) .GE. MID(Q,K))
     &      symbatr_tree_SBINDX = symbatr_tree_SBINDX + 2 ** (3 - K)
 10     CONTINUE
      END
 
C ---------------------------------------------------------
C MKCELL: function to allocate a cell, returning its index.
C ---------------------------------------------------------
 
      INTEGER FUNCTION symbatr_tree_mkcell()
 
        INCLUDE './inc/treedefs.f'
        INTEGER I
 
C       ----------------------------------------------------------
C       Terminate simulation if no remaining space for a new cell.
C       ----------------------------------------------------------
        IF (NCELL .GE. MXCELL)
     &    CALL symbatr_io_TERROR(' MKCELL: NO MORE MEMORY')
C       ----------------------------------------------------
C       Increment cell counter, initialize new cell pointer.
C       ----------------------------------------------------
        NCELL = NCELL + 1
        symbatr_tree_mkcell = NCELL + MXBODY
C       --------------------------------------
C       Zero pointers to subcells of new cell.
C       --------------------------------------
        DO 10 I = 1, NSUBC
          SUBP(symbatr_tree_mkcell,I) = NULL
 10     CONTINUE
      END
 
C ----------------------------------------------------------------
C HACKCM: compute cell masses, c.m. positions, check tree structure,
C assign critical radii, and compute quadrupole moments.  Whew!
C ----------------------------------------------------------------
 
      SUBROUTINE symbatr_tree_HACKCM
 
        INCLUDE './inc/treedefs.f'
        INTEGER IND(MXCELL), P, Q, I, J, K, L, M, N
	REAL*8 POS0(NDIM), DIST2
 
C       ---------------------------------------
C       List cells in order of decreasing size.
C       ---------------------------------------
        CALL symbatr_tree_BFLIST(IND)
C       --------------------------------------------
C       Loop processing cells from smallest to root.
C       --------------------------------------------
        DO 10 I = NCELL, 1, -1
          P = IND(I)
C         --------------------------------------------------------------
C         Zero accumulators for this cell.  A temporary variable is used
C	  for the c. of m. so as to preserve the stored midpoints.
C         --------------------------------------------------------------
          MASS_T(P) = 0.0
          DO 20 K = 1, NDIM
            POS0(K) = 0.0
  20      CONTINUE
C         -------------------------------------------------------------
C         Compute cell properties as sum of properties of its subcells.
C         -------------------------------------------------------------
          DO 30 J = 1, NSUBC
            Q = SUBP(P,J)
C           ------------------------------
C           Only access cells which exist.
C           ------------------------------
            IF (Q .NE. NULL) THEN
C             -------------------------------------------------------
C             Sum properties of subcells to obtain values for cell P.
C             -------------------------------------------------------
              MASS_T(P) = MASS_T(P) + MASS_T(Q)
              DO 40 K = 1, NDIM
                POS0(K) = POS0(K) + MASS_T(Q) * POS(Q,K)
  40          CONTINUE
            ENDIF
  30      CONTINUE
C         --------------------------------------------------------
C         Normalize center of mass coordinates by total cell mass.
C         --------------------------------------------------------
          DO 50 K = 1, NDIM
            POS0(K) = POS0(K) / MASS_T(P)
  50      CONTINUE 
C         -----------------------------------------------------------------
C         Check tree, compute cm-to-mid distance, and assign cell position.
C         -----------------------------------------------------------------
          DIST2 = 0.0
          DO 60 K = 1, NDIM
            IF (POS0(K) .LT. MID(P,K) - CLSIZE(P)/2.0 .OR.
     &		  POS0(K) .GE. MID(P,K) + CLSIZE(P)/2.0) THEN
              WRITE(ULOG, '(/,1X,''TREE ERROR'',2I6,3E14.6)')
     &              P, K, POS0(K), MID(P,K), CLSIZE(P)
              CALL symbatr_io_TERROR(' HACKCM: TREE STRUCTURE ERROR')
            ENDIF
            DIST2 = DIST2 + (POS0(K) - MID(P,K))**2
C	    --------------------------------------------------------
C	    Copy cm position to cell.  This overwrites the midpoint.
C	    --------------------------------------------------------
	    POS(P,K) = POS0(K)
  60      CONTINUE

C	  ------------------------------------------------------------
C	  Assign critical radius for cell, adding offset from midpoint
C	  for more accurate forces.  This overwrites the cell size.
C	  ------------------------------------------------------------

          RCRIT2(P) = (CLSIZE(P) / THETA + SQRT(DIST2))**2

  10    CONTINUE
C	----------------------------------------
C	Compute quadrupole moments, if required.
C	----------------------------------------
	IF (USQUAD) THEN
C         -------------------------------
C         Loop processing cells as above.
C         -------------------------------
          DO 70 I = NCELL, 1, -1
            P = IND(I)
C           --------------------------------------------
C           Zero accumulator for quad moments of cell P.
C           --------------------------------------------
            DO 80 K = 1, NQUAD
              QUAD(P,K) = 0.0
  80        CONTINUE
C           --------------------------------
C           Loop over descendents of cell P.
C           --------------------------------
	    DO 90 J = 1, NSUBC
	      Q = SUBP(P,J)
	      IF (Q .NE. NULL) THEN
C               --------------------------------------------------------
C               Sum properties of subcell Q to obtain values for cell P.
C               --------------------------------------------------------
                DO 100 M = 1, MIN(2,NDIM)
                  DO 110 N = M, NDIM
                    L = (M-1) * (NDIM-1) + N
                    QUAD(P,L) = QUAD(P,L) + 3.0 * MASS_T(Q) * 
     &                  (POS(Q,M) - POS(P,M)) * (POS(Q,N) - POS(P,N))
                    IF (M .EQ. N) THEN
                      DO 120 K = 1, NDIM
                        QUAD(P,L) = QUAD(P,L) - MASS_T(Q) * 
     &                      (POS(Q,K) - POS(P,K))**2
 120                  CONTINUE
                    ENDIF
C                   -------------------------------------------
C                   If Q itself is a cell, add its moments too.
C                   -------------------------------------------
                    IF (Q .GE. INCELL)
     &                QUAD(P,L) = QUAD(P,L) + QUAD(Q,L)
 110              CONTINUE
 100            CONTINUE
	      ENDIF
  90        CONTINUE
  70     CONTINUE
	ENDIF
      END

C -----------------------------------------------------------------
C BFLIST: list cells in breadth-first order, from largest (root) to
C smallest.  Thanks to Jun Makino for this elegant routine.
C -----------------------------------------------------------------
 
      SUBROUTINE symbatr_tree_BFLIST(IND)
      INTEGER IND(*)
 
        INCLUDE './inc/treedefs.f'
	INTEGER FACELL, LACELL, NACELL, K, I
 
C	-----------------------------------------
C	Start scan with root as only active cell.
C	-----------------------------------------
	IND(1) = ROOT
	FACELL = 1
	LACELL = 1
C	-----------------------------------
C	Loop while active cells to process.
C	-----------------------------------
  10	IF (FACELL .LE. LACELL) THEN
C	  ----------------------------------------------
C	  Start counting active cells in next iteration.
C	  ---------------------------------------------- 
	  NACELL = LACELL
C	  ---------------------------------------
C	  Loop over subcells of each active cell.
C	  ---------------------------------------
	  DO 20 K = 1, NSUBC
	    DO 30 I = FACELL, LACELL
C	      -------------------------------------------
C	      Add all cells on next level to active list.
C	      -------------------------------------------
	      IF (SUBP(IND(I),K) .GE. INCELL) THEN
		NACELL = NACELL + 1
		IND(NACELL) = SUBP(IND(I),K)
	      ENDIF
  30	    CONTINUE
  20	  CONTINUE
C	  ------------------------------------------------------
C	  Advance first and last active cell indicies, and loop.
C	  ------------------------------------------------------
	  FACELL = LACELL + 1
	  LACELL = NACELL
	  GOTO 10
	ENDIF
C	--------------------------------------------------
C	Above loop should list all cells; check the count.
C	--------------------------------------------------
	IF (NACELL .NE. NCELL)
     &    CALL symbatr_io_TERROR('  BFLIST: INCONSISTENT CELL COUNT')
      END
















