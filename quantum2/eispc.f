C     ***************************************************************
      SUBROUTINE HTRID3(NM,N,A,D,E,E2,TAU,GA,GIA)
C     ***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NM,NM)
C
      INTEGER I,J,K,L,N,II,NM,JM1,JP1
      DIMENSION D(N),E(N),E2(N),TAU(2,N),GA(NM),GIA(NM)
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX, STORED AS
C     A SINGLE SQUARE ARRAY, TO A REAL SYMMETRIC TRIDIAGONAL MATRIX
C     USING UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE COMPLEX HERMITIAN INPUT
C          MATRIX.  THE REAL PARTS OF THE MATRIX ELEMENTS ARE STORED
C          IN THE FULL LOWER TRIANGLE OF A, AND THE IMAGINARY PARTS
C          ARE STORED IN THE TRANSPOSED POSITIONS OF THE STRICT UPPER
C          TRIANGLE OF A.  NO STORAGE IS REQUIRED FOR THE ZERO
C          IMAGINARY PARTS OF THE DIAGONAL ELEMENTS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED,
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE USE OF THE SUBROUTINES
C     CABS AND CMPLX IN COMPUTING COMPLEX ABSOLUTE VALUES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0000D0
      TAU(2,N) = 0.0000D0
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0000D0
         SCALE = 0.0000D0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE +  DABS(A(I,K)) +  DABS(A(K,I))
C
         IF (SCALE .NE. 0.0000D0) GO TO 140
         TAU(1,L) = 1.0000D0
         TAU(2,L) = 0.0000D0
  130    E(I) = 0.0000D0
         E2(I) = 0.0000D0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            A(K,I) = A(K,I) / SCALE
            H = H + A(I,K) * A(I,K) + A(K,I) * A(K,I)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G =  DSQRT(H)
         E(I) = SCALE * G
c
c  double precision functions not recognised by cray:
c
         F =  CDABS( DCMPLX(A(I,L),A(L,I)))
c         F =  ABS( CMPLX(A(I,L),A(L,I)))
C     ********** FORM NEXT DIAGONAL ELEMENT OF MATRIX T **********
         IF (F .EQ. 0.0000D0) GO TO 160
         TAU(1,L) = (A(L,I) * TAU(2,I) - A(I,L) * TAU(1,I)) / F
         SI = (A(I,L) * TAU(2,I) + A(L,I) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0000D0 + G / F
         A(I,L) = G * A(I,L)
         A(L,I) = G * A(L,I)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         A(I,L) = G
  170    F = 0.0000D0
C
         DO 172 J = 1, L
            GA(J) = 0.0000D0
 172        GIA(J) = 0.0000D0
C     ********** FORM ELEMENT OF A*U **********
            DO 220 K = 1, L
            IF(K .EQ. L) GO TO 190
            JRC3=K+1
            DO 180 J=JRC3,L
               GA(J) = GA(J) + A(J,K) * A(I,K) + A(K,J) * A(K,I)
               GIA(J) = GIA(J) - A(J,K) * A(K,I) + A(K,J) * A(I,K)
  180       CONTINUE
C
  190       GA(K) = GA(K) + A(K,K) * A(I,K)
            GIA(K) = GIA(K) - A(K,K) * A(K,I)
C
            IF(K .EQ. 1) GO TO 220
            JRC4=K-1
            DO 200 J=1,JRC4
               GA(J) = GA(J) + A(K,J) * A(I,K) - A(J,K) * A(K,I)
               GIA(J) = GIA(J) - A(K,J) * A(K,I) - A(J,K) * A(I,K)
  200       CONTINUE
  220    CONTINUE
C     ********** FORM ELEMENT OF P **********
            DO 240 J=1,L
            E(J) = GA(J) / H
            TAU(2,J) = GIA(J) / H
            F = F + E(J) * A(I,J) - TAU(2,J) * A(J,I)
  240    CONTINUE
C
         HH = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -A(J,I)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
            A(J,J) = A(J,J) - 2.0000D0 * (F * G + FI * GI)
            IF (J .EQ. 1) GO TO 260
            JM1 = J - 1
C
CDIR$ IVDEP
            DO 250 K = 1, JM1
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
     &                         + FI * TAU(2,K) + GI * A(K,I)
               A(K,J) = A(K,J) - F * TAU(2,K) - G * A(K,I)
     &                         - FI * E(K) - GI * A(I,K)
  250       CONTINUE
C
  260    CONTINUE
C
  270    DO 280 K = 1, L
            A(I,K) = SCALE * A(I,K)
            A(K,I) = SCALE * A(K,I)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    D(I) = A(I,I)
         A(I,I) = SCALE *  DSQRT(H)
  300 CONTINUE
C
      RETURN
      END
C     ***************************************************************
      SUBROUTINE TQL2(NM,N,KEY,Z,D,E,HA,HIA,IERR)
C     ***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(NM,NM)
C
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      REAL*8    D(N),E(N),HA(NM),HIA(NM)
      REAL*8    B,C,F,G,H,P,R,S,MACHEP
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
c for single precision (real*4)
c      MACHEP = 2.000D0**(-26)
c
      MACHEP = 2.000D0**(-48)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0000D0
      B = 0.0000D0
      E(N) = 0.0000D0
C
      DO 240 L = 1, N
         J = 0
         H = MACHEP * ( DABS(D(L)) +  DABS(E(L)))
         IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 110 M = L, N
            IF ( DABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0000D0 * E(L))
         R =  DSQRT(P*P+1.0000D0)
         D(L) = E(L) / (P +  SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C = 1.0000D0
         S = 0.0000D0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF ( DABS(P) .LT.  DABS(E(I))) GO TO 150
            C = E(I) / P
            R =  DSQRT(C*C+1.0000D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0000D0 / R
            GO TO 160
  150       C = P / E(I)
            R =  DSQRT(C*C+1.0000D0)
            E(I+1) = S * E(I) * R
            S = 1.0000D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
            DO 170 K=1,N
               HA(K)=Z(K,I)
               HIA(K)=Z(K,I+1)
  170       CONTINUE
            DO 180 K = 1, N
               Z(K,I+1) = S * HA(K) + C * HIA(K)
               Z(K,I) = C * HA(K) - S * HIA(K)
  180       CONTINUE
C
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF ( DABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         IF (KEY.GT.0) GO TO 262
         DO 260 J = II, N
            IF ( D(J).LE.P ) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         GO TO 266
C
  262    CONTINUE
         DO 264 J = II, N
            IF ( D(J).GE.P ) GO TO 264
            K = J
            P = D(J)
  264    CONTINUE
C
  266    IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
      END
C     ***************************************************************
      SUBROUTINE HTRIB3(NM,N,TAU,M,A,ZR,ZI,SA,SIA)
C     ***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NM,NM),ZR(NM,NM),ZI(NM,NM)
      INTEGER I,J,K,L,M,N,NM
      REAL*8    TAU(2,N),SA(NM),SIA(NM)
      REAL*8    H,S,SI
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK3, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRID3.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        A CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
C          USED IN THE REDUCTION BY  HTRID3,
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
C     ********** TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. **********
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     ********** RECOVER AND APPLY THE HOUSEHOLDER MATRICES **********
      DO 140 I = 2, N
         L = I - 1
         H = A(I,I)
         IF (H .EQ. 0.0000D0) GO TO 140
C
         DO 100 J = 1, M
            SA(J) = 0.0000D0
  100       SIA(J) = 0.0000D0
C
            DO 110 K = 1, L
            DO 110 J=1,M
               SA(J) = SA(J) + A(I,K) * ZR(K,J) - A(K,I) * ZI(K,J)
               SIA(J) = SIA(J) + A(I,K) * ZI(K,J) + A(K,I) * ZR(K,J)
  110       CONTINUE
C     ********** DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW **********
            DO 112 J=1,M
            SA(J) = (SA(J) / H)
  112       SIA(J) = (SIA(J) / H)
            DO 114 J=1,M
            SA(J) = SA(J)/H
  114       SIA(J)=SIA(J)/H
C
            DO 120 K = 1, L
            DO 120 J=1,M
               ZR(K,J) = ZR(K,J) - SA(J) * A(I,K) - SIA(J) * A(K,I)
               ZI(K,J) = ZI(K,J) - SIA(J) * A(I,K) + SA(J) * A(K,I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
