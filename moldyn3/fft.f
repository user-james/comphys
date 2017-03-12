C     ****************************************************************         
      SUBROUTINE fft1d(A,WORK,N,ISIGN, MNB, TRIG)          
      implicit double precision (a-h, o-z)
C     ****************************************************************         
c
c  This routine returns the Fourier Transform of a function defined
c  on a grid of points. The array A contains the function and is
c  over-written on output with the Fourier transform.
c
c  A defines the function over one period (T) and the function is assumed
c  to be periodic outside that range. The N points at which the function
c  is given are equally spaced:
c         A(2n-1) = real part of function at t = n*T/N.
c         A(2n)   = imag part of function at t = n*T/N.
c
c  The total number of points should be even and of the form 
c                N = 2^p * 3^q * 5^r 
c  i.e. it only has prime factors 2,3 and 5.
c
c  NOTE: For these applications N should not be larger than 10^4.
c
c  On return, A(2n-1) and A(2n) contain the real and imaginary parts
c  of the Fourier Coefficients of the periodic function at values 
c              omega = n * 2*pi/(N*T), 
c  where n=1,...,N, and T is the (assumed) period of the function  
c
c Input Parameters:
c
c    A     =  a double precision array of dimension at least 2*N
c             (defines function as above)
c    WORK  =  a double precision array of dimension at least 2*N
c             (a scratch array for transform routines)
c    N     =  number of points at which function is given (see above).
c    ISIGN =  +1 for forward transform
c          =  -1, for inverse transform
c    MNB   =  (dimension of arrays)/2   (i.e. maximum allowed N)
c    TRIG  =  a double precision array of dimension at least 2*N
c             (a scratch array for transform routines)
c  
c **********************************************************************
c
c The following subroutines form the package for the Fast-Fourier 
c Transforms:
c
c (1) fft1d:
c     This is the main driver subroutine and is the only routine called
c     from outside the package.
c
c (5) CFFT99:
c     Main driver of package routines.
c
c (6) CFTFAX:
c     Set-up routine for transforms.
c     A package routine.
c 
c (7) FACT:
c     Called from CFTFAX. It extracts the prime factorization of the
c     lengths of the transforms in the various directions.
c     A package routine.
c
c (8) CFTRIG:
c     Called from CFTFAX. Calculates some sines and cosines.
c     A package routine.
c
c (9) VPASSM:
c     Called from CFFT99. The core of the Fourier transform process.
c     A package routine.
c
c *********************************************************************
      REAL*8 A(2*MNB),WORK(2*MNB)       
      REAL*8 TRIG(2*MNB)
      INTEGER N, ISIGN, MNB
      INTEGER IFAC(13)

      DATA KEY/0/  
C        
C     1-D FAST FOURIER TRANSFORM       
C        
C     ISIGN = +1 : FORWARD : REAL --> RECIP : EXP(-IQR)    
C     ISIGN = -1 : BACKWARD: RECIP --> REAL : EXP(+IQR)    
C        
c  The initial set-up for the transform routines.
c
      IF (N.GT.MNB) GO TO 900      
c
      CALL CFTFAX(N,IFAC,TRIG)      
      IF (IFAC(1).LE.0) STOP 44       
c
c end of initial set-up
C        
      IF (ISIGN.LT.0) GO TO 100        
      PIN=1.d0/dFLOAT(N) 
      DO 50 I=1,N 
      A(2*I-1)=A(2*I-1) * PIN  
   50 A(2*I)=A(2*I) * PIN  
  100 CONTINUE     
C        
      MI=-ISIGN    
C
C     CFFT99 HAS DIFFERENT CONVENTION : ISIGN=-1 --> FORWARD         
C        
      CALL CFFT99(A,WORK,TRIG,IFAC,1,1,N,1,MI)     
C        
      RETURN       
C        
 900  WRITE (6,5) N,MNB
    5 FORMAT (' === ARRAY DIMENSIONS FOR FFT TOO SMALL ==='      
     #    /' N   =',I5/          
     #    /' MNB =',I5)          
      STOP 43      
      END          
C CFFT99     FROM XLIB                                     05/13/82  
      SUBROUTINE CFFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)      
      implicit double precision (a-h, o-z)
C        
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE        
C              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC          
C              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.         
C              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE          
C              PACKAGE RETURNS A SET OF COMPLEX FOURIER    
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE          
C              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS   
C              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.    
C        
C              THE PACKAGE CFFT99 CONTAINS SEVERAL USER-LEVEL ROUTINES:        
C        
C            SUBROUTINE CFTFAX         
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE  
C                BEFORE A SEQUENCE OF CALLS TO CFFT99      
C                (PROVIDED THAT N IS NOT CHANGED).         
C        
C            SUBROUTINE CFFT99         
C                THE ACTUAL TRANSFORM ROUTINE ROUTINE, CABABLE OF    
C                PERFORMING BOTH THE TRANSFORM AND ITS INVERSE.      
C                HOWEVER, AS THE TRANSFORMS ARE NOT NORMALIZED,      
C                THE APPLICATION OF A TRANSFORM FOLLOWED BY ITS      
C                INVERSE WILL YIELD THE ORIGINAL VALUES MULTIPLIED   
C                BY N.       
C        
C        
C ACCESS       *FORTRAN,P=XLIB,SN=CFFT99         
C        
C        
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0,        
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF   
C              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF  
C              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER       
C              COEFFICIENT VECTORS OF LENGTH N IS          
C        
C                   DIMENSION IFAX(13),TRIGS(2*N)          
C                   COMPLEX A(...), WORK(...)    
C        
C                   CALL CFTFAX (N, IFAX, TRIGS) 
C                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)       
C        
C              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND   
C              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR  
C              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED  
C              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.        
C              SEE THE INDIVIDUAL WRITE-UPS FOR CFTFAX AND 
C              CFFT99 BELOW, FOR A DETAILED DESCRIPTION OF THE       
C              ARGUMENTS.    
C        
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN          
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED        
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.  IT WAS      
C              FURTHER MODIFIED FOR THE FULLY COMPLEX CASE BY DAVE   
C              FULKER IN NOVEMBER, 1980.         
C        
C-----------------------------------------------------------------------       
C        
C SUBROUTINE CFTFAX (N,IFAX,TRIGS)     
C        
C PURPOSE      A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE         
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,     
C              PROVIDED THAT N IS NOT CHANGED.   
C        
C ARGUMENT     IFAX(13),TRIGS(2*N)     
C DIMENSIONS       
C        
C ARGUMENTS        
C        
C ON INPUT     N   
C               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR         
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE        
C               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF   
C               THE TRANSFORMS).       
C        
C              IFAX          
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED        
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING 
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.  
C        
C              TRIGS         
C               A REAL ARRAY OF DIMENSION 2*N    
C        
C ON OUTPUT    IFAX          
C               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE     
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED       
C               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS   
C               GREATER THAN 5, IFAX(1) IS SET TO -99.     
C        
C              TRIGS         
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY         
C               USED BY THE CFT ROUTINES.        
C        
C-----------------------------------------------------------------------       
C        
C SUBROUTINE CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)         
C        
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS (UNNORMALIZED) COMPLEX         
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE  
C              TRANSFORMS.  GIVEN A SET OF COMPLEX GRIDPOINT         
C              VECTORS, THE PACKAGE RETURNS A SET OF       
C              COMPLEX FOURIER COEFFICIENT VECTORS, OR VICE          
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE A        
C              NUMBER HAVING NO PRIME FACTORS OTHER THAN   
C              2, 3, AND 5.  THIS ROUTINE IS     
C              OPTIMIZED FOR USE ON THE CRAY-1.  
C        
C ARGUMENT     COMPLEX A(N*INC+(LOT-1)*JUMP), WORK(N*LOT)  
C DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13) 
C        
C ARGUMENTS        
C        
C ON INPUT     A   
C               A COMPLEX ARRAY OF LENGTH N*INC+(LOT-1)*JUMP CONTAINING        
C               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY IS     
C               OVERWRITTEN BY THE RESULTS.      
C        
C               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF       
C               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS          
C               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT         
C               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER'S      
C               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS 
C               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE         
C               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A    
C               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF         
C               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE          
C               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND 
C               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND      
C               IMAGINARY PORTIONS EXPLICITLY.   
C        
C              WORK          
C               A COMPLEX WORK ARRAY OF LENGTH N*LOT OR A REAL ARRAY 
C               OF LENGTH 2*N*LOT.  SEE N.B. ABOVE.        
C        
C              TRIGS         
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.         
C        
C              IFAX          
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.         
C        
C        
C               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED       
C               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED         
C               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.       
C        
C              INC 
C               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENTS      
C               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR    
C               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).         
C        
C              JUMP          
C               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS       
C               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,     
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8        
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF        
C               INC AND JUMP, SEE THE EXAMPLES BELOW.      
C        
C              N   
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF      
C               TRANSFORMS, BELOW).    
C        
C              LOT 
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.  
C        
C              ISIGN         
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER          
C                    COEFFICIENTS.     
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO    
C                    GRIDPOINT VALUES. 
C        
C ON OUTPUT    A   
C               IF ISIGN = -1, AND LOT GRIDPOINT VECTORS ARE SUPPLIED,         
C               EACH CONTAINING THE COMPLEX SEQUENCE:      
C        
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES) 
C        
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH 
C               CONTAINING THE CORRESPONDING N COEFFICIENT VALUES:   
C        
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES) 
C        
C               DEFINED BY:  
C                 C(K) = SUM(J=0,...,N-1)( G(J)*EXP(-2*I*J*K*PI/N) ) 
C                 WHERE I = SQRT(-1)   
C        
C        
C               IF ISIGN = +1, AND LOT COEFFICIENT VECTORS ARE SUPPLIED,       
C               EACH CONTAINING THE COMPLEX SEQUENCE:      
C        
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES) 
C        
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH 
C               CONTAINING THE CORRESPONDING N GRIDPOINT VALUES:     
C        
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES) 
C        
C               DEFINED BY:  
C                 G(J) = SUM(K=0,...,N-1)( G(K)*EXP(+2*I*J*K*PI/N) ) 
C                 WHERE I = SQRT(-1)   
C        
C        
C               A CALL WITH ISIGN=-1 FOLLOWED BY A CALL WITH ISIGN=+1          
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED          
C               BY THE FACTOR N.       
C        
C        
C EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN    
C               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL        
C               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,         
C               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE  
C               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE    
C               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE 
C               ORTHOGONAL DIRECTION.  
C        
C               COMPLEX A(66,9), WORK(64,9)      
C               REAL TRIGS1(128), TRIGS2(18)     
C               INTEGER IFAX1(13), IFAX2(13)     
C        
C               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:           
C        
C               CALL CFTFAX(64, IFAX1, TRIGS1)   
C               CALL CFTFAX( 9, IFAX2, TRIGS2)   
C        
C               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE     
C               STORED IN MEMORY AS FOLLOWS (USING U AND V TO        
C               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND        
C               ASSUMING CONVENTIONAL FORTRAN STORAGE):    
C        
C   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS, 
C        
C   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS, 
C        
C   .       .       .       .         .   .        .        .        
C   .       .       .       .         .   .        .        .        
C   .       .       .       .         .   .        .        .        
C        
C   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS. 
C        
C               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE  
C               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE    
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE       
C               NUMBER OF TRANSFORMS, LOT, TO BE 9, THE INCREMENT    
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,    
C               AND THE INCREMENT BETWEEN THE STARTING POINTS        
C               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST        
C               DIMENSION OF A).       
C        
C               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)         
C        
C               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT,      
C               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFINE      
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE        
C               NUMBER OF TRANSFORMS, LOT, TO BE 64, THE INCREMENT   
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,   
C               AND THE INCREMENT BETWEEN THE STARTING POINTS        
C               FOR EACH TRANSFORM, JUMP, TO BE 1          
C        
C               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)         
C        
C               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONAL      
C               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT      
C               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN          
C               WITH ISIGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF          
C               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).          
C        
C        
C-----------------------------------------------------------------------       
      DIMENSION A(*),WORK(*),TRIGS(*),IFAX(13)    
C        
C     SUBROUTINE 'CFFT99' - MULTIPLE FAST COMPLEX FOURIER TRANSFORM  
C        
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA      
C     WORK IS AN AREA OF SIZE N*LOT    
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES    
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N   
C     INC IS THE INCREMENT WITHIN EACH DATA "VECTOR"       
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)       
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR    
C     N IS THE LENGTH OF THE DATA VECTORS        
C     LOT IS THE NUMBER OF DATA VECTORS          
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT  
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL  
C        
C        
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN   
C     PARALLEL.    
C        
C        
C THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR 
C     CALL Q8QST4 ( 4HXLIB, 6HCFFT99, 6HCFFT99, 10HVERSION 01)       
C        
      NN = N+N     
      INK=INC+INC  
      JUM = JUMP+JUMP        
      NFAX=IFAX(1) 
      JNK = 2      
      JST = 2      
      IF (ISIGN.GE.0) GO TO 30         
C        
C     THE INNERMOST TEMPERTON ROUTINES HAVE NO FACILITY FOR THE      
C     FORWARD (ISIGN = -1) TRANSFORM.  THEREFORE, THE INPUT MUST BE  
C     REARRANGED AS FOLLOWS:           
C        
C     THE ORDER OF EACH INPUT VECTOR,  
C        
C     G(0), G(1), G(2), ... , G(N-2), G(N-1)     
C        
C     IS REVERSED (EXCLUDING G(0)) TO YIELD      
C        
C     G(0), G(N-1), G(N-2), ... , G(2), G(1).    
C        
C     WITHIN THE TRANSFORM, THE CORRESPONDING EXPONENTIAL MULTIPLIER 
C     IS THEN PRECISELY THE CONJUGATE OF THAT FOR THE NORMAL         
C     ORDERING.  THUS THE FORWARD (ISIGN = -1) TRANSFORM IS          
C     ACCOMPLISHED 
C        
C     FOR NFAX ODD, THE INPUT MUST BE TRANSFERRED TO THE WORK ARRAY, 
C     AND THE REARRANGEMENT CAN BE DONE DURING THE MOVE.   
C        
      JNK = -2     
      JST = NN-2   
      IF (MOD(NFAX,2).EQ.1) GOTO 40    
C        
C     FOR NFAX EVEN, THE REARRANGEMENT MUST BE APPLIED DIRECTLY TO   
C     THE INPUT ARRAY.  THIS CAN BE DONE BY SWAPPING ELEMENTS.       
C        
      IBASE = 1    
      ILAST = (N-1)*INK      
      NH = N/2     
      DO 20 L=1,LOT          
      I1 = IBASE+INK         
      I2 = IBASE+ILAST       
CDIR$ IVDEP        
      DO 10 M=1,NH 
C     SWAP REAL AND IMAGINARY PORTIONS 
      HREAL = A(I1)          
      HIMAG = A(I1+1)        
      A(I1) = A(I2)          
      A(I1+1) = A(I2+1)      
      A(I2) = HREAL          
      A(I2+1) = HIMAG        
      I1 = I1+INK  
      I2 = I2-INK  
   10 CONTINUE     
      IBASE = IBASE+JUM      
   20 CONTINUE     
      GOTO 100     
C        
   30 CONTINUE     
      IF (MOD(NFAX,2).EQ.0) GOTO 100   
C        
   40 CONTINUE     
C        
C     DURING THE TRANSFORM PROCESS, NFAX STEPS ARE TAKEN, AND THE    
C     RESULTS ARE STORED ALTERNATELY IN WORK AND IN A.  IF NFAX IS   
C     ODD, THE INPUT DATA ARE FIRST MOVED TO WORK SO THAT THE FINAL  
C     RESULT (AFTER NFAX STEPS) IS STORED IN ARRAY A.      
C        
      IBASE=1      
      JBASE=1      
      DO 60 L=1,LOT          
C     MOVE REAL AND IMAGINARY PORTIONS OF ELEMENT ZERO     
      WORK(JBASE) = A(IBASE) 
      WORK(JBASE+1) = A(IBASE+1)       
      I=IBASE+INK  
      J=JBASE+JST  
CDIR$ IVDEP        
      DO 50 M=2,N  
C     MOVE REAL AND IMAGINARY PORTIONS OF OTHER ELEMENTS (POSSIBLY IN          
C     REVERSE ORDER, DEPENDING ON JST AND JNK)   
      WORK(J) = A(I)         
      WORK(J+1) = A(I+1)     
      I=I+INK      
      J=J+JNK      
   50 CONTINUE     
      IBASE=IBASE+JUM        
      JBASE=JBASE+NN         
   60 CONTINUE     
C        
  100 CONTINUE     
C        
C     PERFORM THE TRANSFORM PASSES, ONE PASS FOR EACH FACTOR.  DURING          
C     EACH PASS THE DATA ARE MOVED FROM A TO WORK OR FROM WORK TO A. 
C        
C     FOR NFAX EVEN, THE FIRST PASS MOVES FROM A TO WORK   
      IGO = 110    
C     FOR NFAX ODD, THE FIRST PASS MOVES FROM WORK TO A    
      IF (MOD(NFAX,2).EQ.1) IGO = 120  
      LA=1         
      DO 140 K=1,NFAX        
      IF (IGO.EQ.120) GO TO 120        
  110 CONTINUE     
      CALL VPASSM(A(1),A(2),WORK(1),WORK(2),TRIGS,         
     *   INK,2,JUM,NN,LOT,N,IFAX(K+1),LA)        
      IGO=120      
      GO TO 130    
  120 CONTINUE     
      CALL VPASSM(WORK(1),WORK(2),A(1),A(2),TRIGS,         
     *    2,INK,NN,JUM,LOT,N,IFAX(K+1),LA)       
      IGO=110      
  130 CONTINUE     
      LA=LA*IFAX(K+1)        
  140 CONTINUE     
C        
C     AT THIS POINT THE FINAL TRANSFORM RESULT IS STORED IN A.       
C        
      RETURN       
      END          
      SUBROUTINE CFTFAX(N,IFAX,TRIGS)  
      implicit double precision (a-h, o-z)
      DIMENSION IFAX(13),TRIGS(*)      
C        
C     THIS ROUTINE WAS MODIFIED FROM TEMPERTON'S ORIGINAL  
C     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING    
C     ORDER, AND THERE ARE NONE OF THE ORIGINAL "MODE" OPTIONS.      
C        
C ON INPUT     N   
C               THE LENGTH OF EACH COMPLEX TRANSFORM TO BE PERFORMED 
C        
C               N MUST BE GREATER THAN 1 AND CONTAIN NO PRIME        
C               FACTORS GREATER THAN 5.          
C        
C ON OUTPUT    IFAX          
C               IFAX(1)      
C                 THE NUMBER OF FACTORS CHOSEN OR -99 IN CASE OF ERROR         
C               IFAX(2) THRU IFAX( IFAX(1)+1 )   
C                 THE FACTORS OF N IN THE FOLLOWIN ORDER:  APPEARING 
C                 FIRST ARE AS MANY FACTORS OF 4 AS CAN BE OBTAINED. 
C                 SUBSEQUENT FACTORS ARE PRIMES, AND APPEAR IN       
C                 ASCENDING ORDER, EXCEPT FOR MULTIPLE FACTORS.      
C        
C              TRIGS         
C               2N SIN AND COS VALUES FOR USE BY THE TRANSFORM ROUTINE         
C        
      CALL FACT(N,IFAX)      
      K = IFAX(1)  
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99    
      IF (IFAX(1) .LE. 0) WRITE(6,991) IFAX(1)   
      IF (IFAX(1) .LE. 0) STOP 'FFTFAX'          
991   FORMAT(' FFTFAX -- INVALID N', I5)         
      CALL CFTRIG (N, TRIGS) 
      RETURN       
      END          
      SUBROUTINE FACT(N,IFAX)          
      implicit double precision (a-h, o-z)
C     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4     
      DIMENSION IFAX(13)     
      IF (N.GT.1) GO TO 10   
      IFAX(1) = 0  
      IF (N.LT.1) IFAX(1) = -99        
      RETURN       
   10 NN=N         
      K=1          
C     TEST FOR FACTORS OF 4  
   20 IF (MOD(NN,4).NE.0) GO TO 30     
      K=K+1        
      IFAX(K)=4    
      NN=NN/4      
      IF (NN.EQ.1) GO TO 80  
      GO TO 20     
C     TEST FOR EXTRA FACTOR OF 2       
   30 IF (MOD(NN,2).NE.0) GO TO 40     
      K=K+1        
      IFAX(K)=2    
      NN=NN/2      
      IF (NN.EQ.1) GO TO 80  
C     TEST FOR FACTORS OF 3  
   40 IF (MOD(NN,3).NE.0) GO TO 50     
      K=K+1        
      IFAX(K)=3    
      NN=NN/3      
      IF (NN.EQ.1) GO TO 80  
      GO TO 40     
C     NOW FIND REMAINING FACTORS       
   50 L=5          
      MAX = SQRT(dFLOAT(NN))  
      INC=2        
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4    
   60 IF (MOD(NN,L).NE.0) GO TO 70     
      K=K+1        
      IFAX(K)=L    
      NN=NN/L      
      IF (NN.EQ.1) GO TO 80  
      GO TO 60     
   70 IF (L.GT.MAX) GO TO 75 
      L=L+INC      
      INC=6-INC    
      GO TO 60     
   75 K = K+1      
      IFAX(K) = NN 
   80 IFAX(1)=K-1  
C     IFAX(1) NOW CONTAINS NUMBER OF FACTORS     
      RETURN       
      END          
      SUBROUTINE CFTRIG(N,TRIGS)       
      implicit double precision (a-h, o-z)
      DIMENSION TRIGS(*)     
      PI=2.0d0*dASIN(1.0d0)       
      DEL=(PI+PI)/dFLOAT(N)   
      L=N+N        
      DO 10 I=1,L,2          
      ANGLE=0.5d0*dFLOAT(I-1)*DEL         
      TRIGS(I)=COS(ANGLE)    
      TRIGS(I+1)=SIN(ANGLE)  
   10 CONTINUE     
      RETURN       
      END          
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)       
      implicit double precision (a-h, o-z)
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)     
C        
C     SUBROUTINE 'VPASSM' - MULTIPLE VERSION OF 'VPASSA'   
C     PERFORMS ONE PASS THROUGH DATA   
C     AS PART OF MULTIPLE COMPLEX (INVERSE) FFT ROUTINE    
C     A IS FIRST REAL INPUT VECTOR     
C     B IS FIRST IMAGINARY INPUT VECTOR          
C     C IS FIRST REAL OUTPUT VECTOR    
C     D IS FIRST IMAGINARY OUTPUT VECTOR         
C     TRIGS IS PRECALCULATED TABLE OF SINES & COSINES      
C     INC1 IS ADDRESSING INCREMENT FOR A AND B   
C     INC2 IS ADDRESSING INCREMENT FOR C AND D   
C     INC3 IS ADDRESSING INCREMENT BETWEEN A'S & B'S       
C     INC4 IS ADDRESSING INCREMENT BETWEEN C'S & D'S       
C     LOT IS THE NUMBER OF VECTORS     
C     N IS LENGTH OF VECTORS 
C     IFAC IS CURRENT FACTOR OF N      
C     LA IS PRODUCT OF PREVIOUS FACTORS          
C        
      DATA SIN36/0.587785252292473d0/,COS36/0.809016994374947d0/,        
     *     SIN72/0.951056516295154d0/,COS72/0.309016994374947d0/,        
     *     SIN60/0.866025403784437d0/    
C        
      M=N/IFAC     
      IINK=M*INC1  
      JINK=LA*INC2 
      JUMP=(IFAC-1)*JINK     
      IBASE=0      
      JBASE=0      
      IGO=IFAC-1   
      IF (IGO.GT.4) RETURN   
      GO TO (10,50,90,130),IGO         
C        
C     CODING FOR FACTOR 2    
C        
   10 IA=1         
      JA=1         
      IB=IA+IINK   
      JB=JA+JINK   
      DO 20 L=1,LA 
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 15 IJK=1,LOT        
      C(JA+J)=A(IA+I)+A(IB+I)          
      D(JA+J)=B(IA+I)+B(IB+I)          
      C(JB+J)=A(IA+I)-A(IB+I)          
      D(JB+J)=B(IA+I)-B(IB+I)          
      I=I+INC3     
      J=J+INC4     
   15 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
   20 CONTINUE     
      IF (LA.EQ.M) RETURN    
      LA1=LA+1     
      JBASE=JBASE+JUMP       
      DO 40 K=LA1,M,LA       
      KB=K+K-2     
      C1=TRIGS(KB+1)         
      S1=TRIGS(KB+2)         
      DO 30 L=1,LA 
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 25 IJK=1,LOT        
      C(JA+J)=A(IA+I)+A(IB+I)          
      D(JA+J)=B(IA+I)+B(IB+I)          
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))    
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))    
      I=I+INC3     
      J=J+INC4     
   25 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
   30 CONTINUE     
      JBASE=JBASE+JUMP       
   40 CONTINUE     
      RETURN       
C        
C     CODING FOR FACTOR 3    
C        
   50 IA=1         
      JA=1         
      IB=IA+IINK   
      JB=JA+JINK   
      IC=IB+IINK   
      JC=JB+JINK   
      DO 60 L=1,LA 
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 55 IJK=1,LOT        
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))          
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))          
      C(JB+J)=(A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))
     *  -(SIN60*(B(IB+I)-B(IC+I)))        
      C(JC+J)=(A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))
     *  +(SIN60*(B(IB+I)-B(IC+I)))        
      D(JB+J)=(B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))
     *  +(SIN60*(A(IB+I)-A(IC+I)))        
      D(JC+J)=(B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))
     *  -(SIN60*(A(IB+I)-A(IC+I)))        
      I=I+INC3     
      J=J+INC4     
   55 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
   60 CONTINUE     
      IF (LA.EQ.M) RETURN    
      LA1=LA+1     
      JBASE=JBASE+JUMP       
      DO 80 K=LA1,M,LA       
      KB=K+K-2     
      KC=KB+KB     
      C1=TRIGS(KB+1)         
      S1=TRIGS(KB+2)         
      C2=TRIGS(KC+1)         
      S2=TRIGS(KC+2)         
      DO 70 L=1,LA 
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 65 IJK=1,LOT        
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))          
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))          
      C(JB+J)=     
     *  C1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))       
     * -S1*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))       
      D(JB+J)=     
     *  S1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))       
     * +C1*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))       
      C(JC+J)=     
     *  C2*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))       
     * -S2*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))       
      D(JC+J)=     
     *  S2*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))       
     * +C2*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))       
      I=I+INC3     
      J=J+INC4     
   65 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
   70 CONTINUE     
      JBASE=JBASE+JUMP       
   80 CONTINUE     
      RETURN       
C        
C     CODING FOR FACTOR 4    
C        
   90 IA=1         
      JA=1         
      IB=IA+IINK   
      JB=JA+JINK   
      IC=IB+IINK   
      JC=JB+JINK   
      ID=IC+IINK   
      JD=JC+JINK   
      DO 100 L=1,LA          
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 95 IJK=1,LOT        
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))          
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))          
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))          
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))          
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))          
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))          
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))          
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))          
      I=I+INC3     
      J=J+INC4     
   95 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
  100 CONTINUE     
      IF (LA.EQ.M) RETURN    
      LA1=LA+1     
      JBASE=JBASE+JUMP       
      DO 120 K=LA1,M,LA      
      KB=K+K-2     
      KC=KB+KB     
      KD=KC+KB     
      C1=TRIGS(KB+1)         
      S1=TRIGS(KB+2)         
      C2=TRIGS(KC+1)         
      S2=TRIGS(KC+2)         
      C3=TRIGS(KD+1)         
      S3=TRIGS(KD+2)         
      DO 110 L=1,LA          
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 105 IJK=1,LOT       
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))          
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))          
      C(JC+J)=     
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))         
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))         
      D(JC+J)=     
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))         
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))         
      C(JB+J)=     
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))         
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))         
      D(JB+J)=     
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))         
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))         
      C(JD+J)=     
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))         
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))         
      D(JD+J)=     
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))         
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))         
      I=I+INC3     
      J=J+INC4     
  105 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
  110 CONTINUE     
      JBASE=JBASE+JUMP       
  120 CONTINUE     
      RETURN       
C        
C     CODING FOR FACTOR 5    
C        
  130 IA=1         
      JA=1         
      IB=IA+IINK   
      JB=JA+JINK   
      IC=IB+IINK   
      JC=JB+JINK   
      ID=IC+IINK   
      JD=JC+JINK   
      IE=ID+IINK   
      JE=JD+JINK   
      DO 140 L=1,LA          
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 135 IJK=1,LOT       
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))  
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))  
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))) 
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))) 
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))) 
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))) 
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))) 
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))) 
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))) 
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))) 
      I=I+INC3     
      J=J+INC4     
  135 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
  140 CONTINUE     
      IF (LA.EQ.M) RETURN    
      LA1=LA+1     
      JBASE=JBASE+JUMP       
      DO 160 K=LA1,M,LA      
      KB=K+K-2     
      KC=KB+KB     
      KD=KC+KB     
      KE=KD+KB     
      C1=TRIGS(KB+1)         
      S1=TRIGS(KB+2)         
      C2=TRIGS(KC+1)         
      S2=TRIGS(KC+2)         
      C3=TRIGS(KD+1)         
      S3=TRIGS(KD+2)         
      C4=TRIGS(KE+1)         
      S4=TRIGS(KE+2)         
      DO 150 L=1,LA          
      I=IBASE      
      J=JBASE      
CDIR$ IVDEP        
      DO 145 IJK=1,LOT       
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))  
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))  
      C(JB+J)=     
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))      
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))      
      D(JB+J)=     
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))      
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))      
      C(JE+J)=     
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))      
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))      
      D(JE+J)=     
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))        
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))      
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))        
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))      
      C(JC+J)=     
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))      
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))      
      D(JC+J)=     
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))      
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))      
      C(JD+J)=     
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))      
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))      
      D(JD+J)=     
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))        
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))      
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))        
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))      
      I=I+INC3     
      J=J+INC4     
  145 CONTINUE     
      IBASE=IBASE+INC1       
      JBASE=JBASE+INC2       
  150 CONTINUE     
      JBASE=JBASE+JUMP       
  160 CONTINUE     
      RETURN       
      END          
      SUBROUTINE EMPTY(I)    
      RETURN       
      END          
