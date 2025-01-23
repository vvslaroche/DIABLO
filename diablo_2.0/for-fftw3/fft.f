C******************************************************************************|
C fft.f, the FFT package for diablo.                               VERSION 0.9
C
C This file isolates all calls to the FFTW package (available at: www.fftw.org)
C These wrapper routines were written by T. Bewley (spring 2001).
C******************************************************************************|
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C The arrangement of the significant real numbers in the arrays (denoted by +)
C in physical space, in Fourier space, and in Fourier space after packing are
C shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
C an identical matter as the Z direction shown here.
C
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C NZ-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo    -NKZ ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo    -NKZ ++++++++++++oooooo
C      ++++++++++++++++oo     NKZ ++++++++++++oooooo     NKZ ++++++++++++oooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
C   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
C   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
C   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
C   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
C      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
C      0123           NX-1        0 1 2     NKX              0 1 2     NKX
C
C       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
C
C After the Real->Fourier transform, the significant coefficients are put next
C to each other in the array, so a loop such as
C
C        DO K=0,TNKZ           [where TNKZ = 2*NKZ = 2*(NZ/3) ]
C          DO I=0,NKX          [where  NKX = NX/3             ]
C            CP(I,K,J)= ...
C          END DO
C        END DO
C
C includes all the Fourier coefficients of interest.  The subsequent loops in
C Fourier space just work on these coefficients in the matrix.
C  
C Before a Fourier->Real transform, the significant coefficients are unpacked
C and the higher wavenumbers are SET TO ZERO before the inverse transform.
C This has the effect of doing the required dealiasing.
C
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_FFT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE 'fftw3.f'
      INTEGER I,J,K

      INTEGER HOWMANY_N(2), HOWMANY_STRIDE_R(2), HOWMANY_STRIDE_C(2)
      REAL*8      V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16  VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16  TMP(0:NX/2,0:NZP+1,0:NY+1)

    !   INTEGER         FFTW_FORWARD,      FFTW_BACKWARD,
    !  *                FFTW_ESTIMATE,     FFTW_MEASURE,
    !  *                FFTW_OUT_OF_PLACE, FFTW_IN_PLACE,
    !  *                FFTW_USE_WISDOM,   FFTW_THREADSAFE,
    !  *                FFTW_PATIENT,      FFTW_DESTROY_INPUT
    !   PARAMETER(      FFTW_FORWARD=-1,      FFTW_BACKWARD=1,
    !  *                FFTW_ESTIMATE=0,      FFTW_MEASURE=1,
    !  *                FFTW_OUT_OF_PLACE=0,  FFTW_IN_PLACE=8,
    !  *                FFTW_USE_WISDOM=16,   FFTW_THREADSAFE=128,
    !  *                FFTW_PATIENT=1,       FFTW_DESTROY_INPUT=1 )

      IF (RANK.EQ.0) 
     &     WRITE(6,*) 'Initializing FFTW package.'

      PI = 4. * ATAN(1.D0)
      CI = CMPLX(0.0,1.0)
      EPS= 0.000000001

      IF (NUM_PER_DIR .GT. 0) THEN
         HOWMANY_N(1) = NZP
         HOWMANY_N(2) = NY + 2
         HOWMANY_STRIDE_R(1) = NX + 2
         HOWMANY_STRIDE_R(2) = (NZP + 2)*(NX + 2)
         HOWMANY_STRIDE_C(1) = INT(NX/2 + 1)
         HOWMANY_STRIDE_C(2) = INT(NX/2 + 1)*(NZP + 2)
    !      CALL DFFTW_PLAN_GURU_DFT_R2C(FFTW_X_TO_F_PLAN, 1, NX, 1, 1,
    !  *        2, HOWMANY_N, HOWMANY_STRIDE_R, HOWMANY_STRIDE_C,
    !  *        V(0,0,0), TEMP_FFT(0,0,0),
    !  *        FFTW_PATIENT + FFTW_DESTROY_INPUT)
    !      CALL DFFTW_PLAN_GURU_DFT_C2R(FFTW_X_TO_P_PLAN, 1, NX, 1, 1,
    !  *        2, HOWMANY_N, HOWMANY_STRIDE_C, HOWMANY_STRIDE_R,
    !  *        TEMP_FFT(0,0,0), V(0,0,0),
    !  *        FFTW_PATIENT + FFTW_DESTROY_INPUT)
    !      CALL DFFTW_PLAN_DFT_R2C_1D(FFTW_X_TO_F_PLAN, NX,
    !  *        V, TMP, FFTW_MEASURE)
    !      CALL DFFTW_PLAN_DFT_C2R_1D(FFTW_X_TO_P_PLAN, NX,
    !  *        TMP, V, FFTW_MEASURE)
         CALL DFFTW_PLAN_GURU_DFT_R2C(FFTW_X_TO_F_PLAN, 1, NX, 1, 1,
     *        1, NZP, NX+2, NX/2+1, V, TMP, FFTW_MEASURE)
         CALL DFFTW_PLAN_GURU_DFT_C2R(FFTW_X_TO_P_PLAN, 1, NX, 1, 1,
     *        1, NZP, NX/2+1, NX+2, TMP, V, FFTW_MEASURE)
    !      CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN, 1, NX,      
    !  *        FFTW_FORWARD,  FFTW_MEASURE  ) 
    !      CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX,            
    !  *        FFTW_BACKWARD,  FFTW_MEASURE  ) 
        RNX=1.0*NX
        DO I=0,NXPP-1
          KX(I)=(I+NXPP*RANKZ)*(2.*PI)/LX
          KX2(I)=KX(I)*KX(I)
          CIKX(I)=CI*KX(I)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
        HOWMANY_N(1) = NXP
        HOWMANY_N(2) = NY + 2
        HOWMANY_STRIDE_C(1) = 1
        HOWMANY_STRIDE_C(2) = (NXP + 1)*(NZ + 2)
    !     CALL DFFTW_PLAN_GURU_DFT(FFTW_Z_TO_F_PLAN, 1, NZ, NXP+1, NXP+1
    !  *        2, HOWMANY_N, HOWMANY_STRIDE_C, HOWMANY_STRIDE_C,
    !  *        VV(0,0,0), VV(0,0,0),
    !  *        FFTW_FORWARD, FFTW_PATIENT)
    !     CALL DFFTW_PLAN_GURU_DFT(FFTW_Z_TO_P_PLAN, 1, NZ,
    !  *        1, 1, 2, HOWMANY_N, HOWMANY_STRIDE_C, HOWMANY_STRIDE_C,
    !  *        VV(0,0,0), VV(0,0,0),
    !  *        FFTW_FORWARD, FFTW_PATIENT)
    !     CALL DFFTW_PLAN_DFT_1D(FFTW_Z_TO_F_PLAN, NZ,
    !  *        VV, VV, FFTW_FORWARD, FFTW_MEASURE)
    !     CALL DFFTW_PLAN_DFT_1D(FFTW_Z_TO_P_PLAN, NZ,
    !  *        VV, VV, FFTW_BACKWARD, FFTW_MEASURE)
        CALL DFFTW_PLAN_GURU_DFT(FFTW_Z_TO_F_PLAN, 1, NZ, NXP+1, NXP+1,
     *       1, NXP, 1, 1, VV, VV, FFTW_FORWARD, FFTW_MEASURE)
        CALL DFFTW_PLAN_GURU_DFT(FFTW_Z_TO_P_PLAN, 1, NZ, NXP+1, NXP+1,
     *       1, NXP, 1, 1, VV, VV, FFTW_BACKWARD, FFTW_MEASURE)
    !     CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_F_PLAN, 1, NZ,
    !  *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
    !     CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_P_PLAN, 1, NZ,
    !  *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
        RNZ=1.0*NZ
        DO K=0,NKZ
          KZ(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
           KZ(TNKZ+1-K)=-K*(2.*PI)/LZ
        END DO
        DO K=0,TNKZ
          KZ2(K)=KZ(K)*KZ(K)
          CIKZ(K)=CI*KZ(K)
        END DO
      END IF

      IF (RANK.EQ.0) 
     &     write(*,*) 'In fft: NKX,TNKZ: ',NKX,TNKZ

    !   IF (NUM_PER_DIR .GT. 2) THEN
    !     CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN, 1, NY,
    !  *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
    !     CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN, 1, NY,
    !  *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE )
    !     RNY=1.0*NY
	  ! KY(0) = 0.
	  ! KY2(0) = 0.
	  ! CIKY(0) = (0.0,0.0)
    !     DO J=1,NKY
    !       KY(J)=J*(2.*PI)/LY
    !     END DO
    !     DO J=1,NKY
    !       KY(TNKY+1-J)=-J*(2.*PI)/LY
    !     END DO
    !     DO J=1,TNKY
    !       KY2(J)=KY(J)*KY(J)
    !       CIKY(J)=CI*KY(J)
    !     END DO
    !   END IF

      DO I=0,NKX
        DO K=0,NZ
          CZX_PLANE(K,I)=CMPLX(0.0,0.0)
        END DO
      END DO
      DO K=0,TNKZ
        DO J=0,NY
          CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
        END DO
      END DO

      IF (RANK.EQ.0) 
     &     WRITE(6,*) 'FFTW package initialized.'

      RETURN
      END

C******************************************************************************|
C-------------> The transform routines for the duct flow follow. <-------------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER(U,CU,JMIN,JMAX,KMIN,KMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C Looping over the planes of interest, simply perform a real -> complex
C transform in place in the big storage array, scaling appropriately.

    !   DO J=JMIN,JMAX
    !    CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN,(KMAX-KMIN+1),
    !  *    U(0,KMIN,J), 1, NX+2, CU(0,KMIN,J), 1, NX/2+1)
    !     DO K=KMIN,KMAX
    !       DO I=0,NKX
    !         CU(I,K,J)=CU(I,K,J)/RNX
    !       END DO
    !     END DO
    !   END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL(CU,U,JMIN,JMAX,KMIN,KMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C Looping over the planes of interest, simply set the higher wavenumbers to
C zero and then perform a complex -> real transform in place in the big
C storage array.

    !   DO J=JMIN,JMAX
    !     DO K=KMIN,KMAX
    !       DO I=NKX+1,NX/2
    !         CU(I,K,J)=0.
    !       END DO
    !     END DO
    !    CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN,(KMAX-KMIN+1),
    !  *    CU(0,KMIN,J), 1, NX/2+1, U(0,KMIN,J), 1, NX+2)
    !   END DO

      RETURN
      END

C******************************************************************************|
C-----------> The transform routines for the channel flow follow. <------------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

c
c     THE CHANNEL TRANSFORM FFT_XZ_TO_* HAVE BEEN MOVED TO THE MPI FILE 
c     SINCE THEY REQUIRE 2D PARALLELIZATION
c


C******************************************************************************|
C--------> The transform routines for the fully-periodic box follow. <---------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C First, transform in the X & Z directions using FFT_XZ_TO_FOURIER.  Then,
C looping over the planes of interest, put the data into the CYZ_PLANE
C temporary storage variable, perform a complex -> complex transform in the
C y direction, then put the data back into the big storage array, packing the
C data towards J=0 and scaling appropriately.

    !   CALL FFT_XZ_TO_FOURIER(U,CU,0,NYM)
    !   DO I=0,NKX
    !     DO K=0,TNKZ
    !       DO J=0,NYM       
    !         CYZ_PLANE(J,K)=CU(I,K,J)
    !       END DO
    !     END DO        
    !     CALL FFTWND_F77(FFTW_Y_TO_F_PLAN, TNKZ+1,
    !  *    CYZ_PLANE(0,0), 1, NY+1, CYZ_PLANE(0,0), 1, NY+1)
    !     DO K=0,TNKZ 
    !       DO J=0,NKY
    !         CU(I,K,J)=CYZ_PLANE(J,K)/RNY
    !       END DO
    !       DO J=1,NKY
    !         CU(I,K,NKY+J)=CYZ_PLANE(NYM-NKY+J,K)/RNY
    !       END DO
    !     END DO
    !   END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C Looping over the planes of interest, unpack data into the CYZ_PLANE temporary
C storage variable (setting higher wavenumbers to zero), perform a complex ->
C complex transform in the y direction, then put the data back into the big
C storage array.  Finally, transform in the X & Z directions using
C FFT_XZ_TO_PHYSICAL.

    !   DO I=0,NKX
    !     DO K=0,TNKZ 
    !       DO J=0,NKY
    !         CYZ_PLANE(J,K)=CU(I,K,J)
    !       END DO
    !       DO J=NKY+1,NYM-NKY
    !         CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
    !       END DO
    !       DO J=1,NKY
    !         CYZ_PLANE(NYM-NKY+J,K)=CU(I,K,NKY+J)
    !       END DO
    !     END DO
    !     CALL FFTWND_F77(FFTW_Y_TO_P_PLAN, TNKZ+1,
    !  *    CYZ_PLANE(0,0), 1, NY+1, CYZ_PLANE(0,0), 1, NY+1)
    !     DO K=0,TNKZ
    !       DO J=0,NYM       
    !         CU(I,K,J)=CYZ_PLANE(J,K)
    !       END DO
    !     END DO        
    !   END DO
    !   CALL FFT_XZ_TO_PHYSICAL(CU,U,0,NYM)

      RETURN
      END






