C******************************************************************************|
C duct.f, the duct-flow solvers for diablo.                        VERSION 0.9
C These solvers were written by ? and ? (spring 2001).
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the duct-flow case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the duct-flow case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      INTEGER I,J,K

         IF (RANK.EQ.0)
     &     WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3 .AND. RANK.EQ.0)
     &          WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
         IF (RANK.EQ.0)
     &        WRITE (6,*) 'Finite-difference in Z'
         IF (RANK.EQ.0)
     &        WRITE (6,*) 'Finite-difference in Y'

         IF (RANK.EQ.0)
     &        write(*,*) 'USE_MPI: ',USE_MPI

         FNAME='grid.h5'
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if

         call ReadGridHDF5(FNAME,2)
         DO J=0,NY+1
           IF (VERBOSITY .GT. 3 .AND. RANKZ.EQ.0)
     &          WRITE(6,*) 'RANKY=',RANKY,'GYF(',J,') = ',GYF(J)
         END DO

         call ReadGridHDF5(FNAME,3)
         DO K=0,NZP+1
           IF (VERBOSITY .GT. 3 .AND. RANKY.EQ.0)
     &          WRITE(6,*) 'RANKZ=',RANKZ,'GZFP(',K,') = ',GZFP(K)
         END DO


C Define grid spacing - Y
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DYF(NY+1)=DYF(NY)
C Define grid spacing - Z
         DO K=1,NZP+1
           DZP(K)=(GZFP(K)-GZFP(K-1))
         END DO
         DO K=1,NZP
           DZFP(K)=(GZP(K+1)-GZP(K))
         END DO
         DZFP(NZP+1)=DZFP(NZP)

         RETURN
         END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL    VERSION, CURRENT_VERSION
      INTEGER I,J,K,N

! Read in input parameters specific for duct flow case
      OPEN (11,file='input_duct.dat',form='formatted',status='old')
C Read input file.

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION)
     &         STOP 'Wrong input data format input_duct'
      READ(11,*)
      READ(11,*) TIME_AD_METH
      READ(11,*)
      READ(11,*) LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      READ(11,*)
      READ(11,*) I_RO
      READ(11,*)
      READ(11,*) GRAV_X, GRAV_Y, GRAV_Z
      READ(11,*)
      READ(11,*) F_TYPE, UBULK0, PX0, OMEGA0, AMP_OMEGA0
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*)
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3
      READ(11,*)
      READ(11,*) U_BC_ZMIN, U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) V_BC_ZMIN, V_BC_ZMIN_C1, V_BC_ZMIN_C2, V_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) W_BC_ZMIN, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) U_BC_ZMAX, U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) V_BC_ZMAX, V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) W_BC_ZMAX, W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3
      READ(11,*)
! Read in boundary conditions and background gradients for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N)
     &             ,TH_BC_YMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N)
     &             ,TH_BC_YMAX_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_ZMIN(N),TH_BC_ZMIN_C1(N),TH_BC_ZMIN_C2(N)
     &             ,TH_BC_ZMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_ZMAX(N),TH_BC_ZMAX_C1(N),TH_BC_ZMAX_C2(N)
     &             ,TH_BC_ZMAX_C3(N)
        READ(11,*)
        READ(11,*) DRHODX(N), DRHODZ(N)
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_DUCT_MOVIE
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      INTEGER I,J,K

! Set parameters for writing "movie" files with 2d slices

! First, read in the center of each 2D slice
      OPEN(unit=650,file='MOVIE.dat',status='old',
     &           form='formatted')
      READ(650, *) XcMovie, YcMovie, ZcMovie
      CLOSE(650)

! Get the indices
      NxMovie=int(XcMovie*NX/LX)

      RankYMovie=-1
      IF (GYF(JSTART).LE.YcMovie.and.GYF(JEND+1).GT.YcMovie) THEN
         RankYMovie=RANKY
         I=1
         DO WHILE (.not.(GYF(I).LE.YcMovie .and. GYF(I+1).GT.YcMovie))
            I=I+1
         END DO
         NyMovie=I;
      END IF

      RankZMovie=-1
      IF (GZFP(KSTART).LE.ZcMovie.and.GZFP(KEND+1).GT.ZcMovie) THEN
         RankZMovie=RANKZ
         I=1
         DO WHILE (.not.(GZFP(I).LE.ZcMovie .and. GZFP(I+1).GT.ZcMovie))
            I=I+1
         END DO
         NzMovie=I;
      END IF

      IF (RANKY.eq.RankYMovie .and. RANKZ.eq.RankZMovie) THEN
         write(*,*) 'Movie Parameters, RANK:', RANK
         write(*,*) '    Xc: ', GX(NxMovie),
     &      ' (NxMovie: ', NxMovie, ')'
         write(*,*) '    Yc: ', GYF(NyMovie),
     &      ' (NyMovie: ', NyMovie, ')'
         write(*,*) '    Zc: ', GZFP(NzMovie),
     &      ' (NzMovie: ', NzMovie, ')'
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_DUCT(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      LOGICAL FINAL

      RETURN
      END
            
