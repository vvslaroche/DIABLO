C******************************************************************************|
C duct.f, the duct-flow solvers for diablo.                        VERSION 0.9
C These solvers were written by Vincent Laroche (spring 2025).
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

C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

      INCLUDE 'header'
      INTEGER I,J,K

C First, Initialize the matrix components
      DO J=0,NY+1
        DO K=0,NZ+1
            MATL_C_YZ(J,K)=0.
            MATD_C_YZ(J,K)=1.
            MATU_C_YZ(J,K)=0.
            VEC_C_YZ(J,K)=(0.,0.)
        END DO
      END DO

      IF (RANK.EQ.0) write(*,*) 'CHECK 3.1'
C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO I=0,NXP-1
        DO J=1,NY
          DO K=1,NZ
            MATL_C_YZ(J,K)=1./(DY(J)*DYF(J)) + 1./(DZ(K)*DZF(K))
            MATD_C_YZ(J,K)=-KX2(I)
     &          -1./(DY(J+1)*DYF(J)) - 1./(DY(J)*DYF(J))
     &          -1./(DZ(K+1)*DZF(K)) - 1./(DZ(K)*DZF(K))
            MATU_C_YZ(J,K)=1./(DY(J+1)*DYF(J)) + 1./(DZ(K+1)*DZF(K))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY
          DO K=1,NZ
            VEC_C(J,K)=(CIKX(I)*CU1(I,K,J)
     &           + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
     &           + (CU2(I,K+1,J)-CU2(I,K,J))/DZF(K))
          END DO
        END DO
        IF (RANK.EQ.0) write(*,*) 'CHECK 3.2'
C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
C First, do all forward sweeps
          IF (RANK.EQ.0) write(*,*) 'CHECK 3.3'
          CALL THOMAS_FORWARD_COMPLEX_DUCT_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
C Now, do the backward sweeps
          IF (RANK.EQ.0) write(*,*) 'CHECK 3.4'
          CALL THOMAS_BACKWARD_COMPLEX_DUCT_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)

C Else we are running in serial mode
        ELSE
          IF (RANK.EQ.0) 
     &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

        DO J=1,NY
          DO K=1,NZ
            CR1(I,K,J)=VEC_C(J,K)
          END DO
        END DO
      END DO


C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CU2(I,K,J)=CU2(I,K,J)-(CR1(I,K,J)
     &             -CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      DO J=1,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CU3(I,K,J)=CU3(I,K,J)-(CR1(I,K,J)
     &             -CR1(I,K-1,J))/DZ(K)
          END DO
        END DO
      END DO

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
      INTEGER I,J,K,N
      REAL*8 TMP(0:NPROCZ*(NZP+2)-1)

      IF (RANK.EQ.0)
     &     WRITE (6,*) 'Fourier in X'
      DO I=0,NX
        GX(I)=(I*LX)/NX
        DX(I)=LX/NX
        IF (VERBOSITY .GT. 3 .AND. RANK.EQ.0)
     &       WRITE(6,*) 'GX(',I,') = ',GX(I)
      END DO
      IF (RANK.EQ.0)
     &     WRITE (6,*) 'Finite-difference in Z'
      IF (RANK.EQ.0)
     &     WRITE (6,*) 'Finite-difference in Y'

      IF (RANK.EQ.0)
     &     WRITE(*,*) 'USE_MPI: ',USE_MPI

      FNAME='grid.h5'
      if (USE_MPI) then
        call mpi_barrier(MPI_COMM_WORLD,ierror)
      end if

      call ReadGridHDF5(FNAME,2)
      DO J=0,NY+1
        IF (VERBOSITY .GT. 3 .AND. RANKZ.EQ.0)
     &       WRITE(6,*) 'RANKY=',RANKY,'GYF(',J,') = ',GYF(J)
      END DO

      call ReadGridHDF5(FNAME,3)
      DO K=0,NZP+1
        IF (VERBOSITY .GT. 3 .AND. RANKY.EQ.0)
     &       WRITE(6,*) 'RANKZ=',RANKZ,'GZFP(',K,') = ',GZFP(K)
      END DO


C Reassemble total Z grid on each process
      CALL MPI_GATHER(GZP(0),NZP+2,MPI_DOUBLE_PRECISION,
     &        TMP,NZP+2,MPI_DOUBLE_PRECISION,0, MPI_COMM_Z,IERROR)
      IF (RANK.EQ.0) THEN
      !   DO K=0,NPROCZ*(NZP+2)-1
      !     WRITE(6,*) 'TMP(',K,') = ',TMP(K)
      !   END DO
        ! first Z processor
        DO K=0,NZP
          GZ(K)=TMP(K)
        END DO
        ! last Z processor
        DO K=0,NZP-1
          GZ((NZP-1)*(NPROCZ-1)+2+K)
     &         =TMP((NZP+2)*(NPROCZ-1)+2+K)
        END DO
        ! intermediate Z processors (if NPROCZ>2)
        DO N=1,NPROCZ-2
          DO K=0,NZP-2
            GZ((NZP-1)*N+2+K)=TMP((NZP+2)*N+2+K)
          END DO
        END DO
      !   DO K=0,NZ+1
      !     WRITE(6,*) 'GZ(',K,') = ',GZ(K)
      !   END DO
      END IF
      CALL MPI_BCAST(GZ,NZ+2,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERROR)
      DO K=1,NZ
        GZF(K)=0.5*(GZ(K)+GZ(K+1))
      END DO
      GZF(0)=2.d0*GZF(1)-GZF(2)
      GZF(NZ+1)=2.d0*GZF(NZ)-GZF(NZ-1)


C Define grid spacing - Y
      DO J=1,NY+1
        DY(J)=GYF(J)-GYF(J-1)
      END DO
      DO J=1,NY
        DYF(J)=GY(J+1)-GY(J)
      END DO
      DYF(NY+1)=DYF(NY)
C Define grid spacing - Z (per process)
      DO K=1,NZP+1
        DZP(K)=GZFP(K)-GZFP(K-1)
      END DO
      DO K=1,NZP
        DZFP(K)=GZP(K+1)-GZP(K)
      END DO
      DZFP(NZP+1)=DZFP(NZP)
C Define grid spacing - Z (total)
      DO K=1,NZ+1
        DZ(K)=GZF(K)-GZF(K-1)
      END DO
      DO K=1,NZ
        DZF(K)=GZ(K+1)-GZ(K)
      END DO
      DZF(NZ+1)=DZF(NZ)



!       WRITE(*,*) 'writing grids on rank ', RANK  
!       FNAME='grids_'//char(RANK+48)//'.dat'
!       OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       WRITE(15,151) GX,51.0,DX,52.0,GXF,53.0,DXF,54.0,
!      &              GY,61.0,DY,62.0,GYF,63.0,DYF,64.0,
!      &              GZ,71.0,DZ,72.0,GZF,73.0,DZF,74.0,
!      &              GZP,81.0,DZP,82.0,GZFP,83.0,DZFP
!      &              ,999.0,TMP
!       CLOSE(15)
! 151   format(F30.20)

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

            


C******************************************************************************|
C--------> The boundary condition routines for the duct flow follow. <---------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_VEL_YLOWER_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K

C Dirichlet
      IF (U_BC_YMIN.EQ.0) THEN
C Start with zero
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU1(I,K,0)=0.d0
             CU1(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,K,1)=U_BC_YMIN_C1
           CU1(0,K,0)=U_BC_YMIN_C1
         END IF
C Neumann
      ELSE IF (U_BC_YMIN.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU1(I,K,1)=CU1(I,K,2)
             CU1(I,K,0)=CU1(I,K,1)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,K,1)=CU1(0,K,2)-DY(2)*U_BC_YMIN_C1
           CU1(0,K,0)=CU1(0,K,1)-DY(1)*U_BC_YMIN_C1
         END IF
      ELSE
         STOP 'Error: U_BC_YMIN must be 0, or 1'
      END IF


C Dirichlet
      IF (W_BC_YMIN.EQ.0) THEN
C Start with zero
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU3(I,K,0)=0.d0
             CU3(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,K,1)=W_BC_YMIN_C1
           CU3(0,K,0)=W_BC_YMIN_C1
         END IF
C Neumann
      ELSE IF (W_BC_YMIN.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU3(I,K,1)=CU3(I,K,2)
             CU3(I,K,0)=CU3(I,K,1)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,K,1)=CU3(0,K,2)-DY(2)*W_BC_YMIN_C1
           CU3(0,K,0)=CU3(0,K,1)-DY(1)*W_BC_YMIN_C1
         END IF
      ELSE
         STOP 'Error: W_BC_YMIN must be 0, or 1'
      END IF


C Dirichlet
      IF (V_BC_YMIN.EQ.0) THEN
C Set the vertical velocity at GYF(1) (halfway between GY(1) and GY(2))
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU2(I,K,1)=2.d0*V_BC_YMIN_C1-CU2(I,K,2)
             CU2(I,K,0)=CU2(I,K,1)
           END DO
         END DO
C Neumann
      ELSE IF (V_BC_YMIN.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU2(I,K,1)=CU2(I,K,2)
             CU2(I,K,0)=CU2(I,K,1)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,K,1)=CU2(0,K,2)-DYF(1)*V_BC_YMIN_C1
           CU2(0,K,0)=CU2(0,K,1)-DYF(1)*V_BC_YMIN_C1
         END IF
      ELSE
         STOP 'Error: V_BC_YMIN must be 0, or 1'
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_VEL_YUPPER_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K

C Dirichlet
      IF (U_BC_YMAX.EQ.0) THEN
C Start with zero
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU1(I,K,NY)=0.d0
             CU1(I,K,NY+1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,K,NY)=U_BC_YMAX_C1
           CU1(0,K,NY+1)=U_BC_YMAX_C1
         END IF
C Neumann
      ELSE IF (U_BC_YMAX.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU1(I,K,NY)=CU1(I,K,NY-1)
             CU1(I,K,NY+1)=CU1(I,K,NY)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,K,NY)=CU1(0,K,NY-1)+DY(NY)*U_BC_YMAX_C1
           CU1(0,K,NY+1)=CU1(0,K,NY)+DY(NY)*U_BC_YMAX_C1
         END IF
      ELSE
         STOP 'Error: U_BC_YMAX must be 0, or 1'
      END IF


C Dirichlet
      IF (W_BC_YMAX.EQ.0) THEN
C Start with zero
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU3(I,K,NY)=0.d0
             CU3(I,K,NY+1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,K,NY)=W_BC_YMAX_C1
           CU3(0,K,NY+1)=W_BC_YMAX_C1
         END IF
C Neumann
      ELSE IF (W_BC_YMAX.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU3(I,K,NY)=CU3(I,K,NY-1)
             CU3(I,K,NY+1)=CU3(I,K,NY)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,K,NY)=CU3(0,K,NY-1)+DY(NY)*W_BC_YMAX_C1
           CU3(0,K,NY+1)=CU3(0,K,NY)+DY(NY)*W_BC_YMAX_C1
         END IF
      ELSE
         STOP 'Error: W_BC_YMAX must be 0, or 1'
      END IF


C Dirichlet
      IF (V_BC_YMAX.EQ.0) THEN
C Set the vertical velocity at GYF(NY) (halfway between GY(NY) and GY(NY+1))
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU2(0,K,NY+1)=2.d0*V_BC_YMAX_C1-CU2(0,K,NY)
           END DO
         END DO
C Neumann
      ELSE IF (V_BC_YMAX.EQ.1) THEN
         DO K=0,NZ+1
           DO I=0,NXP-1
             CU2(I,K,NY)=CU2(I,K,NY-1)
             CU2(I,K,NY+1)=CU2(I,K,NY)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,K,NY+1)=CU2(0,K,NY)+DY(NY)*V_BC_YMAX_C1
         END IF
      ELSE
         STOP 'Error: V_BC_YMAX must be 0, or 1'
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_VEL_ZLOWER_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,J

C Dirichlet
      IF (U_BC_ZMIN.EQ.0) THEN
C Start with zero
         DO J=0,NY+1
           DO I=0,NXP-1
             CU1(I,0,J)=0.d0
             CU1(I,1,J)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,1,J)=U_BC_ZMIN_C1
           CU1(0,0,J)=U_BC_ZMIN_C1
         END IF
C Neumann
      ELSE IF (U_BC_ZMIN.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU1(I,1,J)=CU1(I,2,J)
             CU1(I,0,J)=CU1(I,1,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,1,J)=CU1(0,2,J)-DZ(2)*U_BC_ZMIN_C1
           CU1(0,0,J)=CU1(0,1,J)-DZ(1)*U_BC_ZMIN_C1
         END IF
      ELSE
         STOP 'Error: U_BC_ZMIN must be 0, or 1'
      END IF


C Dirichlet
      IF (V_BC_ZMIN.EQ.0) THEN
C Start with zero
         DO J=0,NY+1
           DO I=0,NXP-1
             CU2(I,0,J)=0.d0
             CU2(I,1,J)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,1,J)=V_BC_ZMIN_C1
           CU2(0,0,J)=V_BC_ZMIN_C1
         END IF
C Neumann
      ELSE IF (V_BC_ZMIN.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU2(I,1,J)=CU2(I,2,J)
             CU2(I,0,J)=CU2(I,1,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,1,J)=CU2(0,2,J)-DZ(2)*V_BC_ZMIN_C1
           CU2(0,0,J)=CU2(0,1,J)-DZ(1)*V_BC_ZMIN_C1
         END IF
      ELSE
         STOP 'Error: V_BC_ZMIN must be 0, or 1'
      END IF


C Dirichlet
      IF (W_BC_ZMIN.EQ.0) THEN
C Set the vertical velocity at GZF(1) (halfway between GZ(1) and GZ(2))
         DO J=0,NY+1
           DO I=0,NXP-1
             CU3(I,1,J)=2.d0*W_BC_ZMIN_C1-CU3(I,2,J)
             CU3(I,0,J)=CU2(I,1,J)
           END DO
         END DO
C Neumann
      ELSE IF (W_BC_ZMIN.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU3(I,1,J)=CU3(I,2,J)
             CU3(I,0,J)=CU3(I,1,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,1,J)=CU3(0,2,J)-DZF(1)*W_BC_ZMIN_C1
           CU3(0,0,J)=CU3(0,1,J)-DZF(1)*W_BC_ZMIN_C1
         END IF
      ELSE
         STOP 'Error: W_BC_ZMIN must be 0, or 1'
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_VEL_ZUPPER_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,J

C Dirichlet
      IF (U_BC_ZMAX.EQ.0) THEN
C Start with zero
         DO J=0,NY+1
           DO I=0,NXP-1
             CU1(I,NZ,J)=0.d0
             CU1(I,NZ+1,J)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,NZ,J)=U_BC_ZMAX_C1
           CU1(0,NZ+1,J)=U_BC_ZMAX_C1
         END IF
C Neumann
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU1(I,NZ,J)=CU1(I,NZ-1,J)
             CU1(I,NZ+1,J)=CU1(I,NZ,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,NZ,J)=CU1(0,NZ-1,J)+DZ(NZ)*U_BC_ZMAX_C1
           CU1(0,NZ+1,J)=CU1(0,NZ,J)+DZ(NZ)*U_BC_ZMAX_C1
         END IF
      ELSE
         STOP 'Error: U_BC_ZMAX must be 0, or 1'
      END IF


C Dirichlet
      IF (V_BC_ZMAX.EQ.0) THEN
C Start with zero
         DO J=0,NY+1
           DO I=0,NXP-1
             CU2(I,NZ,J)=0.d0
             CU2(I,NZ+1,J)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,NZ,J)=V_BC_ZMAX_C1
           CU2(0,NZ+1,J)=V_BC_ZMAX_C1
         END IF
C Neumann
      ELSE IF (V_BC_ZMAX.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU2(I,NZ,J)=CU2(I,NZ-1,J)
             CU2(I,NZ+1,J)=CU2(I,NZ,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,NZ,J)=CU2(0,NZ-1,J)+DZ(NZ)*V_BC_ZMAX_C1
           CU2(0,NZ+1,J)=CU2(0,NZ,J)+DZ(NZ)*V_BC_ZMAX_C1
         END IF
      ELSE
         STOP 'Error: V_BC_ZMAX must be 0, or 1'
      END IF


C Dirichlet
      IF (W_BC_ZMAX.EQ.0) THEN
C Set the vertical velocity at GZF(1) (halfway between GZ(1) and GZ(2))
         DO J=0,NY+1
           DO I=0,NXP-1
             CU3(I,NZ+1,J)=2.d0*W_BC_ZMAX_C1-CU3(I,NZ,J)
           END DO
         END DO
C Neumann
      ELSE IF (W_BC_ZMAX.EQ.1) THEN
         DO J=0,NY+1
           DO I=0,NXP-1
             CU3(I,NZ,J)=CU3(I,NZ-1,J)
             CU3(I,NZ+1,J)=CU3(I,NZ,J)
           END DO
         END DO
C Now, apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,NZ+1,J)=CU3(0,NZ,J)+DZ(NZ)*W_BC_ZMAX_C1
         END IF
      ELSE
         STOP 'Error: W_BC_ZMAX must be 0, or 1'
      END IF

      RETURN
      END