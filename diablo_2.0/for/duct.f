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
C This algorithm uses Crank-Nicolson for viscous/diffusive terms and
C and 3rd order Runge-Kutta for the rest of the terms
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N,ISTART
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK
      REAL*8 HAVG_TH
      CHARACTER*55 FNAME
      REAL*8 TIME1, TIME2
      REAL*8 TIME3, TIME4
      REAL*8 TIME5, TIME6


      ! CALL WALL_TIME(TIME1)

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.d0
      TEMP2=H_BAR(RK_STEP) / 2.d0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

C First, we will compute the explicit RHS terms and store in Ri
C Note, Momentum equation and hence the RHS is evaluated at the
C corresponding velocity points.

C Store the old velocity in the RHS vector
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CR1(I,K,J)=CU1(I,K,J)
            CR2(I,K,J)=CU2(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO


      ! WRITE(*,*) 'writing CU CHECK1_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CU1_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU1)
      ! CLOSE(15)
      ! FNAME='CU2_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU2)
      ! CLOSE(15)
      ! FNAME='CU3_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU3)
      ! CLOSE(15)
      ! FNAME='CTH_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CTH)
      ! CLOSE(15)
      ! FNAME='CP_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CP)
      ! CLOSE(15)

      ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ! IF ((RK_STEP.eq.1).and.(TIME_STEP.eq.2)) STOP

C Add the R-K term from the rk-1 step
      IF (RK_STEP .GT. 1) THEN
        DO J=1,NY
          DO K=1,NZ
            DO I=0,NXP-1
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
      END IF


      ! WRITE(*,*) 'writing CR CHECK1_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CR1_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)
      ! FNAME='CR2_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR2)
      ! CLOSE(15)
      ! FNAME='CR3_CHK1_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR3)
      ! CLOSE(15)

C Add the pressure gradient to the RHS for X as explicit Euler
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
          END DO
        END DO
      END DO

C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO
C Add the pressure gradient to the RHS for Y as explicit Euler
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

C Take the z-derivative of the pressure at GZ points in Fourier space
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K-1,J)) / DZ(K)
          END DO
        END DO
      END DO
C Add the pressure gradient to the RHS for Z as explicit Euler
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO


      ! WRITE(*,*) 'writing CR CHECK2_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CR1_CHK2_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)
      ! FNAME='CR2_CHK2_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR2)
      ! CLOSE(15)
      ! FNAME='CR3_CHK2_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR3)
      ! CLOSE(15)

C Here, add the constant, forcing pressure gradient
C There are two built-in ways of doing this
C F_TYPE=1 -> Constant pressure gradient in the x-direction
C F_TYPE=2 -> Oscillatory pressure gradient in the x-direction
C OTHER VALUES -> No forcing added
      IF (F_TYPE.EQ.1) THEN
C Add forcing for a constant pressure gradient
        DO J=2,NY
          DO K=2,NZ
            IF (RANKZ.eq.0) CR1(0,K,J)=CR1(0,K,J)-TEMP4*PX0
          END DO
        END DO
      ELSE IF (F_TYPE.EQ.2) THEN
C Oscillatory pressure gradient
        DO J=2,NY
          DO K=2,NZ
            IF (RANKZ.eq.0) CR1(0,K,J)=CR1(0,K,J)-
     &          TEMP4*(PX0+AMP_OMEGA0*cos(OMEGA0*TIME))
          END DO
        END DO
C End if forcing type
      END IF

C Now compute the term R-K term Ai
C Compile terms of Ai in CFi which will be saved for next time step
C First, store the x-derivative viscous terms in CFi
C Note that BETA is the order of the Laplacian operator
C e.g. BETA=1 for second order viscosity/diffusivity
C or   BETA=2, etc. for hyper viscosity/diffusivity
      IF (BETA.NE.1) THEN
        IF (RANK.EQ.0) write(*,*) 
     &    'Hyperviscosity (BETA>1) not supported in duct flow'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
        STOP
      END IF
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=-NU * (KX2(I))**BETA * CU1(I,K,J)
            CF2(I,K,J)=-NU * (KX2(I))**BETA * CU2(I,K,J)
            CF3(I,K,J)=-NU * (KX2(I))**BETA * CU3(I,K,J)
          END DO
        END DO
      END DO

! Add the terms owing to the system rotation (Coriolis terms)
! Assume that the flow is on an f-plane
      DO K=2,NZ
        DO I=0,NXP-1
          DO J=2,NY
            CF1(I,K,J)=CF1(I,K,J)+I_RO*
     &      ((CU3(I,K,J)*DZ(K-1)+CU3(I,K-1,J)*DZ(K))/(2.d0*DZF(K)))
            CF3(I,K,J)=CF3(I,K,J)-I_RO*
     &      ((CU1(I,K,J)*DZF(K-1)+CU1(I,K-1,J)*DZF(K))/(2.d0*DZ(K)))
          END DO
        END DO
      END DO

! Do for each scalar
      DO N=1,N_TH

C !!REVISIT THIS!! do I include the 0 mode in X?
C its not a plane average - but still average in X
! If a scalar contributes to the density, RI is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode in the
! y-direction, which corresponds to the plane-average.
! The plane averaged density balances the hydrostatic pressure
!       IF (RANKZ.eq.0) THEN
!         ISTART=1
!       ELSE
!         ISTART=0
!       END IF
!       DO J=2,NY
!         DO K=2,NZ
!           DO I=ISTART,NXP-1
! ! Use second order interpolation
!             CF2(I,K,J)=CF2(I,K,J)+RI(N)*GRAV_Y*
!      &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
!           END DO
!         END DO
!       END DO
!       DO J=2,NY
!         DO K=2,NZ
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)+RI(N)*GRAV_X*CTH(I,K,J,N)
!             CF3(I,K,J)=CF3(I,K,J)+RI(N)*GRAV_Z*CTH(I,K,J,N)
!           END DO
!         END DO
!       END DO

! alternate way of doing the buoyancy term
! use second order interpolation to get contribution of CTH from GYF and GZF
! onto GY (for U2 contribution via CF2) and GZ (for U3 contribution via CF3)
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J)+RI(N)*GRAV_X*CTH(I,K,J,N)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J)+RI(N)*GRAV_Y*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
      END DO
      DO J=1,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J)+RI(N)*GRAV_Z*
     &      (CTH(I,K,J,N)*DZF(K-1)+CTH(I,K-1,J,N)*DZF(K))/(2.d0*DZ(K))
          END DO
        END DO
      END DO

! calculate horizontal average of buoyancy contribution and subtract off
! interpolate CTH from GYF and GZF onto GY and GZ, then average in Z
      IF (RANKZ.EQ.0) THEN
        I=0
        DO J=2,NY
          HAVG_TH=0.d0
          DO K=2,NZ
            HAVG_TH=HAVG_TH+ RI(N)*GRAV_Y*DZ(K)* (
     &         ( (CTH(I,K  ,J,N)*DYF(J-1)+CTH(I,K  ,J-1,N)*DYF(J))
     &            /(2.d0*DY(J)) ) * DZF(K-1)
     &       + ( (CTH(I,K-1,J,N)*DYF(J-1)+CTH(I,K-1,J-1,N)*DYF(J))
     &            /(2.d0*DY(J)) ) * DZF(K)  ) / (2.d0*DZ(K))
          END DO
          HAVG_TH=HAVG_TH/LZ
          DO K=1,NZ
            CF2(I,K,J)=CF2(I,K,J)-HAVG_TH
          END DO
        END DO
      END IF



! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTH
! The k-1 term for the R-K stepping is saved in FTH, CFTH

! First, build the RHS vector, use CRTH
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CRTH(I,K,J,N)=CTH(I,K,J,N)
          END DO
        END DO
      END DO
! Add term from k-2 step to free up CFTH variable
      IF (RK_STEP .GT. 1) THEN
        DO J=2,NY
          DO K=2,NZ
            DO I=0,NXP-1
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
      END IF

! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=-(NU/PR(N))*(KX2(I))**BETA*CTH(I,K,J,N)
          END DO
        END DO
      END DO

! Add advection acting on the background scalar gradient if present
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CU1(I,K,J)*DRHODX(N)
     &      -((CU3(I,K,J)*DZ(K-1)+CU3(I,K-1,J)*DZ(K))/(2.d0*DZF(K)))
     &           *DRHODZ(N)
          END DO
        END DO
      END DO

C End do number of passive scalars (N_TH)
      END DO

C Optionally, add user forcing to the right hand side
C Here, we have U1, U2, U3, and TH in Fourier space
      CALL USER_RHS_DUCT_FOURIER

C If we are considering an LES, then add the subgrid scale stress:
C Here, velocity and CFi should be in Fourier space
C The subgrid scale stress is added to CFi:   CFi=CFi - d/dx_i tau_ij

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.100))) THEN
C If we have created new flow with random perturbations, wait for a
C spinup before applying the subgrid model for stability purposes
C In the process, Ui is converted to physical space
        IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
        STOP
        CALL LES_DUCT


      ELSE
C If the subgrid model hasn't been called, then it is necessary to
C convert to physical space.
        CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1)
        CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1)
        CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1)

! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
        DO N=1,N_TH
          CS1(:,:,:)=CTH(:,:,:,N)
          CALL FFT_X_TO_PHYSICAL(CS1,S1,0,NY+1)
          TH(:,:,:,N)=S1(:,:,:)
        END DO
      END IF


C Compute the nonlinear products in physical space, then transform
C back to Fourier space to compute the derivative.
C Here, we compute the horizontal derivatives of the nonlinear terms
C which will be treated with RKW3.
C Do terms one at a time to save on memory

C  d( U1*U1 ) /dX
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U2*U1 ) /dY
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(U2(I,K,J+1)
     &             *(U1(I,K,J+1)+U1(I,K,J))/(2.d0)
     &                -U2(I,K,J)
     &             *(U1(I,K,J)+U1(I,K,J-1))/(2.d0)  )/DYF(J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U3*U1 ) /dZ
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(U3(I,K+1,J)
     &             *(U1(I,K+1,J)+U1(I,K,J))/(2.d0)
     &                -U3(I,K,J)
     &             *(U1(I,K,J)+U1(I,K-1,J))/(2.d0)  )/DZFP(K)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U1*U2 ) /dX
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=((U1(I,K,J)*DYF(J)
     &                 +U1(I,K,J-1)*DYF(J-1))/(2.d0*DY(J)))
     &                 *U2(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U2*U2 ) /dY
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=
     &           (0.25d0*(U2(I,K,J+1)+U2(I,K,J))**2.d0
     &           -0.25d0*(U2(I,K,J)+U2(I,K,J-1))**2.d0)/DY(J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U3*U2 ) /dZ
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(
     &         (U3(I,K+1,J)*DYF(J)+U3(I,K+1,J-1)*DYF(J-1))/(2.d0*DY(J))
     &        *(U2(I,K+1,J)+U2(I,K,J))/(2.d0)
     &        -(U3(I,K  ,J)*DYF(J)+U3(I,K  ,J-1)*DYF(J-1))/(2.d0*DY(J))
     &        *(U2(I,K,J)+U2(I,K-1,J))/(2.d0)  )/DZFP(K)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U1*U3 ) /dX
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=((U1(I,K,J)*DZFP(K)
     &                 +U1(I,K-1,J)*DZFP(K-1))/(2.d0*DZP(K)))
     &                 *U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U2*U3 ) /dY
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(
     &       (U2(I,K,J+1)*DZFP(K)+U2(I,K-1,J+1)*DZFP(K-1))/(2.d0*DZP(K))
     &      *(U3(I,K,J+1)+U3(I,K,J))/(2.d0)
     &      -(U2(I,K,J  )*DZFP(K)+U2(I,K-1,J  )*DZFP(K-1))/(2.d0*DZP(K))
     &      *(U3(I,K,J)+U3(I,K,J-1))/(2.d0)  )/DYF(J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U3*U3 ) /dZ
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=
     &           (0.25d0*(U3(I,K+1,J)+U3(I,K,J))**2.d0
     &           -0.25d0*(U3(I,K,J)+U3(I,K-1,J))**2.d0)/DZP(K)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO


C -- At this point, we are done computing the nonlinear terms --

C Optionally, add user forcing to the right hand side
C Here, we have U1, U2, U3, and TH in physical space
      CALL USER_RHS_DUCT_PHYSICAL


      ! WRITE(*,*) 'writing CR CHECK3_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CR1_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)
      ! FNAME='CR2_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR2)
      ! CLOSE(15)
      ! FNAME='CR3_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR3)
      ! CLOSE(15)

      ! WRITE(*,*) 'writing CF CHECK3_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CF1_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CF1)
      ! CLOSE(15)
      ! FNAME='CF2_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CF2)
      ! CLOSE(15)
      ! FNAME='CF3_CHK3_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CF3)
      ! CLOSE(15)


C Finally, Add CFi to CRi
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CR1(I,K,J)=CR1(I,K,J) + TEMP5 * CF1(I,K,J)
            CR2(I,K,J)=CR2(I,K,J) + TEMP5 * CF2(I,K,J)
            CR3(I,K,J)=CR3(I,K,J) + TEMP5 * CF3(I,K,J)
          END DO
        END DO
      END DO

      ! WRITE(*,*) 'writing CR CHECK4_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CR1_CHK4_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)
      ! FNAME='CR2_CHK4_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR2)
      ! CLOSE(15)
      ! FNAME='CR3_CHK4_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR3)
      ! CLOSE(15)

C Convert RHS terms to physical space
      CALL FFT_X_TO_PHYSICAL(CR1,R1,0,NY+1)
      CALL FFT_X_TO_PHYSICAL(CR2,R2,0,NY+1)
      CALL FFT_X_TO_PHYSICAL(CR3,R3,0,NY+1)


      ! WRITE(*,*) 'writing CHECK5_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='R1_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R1)
      ! CLOSE(15)
      ! FNAME='R2_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R2)
      ! CLOSE(15)
      ! FNAME='R3_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R3)
      ! CLOSE(15)
      ! FNAME='RTH_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(RTH)
      ! CLOSE(15)

C Compute the y- and z-derivative viscous terms in physical space
C and add to RHS
C This is the explicit part of the Crank-Nicolson term
C NU_V_SCALE is a coefficient defined as the ratio of the vertical
C to horizontal viscosity and can be used to add anisotropic viscosity
      DO J=1,NY
        DO K=1,NZP
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)
     &           -(U1(I,K,J)   - U1(I,K,J-1)) / DY(J)) /DYF(J)  )
     &                         +TEMP1*
     &        (  ((U1(I,K+1,J) - U1(I,K,J)) / DZP(K+1)
     &           -(U1(I,K,J)   - U1(I,K-1,J)) / DZP(K)) /DZFP(K)  )
            R2(I,K,J)=R2(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U2(I,K,J+1) - U2(I,K,J))  / DYF(J)
     &           -(U2(I,K,J)   - U2(I,K,J-1))/ DYF(J-1)) /DY(J)  )
     &                         +TEMP1*
     &        (  ((U2(I,K+1,J) - U2(I,K,J)) / DZP(K+1)
     &           -(U2(I,K,J)   - U2(I,K-1,J)) / DZP(K)) /DZFP(K)  )
            R3(I,K,J)=R3(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U3(I,K,J+1) - U3(I,K,J)) / DY(J+1)
     &           -(U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) /DYF(J)  )
     &                         +TEMP1*
     &        (  ((U3(I,K+1,J) - U3(I,K,J)) / DZFP(K+1)
     &           -(U3(I,K,J)   - U3(I,K-1,J)) / DZFP(K)) /DZP(K)  )
          END DO
        END DO
      END DO



      ! WRITE(*,*) 'writing CHECK6_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='R1_CHK6_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R1)
      ! CLOSE(15)
      ! FNAME='R2_CHK6_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R2)
      ! CLOSE(15)
      ! FNAME='R3_CHK6_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R3)
      ! CLOSE(15)
      ! FNAME='RTH_CHK6_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(RTH)
      ! CLOSE(15)

C If we are using a subgrid model, add the eddy viscosity term
C This is an added viscosity that will be treated just like the
C molecular viscosity with Crank-Nicolson for the vertical derivatives
      IF (LES) then
        IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
        STOP
      END IF

C -- Here, we are done with computation of Velocity RHS, explicit terms --

C Now, build the explicit RHS terms for the passive scalar(s)

      DO N=1,N_TH
! Do for each scalar:

! Compute the nonlinear terms that are present in the explicit term A
C  d( U1*TH ) /dX
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=U1(I,K,J)*TH(I,K,J,N)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U2*TH ) /dY
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(U2(I,K,J+1)
     &             *(TH(I,K,J+1,N)+TH(I,K,J,N))/(2.d0)
     &                -U2(I,K,J)
     &             *(TH(I,K,J,N)+TH(I,K,J-1,N))/(2.d0)  )/DYF(J)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
          END DO
        END DO
      END DO

C  d( U3*TH ) /dZ
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            S1(I,K,J)=(U3(I,K+1,J)
     &             *(TH(I,K+1,J,N)+TH(I,K,J,N))/(2.d0)
     &                -U3(I,K,J)
     &             *(TH(I,K,J,N)+TH(I,K-1,J,N))/(2.d0)  )/DZFP(K)
          END DO
        END DO
      END DO
      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
          END DO
        END DO
      END DO
      

! Add CFTH to the RHS vector CRTH
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CRTH(I,K,J,N)=CRTH(I,K,J,N) + TEMP5 * CFTH(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space

      CS1(:,:,:)=CRTH(:,:,:,N)
      CALL FFT_X_TO_PHYSICAL(CS1,S1,0,NY+1)
      RTH(:,:,:,N)=S1(:,:,:)

! Compute the Explicit part of the Crank-Nicolson terms for the TH equation
! First, the y- and z-derivative viscous terms
      DO J=2,NY
        DO K=2,NZP
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+(TEMP1/PR(N))*NU_V_SCALE*
     &        (  ((TH(I,K,J+1,N) - TH(I,K,J,N)) / DY(J+1)
     &           -(TH(I,K,J,N) - TH(I,K,J-1,N)) / DY(J)) /DYF(J)  )
     &                         +(TEMP1/PR(N))*
     &        (  ((TH(I,K+1,J,N) - TH(I,K,J,N)) / DZP(K+1)
     &           -(TH(I,K,J,N) - TH(I,K-1,J,N)) / DZP(K)) /DZFP(K)  )
          END DO
        END DO
      END DO
! If we are using a subgrid model (LES) then add the eddy diffusivity here
! Note, KAPPA_T is defined at GY points
      IF (LES) THEN
        IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
        STOP
      END IF




      ! WRITE(*,*) 'writing CHECK7 on rank ', RANK  
      ! FNAME='R1_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R1)
      ! CLOSE(15)
      ! FNAME='R2_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R2)
      ! CLOSE(15)
      ! FNAME='R3_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(R3)
      ! CLOSE(15)
      ! FNAME='RTH_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(RTH)
      ! CLOSE(15)

      ! WRITE(*,*) 'writing CHECK7 on rank ', RANK  
      ! FNAME='U1_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U1)
      ! CLOSE(15)
      ! FNAME='U2_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U2)
      ! CLOSE(15)
      ! FNAME='U3_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U3)
      ! CLOSE(15)
      ! FNAME='TH_CHECK7_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(TH)
      ! CLOSE(15)


      ! WRITE(*,*) 'reached start of TH advancement RANK=', RANK  

      ! CALL WALL_TIME(TIME5)

      ! CALL ADI_YZ_TH_DUCT(N)
      IF (MOD(TIME_STEP,2).EQ.0) THEN
        CALL ADI_ZY_TH_DUCT(N)
      ELSE
        CALL ADI_YZ_TH_DUCT(N)
      END IF

! End do number of passive scalars
        END DO


      ! CALL ADI_YZ_U2_DUCT
      IF (MOD(TIME_STEP,2).EQ.0) THEN
        CALL ADI_ZY_U2_DUCT
      ELSE
        CALL ADI_YZ_U2_DUCT
      END IF

      ! CALL ADI_YZ_U3_DUCT
      IF (MOD(TIME_STEP,2).EQ.0) THEN
        CALL ADI_ZY_U3_DUCT
      ELSE
        CALL ADI_YZ_U3_DUCT
      END IF

      ! CALL ADI_YZ_U1_DUCT
      IF (MOD(TIME_STEP,2).EQ.0) THEN
        CALL ADI_ZY_U1_DUCT
      ELSE
        CALL ADI_YZ_U1_DUCT
      END IF


    !   CALL WALL_TIME(TIME6)
    !   WRITE(*,*) 'RK', RK_STEP, ' on RANK', RANK,
    !  &     ' ADI: ', TIME6-TIME5, ' seconds'


      ! WRITE(*,*) 'writing CHECK9_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='U1_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U1)
      ! CLOSE(15)
      ! FNAME='U2_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U2)
      ! CLOSE(15)
      ! FNAME='U3_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(U3)
      ! CLOSE(15)
      ! FNAME='TH_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(TH)
      ! CLOSE(15)

      ! WRITE(*,*) 'finished updating all U on RANK=', RANK  
      ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ! IF (RK_STEP.eq.3) STOP
      ! STOP

C -- Done getting U1hat, U2hat, U3hat at new RK Step --

C If we are on the final RK step, optionally update the timestep
C based on the CFL criteria. This won't affect the current timestep
C since the TEMP1, etc. variables have already been set using
C the current timestep
      IF (VARIABLE_DT.and.(RK_STEP.eq.3)
     &      .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
        CALL COURANT
      END IF

! Transform TH and U to Fourier Space
      CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1)
      DO N=1,N_TH
        S1(:,:,:)=TH(:,:,:,N)
        CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1)
        CTH(:,:,:,N)=CS1(:,:,:)
      END DO

      ! Fix disparities at the boundary due to the thoms algorithm in parallel
      IF (USE_MPI) THEN
         CALL GHOST_DUCT_MPI
      END IF


C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space


      ! WRITE(*,*) 'writing CU CHECK5_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CU1_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU1)
      ! CLOSE(15)
      ! FNAME='CU2_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU2)
      ! CLOSE(15)
      ! FNAME='CU3_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU3)
      ! CLOSE(15)
      ! FNAME='CTH_CHK5_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CTH)
      ! CLOSE(15)

    !   CALL WALL_TIME(TIME3)
      CALL REM_DIV_DUCT6
    !   CALL WALL_TIME(TIME4)
    !   WRITE(*,*) 'RK', RK_STEP, ' on RANK', RANK,
    !  &     ' REM_DIV: ', TIME4-TIME3, ' seconds'

C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      ! WRITE(*,*) 'writing CU CHECK9_RK', RK_STEP, ' on rank ', RANK  
      ! FNAME='CU1_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU1)
      ! CLOSE(15)
      ! FNAME='CU2_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU2)
      ! CLOSE(15)
      ! FNAME='CU3_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU3)
      ! CLOSE(15)
      ! FNAME='CTH_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CTH)
      ! CLOSE(15)

      ! FNAME='CR1_CHK9_RK'//char(RK_STEP+48)//'_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NXP-1
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

      ! Fix disparities at the boundary due to the thoms algorithm in parallel
      IF (USE_MPI) THEN
         CALL GHOST_DUCT_MPI
      END IF


    !   CALL WALL_TIME(TIME2)
    !   WRITE(*,*) 'RK', RK_STEP, ' on RANK', RANK,
    !  &     ' total RKSTEP: ', TIME2-TIME1, ' seconds'


      ! call mpi_barrier(MPI_COMM_WORLD,ierror)
      ! WRITE(*,*) 'Ending RKSTEP ',RK_STEP, ' on RANK=', RANK

      ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
      ! IF ((RK_STEP.eq.1).and.(TIME_STEP.eq.2)) STOP
      ! IF (RK_STEP.eq.3) STOP

! 151   format(F40.20)  
151   format(E40.30)    
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the duct-flow case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      IF (RANK.EQ.0)
     &     WRITE(*,*) ' RK_DUCT_2 not supported yet'
      CALL MPI_FINALIZE(ierror)
      STOP

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
      INTEGER NI
      CHARACTER*55 FNAME
      INTEGER JLIM1,JLIM2

      JLIM1=1
      JLIM2=NY
      IF (RANKY.EQ.0) THEN
        JLIM1=2
      ELSE IF (RANKY.EQ.NPROCY-1) THEN
        JLIM2=NY-1
      END IF

! C Alternating Direction Implicit method, first we implicitly solve the
! C system with Y derivatives on the LHS and Z derivatives on the RHS
    
! C First, Initialize the matrix components
!       DO J=0,NY+1
!         DO K=0,NZ+1
!             MATL_C_YZ(J,K)=0.
!             MATD_C_YZ(J,K)=1.
!             MATU_C_YZ(J,K)=0.
!             VEC_C_YZ(J,K)=(0.,0.)
!         END DO
!       END DO

! C The 2d FFT of Ui should have been taken and stored in CUi
! C Solving for phi amounts to solving a tridiagonal system
! C First, construct the system to be solved
!       DO I=0,NXP-1
!         DO J=1,NY
!           DO K=1,NZ
!             MATL_C_YZ(J,K)=1./(DY(J)*DYF(J))
!             MATD_C_YZ(J,K)=-KX2(I)
!      &          -1./(DY(J+1)*DYF(J)) - 1./(DY(J)*DYF(J))
!             MATU_C_YZ(J,K)=1./(DY(J+1)*DYF(J))
!           END DO
!         END DO

! C Now, create the RHS vector
!         DO J=1,NY
!           DO K=1,NZ
!             VEC_C_YZ(J,K)=CIKX(I)*CU1(I,K,J)
!      &           + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
!      &           + (CU3(I,K+1,J)-CU3(I,K,J))/DZF(K)
!     !  &           - 1./(DZ(K)*DZF(K)) CU
!     !  &           - 1./(DZ(K+1)*DZF(K))
!     !  &           + 1./(DZ(K+1)*DZF(K)) + 1./(DZ(K)*DZF(K))
!           END DO
!         END DO


!         ! IF (I.EQ.0) THEN
!       ! WRITE(*,*) 'writing MATL_C_YZ on rank ', RANK  
!       ! FNAME='MATL_C_YZ_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(MATL_C_YZ)
!       ! CLOSE(15)

!       ! WRITE(*,*) 'writing MATD_C_YZ on rank ', RANK  
!       ! FNAME='MATD_C_YZ_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(MATD_C_YZ)
!       ! CLOSE(15)

!       ! WRITE(*,*) 'writing MATU_C_YZ on rank ', RANK  
!       ! FNAME='MATU_C_YZ_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(MATU_C_YZ)
!       ! CLOSE(15)

!       ! WRITE(*,*) 'writing VEC_C_YZ_CHECK1 on rank ', RANK  
!       ! FNAME='VEC_C_YZ_CHECK1_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(VEC_C_YZ)
!       ! CLOSE(15)
!       !   END IF


! C If we are using the MPI package...
!         IF (USE_MPI) THEN
!           CALL APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
!      &         MATU_C_YZ,VEC_C_YZ,I)

!       !   IF (I.EQ.0) THEN
!       ! WRITE(*,*) 'writing VEC_C_YZ_CHECK2 on rank ', RANK  
!       ! FNAME='VEC_C_YZ_CHECK2_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(VEC_C_YZ)
!       ! CLOSE(15)
!       !   END IF

!           CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
!      &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
!           CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
!      &         MATU_C_YZ,VEC_C_YZ,NY,NZ)

! C Else we are running in serial mode
!         ELSE
!           IF (RANK.EQ.0) 
!      &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
!           STOP
!         END IF

!         DO J=1,NY
!           DO K=1,NZ
!             CR1(I,K,J)=VEC_C_YZ(J,K)
!           END DO
!         END DO

!       !   IF (I.EQ.0) THEN
!       ! WRITE(*,*) 'writing VEC_C_YZ_CHECK3 on rank ', RANK  
!       ! FNAME='VEC_C_YZ_CHECK3_'//char(RANK+48)//'.dat'
!       ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       ! WRITE(15,151) real(VEC_C_YZ)
!       ! CLOSE(15)
!       !   END IF

!       END DO

! C End ADI Y sweep

      ! WRITE(*,*) 'writing CR1 on rank ', RANK  
      ! FNAME='CR1_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)


C Alternating Direction Implicit method, now we implicitly solve the
C system with Z derivatives on the LHS and Y derivatives on the RHS

C First, Initialize the matrix components
      DO J=0,NY+1
        DO K=0,NZ+1
            MATL_C_YZ(J,K)=0.
            MATD_C_YZ(J,K)=1.
            MATU_C_YZ(J,K)=0.
            VEC_C_YZ(J,K)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO I=0,NXP-1
        DO J=1,NY
          DO K=1,NZ
            MATL_C_YZ(J,K)=1./(DZ(K)*DZF(K))
            MATD_C_YZ(J,K)=-KX2(I)
     &          -1./(DZ(K+1)*DZF(K)) - 1./(DZ(K)*DZF(K))
            MATU_C_YZ(J,K)=1./(DZ(K+1)*DZF(K))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY
          DO K=1,NZ
            VEC_C_YZ(J,K)=CIKX(I)*CU1(I,K,J)
     &           + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
     &           + (CU3(I,K+1,J)-CU3(I,K,J))/DZF(K)
    !  &           - CR1(I,K,J-1)/(DY(J)*DYF(J))
    !  &           - CR1(I,K,J+1)/(DY(J+1)*DYF(J))
    !  &           + CR1(I,K,J)/(DY(J)*DYF(J))
    !  &           + CR1(I,K,J)/(DY(J+1)*DYF(J))
          END DO
        END DO

      !   IF (I.EQ.0) THEN
      ! WRITE(*,*) 'writing MATL_C_YZ on rank ', RANK  
      ! FNAME='MATL_C_YZ_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(MATL_C_YZ)
      ! CLOSE(15)

      ! WRITE(*,*) 'writing MATD_C_YZ on rank ', RANK  
      ! FNAME='MATD_C_YZ_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(MATD_C_YZ)
      ! CLOSE(15)

      ! WRITE(*,*) 'writing MATU_C_YZ on rank ', RANK  
      ! FNAME='MATU_C_YZ_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(MATU_C_YZ)
      ! CLOSE(15)

      ! WRITE(*,*) 'writing VEC_C_YZ on rank ', RANK  
      ! FNAME='VEC_C_YZ_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(VEC_C_YZ)
      ! CLOSE(15)
      !   END IF
        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK4 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK4_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing MATL_C_YZ on rank ', RANK  
      FNAME='MATL_C_YZ_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(MATL_C_YZ)
      CLOSE(15)

      WRITE(*,*) 'writing MATD_C_YZ on rank ', RANK  
      FNAME='MATD_C_YZ_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(MATD_C_YZ)
      CLOSE(15)

      WRITE(*,*) 'writing MATU_C_YZ on rank ', RANK  
      FNAME='MATU_C_YZ_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(MATU_C_YZ)
      CLOSE(15)

      ! WRITE(*,*) 'writing VEC_C_YZ on rank ', RANK  
      ! FNAME='VEC_C_YZ_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(VEC_C_YZ)
      ! CLOSE(15)
        END IF
        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK5 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK5_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

          CALL TDMA_F_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK6 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK6_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF
C Else we are running in serial mode
        ELSE
          IF (RANK.EQ.0) 
     &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

        ! DO J=1,NY
        DO J=JLIM1,JLIM2
          DO K=1,NZ
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO
      END DO


C End ADI Z sweep

      ! WRITE(*,*) 'writing CU2 on rank ', RANK  
      ! FNAME='CU2_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CU2)
      ! CLOSE(15)

      WRITE(*,*) 'writing CR1_CHECK1 on rank ', RANK  
      FNAME='CR1_CHECK1_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)

      CALL APPLY_BC_CR1_CHI

      WRITE(*,*) 'writing CR1_CHECK2 on rank ', RANK  
      FNAME='CR1_CHECK2_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)


151   format(E40.30)



ction Implicit method, first we implicitly solve the
C system with Y derivatives on the LHS and Z derivatives on the RHS
    
C First, Initialize the matrix components
      DO J=0,NY+1
        DO K=0,NZ+1
            MATL_C_YZ(J,K)=0.
            MATD_C_YZ(J,K)=1.
            MATU_C_YZ(J,K)=0.
            VEC_C_YZ(J,K)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO I=0,NXP-1
        DO J=1,NY
          DO K=1,NZ
            MATL_C_YZ(J,K)=1./(DY(J)*DYF(J))
            MATD_C_YZ(J,K)=-KX2(I)
     &          -1./(DY(J+1)*DYF(J)) - 1./(DY(J)*DYF(J))
            MATU_C_YZ(J,K)=1./(DY(J+1)*DYF(J))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY
          DO K=1,NZ
            VEC_C_YZ(J,K)=CIKX(I)*CU1(I,K,J)
     &           + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
     &           + (CU3(I,K+1,J)-CU3(I,K,J))/DZF(K)
     &           - CR1(I,K-1,J)/(DZ(K)*DZF(K))
     &           - CR1(I,K+1,J)/(DZ(K+1)*DZF(K))
     &           + CR1(I,K,J)/(DZ(K)*DZF(K))
     &           + CR1(I,K,J)/(DZ(K+1)*DZF(K))
          END DO
        END DO

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK7 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK7_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK8 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK8_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

          CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK9 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK9_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

C Else we are running in serial mode
        ELSE
          IF (RANK.EQ.0) 
     &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

        DO J=1,NY
          ! DO K=1,NZ
          DO K=2,NZ-1
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO

      END DO

C End ADI Y sweep

      WRITE(*,*) 'writing CR1_CHECK3 on rank ', RANK  
      FNAME='CR1_CHECK3_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)

      CALL APPLY_BC_CR1_CHI

      WRITE(*,*) 'writing CR1_CHECK4 on rank ', RANK  
      FNAME='CR1_CHECK4_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)



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
      SUBROUTINE REM_DIV_DUCT4
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C
C Direct block-Thomas pentadiagonal solver for the duct divergence-removal
C Poisson problem. At each kx mode I, solves
C   d2(phi)/dy2 + d2(phi)/dz2 - kx^2 phi = div(u)
C on the (Y, Z) plane with Neumann walls (Y and Z) plus a Dirichlet anchor
C at the kx=0 corner. The discretization gives a block-tridiagonal-in-J
C system with NZ x NZ blocks; we solve directly via block-Thomas
C parallelized in Y, mirroring the existing TDMA_*_DUCT_Y_MPI pipeline
C structure but carrying NZ x NZ blocks and NZ vectors per neighbor exchange.
C
C Replaces REM_DIV_DUCT3's iterative ADI scheme with a direct solve.
C Stores phi in CR1, then corrects CU1, CU2, CU3 to be divergence-free.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      INTEGER I, J, K
      CHARACTER*55 FNAME

      DO I = 0, NXP-1

C Initialize all entries (including halos) to identity-block / zero so
C that any out-of-range access during BC application or block sweeps
C produces a benign default.
        DO J = 0, NY+1
          DO K = 0, NZ+1
            A_PENTA(J, K)    = 0.d0
            C_PENTA(J, K)    = 0.d0
            B_PENTA(J, K, 1) = 0.d0
            B_PENTA(J, K, 2) = 1.d0
            B_PENTA(J, K, 3) = 0.d0
            G_PENTA(J, K)    = (0.d0, 0.d0)
          END DO
        END DO

C Build the discrete operator and divergence RHS on the interior
C J=1..NY, K=1..NZ. BC rows (Y- and Z-walls) get overwritten in
C APPLY_BC_REM_DIV_DUCT4_MPI below.
        DO J = 1, NY
          DO K = 1, NZ
C Y-couplings (diagonal in K)
            A_PENTA(J, K) = 1.d0/(DY(J)  *DYF(J))
            C_PENTA(J, K) = 1.d0/(DY(J+1)*DYF(J))
C Z-couplings (tridiag in K within block J)
            B_PENTA(J, K, 1) = 1.d0/(DZ(K)  *DZF(K))
            B_PENTA(J, K, 3) = 1.d0/(DZ(K+1)*DZF(K))
C Main diagonal: -kx^2 - (Z curvatures) - (Y curvatures)
            B_PENTA(J, K, 2) = -KX2(I)
     &        - 1.d0/(DZ(K)  *DZF(K)) - 1.d0/(DZ(K+1)*DZF(K))
     &        - 1.d0/(DY(J)  *DYF(J)) - 1.d0/(DY(J+1)*DYF(J))
C RHS: divergence of intermediate velocity at (J, K)
            G_PENTA(J, K) = CIKX(I)*CU1(I,K,J)
     &        + (CU2(I,K,J+1) - CU2(I,K,J))/DYF(J)
     &        + (CU3(I,K+1,J) - CU3(I,K,J))/DZF(K)
          END DO
        END DO

        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT4_MPI(I)

C Diagnostic: dump G_PENTA after BCs, before forward sweep, on kx=0 mode
          ! IF (I.EQ.0) THEN
          !   WRITE(*,*) 'writing G_PENTA_PRE on rank ', RANK
          !   FNAME='G_PENTA_PRE_'//char(RANK+48)//'.dat'
          !   OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
          !   WRITE(15,151) real(G_PENTA)
          !   CLOSE(15)
          ! END IF

          CALL BLOCK_TDMA_F_COMP_DUCT_Y_MPI(NY, NZ)

C Diagnostic: dump G_PENTA after forward sweep (now holds G'_J)
          ! IF (I.EQ.0) THEN
          !   WRITE(*,*) 'writing G_PENTA_FWD on rank ', RANK
          !   FNAME='G_PENTA_FWD_'//char(RANK+48)//'.dat'
          !   OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
          !   WRITE(15,151) real(G_PENTA)
          !   CLOSE(15)
          ! END IF

          CALL BLOCK_TDMA_B_COMP_DUCT_Y_MPI(NY, NZ)

C Diagnostic: dump G_PENTA after backward sweep (now holds the solution)
          ! IF (I.EQ.0) THEN
          !   WRITE(*,*) 'writing G_PENTA_BWD on rank ', RANK
          !   FNAME='G_PENTA_BWD_'//char(RANK+48)//'.dat'
          !   OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
          !   WRITE(15,151) real(G_PENTA)
          !   CLOSE(15)
          ! END IF
        ELSE
          IF (RANK.EQ.0) WRITE(*,*)
     &      'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

C Copy solution into CR1
        DO J = 1, NY
          DO K = 1, NZ
            CR1(I, K, J) = G_PENTA(J, K)
          END DO
        END DO

      END DO

C Diagnostic: dump CR1 after the full kx loop
      ! WRITE(*,*) 'writing CR1_DUCT4 on rank ', RANK
      ! FNAME='CR1_DUCT4_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(CR1)
      ! CLOSE(15)

C ----------------------------------------------------------------------
C Velocity correction: subtract grad(phi) from the intermediate velocity
C (mirrors the post-solve correction in REM_DIV_DUCT3)
C ----------------------------------------------------------------------
      DO J = 1, NY
        DO K = 1, NZ
          DO I = 0, NXP-1
            CU1(I, K, J) = CU1(I, K, J) - CIKX(I)*CR1(I, K, J)
          END DO
        END DO
      END DO
      DO J = 2, NY
        DO K = 1, NZ
          DO I = 0, NXP-1
            CU2(I, K, J) = CU2(I, K, J)
     &        - (CR1(I, K, J) - CR1(I, K, J-1))/DY(J)
          END DO
        END DO
      END DO
      DO J = 1, NY
        DO K = 2, NZ
          DO I = 0, NXP-1
            CU3(I, K, J) = CU3(I, K, J)
     &        - (CR1(I, K, J) - CR1(I, K-1, J))/DZ(K)
          END DO
        END DO
      END DO

151   format(E40.30)
      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_DUCT6
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C DCT-I based direct solver for the duct divergence-removal Poisson problem.
C At each kx mode I, diagonalizes the uniform-Z endpoint-Neumann operator
C using DCT-I, then solves NZ scalar tridiagonal systems in Y.
C This is additive and does not modify REM_DIV_DUCT4 or its block solver.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      INTEGER I,J,K,M

      DO I=0,NXP-1

C Initialize scalar tridiagonal work arrays.
        DO J=0,NY+1
          DO K=0,NZ+1
            MATL_C_YZ(J,K)=0.d0
            MATD_C_YZ(J,K)=1.d0
            MATU_C_YZ(J,K)=0.d0
            VEC_C_YZ(J,K)=(0.d0,0.d0)
          END DO
        END DO

C Build divergence RHS on the physical Z endpoint grid.
        DO J=1,NY
          DO K=1,NZ
            VEC_C_YZ(J,K)=CIKX(I)*CU1(I,K,J)
     &        + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
     &        + (CU3(I,K+1,J)-CU3(I,K,J))/DZF(K)
          END DO
        END DO

C Transform each local Y row from physical Z to DCT-I Z modes.
        DO J=1,NY
          CALL DCT1_Z_FORWARD_COMPLEX(VEC_C_YZ,J)
        END DO

C Build modal Y tridiagonal systems.  The second index is now mode M.
        DO M=1,NZ
          DO J=1,NY
            MATL_C_YZ(J,M)=1.d0/(DY(J)*DYF(J))
            MATU_C_YZ(J,M)=1.d0/(DY(J+1)*DYF(J))
            MATD_C_YZ(J,M)=-KX2(I)+DCT1_Z_LAM(M)
     &        -1.d0/(DY(J)*DYF(J))
     &        -1.d0/(DY(J+1)*DYF(J))
          END DO
        END DO

        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT6_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
          CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
        ELSE
          IF (RANK.EQ.0) WRITE(*,*)
     &      'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

C Inverse DCT-I each Y row and copy phi into CR1.
        DO J=1,NY
          CALL DCT1_Z_BACKWARD_COMPLEX(VEC_C_YZ,J)
          DO K=1,NZ
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO

      END DO

C Velocity correction: subtract grad(phi) from the intermediate velocity.
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
            CU2(I,K,J)=CU2(I,K,J)
     &        -(CR1(I,K,J)-CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      DO J=1,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CU3(I,K,J)=CU3(I,K,J)
     &        -(CR1(I,K,J)-CR1(I,K-1,J))/DZ(K)
          END DO
        END DO
      END DO

      RETURN
      END








C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_YZ_U1_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U1
C does Y first, then Z
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U1 is cell-centered in Y (lives on GYF) -- standard Laplacian stencil
! U1 lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=R1(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U1_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U1(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U1 is cell-centered in Z (lives on GZP) -- standard Laplacian stencil
! U1 lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZP(K+1)*DZFP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=U1(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U1_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U1(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_ZY_U1_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U1
C does Z first, then Y
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U1 is cell-centered in Z (lives on GZP) -- standard Laplacian stencil
! U1 lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZP(K+1)*DZFP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=R1(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U1_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U1(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|

C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U1 is cell-centered in Y (lives on GYF) -- standard Laplacian stencil
! U1 lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=U1(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U1_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U1(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_YZ_U2_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U2
C does Y first, then Z
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U2 is face-centered in Y (lives on GY) -- staggered Laplacian stencil
! U2 lives on GY so we define 2:NY -> 1,NY+1 set by BCs (0 placeholder)
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DYF(J-1)*DY(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DYF(J)*DY(J))
     &           +(TEMP1*NU_V_SCALE) / (DYF(J-1)*DY(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DYF(J)*DY(J))
! Define RHS vector
            VEC_YZ(J,K)=R2(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U2_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U2(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U2 is cell-centered in Z (lives on GZP) -- standard Laplacian stencil
! U2 lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZP(K+1)*DZFP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=U2(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U2_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U2(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_ZY_U2_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U2
C does Z first, then Y
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U2 is cell-centered in Z (lives on GZP) -- standard Laplacian stencil
! U2 lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZP(K+1)*DZFP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=R2(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U2_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U2(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|

C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U2 is face-centered in Y (lives on GY) -- staggered Laplacian stencil
! U2 lives on GY so we define 2:NY -> 1,NY+1 set by BCs (0 placeholder)
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DYF(J-1)*DY(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DYF(J)*DY(J))
     &           +(TEMP1*NU_V_SCALE) / (DYF(J-1)*DY(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DYF(J)*DY(J))
! Define RHS vector
            VEC_YZ(J,K)=U2(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U2_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U2(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_YZ_U3_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U3
C does Y first, then Z
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U3 is cell-centered in Y (lives on GYF) -- standard Laplacian stencil
! U3 lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=R3(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U3_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U3(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U3 is face-centered in Z (lives on GZ) -- staggered Laplacian stencil
! U3 lives on GZ so we define 2:NZP -> 1,NZP+1 set by BCs (0 placeholder)
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZFP(K-1)*DZP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZFP(K)*DZP(K))
     &           +(TEMP1) / (DZFP(K-1)*DZP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZFP(K)*DZP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=U3(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U3_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U3(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_ZY_U3_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving U3
C does Z first, then Y
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U3 is face-centered in Z (lives on GZ) -- staggered Laplacian stencil
! U3 lives on GZ so we define 2:NZP -> 1,NZP+1 set by BCs (0 placeholder)
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1) / (DZFP(K-1)*DZP(K))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1) / (DZFP(K)*DZP(K))
     &           +(TEMP1) / (DZFP(K-1)*DZP(K))
            MATU_YZ(J,K)= -(TEMP1) / (DZFP(K)*DZP(K))
! Define RHS vector (D'Yakonov-style: feed Y-sweep result)
            VEC_YZ(J,K)=R3(I,K,J)
          END DO
        END DO

        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

          IF (USE_MPI) THEN
             CALL APPLY_BC_U3_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U3(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|

C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! U3 is cell-centered in Y (lives on GYF) -- standard Laplacian stencil
! U3 lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1.
     &           +(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=U3(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*)
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_U3_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            U3(I,K,J)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_YZ_TH_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving TH
C does Y first, then Z
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! Use quasi-second order interpolation for TH on GY points
! TH lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1. 
     &           +(TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=RTH(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_TH_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            TH(I,K,J,N)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! Use quasi-second order interpolation for TH on GZ points
! TH lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1/PR(N)) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1. 
     &           +(TEMP1/PR(N)) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1/PR(N)) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1/PR(N)) / (DZP(K+1)*DZFP(K))
! Define RHS vector
            VEC_YZ(J,K)=TH(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_TH_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,N)
      ! WRITE(*,*) 'finished apply BC Z RANK=', RANK , 'I=', I
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
      ! WRITE(*,*) 'finished TDMA F Z RANK=', RANK , 'I=', I
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
      ! WRITE(*,*) 'finished TDMA B Z RANK=', RANK , 'I=', I
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            TH(I,K,J,N)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ADI_ZY_TH_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternating Direction Implicit method for solving TH
C does Z first, then Y
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N
      REAL*8 TEMP1

      TEMP1=NU * H_BAR(RK_STEP) / 2.d0


C------------------------   BEGIN Z sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! Use quasi-second order interpolation for TH on GZ points
! TH lives on GZF so we define 2:NZP-1 -> 0,1,NZP,NZP+1 set by BCs
      DO I=0,NXM
        DO J=1,NY
          DO K=2,NZP
            MATL_YZ(J,K)= -(TEMP1/PR(N)) / (DZP(K)*DZFP(K))
            MATD_YZ(J,K)= 1. 
     &           +(TEMP1/PR(N)) / (DZP(K+1)*DZFP(K))
     &           +(TEMP1/PR(N)) / (DZP(K)*DZFP(K))
            MATU_YZ(J,K)= -(TEMP1/PR(N)) / (DZP(K+1)*DZFP(K))
! Define RHS vector
            VEC_YZ(J,K)=RTH(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_TH_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,N)
      ! WRITE(*,*) 'finished apply BC Z RANK=', RANK , 'I=', I
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
      ! WRITE(*,*) 'finished TDMA F Z RANK=', RANK , 'I=', I
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
      ! WRITE(*,*) 'finished TDMA B Z RANK=', RANK , 'I=', I
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            TH(I,K,J,N)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|

C------------------------   BEGIN Y sweep      ------------------------|
! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO K=0,NZP+1
          MATL_YZ(J,K)=0.
          MATD_YZ(J,K)=1.
          MATU_YZ(J,K)=0.
          VEC_YZ(J,K)=0.
        END DO
      END DO

! Build implicit matrix
! Use quasi-second order interpolation for TH on GY points
! TH lives on GYF so we define 2:NY-1 -> 0,1,NY,NY+1 set by BCs
      DO I=0,NXM
        DO J=2,NY
          DO K=1,NZP
            MATL_YZ(J,K)= -(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD_YZ(J,K)= 1. 
     &           +(TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU_YZ(J,K)= -(TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC_YZ(J,K)=TH(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
        IF (LES) THEN
          IF (RANK.EQ.0) write(*,*) 
     &    'LES not supported in duct flow'
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
          STOP
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN
             CALL APPLY_BC_TH_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            TH(I,K,J,N)=VEC_YZ(J,K)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Y sweep        ------------------------|
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
!     !  &              ,999.0,TMP
!       CLOSE(15)
! 151   format(E40.30)

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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,K,1)=U_BC_YMIN_C1
            CU1(0,K,0)=U_BC_YMIN_C1
          END IF
        END DO
C Neumann
      ELSE IF (U_BC_YMIN.EQ.1) THEN
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,1)=CU1(I,K,2)
            CU1(I,K,0)=CU1(I,K,1)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,K,1)=CU1(0,K,2)-DY(2)*U_BC_YMIN_C1
            CU1(0,K,0)=CU1(0,K,1)-DY(1)*U_BC_YMIN_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,K,1)=W_BC_YMIN_C1
            CU3(0,K,0)=W_BC_YMIN_C1
          END IF
        END DO
C Neumann
      ELSE IF (W_BC_YMIN.EQ.1) THEN
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU3(I,K,1)=CU3(I,K,2)
            CU3(I,K,0)=CU3(I,K,1)
            ! CU3(I,K,0)=CU3(I,K,2)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,K,1)=CU3(0,K,2)-DY(2)*W_BC_YMIN_C1
            CU3(0,K,0)=CU3(0,K,1)-DY(1)*W_BC_YMIN_C1
          END IF
        END DO
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
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,K,1)=CU2(0,K,2)-DYF(1)*V_BC_YMIN_C1
            CU2(0,K,0)=CU2(0,K,1)-DYF(1)*V_BC_YMIN_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,K,NY)=U_BC_YMAX_C1
            CU1(0,K,NY+1)=U_BC_YMAX_C1
          END IF
        END DO
C Neumann
      ELSE IF (U_BC_YMAX.EQ.1) THEN
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,NY)=CU1(I,K,NY-1)
            CU1(I,K,NY+1)=CU1(I,K,NY)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,K,NY)=CU1(0,K,NY-1)+DY(NY)*U_BC_YMAX_C1
            CU1(0,K,NY+1)=CU1(0,K,NY)+DY(NY)*U_BC_YMAX_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,K,NY)=W_BC_YMAX_C1
            CU3(0,K,NY+1)=W_BC_YMAX_C1
          END IF
        END DO
C Neumann
      ELSE IF (W_BC_YMAX.EQ.1) THEN
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU3(I,K,NY)=CU3(I,K,NY-1)
            CU3(I,K,NY+1)=CU3(I,K,NY)
            ! CU3(I,K,NY+1)=CU3(I,K,NY-1)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,K,NY)=CU3(0,K,NY-1)+DY(NY)*W_BC_YMAX_C1
            CU3(0,K,NY+1)=CU3(0,K,NY)+DY(NY)*W_BC_YMAX_C1
          END IF
        END DO
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
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,K,NY+1)=CU2(0,K,NY)+DY(NY)*V_BC_YMAX_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,1,J)=U_BC_ZMIN_C1
            CU1(0,0,J)=U_BC_ZMIN_C1
          END IF
        END DO
C Neumann
      ELSE IF (U_BC_ZMIN.EQ.1) THEN
        DO J=0,NY+1
          DO I=0,NXP-1
            CU1(I,1,J)=CU1(I,2,J)
            CU1(I,0,J)=CU1(I,1,J)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,1,J)=CU1(0,2,J)-DZ(2)*U_BC_ZMIN_C1
            CU1(0,0,J)=CU1(0,1,J)-DZ(1)*U_BC_ZMIN_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,1,J)=V_BC_ZMIN_C1
            CU2(0,0,J)=V_BC_ZMIN_C1
          END IF
        END DO
C Neumann
      ELSE IF (V_BC_ZMIN.EQ.1) THEN
        DO J=0,NY+1
          DO I=0,NXP-1
            CU2(I,1,J)=CU2(I,2,J)
            CU2(I,0,J)=CU2(I,1,J)
            ! CU2(I,0,J)=CU2(I,2,J)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,1,J)=CU2(0,2,J)-DZ(2)*V_BC_ZMIN_C1
            CU2(0,0,J)=CU2(0,1,J)-DZ(1)*V_BC_ZMIN_C1
          END IF
        END DO
      ELSE
        STOP 'Error: V_BC_ZMIN must be 0, or 1'
      END IF


C Dirichlet
      IF (W_BC_ZMIN.EQ.0) THEN
C Set the vertical velocity at GZF(1) (halfway between GZ(1) and GZ(2))
        DO J=0,NY+1
          DO I=0,NXP-1
            CU3(I,1,J)=2.d0*W_BC_ZMIN_C1-CU3(I,2,J)
            CU3(I,0,J)=CU3(I,1,J)
          END DO
        END DO
C Neumann
      ELSE IF (W_BC_ZMIN.EQ.1) THEN
        DO J=0,NY+1
          DO I=0,NXP-1
            CU3(I,1,J)=CU3(I,2,J)
            CU3(I,0,J)=CU3(I,1,J)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,1,J)=CU3(0,2,J)-DZF(1)*W_BC_ZMIN_C1
            CU3(0,0,J)=CU3(0,1,J)-DZF(1)*W_BC_ZMIN_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,NZ,J)=U_BC_ZMAX_C1
            CU1(0,NZ+1,J)=U_BC_ZMAX_C1
          END IF
        END DO
C Neumann
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
        DO J=0,NY+1
          DO I=0,NXP-1
            CU1(I,NZ,J)=CU1(I,NZ-1,J)
            CU1(I,NZ+1,J)=CU1(I,NZ,J)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU1(0,NZ,J)=CU1(0,NZ-1,J)+DZ(NZ)*U_BC_ZMAX_C1
            CU1(0,NZ+1,J)=CU1(0,NZ,J)+DZ(NZ)*U_BC_ZMAX_C1
          END IF
        END DO
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
C Now, set only the mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,NZ,J)=V_BC_ZMAX_C1
            CU2(0,NZ+1,J)=V_BC_ZMAX_C1
          END IF
        END DO
C Neumann
      ELSE IF (V_BC_ZMAX.EQ.1) THEN
        DO J=0,NY+1
          DO I=0,NXP-1
            CU2(I,NZ,J)=CU2(I,NZ-1,J)
            CU2(I,NZ+1,J)=CU2(I,NZ,J)
            ! CU2(I,NZ+1,J)=CU2(I,NZ-1,J)
          END DO
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU2(0,NZ,J)=CU2(0,NZ-1,J)+DZ(NZ)*V_BC_ZMAX_C1
            CU2(0,NZ+1,J)=CU2(0,NZ,J)+DZ(NZ)*V_BC_ZMAX_C1
          END IF
        END DO
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
C Now, apply BC to mean
          IF (RANKZ.EQ.0) THEN
            CU3(0,NZ+1,J)=CU3(0,NZ,J)+DZ(NZ)*W_BC_ZMAX_C1
          END IF
        END DO
      ELSE
        STOP 'Error: W_BC_ZMAX must be 0, or 1'
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_CR1_CHI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K,J

C Y Direction Lower
      IF (RANKY.EQ.0) THEN
        DO I=0,NXP-1
          DO K=0,NZ+1
            CR1(I,K,1)=CR1(I,K,2)
            CR1(I,K,0)=CR1(I,K,1)
          END DO
        END DO
      END IF
C Y Direction Upper
      IF (RANKY.EQ.NPROCY-1) THEN
        DO I=0,NXP-1
          DO K=0,NZ+1
            CR1(I,K,NY)=CR1(I,K,NY-1)
            CR1(I,K,NY+1)=CR1(I,K,NY)
          END DO
        END DO
      END IF

C Z Direction Lower
      DO I=0,NXP-1
        DO J=0,NY+1
          CR1(I,1,J)=CR1(I,2,J)
          CR1(I,0,J)=CR1(I,1,J)
        END DO
      END DO
C Z Direction Upper
      DO I=0,NXP-1
        DO J=0,NY+1
          CR1(I,NZ,J)=CR1(I,NZ-1,J)
          CR1(I,NZ+1,J)=CR1(I,NZ,J)
        END DO
      END DO

      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-floor BC applied during the Y-sweep of the duct TH ADI step.
! Sets the J=0 (ghost) and J=1 (wall) rows of MATL_YZ/MATD_YZ/
! MATU_YZ/VEC_YZ for every K.
! Matrix arrays are passed by reference; their shapes come from
! the declarations in 'header'.
      INCLUDE 'header'
      INTEGER N, K

      IF (TH_BC_YMIN(N).EQ.0) THEN
! Dirichlet at the lower Y wall (TH = TH_BC_YMIN_C1 at J=1)
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = 1.d0
          MATU_YZ(1,K) = 0.d0
          VEC_YZ (1,K) = TH_BC_YMIN_C1(N)
        END DO
      ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
! Neumann at the lower Y wall (dTH/dy = TH_BC_YMIN_C1 enforced at GY(2))
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = -1.d0
          MATU_YZ(1,K) = 1.d0
          VEC_YZ (1,K) = DY(2)*TH_BC_YMIN_C1(N)
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: TH_BC_YMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-ceiling BC applied during the Y-sweep of the duct TH ADI step.
! Sets the J=NY (wall) and J=NY+1 (ghost) rows of MATL_YZ/MATD_YZ/
! MATU_YZ/VEC_YZ for every K.
! Matrix arrays are passed by reference; their shapes come from
! the declarations in 'header'.
      INCLUDE 'header'
      INTEGER N, K

      IF (TH_BC_YMAX(N).EQ.0) THEN
! Dirichlet at the upper Y wall (TH = TH_BC_YMAX_C1 at J=NY)
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = 0.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = TH_BC_YMAX_C1(N)

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
! Neumann at the upper Y wall (dTH/dy = TH_BC_YMAX_C1 enforced at GY(NY))
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = -1.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = DY(NY)*TH_BC_YMAX_C1(N)

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: TH_BC_YMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZLOWER_DUCT_Y(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZUPPER_DUCT_Y(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YLOWER_DUCT_Z(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YUPPER_DUCT_Z(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-floor BC applied during the Z-sweep of the duct TH ADI step.
! Sets the K=0 (ghost) and K=1 (wall) columns of MATL_YZ/MATD_YZ/
! MATU_YZ/VEC_YZ for every J.
! Matrix arrays are passed by reference; their shapes come from
! the declarations in 'header'.
      INCLUDE 'header'
      INTEGER N, J

      IF (TH_BC_ZMIN(N).EQ.0) THEN
! Dirichlet at the lower Z wall (TH = TH_BC_ZMIN_C1 at K=1)
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = 1.d0
          MATU_YZ(J,1) = 0.d0
          VEC_YZ (J,1) = TH_BC_ZMIN_C1(N)
        END DO
      ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
! Neumann at the lower Z wall (dTH/dz = TH_BC_ZMIN_C1 enforced at GZ(2))
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = -1.d0
          MATU_YZ(J,1) = 1.d0
          VEC_YZ (J,1) = DZP(2)*TH_BC_ZMIN_C1(N)
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: TH_BC_ZMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-ceiling BC applied during the Z-sweep of the duct TH ADI step.
! Sets the K=NZP (wall) and K=NZP+1 (ghost) columns of MATL_YZ/
! MATD_YZ/MATU_YZ/VEC_YZ for every J.
! Matrix arrays are passed by reference; their shapes come from
! the declarations in 'header'.
      INCLUDE 'header'
      INTEGER N, J

      IF (TH_BC_ZMAX(N).EQ.0) THEN
! Dirichlet at the upper Z wall (TH = TH_BC_ZMAX_C1 at K=NZP)
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = 0.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = TH_BC_ZMAX_C1(N)

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
! Neumann at the upper Z wall (dTH/dz = TH_BC_ZMAX_C1 enforced at GZ(NZP))
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = -1.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = DZP(NZP)*TH_BC_ZMAX_C1(N)

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: TH_BC_ZMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-floor BC for the Y-sweep of the duct U1 ADI step.
! U1 is cell-centered in Y -- TH-style convention.
      INCLUDE 'header'
      INTEGER K

      IF (U_BC_YMIN.EQ.0) THEN
! Dirichlet at the lower Y wall (U1 = U_BC_YMIN_C1 at J=1)
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = 1.d0
          MATU_YZ(1,K) = 0.d0
          VEC_YZ (1,K) = U_BC_YMIN_C1
        END DO
      ELSE IF (U_BC_YMIN.EQ.1) THEN
! Neumann at the lower Y wall (dU1/dy = U_BC_YMIN_C1 enforced at GY(2))
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = -1.d0
          MATU_YZ(1,K) = 1.d0
          VEC_YZ (1,K) = DY(2)*U_BC_YMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: U_BC_YMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-ceiling BC for the Y-sweep of the duct U1 ADI step.
! U1 is cell-centered in Y -- TH-style convention.
      INCLUDE 'header'
      INTEGER K

      IF (U_BC_YMAX.EQ.0) THEN
! Dirichlet at the upper Y wall (U1 = U_BC_YMAX_C1 at J=NY)
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = 0.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = U_BC_YMAX_C1

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE IF (U_BC_YMAX.EQ.1) THEN
! Neumann at the upper Y wall (dU1/dy = U_BC_YMAX_C1 enforced at GY(NY))
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = -1.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = DY(NY)*U_BC_YMAX_C1

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: U_BC_YMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-floor BC for the Z-sweep of the duct U1 ADI step.
! U1 is cell-centered in Z -- TH-style convention.
      INCLUDE 'header'
      INTEGER J

      IF (U_BC_ZMIN.EQ.0) THEN
! Dirichlet at the lower Z wall (U1 = U_BC_ZMIN_C1 at K=1)
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = 1.d0
          MATU_YZ(J,1) = 0.d0
          VEC_YZ (J,1) = U_BC_ZMIN_C1
        END DO
      ELSE IF (U_BC_ZMIN.EQ.1) THEN
! Neumann at the lower Z wall (dU1/dz = U_BC_ZMIN_C1 enforced at GZ(2))
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = -1.d0
          MATU_YZ(J,1) = 1.d0
          VEC_YZ (J,1) = DZP(2)*U_BC_ZMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: U_BC_ZMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-ceiling BC for the Z-sweep of the duct U1 ADI step.
! U1 is cell-centered in Z -- TH-style convention.
      INCLUDE 'header'
      INTEGER J

      IF (U_BC_ZMAX.EQ.0) THEN
! Dirichlet at the upper Z wall (U1 = U_BC_ZMAX_C1 at K=NZP)
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = 0.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = U_BC_ZMAX_C1

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
! Neumann at the upper Z wall (dU1/dz = U_BC_ZMAX_C1 enforced at GZ(NZP))
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = -1.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = DZP(NZP)*U_BC_ZMAX_C1

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: U_BC_ZMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-floor BC for the Y-sweep of the duct U2 ADI step.
! U2 is face-centered in Y (lives on GY) -- channel-style staggered
! convention. Dirichlet writes rows J=1 AND J=2; Neumann writes row
! J=1 only. Row J=0 is an unused identity placeholder.
      INCLUDE 'header'
      INTEGER K

      IF (V_BC_YMIN.EQ.0) THEN
! Dirichlet at the lower Y wall (U2 = V_BC_YMIN_C1 at J=2)
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = 1.d0
          MATU_YZ(1,K) = 0.d0
          VEC_YZ (1,K) = V_BC_YMIN_C1

          MATL_YZ(2,K) = 0.d0
          MATD_YZ(2,K) = 1.d0
          MATU_YZ(2,K) = 0.d0
          VEC_YZ (2,K) = V_BC_YMIN_C1
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
! Neumann at the lower Y wall (dU2/dy = V_BC_YMIN_C1 enforced at GYF(1))
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = -1.d0
          MATU_YZ(1,K) = 1.d0
          VEC_YZ (1,K) = DYF(1)*V_BC_YMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: V_BC_YMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-ceiling BC for the Y-sweep of the duct U2 ADI step.
! U2 is face-centered in Y (lives on GY) -- channel-style staggered
! convention. Dirichlet writes rows J=NY AND J=NY+1; Neumann writes
! row J=NY+1 only.
      INCLUDE 'header'
      INTEGER K

      IF (V_BC_YMAX.EQ.0) THEN
! Dirichlet at the upper Y wall (U2 = V_BC_YMAX_C1 at J=NY+1)
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = 0.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = V_BC_YMAX_C1

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
! Neumann at the upper Y wall (dU2/dy = V_BC_YMAX_C1 enforced at GYF(NY))
        DO K = 0, NZP+1
          MATL_YZ(NY+1,K) = -1.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = DYF(NY)*V_BC_YMAX_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: V_BC_YMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-floor BC for the Z-sweep of the duct U2 ADI step.
! U2 is cell-centered in Z -- TH-style convention.
      INCLUDE 'header'
      INTEGER J

      IF (V_BC_ZMIN.EQ.0) THEN
! Dirichlet at the lower Z wall (U2 = V_BC_ZMIN_C1 at K=1)
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = 1.d0
          MATU_YZ(J,1) = 0.d0
          VEC_YZ (J,1) = V_BC_ZMIN_C1
        END DO
      ELSE IF (V_BC_ZMIN.EQ.1) THEN
! Neumann at the lower Z wall (dU2/dz = V_BC_ZMIN_C1 enforced at GZ(2))
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = -1.d0
          MATU_YZ(J,1) = 1.d0
          VEC_YZ (J,1) = DZP(2)*V_BC_ZMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: V_BC_ZMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-ceiling BC for the Z-sweep of the duct U2 ADI step.
! U2 is cell-centered in Z -- TH-style convention.
      INCLUDE 'header'
      INTEGER J

      IF (V_BC_ZMAX.EQ.0) THEN
! Dirichlet at the upper Z wall (U2 = V_BC_ZMAX_C1 at K=NZP)
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = 0.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = V_BC_ZMAX_C1

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE IF (V_BC_ZMAX.EQ.1) THEN
! Neumann at the upper Z wall (dU2/dz = V_BC_ZMAX_C1 enforced at GZ(NZP))
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = -1.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = DZP(NZP)*V_BC_ZMAX_C1

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: V_BC_ZMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-floor BC for the Y-sweep of the duct U3 ADI step.
! U3 is cell-centered in Y -- TH-style convention.
      INCLUDE 'header'
      INTEGER K

      IF (W_BC_YMIN.EQ.0) THEN
! Dirichlet at the lower Y wall (U3 = W_BC_YMIN_C1 at J=1)
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = 1.d0
          MATU_YZ(1,K) = 0.d0
          VEC_YZ (1,K) = W_BC_YMIN_C1
        END DO
      ELSE IF (W_BC_YMIN.EQ.1) THEN
! Neumann at the lower Y wall (dU3/dy = W_BC_YMIN_C1 enforced at GY(2))
        DO K = 0, NZP+1
          MATL_YZ(0,K) = 0.d0
          MATD_YZ(0,K) = 1.d0
          MATU_YZ(0,K) = 0.d0
          VEC_YZ (0,K) = 0.d0

          MATL_YZ(1,K) = 0.d0
          MATD_YZ(1,K) = -1.d0
          MATU_YZ(1,K) = 1.d0
          VEC_YZ (1,K) = DY(2)*W_BC_YMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: W_BC_YMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Y-ceiling BC for the Y-sweep of the duct U3 ADI step.
! U3 is cell-centered in Y -- TH-style convention.
      INCLUDE 'header'
      INTEGER K

      IF (W_BC_YMAX.EQ.0) THEN
! Dirichlet at the upper Y wall (U3 = W_BC_YMAX_C1 at J=NY)
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = 0.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = W_BC_YMAX_C1

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE IF (W_BC_YMAX.EQ.1) THEN
! Neumann at the upper Y wall (dU3/dy = W_BC_YMAX_C1 enforced at GY(NY))
        DO K = 0, NZP+1
          MATL_YZ(NY,K)   = -1.d0
          MATD_YZ(NY,K)   = 1.d0
          MATU_YZ(NY,K)   = 0.d0
          VEC_YZ (NY,K)   = DY(NY)*W_BC_YMAX_C1

          MATL_YZ(NY+1,K) = 0.d0
          MATD_YZ(NY+1,K) = 1.d0
          MATU_YZ(NY+1,K) = 0.d0
          VEC_YZ (NY+1,K) = 0.d0
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: W_BC_YMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-floor BC for the Z-sweep of the duct U3 ADI step.
! U3 is face-centered in Z (lives on GZ) -- channel-style staggered
! convention. Dirichlet writes columns K=1 AND K=2; Neumann writes
! column K=1 only. Column K=0 is an unused identity placeholder.
      INCLUDE 'header'
      INTEGER J

      IF (W_BC_ZMIN.EQ.0) THEN
! Dirichlet at the lower Z wall (U3 = W_BC_ZMIN_C1 at K=2)
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = 1.d0
          MATU_YZ(J,1) = 0.d0
          VEC_YZ (J,1) = W_BC_ZMIN_C1

          MATL_YZ(J,2) = 0.d0
          MATD_YZ(J,2) = 1.d0
          MATU_YZ(J,2) = 0.d0
          VEC_YZ (J,2) = W_BC_ZMIN_C1
        END DO
      ELSE IF (W_BC_ZMIN.EQ.1) THEN
! Neumann at the lower Z wall (dU3/dz = W_BC_ZMIN_C1 enforced at GZP(1))
        DO J = 0, NY+1
          MATL_YZ(J,0) = 0.d0
          MATD_YZ(J,0) = 1.d0
          MATU_YZ(J,0) = 0.d0
          VEC_YZ (J,0) = 0.d0

          MATL_YZ(J,1) = 0.d0
          MATD_YZ(J,1) = -1.d0
          MATU_YZ(J,1) = 1.d0
          VEC_YZ (J,1) = DZFP(1)*W_BC_ZMIN_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: W_BC_ZMIN must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Z-ceiling BC for the Z-sweep of the duct U3 ADI step.
! U3 is face-centered in Z (lives on GZ) -- channel-style staggered
! convention. Dirichlet writes columns K=NZP AND K=NZP+1; Neumann
! writes column K=NZP+1 only.
      INCLUDE 'header'
      INTEGER J

      IF (W_BC_ZMAX.EQ.0) THEN
! Dirichlet at the upper Z wall (U3 = W_BC_ZMAX_C1 at K=NZP+1)
        DO J = 0, NY+1
          MATL_YZ(J,NZP)   = 0.d0
          MATD_YZ(J,NZP)   = 1.d0
          MATU_YZ(J,NZP)   = 0.d0
          VEC_YZ (J,NZP)   = W_BC_ZMAX_C1

          MATL_YZ(J,NZP+1) = 0.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = W_BC_ZMAX_C1
        END DO
      ELSE IF (W_BC_ZMAX.EQ.1) THEN
! Neumann at the upper Z wall (dU3/dz = W_BC_ZMAX_C1 enforced at GZP(NZP))
        DO J = 0, NY+1
          MATL_YZ(J,NZP+1) = -1.d0
          MATD_YZ(J,NZP+1) = 1.d0
          MATU_YZ(J,NZP+1) = 0.d0
          VEC_YZ (J,NZP+1) = DZFP(NZP)*W_BC_ZMAX_C1
        END DO
      ELSE
        IF (RANK.EQ.0) WRITE(*,*)
     &    'ERROR: W_BC_ZMAX must be 0 (Dirichlet) or 1 (Neumann)'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        STOP
      END IF

      RETURN
      END




C depreciated routines below

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YLOWER_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YUPPER_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZLOWER_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZUPPER_DUCT(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END



















C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE BUOYANCY_PDF_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Here the buoyancy pdf is computed using a rescaled histogram for a
C single snapshot of the complete 3D flow field and then written to file
C NOTE THIS REQUIRES UNIFORM Z with how its currently programmed

      include 'header'

      integer i,j,k,n
      real*8 bins(nbins+1),dbins(nbins)
      integer testhist(nbins)
      real*8 testhist_rescaled(nbins),counts(nbins)
      character*55 fname
      character*20 gname
      real*8 Diagbins(1:nbins)

 ! bin variables
      deltabin=(binmax-binmin)/real(nbins)
      
! generate bin edge values
      do i=1,nbins+1
         bins(i)=binmin+real(i-1)*deltabin
      end do
      
! generate bin edge values
      do i=1,nbins+1
         bins(i)=binmin+real(i-1)*deltabin
      end do
      
! calculate bin widths
      do i=1,nbins
         dbins(i)=bins(i+1)-bins(i)
      end do
      
! zero initial histogram matrices
      do i=1,nbins
         testhist(i)=0
         testhist_rescaled(i)=0.d0
         counts(i)=0.d0
      end do

! calculate buoyancy histogram over (x,z) directions at each height y
! and rescale by delta(y), the height of the layer

! first do edge values
      do i=0,nxm
         do k=2,nzp
            binlabel=floor((
     &       ( (DZFP(K-1)*TH(i,k,1,1) + DZFP(K)*TH(i,k-1,1,1))
     &           /(2.0*DZP(K)) )
     &          -binmin)/deltabin)+1
            IF (BINLABEL.LT.1) BINLABEL=1
            IF (BINLABEL.GT.NBINS) BINLABEL=NBINS
            testhist(binlabel)=testhist(binlabel)+1
         end do
      end do

      do i=1,nbins
         testhist_rescaled(i)=real(testhist(i))/dbins(i)
         counts(i)=counts(i)+testhist_rescaled(i)*dyf(1)/2.d0
         testhist(i)=0
      end do

      do i=0,nxm
         do k=2,nzp
            binlabel=floor((
     &       ( (DZFP(K-1)*TH(i,k,ny,1) + DZFP(K)*TH(i,k-1,ny,1))
     &           /(2.0*DZP(K)) )
     &          -binmin)/deltabin)+1
            IF (BINLABEL.LT.1) BINLABEL=1
            IF (BINLABEL.GT.NBINS) BINLABEL=NBINS
            testhist(binlabel)=testhist(binlabel)+1
         end do
      end do

      do i=1,nbins
         testhist_rescaled(i)=real(testhist(i))/dbins(i)
         counts(i)=counts(i)+testhist_rescaled(i)*dyf(ny)/2.d0
         testhist(i)=0
      end do

! do interior points      
      do j=2,ny-1
         ! bin data
         do i=0,nxm
            do k=2,nzp
               binlabel=floor((
     &       ( (DZFP(K-1)*TH(i,k,j,1) + DZFP(K)*TH(i,k-1,j,1))
     &           /(2.0*DZP(K)) )
     &          -binmin)/deltabin)+1
               IF (BINLABEL.LT.1) BINLABEL=1
               IF (BINLABEL.GT.NBINS) BINLABEL=NBINS
               testhist(binlabel)=testhist(binlabel)+1
            end do
         end do
         ! rescale binned data to account for varying bin widths
         ! and store in counts for the layer contribution
         ! also rezero testhist for next layer  
         do i=1,nbins
            testhist_rescaled(i)=real(testhist(i))/dbins(i)
            counts(i)=counts(i)+testhist_rescaled(i)*dy(j)
            testhist(i)=0
         end do
      end do

!     sum resulting counts from all mpi processors
      if (use_mpi) then
      call mpi_allreduce(mpi_in_place,counts,nbins,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      end if
!      write(*,*) 'gathered all data from all processors'

#ifdef HDF5
! write out resulting rescaled histogram to hdf5 file
      fname='buoypdf.h5'

      gname='time'
      call WriteHDF5_real(fname,gname,TIME)
!      write(*,*) 'time written to buoypdf.h5'
      
      ! only need to write from one process  
      IF ((RANKZ.eq.0).and.(RANKY.eq.0)) then
!         write(*,*) 'starting write to buoypdf.h5'
!         do j=1,nbins
!            write(*,*) 'bins,dbins,counts=',bins(j),dbins(j),
!     &           counts(j)
!         end do
         gname='bins'
         Diagbins=bins(1:nbins)
         call WriteBinsH5(fname,bins,dbins,counts)
!         write(*,*) 'bin locations written to buoypdf.h5'

!         write(*,*) 'binned data written to file'
      END IF

      if (use_mpi) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
      end if
         
#else      

#endif

      do i=1,nbins
         counts(i)=0.d0
         testhist_rescaled(i)=0.d0
         testhist(i)=0
      end do

      RETURN
      END
