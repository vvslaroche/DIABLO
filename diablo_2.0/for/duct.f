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
      CHARACTER*55 TEMPFNAME

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
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CR1(I,K,J)=CU1(I,K,J)
            CR2(I,K,J)=CU2(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO

C Add the R-K term from the rk-1 step
      IF (RK_STEP .GT. 1) THEN
        DO J=2,NY
          DO K=2,NZ
            DO I=0,NXP-1
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
      END IF

C Add the pressure gradient to the RHS for X as explicit Euler
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
          END DO
        END DO
      END DO

C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO
C Add the pressure gradient to the RHS for Y as explicit Euler
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

C Take the z-derivative of the pressure at GZ points in Fourier space
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K-1,J)) / DZ(J)
          END DO
        END DO
      END DO
C Add the pressure gradient to the RHS for Z as explicit Euler
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

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
! If a scalar contributes to the denisty, RI is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode in the
! y-direction, which corresponds to the plane-average.
! The plane averaged density balances the hydrostatic pressure
      IF (RANKZ.eq.0) THEN
        ISTART=1
      ELSE
        ISTART=0
      END IF
      DO J=2,NY
        DO K=2,NZ
          DO I=ISTART,NXP-1
! Use second order interpolation
            CF2(I,K,J)=CF2(I,K,J)+RI(N)*GRAV_Y*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J)+RI(N)*GRAV_X*CTH(I,K,J,N)
            CF3(I,K,J)=CF3(I,K,J)+RI(N)*GRAV_Z*CTH(I,K,J,N)
          END DO
        END DO
      END DO

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

C Convert RHS terms to physical space
      CALL FFT_X_TO_PHYSICAL(CR1,R1,0,NY+1)
      CALL FFT_X_TO_PHYSICAL(CR2,R2,0,NY+1)
      CALL FFT_X_TO_PHYSICAL(CR3,R3,0,NY+1)

C Compute the y- and z-derivative viscous terms in physical space
C and add to RHS
C This is the explicit part of the Crank-Nicolson term
C NU_V_SCALE is a coefficient defined as the ratio of the vertical
C to horizontal viscosity and can be used to add anisotropic viscosity
      DO J=2,NY
        DO K=2,NZP
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


C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    BEGIN timestep scalar equation       X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
C      We solve the the passive scalar before the velocity so that
C      it is advected with the velocity from the previous R-K step
C      which we have already made divergence free

! Solve the implicit equation for THETA
! Note that the system size is NY+1, but only 1..NY are used

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
     &           - TH(I,K,J,N)
     &           + TH(I,K-1,J,N) / (DZP(K)*DZFP(K))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           + TH(I,K+1,J,N) / (DZP(K+1)*DZFP(K))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           - TH(I,K,J,N) / (DZP(K)*DZFP(K))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           - TH(I,K,J,N) / (DZP(K+1)*DZFP(K))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
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
             CALL APPLY_BC_TH_DUCT_MPI(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
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
            TH(I,K,J,N)=VEC(I,J)
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
     &           - TH(I,K,J,N)
     &           + TH(I,K,J-1,N) / (DY(J)*DYF(J))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           + TH(I,K,J+1,N) / (DY(J+1)*DYF(J))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           - TH(I,K,J,N) / (DY(J)*DYF(J))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
     &           - TH(I,K,J,N) / (DY(J+1)*DYF(J))
     &                 *(TEMP1/PR(N)*NU_V_SCALE)
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
             CALL APPLY_BC_TH_DUCT_MPI(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL TDMA_F_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
             CALL TDMA_B_REAL_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &            MATU_YZ,VEC_YZ,NY,NZP)
          ELSE
! Else we are running in serial mode
          END IF

        DO J=0,NY+1
          DO K=0,NZP+1
            TH(I,K,J,N)=VEC(I,J)
          END DO
        END DO

! END DO I
      END DO
C------------------------   END Z sweep        ------------------------|

! End do number of passive scalars
        END DO

C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    END timestep scalar equation         X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|






C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    BEGIN timestep scalar equation       X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|


C------------------------   BEGIN Y sweep      ------------------------|

C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|

C------------------------   END Z sweep        ------------------------|

C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    END timestep scalar equation         X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|


C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    BEGIN timestep U1                    X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|


C------------------------   BEGIN Y sweep      ------------------------|

C------------------------   END Y sweep        ------------------------|

C------------------------   BEGIN Z sweep      ------------------------|

C------------------------   END Z sweep        ------------------------|

C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|
CX*X*X*X*X*X*X*X    END timestep U1                      X*X*X*X*X*X*X*|
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X|


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

C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space

      CALL REM_DIV_DUCT

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=2,NY
        DO K=2,NZ
          DO I=0,NXP-1
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

      ! Fix disparities at the boundary due to the thoms algorithm in parallel
      IF (USE_MPI) THEN
         CALL GHOST_DUCT_MPI
      END IF


      call mpi_barrier(MPI_COMM_WORLD,ierror)
      IF (RANK.EQ.0) write(*,*) 'Ending RKSTEP ',RK_STEP

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


151   format(F30.20)



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
      SUBROUTINE REM_DIV_DUCT2
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
        ! DO J=JLIM1,JLIM2
          DO K=1,NZ
            MATL_C_YZ(J,K)=1./(DZ(K)*DZF(K))
            MATD_C_YZ(J,K)=-KX2(I)
     &          -1./(DZ(K+1)*DZF(K)) - 1./(DZ(K)*DZF(K))
            MATU_C_YZ(J,K)=1./(DZ(K+1)*DZF(K))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY
        ! DO J=JLIM1,JLIM2
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


      ! CALL APPLY_BC_CR1_CHI

      WRITE(*,*) 'writing CR1_CHECK1 on rank ', RANK  
      FNAME='CR1_CHECK1_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)



C START LOOP TO ITERATE REM_DIV BOTH DIRECTIONS

      DO NI=1,10


C ----------------- BEGIN Y
    
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

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
          CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
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

      ! CALL APPLY_BC_CR1_CHI
C ----------------- END Y


C ----------------- BEGIN Z

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
     &           - CR1(I,K,J-1)/(DY(J)*DYF(J))
     &           - CR1(I,K,J+1)/(DY(J+1)*DYF(J))
     &           + CR1(I,K,J)/(DY(J)*DYF(J))
     &           + CR1(I,K,J)/(DY(J+1)*DYF(J))
          END DO
        END DO

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
          CALL TDMA_F_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
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

      ! CALL APPLY_BC_CR1_CHI
C ----------------- END Z



C END LOOP ITERATE
      END DO

      WRITE(*,*) 'writing CR1_CHECK2 on rank ', RANK  
      FNAME='CR1_CHECK2_'//char(RANK+48)//'.dat'
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
151   format(F30.20)
      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_DUCT3
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

      INCLUDE 'header'
      INTEGER I,J,K
      INTEGER NI
      CHARACTER*55 FNAME


C Alternating Direction Implicit method, now we implicitly solve the
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
    !  &           - CR1(I,K-1,J)/(DZ(K)*DZF(K))
    !  &           - CR1(I,K+1,J)/(DZ(K+1)*DZF(K))
    !  &           + CR1(I,K,J)/(DZ(K)*DZF(K))
    !  &           + CR1(I,K,J)/(DZ(K+1)*DZF(K))
          END DO
        END DO

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK4 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK4_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)

        IF (I.EQ.0) THEN
      WRITE(*,*) 'writing VEC_C_YZ_CHECK5 on rank ', RANK  
      FNAME='VEC_C_YZ_CHECK5_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(VEC_C_YZ)
      CLOSE(15)
        END IF

          CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
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

        DO J=1,NY
          DO K=1,NZ
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO

      END DO

      ! CALL APPLY_BC_CR1_CHI

      WRITE(*,*) 'writing CR1 on rank ', RANK  
      FNAME='CR1_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) real(CR1)
      CLOSE(15)


C START LOOP TO ITERATE REM_DIV BOTH DIRECTIONS

      DO NI=1,5
C ----------------- BEGIN Z

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
     &           - CR1(I,K,J-1)/(DY(J)*DYF(J))
     &           - CR1(I,K,J+1)/(DY(J+1)*DYF(J))
     &           + CR1(I,K,J)/(DY(J)*DYF(J))
     &           + CR1(I,K,J)/(DY(J+1)*DYF(J))
          END DO
        END DO

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
          CALL TDMA_F_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
C Else we are running in serial mode
        ELSE
          IF (RANK.EQ.0) 
     &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

        DO J=1,NY
          DO K=1,NZ
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO
      END DO

      ! CALL APPLY_BC_CR1_CHI
C ----------------- END Z

C ----------------- BEGIN Y
    
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

C If we are using the MPI package...
        IF (USE_MPI) THEN
          CALL APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
          CALL TDMA_F_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
          CALL TDMA_B_COMP_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,NY,NZ)
C Else we are running in serial mode
        ELSE
          IF (RANK.EQ.0) 
     &      WRITE(*,*) 'ERROR: serial duct solver not implemented yet'
          STOP
        END IF

        DO J=1,NY
          DO K=1,NZ
            CR1(I,K,J)=VEC_C_YZ(J,K)
          END DO
        END DO

      END DO

      ! CALL APPLY_BC_CR1_CHI
C ----------------- END Y




C END LOOP ITERATE
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
151   format(F30.20)
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



      WRITE(*,*) 'writing grids on rank ', RANK  
      FNAME='grids_'//char(RANK+48)//'.dat'
      OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(15,151) GX,51.0,DX,52.0,GXF,53.0,DXF,54.0,
     &              GY,61.0,DY,62.0,GYF,63.0,DYF,64.0,
     &              GZ,71.0,DZ,72.0,GZF,73.0,DZF,74.0,
     &              GZP,81.0,DZP,82.0,GZFP,83.0,DZFP
    !  &              ,999.0,TMP
      CLOSE(15)
151   format(F10.5)

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
      SUBROUTINE APPLY_BC_TH_YLOWER_DUCT_Y(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_YUPPER_DUCT_Y(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

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
      SUBROUTINE APPLY_BC_TH_ZLOWER_DUCT_Z(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_ZUPPER_DUCT_Z(N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER N

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


      RETURN
      END