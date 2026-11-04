      SUBROUTINE GHOST_CHAN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NXP-1,0:TNKZ,4+N_TH)
      COMPLEX*16 ICPACK(0:NXP-1,0:TNKZ,4+N_TH)

! If we are using more than one processor, then we need to pass data

      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain
      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END


 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GHOST_DUCT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This subroutine is part of the MPI package for the duct flow
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep
! Because of how the data is stored in complex space, only need to 
! communicate in Y, not Z

      include 'header'

      integer I,J,K,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is (0:NXP-1,0:NZ+1,# of variables)
      COMPLEX*16 OCPACK(0:NXP-1,0:NZ+1,4+N_TH)
      COMPLEX*16 ICPACK(0:NXP-1,0:NZ+1,4+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,NZ+1
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,NZ+1
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain
      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,NZ+1
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,NZ+1
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(NZ+2)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZ+1
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF


      RETURN
      END




      SUBROUTINE GHOST_LES_MPI_KAPPA_T
! This subroutine is part of the MPI package for the LES subroutine
! Here, after calculating the SGS viscosity, NU_T on each core,
! We need to share the ghost cells between neighboring processors

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is 0:NXM,0:NZP-1
      REAL*8 OCPACK(0:NXM,0:NZP-1)
      REAL*8 ICPACK(0:NXM,0:NZP-1)

      DO N=1,N_TH

! If we are using more than one processor, then we need to pass data

      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells. Instead, set NU_T=0 at the lower wall
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,1,N)=0.d0
            KAPPA_T(I,K,2,N)=0.d0
          END DO
        END DO

! Pass data up to the next process from GY(NY)
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=KAPPA_T(I,K,NY,N)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=KAPPA_T(I,K,NY,N)
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,1,N)=ICPACK(I,K)
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,1,N)=ICPACK(I,K)
          END DO
        END DO
      END IF

! Now, we have hit the top process.  Set the BCs and pass data down

      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells.
! Set KAPPA_T=0 at the upper wall
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,NY,N)=0.d0
            KAPPA_T(I,K,NY+1,N)=0.d0
          END DO
        END DO

! Now, send data down the chain
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=KAPPA_T(I,K,2,N)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=KAPPA_T(I,K,2,N)
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,NY+1,N)=ICPACK(I,K)
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,NY+1,N)=ICPACK(I,K)
          END DO
        END DO
      END IF

      ELSE
! Here, NPROCY=1, so we just need to set the boundary values
! Set NU_T=0 at the lower wall
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,1,N)=0.d0
            KAPPA_T(I,K,2,N)=0.d0
          END DO
        END DO
! Set NU_T=0 at the upper wall
        DO K=0,NZP-1
          DO I=0,NXM
            KAPPA_T(I,K,NY,N)=0.d0
            KAPPA_T(I,K,NY+1,N)=0.d0
          END DO
        END DO

      END IF

      END DO

      RETURN
      END

      SUBROUTINE GHOST_LES_MPI
! This subroutine is part of the MPI package for the LES subroutine
! Here, after calculating the SGS viscosity, NU_T on each core,
! We need to share the ghost cells between neighboring processors

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is 0:NXM,0:NZP-1
      REAL*8 OCPACK(0:NXM,0:NZP-1)
      REAL*8 ICPACK(0:NXM,0:NZP-1)

! If we are using more than one processor, then we need to pass data

      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells. Instead, set NU_T=0 at the lower wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=0.d0
            NU_T(I,K,2)=0.d0
          END DO
        END DO

! Pass data up to the next process from GY(NY)
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,NY)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,NY)
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=ICPACK(I,K)
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=ICPACK(I,K)
          END DO
        END DO
      END IF

! Now, we have hit the top process.  Set the BCs and pass data down

      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells.
! Set NU_T=0 at the upper wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY)=0.d0
            NU_T(I,K,NY+1)=0.d0
          END DO
        END DO

! Now, send data down the chain
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,2)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,2)
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY+1)=ICPACK(I,K)
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY+1)=ICPACK(I,K)
          END DO
        END DO
      END IF

      ELSE
! Here, NPROCY=1, so we just need to set the boundary values
! Set NU_T=0 at the lower wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=0.d0
            NU_T(I,K,2)=0.d0
          END DO
        END DO
! Set NU_T=0 at the upper wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY)=0.d0
            NU_T(I,K,NY+1)=0.d0
          END DO
        END DO

      END IF

      RETURN
      END


      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

      real*8 OCPACK(3),ICPACK(3)

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN

! Set the lower ghost cells
        GYF(0)=2.d0*GYF(1)-GYF(2)
        GY(0)=2.d0*GY(1)-GY(2)

        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Set the top ghost cell
        GYF(NY+1)=2.*GYF(NY)-GYF(NYM)

! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      END IF
! AT this point we have passed data up the chain

      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)


! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX-1,0:INY+1), B(0:INX-1,0:INY+1), C(0:INX-1,0:INY+1)
     &         , G(0:INX-1,0:INY+1)

      REAL*8 OCPACK(0:INX-1,4),ICPACK(0:INX-1,4)


      IF (RANKY.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*INX,MPI_DOUBLE_PRECISION,RANKY-1,12
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
        DO I=0,INX-1
        A(I,1)=OCPACK(I,1)
        B(I,1)=OCPACK(I,2)
        C(I,1)=OCPACK(I,3)
        G(I,1)=OCPACK(I,4)
        END DO
! If we aren't the lowest process, start at J=2
        DO J=2,INY
          DO I=0,INX-1
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO
      ELSE
! Here, we are the lowest process, start solving at J=1
      DO J=1,INY
        DO I=0,INX-1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
      END DO
      END IF

      IF (RANKY.NE.NPROCY-1) THEN
        DO I=0,INX-1
        ICPACK(I,1)=A(I,INY)
        ICPACK(I,2)=B(I,INY)
        ICPACK(I,3)=C(I,INY)
        ICPACK(I,4)=G(I,INY)
        END DO
        CALL MPI_SEND(ICPACK,4*INX,MPI_DOUBLE_PRECISION,RANKY+1,12
     &               ,MPI_COMM_Y,ierror)
      ELSE
! Here, we are at the upper process, so solve one more row containing
! the boundary conditions
        J=INY+1
        DO I=0,INX-1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX,0:INY+1), B(0:INX,0:INY+1), C(0:INX,0:INY+1)
      COMPLEX*16 G(0:INX,0:INY+1)

      COMPLEX*16 OCPACK(4),ICPACK(4)

      DO I=0,INX

      IF (RANKY.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_COMPLEX,RANKY-1,13
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
        A(I,1)=REAL(OCPACK(1))
        B(I,1)=REAL(OCPACK(2))
        C(I,1)=REAL(OCPACK(3))
        G(I,1)=OCPACK(4)
      END IF

      DO J=2,INY
        IF (B(I,J-1).ne.0.d0) THEN
          A(I,J)=-A(I,J)/B(I,J-1)
        ELSE
          A(I,J)=0.d0
        END IF
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
      END DO

      IF (RANKY.NE.NPROCY-1) THEN
        ICPACK(1)=A(I,INY)
        ICPACK(2)=B(I,INY)
        ICPACK(3)=C(I,INY)
        ICPACK(4)=G(I,INY)
        CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_COMPLEX,RANKY+1,13
     &               ,MPI_COMM_Y,ierror)
      END IF


      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
c ￼￼￼C[b1c1 0 C[a2 b2 c2 C[0a3b3
c 0 0... 0 0... c30...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX-1,0:INY+1), B(0:INX-1,0:INY+1), C(0:INX-1,0:INY+1)
     &     , G(0:INX-1,0:INY+1)
      REAL*8 ICPACK(1),OCPACK(1)
      DO I = 0,INX-1
       IF (RANKY.NE.NPROCY-1) THEN
! If we aren’t the highest process, then wait for data
        CALL MPI_RECV(OCPACK,1,MPI_DOUBLE_PRECISION,RANKY+1,10
     &               ,MPI_COMM_Y,status,ierror)
          G(I,INY+1)=OCPACK(1)
       ELSE
! Else, if we are the highest process, compute the solution at j=INY
          G(I,INY+1)=G(I,INY+1)/B(I,INY+1)
       END IF
! All processes solve from INY..1
      DO J=INY,0,-1
        IF (B(I,J).ne.0.d0) THEN
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        ELSE
          G(I,J)=0.d0
        END IF
      END DO
      IF (RANKY.NE.0) THEN
          ICPACK(1)=G(I,2)
        CALL MPI_SEND(ICPACK,1,MPI_DOUBLE_PRECISION,RANKY-1,10
     &               ,MPI_COMM_Y,ierror)
      END IF
      END DO
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX,0:INY+1), B(0:INX,0:INY+1), C(0:INX,0:INY+1)
      COMPLEX*16 G(0:INX,0:INY+1)

      DO I=0,INX

      IF (RANKY.NE.NPROCY-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(G(I,INY+1),1,MPI_DOUBLE_COMPLEX,RANKY+1,11
     &               ,MPI_COMM_Y,status,ierror)
        J=INY
        if (B(I,J).ne.0.d0) THEN
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END IF
      ELSE
C Else, if we are the highest process, then compute the solution at j=INY
        IF (B(I,INY).ne.0.d0) THEN
          G(I,INY)=G(I,INY)/B(I,INY)
        ELSE
          G(I,INY)=0.d0
        END IF
      END IF

      DO J=INY-1,1,-1
        IF (B(I,J).ne.0.d0) THEN
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        ELSE
          G(I,J)=0.d0
        END IF
      END DO

      IF (RANKY.NE.0) THEN
        CALL MPI_SEND(G(I,2),1,MPI_DOUBLE_COMPLEX,RANKY-1,11
     &               ,MPI_COMM_Y,ierror)
      END IF

      END DO

      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_F_REAL_DUCT_Y_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm FORWARD sweep in Y-DIRECTION. Input data is in
C physical (REAL) space and formatted for DUCT geometry in YZ.
C Requires communication in Y.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
     &     , G(0:INY+1,0:INZ+1)
      REAL*8 OCPACK(1:INZ,4),ICPACK(1:INZ,4)

      IF (RANKY.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*INZ,MPI_DOUBLE_PRECISION,RANKY-1,10
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
        DO K=1,INZ
          A(1,K)=OCPACK(K,1)
          B(1,K)=OCPACK(K,2)
          C(1,K)=OCPACK(K,3)
          G(1,K)=OCPACK(K,4)
        END DO
! If we aren't the lowest process, start at J=2
        DO J=2,INY
          DO K=1,INZ
            A(J,K)=-A(J,K)/B(J-1,K)
            B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
            G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
          END DO
        END DO
      ELSE
! Here, we are the lowest process, start solving at J=1
        DO J=1,INY
          DO K=1,INZ
            A(J,K)=-A(J,K)/B(J-1,K)
            B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
            G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
          END DO
        END DO
      END IF

      IF (RANKY.NE.NPROCY-1) THEN
        DO K=1,INZ
          ICPACK(K,1)=A(INY,K)
          ICPACK(K,2)=B(INY,K)
          ICPACK(K,3)=C(INY,K)
          ICPACK(K,4)=G(INY,K)
        END DO
        CALL MPI_SEND(ICPACK,4*INZ,MPI_DOUBLE_PRECISION,RANKY+1,10
     &               ,MPI_COMM_Y,ierror)
      ELSE
! Here, we are at the upper process, so solve one more row containing
! the boundary conditions
        J=INY+1
        DO K=1,INZ
          A(J,K)=-A(J,K)/B(J-1,K)
          B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
          G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_B_REAL_DUCT_Y_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm BACKWARD sweep in Y-DIRECTION. Input data is in
C physical (REAL) space and formatted for DUCT geometry in YZ.
C Requires communication in Y.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
     &     , G(0:INY+1,0:INZ+1)
      REAL*8 OCPACK(1),ICPACK(1)

      DO K=1,INZ
        IF (RANKY.NE.NPROCY-1) THEN
! If we aren’t the highest process, then wait for data
          CALL MPI_RECV(OCPACK,1,MPI_DOUBLE_PRECISION,RANKY+1,11
     &                 ,MPI_COMM_Y,status,ierror)
          G(INY+1,K)=OCPACK(1)
        ELSE
! Else, if we are the highest process, compute the solution at J=INY
          G(INY+1,K)=G(INY+1,K)/B(INY+1,K)
        END IF

! All processes solve from INY..1
        DO J=INY,0,-1
          IF (B(J,K).NE.0.d0) THEN
            G(J,K)=(G(J,K)-C(J,K)*G(J+1,K))/B(J,K)
          ELSE
            G(J,K)=0.d0
          END IF
        END DO

        IF (RANKY.NE.0) THEN
          ICPACK(1)=G(2,K)
          CALL MPI_SEND(ICPACK,1,MPI_DOUBLE_PRECISION,RANKY-1,11
     &                 ,MPI_COMM_Y,ierror)
        END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_F_COMP_DUCT_Y_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm FORWARD sweep in Y-DIRECTION. Input data is in
C Fourier (COMPLEX) space and formatted for DUCT geometry in YZ.
C Requires communication in Y.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)
      COMPLEX*16 OCPACK(4),ICPACK(4)

      DO K=1,INZ
C If we aren't the lowest process, then wait for data from below
        IF (RANKY.NE.0) THEN
          CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_COMPLEX,RANKY-1,12
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
          A(1,K)=REAL(OCPACK(1))
          B(1,K)=REAL(OCPACK(2))
          C(1,K)=REAL(OCPACK(3))
          G(1,K)=OCPACK(4)
        END IF

C Do interior points on each process
        DO J=2,INY
          IF (B(J-1,K).NE.0.d0) THEN
            A(J,K)=-A(J,K)/B(J-1,K)
          ELSE
            A(J,K)=0.d0
          END IF
          B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
          G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
        END DO

C If we aren't the highest process, then send data upwards
        IF (RANKY.NE.NPROCY-1) THEN
          ICPACK(1)=A(INY,K)
          ICPACK(2)=B(INY,K)
          ICPACK(3)=C(INY,K)
          ICPACK(4)=G(INY,K)
          CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_COMPLEX,RANKY+1,12
     &               ,MPI_COMM_Y,ierror)
        END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_B_COMP_DUCT_Y_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm BACKWARD sweep in Y-DIRECTION. Input data is in
C Fourier (COMPLEX) space and formatted for duct geometry in YZ.
C Requires communication in Y.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)

      DO K=1,INZ
C If we aren't the highest process, then wait for data from above
        IF (RANKY.NE.NPROCY-1) THEN
          CALL MPI_RECV(G(INY+1,K),1,MPI_DOUBLE_COMPLEX,RANKY+1,13
     &               ,MPI_COMM_Y,status,ierror)
          IF (B(INY,K).NE.0.d0) THEN
            G(INY,K)=(G(INY,K)-C(INY,K)*G(INY+1,K))/B(INY,K)
          END IF
        ELSE
C If we are the highest process, then compute the solution at J=INY
          IF (B(INY,K).NE.0.d0) THEN
            G(INY,K)=G(INY,K)/B(INY,K)
          ELSE
            G(INY,K)=0.d0
          END IF
        END IF

C Do interior points on each process
        DO J=INY-1,1,-1
          IF (B(J,K).NE.0.d0) THEN
            G(J,K)=(G(J,K)-C(J,K)*G(J+1,K))/B(J,K)
          ELSE
            G(J,K)=0.d0
          END IF
        END DO

C If we aren't the lowest process, then send data down
        IF (RANKY.NE.0) THEN
          CALL MPI_SEND(G(2,K),1,MPI_DOUBLE_COMPLEX,RANKY-1,13
     &               ,MPI_COMM_Y,ierror)
        END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GAUSS_NOPIVOT(M, N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C In-place LU decomposition of a real N x N matrix M, no pivoting.
C On entry: M holds the matrix to factor.
C On exit:  M holds the combined L (strict lower triangle, unit diagonal
C           implied) and U (upper triangle including diagonal) factors.
C Used by REM_DIV_DUCT4 block-Thomas inner solve (NZ x NZ dense block).
C The diagonal blocks D_J in that solve are diagonally dominant for
C kx /= 0 modes; the kx=0 anchor BC keeps the kx=0 system well-posed.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER N
      REAL*8 M(0:NZ+1, 0:NZ+1)
      INTEGER I, J, K
      REAL*8 PIV

      DO K = 1, N-1
        PIV = M(K, K)
        DO I = K+1, N
          M(I, K) = M(I, K) / PIV
          DO J = K+1, N
            M(I, J) = M(I, J) - M(I, K) * M(K, J)
          END DO
        END DO
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GAUSS_SOLVE_R(M, RHS, NRHS, N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Forward and back substitution against the LU factors stored in M
C (output of GAUSS_NOPIVOT). Real RHS variant.
C RHS has shape (0:NZ+1, 0:NZ+1) with active entries (1..N, 1..NRHS),
C overwritten in-place with the solution.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER N, NRHS
      REAL*8 M(0:NZ+1, 0:NZ+1)
      REAL*8 RHS(0:NZ+1, 0:NZ+1)
      INTEGER I, J, R

      DO R = 1, NRHS
C Forward solve L y = RHS (L has unit diagonal, stored in strict lower)
        DO I = 2, N
          DO J = 1, I-1
            RHS(I, R) = RHS(I, R) - M(I, J) * RHS(J, R)
          END DO
        END DO
C Back solve U x = y (U stored in upper triangle including diagonal)
        DO I = N, 1, -1
          DO J = I+1, N
            RHS(I, R) = RHS(I, R) - M(I, J) * RHS(J, R)
          END DO
          RHS(I, R) = RHS(I, R) / M(I, I)
        END DO
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE GAUSS_SOLVE_C(M, RHS, N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Forward and back substitution against the LU factors stored in M
C (output of GAUSS_NOPIVOT). Single-RHS complex vector variant.
C M is real LU; RHS is complex; mixed-type arithmetic is handled by Fortran.
C RHS has shape (0:NZ+1) with active entries (1..N), overwritten in-place
C with the solution.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER N
      REAL*8 M(0:NZ+1, 0:NZ+1)
      COMPLEX*16 RHS(0:NZ+1)
      INTEGER I, J

      DO I = 2, N
        DO J = 1, I-1
          RHS(I) = RHS(I) - M(I, J) * RHS(J)
        END DO
      END DO
      DO I = N, 1, -1
        DO J = I+1, N
          RHS(I) = RHS(I) - M(I, J) * RHS(J)
        END DO
        RHS(I) = RHS(I) / M(I, I)
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE BLOCK_TDMA_F_COMP_DUCT_Y_MPI(INY, INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Block-Thomas FORWARD sweep in Y for the REM_DIV_DUCT4 pentadiagonal solver.
C Operates on /PENTA_VARS/ arrays from header. Pipeline mirrors
C TDMA_F_COMP_DUCT_Y_MPI: each rank waits for the rank below to finish,
C does its NY block-row eliminations, then sends to the rank above.
C Per neighbor exchange: INZ*INZ + INZ complex doubles (F block + G' vector).
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER INY, INZ
      INTEGER J, K, K1, K2, IDX
      INTEGER MSGSIZE
      COMPLEX*16 SBUF(0:NZ*NZ + NZ)
      COMPLEX*16 GTMP(0:NZ+1)

      MSGSIZE = INZ*INZ + INZ

C Initialize J=1 state on this rank
      IF (RANKY.NE.0) THEN
C Receive packed (F_block, G'_vec) from rank below; place at J=1 slot
        CALL MPI_RECV(SBUF, MSGSIZE, MPI_DOUBLE_COMPLEX, RANKY-1, 15
     &               ,MPI_COMM_Y, status, ierror)
        IDX = 0
        DO K2 = 1, INZ
          DO K1 = 1, INZ
            IDX = IDX + 1
            F_WORK(K1, K2, 1) = REAL(SBUF(IDX))
          END DO
        END DO
        DO K = 1, INZ
          IDX = IDX + 1
          G_PENTA(1, K) = SBUF(IDX)
        END DO
      ELSE
C Rank 0: factor B_1 to obtain F_1 = B_1^{-1} C_1 and G'_1 = B_1^{-1} G_1.
C With Y-floor BCs applied, B_1 is identity, but we use the general path
C so any future BC tweak still produces correct F_1 / G'_1.
        DO K2 = 0, INZ+1
          DO K1 = 0, INZ+1
            D_WORK(K1, K2) = 0.d0
          END DO
        END DO
        DO K = 1, INZ
          D_WORK(K, K) = B_PENTA(1, K, 2)
          IF (K.GT.1)   D_WORK(K, K-1) = B_PENTA(1, K, 1)
          IF (K.LT.INZ) D_WORK(K, K+1) = B_PENTA(1, K, 3)
        END DO
        CALL GAUSS_NOPIVOT(D_WORK, INZ)

C F_WORK(:,:,1) = D_WORK^{-1} * diag(C_PENTA(1,:))
        DO K2 = 0, INZ+1
          DO K1 = 0, INZ+1
            F_WORK(K1, K2, 1) = 0.d0
          END DO
        END DO
        DO K = 1, INZ
          F_WORK(K, K, 1) = C_PENTA(1, K)
        END DO
        CALL GAUSS_SOLVE_R(D_WORK, F_WORK(0,0,1), INZ, INZ)

C G_PENTA(1,:) <- D_WORK^{-1} * G_PENTA(1,:)
        DO K = 0, INZ+1
          GTMP(K) = G_PENTA(1, K)
        END DO
        CALL GAUSS_SOLVE_C(D_WORK, GTMP, INZ)
        DO K = 0, INZ+1
          G_PENTA(1, K) = GTMP(K)
        END DO
      END IF

C Block elimination over J=2..INY
      DO J = 2, INY
C Inflate B_PENTA(J,:,:) (tridiagonal in K) into dense D_WORK
        DO K2 = 0, INZ+1
          DO K1 = 0, INZ+1
            D_WORK(K1, K2) = 0.d0
          END DO
        END DO
        DO K = 1, INZ
          D_WORK(K, K) = B_PENTA(J, K, 2)
          IF (K.GT.1)   D_WORK(K, K-1) = B_PENTA(J, K, 1)
          IF (K.LT.INZ) D_WORK(K, K+1) = B_PENTA(J, K, 3)
        END DO

C D_WORK <- D_WORK - diag(A_PENTA(J,:)) * F_WORK(:,:,J-1)
C   (A_J is diagonal -> row-scaling of F_{J-1})
        DO K2 = 1, INZ
          DO K1 = 1, INZ
            D_WORK(K1, K2) = D_WORK(K1, K2)
     &        - A_PENTA(J, K1) * F_WORK(K1, K2, J-1)
          END DO
        END DO

        CALL GAUSS_NOPIVOT(D_WORK, INZ)

C F_WORK(:,:,J) = D_WORK^{-1} * diag(C_PENTA(J,:))
        DO K2 = 0, INZ+1
          DO K1 = 0, INZ+1
            F_WORK(K1, K2, J) = 0.d0
          END DO
        END DO
        DO K = 1, INZ
          F_WORK(K, K, J) = C_PENTA(J, K)
        END DO
        CALL GAUSS_SOLVE_R(D_WORK, F_WORK(0,0,J), INZ, INZ)

C G_PENTA(J,:) <- D_WORK^{-1} * (G_PENTA(J,:) - A_PENTA(J,:) * G_PENTA(J-1,:))
        DO K = 0, INZ+1
          GTMP(K) = G_PENTA(J, K) - A_PENTA(J, K) * G_PENTA(J-1, K)
        END DO
        CALL GAUSS_SOLVE_C(D_WORK, GTMP, INZ)
        DO K = 0, INZ+1
          G_PENTA(J, K) = GTMP(K)
        END DO
      END DO

C Send (F_block, G'_vec) at J=INY upward to RANKY+1
      IF (RANKY.NE.NPROCY-1) THEN
        IDX = 0
        DO K2 = 1, INZ
          DO K1 = 1, INZ
            IDX = IDX + 1
            SBUF(IDX) = F_WORK(K1, K2, INY)
          END DO
        END DO
        DO K = 1, INZ
          IDX = IDX + 1
          SBUF(IDX) = G_PENTA(INY, K)
        END DO
        CALL MPI_SEND(SBUF, MSGSIZE, MPI_DOUBLE_COMPLEX, RANKY+1, 15
     &               ,MPI_COMM_Y, ierror)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE BLOCK_TDMA_B_COMP_DUCT_Y_MPI(INY, INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Block-Thomas BACKWARD sweep in Y for the REM_DIV_DUCT4 pentadiagonal solver.
C Operates on /PENTA_VARS/ arrays from header. Pipeline mirrors
C TDMA_B_COMP_DUCT_Y_MPI: receive x_{INY+1} from rank above, back-substitute
C from J=INY down to J=1, then send x_2 to rank below.
C Per neighbor exchange: INZ complex doubles (the solution vector).
C On entry, G_PENTA(J,:) holds G'_J (output of forward sweep), F_WORK(:,:,J)
C holds F_J. On exit, G_PENTA(J,:) holds the solution x_J.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER INY, INZ
      INTEGER J, K, K1, K2
      COMPLEX*16 XNEXT(0:NZ+1)
      COMPLEX*16 GACC(0:NZ+1)

C Receive x_{INY+1} from RANKY+1 if it exists; otherwise INY is the last
C global block-row, so x_{INY} = G'_{INY} already in G_PENTA(INY,:) and
C no F * x_{INY+1} correction is needed.
      IF (RANKY.NE.NPROCY-1) THEN
        CALL MPI_RECV(XNEXT, INZ, MPI_DOUBLE_COMPLEX, RANKY+1, 16
     &               ,MPI_COMM_Y, status, ierror)
C Apply at J=INY: G_PENTA(INY,:) <- G_PENTA(INY,:) - F_WORK(:,:,INY) * XNEXT
        DO K = 1, INZ
          GACC(K) = G_PENTA(INY, K)
        END DO
        DO K2 = 1, INZ
          DO K1 = 1, INZ
            GACC(K1) = GACC(K1) - F_WORK(K1, K2, INY) * XNEXT(K2)
          END DO
        END DO
        DO K = 1, INZ
          G_PENTA(INY, K) = GACC(K)
        END DO
      END IF

C For J = INY-1, ..., 1 do x_J = G'_J - F_J * x_{J+1}
      DO J = INY-1, 1, -1
        DO K = 1, INZ
          XNEXT(K) = G_PENTA(J+1, K)
        END DO
        DO K = 1, INZ
          GACC(K) = G_PENTA(J, K)
        END DO
        DO K2 = 1, INZ
          DO K1 = 1, INZ
            GACC(K1) = GACC(K1) - F_WORK(K1, K2, J) * XNEXT(K2)
          END DO
        END DO
        DO K = 1, INZ
          G_PENTA(J, K) = GACC(K)
        END DO
      END DO

C Send x_2 down to RANKY-1 (mirrors scalar TDMA_B which sends G(2,K))
      IF (RANKY.NE.0) THEN
        DO K = 1, INZ
          XNEXT(K) = G_PENTA(2, K)
        END DO
        CALL MPI_SEND(XNEXT, INZ, MPI_DOUBLE_COMPLEX, RANKY-1, 16
     &               ,MPI_COMM_Y, ierror)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_F_REAL_DUCT_Z_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm FORWARD sweep in Z-DIRECTION. Input data is in
C physical (REAL) space and formatted for DUCT geometry in YZ.
C Requires communication in Z.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
     &     , G(0:INY+1,0:INZ+1)
      REAL*8 OCPACK(1:INY,4),ICPACK(1:INY,4)

      IF (RANKZ.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*INY,MPI_DOUBLE_PRECISION,RANKZ-1,14
     &               ,MPI_COMM_Z,status,ierror)
C Unpack the data
        DO J=1,INY
          A(J,1)=OCPACK(J,1)
          B(J,1)=OCPACK(J,2)
          C(J,1)=OCPACK(J,3)
          G(J,1)=OCPACK(J,4)
        END DO
! If we aren't the lowest process, start at K=2
        DO K=2,INZ
          DO J=1,INY
            A(J,K)=-A(J,K)/B(J,K-1)
            B(J,K)=B(J,K)+A(J,K)*C(J,K-1)
            G(J,K)=G(J,K)+A(J,K)*G(J,K-1)
          END DO
        END DO
      ELSE
! Here, we are the lowest process, start solving at K=1
        DO K=1,INZ
          DO J=1,INY
            A(J,K)=-A(J,K)/B(J,K-1)
            B(J,K)=B(J,K)+A(J,K)*C(J,K-1)
            G(J,K)=G(J,K)+A(J,K)*G(J,K-1)
          END DO
        END DO
      END IF

      IF (RANKZ.NE.NPROCZ-1) THEN
        DO J=1,INY
          ICPACK(J,1)=A(J,INZ)
          ICPACK(J,2)=B(J,INZ)
          ICPACK(J,3)=C(J,INZ)
          ICPACK(J,4)=G(J,INZ)
        END DO
        CALL MPI_SEND(ICPACK,4*INY,MPI_DOUBLE_PRECISION,RANKZ+1,14
     &               ,MPI_COMM_Z,ierror)
      ELSE
! Here, we are at the upper process, so solve one more row containing
! the boundary conditions
        K=INZ+1
        DO J=1,INY
          A(J,K)=-A(J,K)/B(J,K-1)
          B(J,K)=B(J,K)+A(J,K)*C(J,K-1)
          G(J,K)=G(J,K)+A(J,K)*G(J,K-1)
        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_B_REAL_DUCT_Z_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm BACKWARD sweep in Z-DIRECTION. Input data is in
C physical (REAL) space and formatted for DUCT geometry in YZ.
C Requires communication in Z.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
     &     , G(0:INY+1,0:INZ+1)
      REAL*8 OCPACK(1),ICPACK(1)

      DO J=1,INY
        IF (RANKZ.NE.NPROCZ-1) THEN
! If we aren’t the highest process, then wait for data
          CALL MPI_RECV(OCPACK,1,MPI_DOUBLE_PRECISION,RANKZ+1,15
     &                 ,MPI_COMM_Z,status,ierror)
          G(J,INZ+1)=OCPACK(1)
        ELSE
! Else, if we are the highest process, compute the solution at K=INZ
          G(J,INZ+1)=G(J,INZ+1)/B(J,INZ+1)
        END IF

! All processes solve from INZ..1
        DO K=INZ,0,-1
          IF (B(J,K).NE.0.d0) THEN
            G(J,K)=(G(J,K)-C(J,K)*G(J,K+1))/B(J,K)
          ELSE
            G(J,K)=0.d0
          END IF
        END DO

        IF (RANKZ.NE.0) THEN
          ICPACK(1)=G(J,2)
          CALL MPI_SEND(ICPACK,1,MPI_DOUBLE_PRECISION,RANKZ-1,15
     &                 ,MPI_COMM_Z,ierror)
        END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_F_COMP_DUCT_Z_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm FORWARD sweep in Z-DIRECTION. Input data is in
C Fourier (COMPLEX) space and formatted for DUCT geometry in YZ.
C Requires no communication.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)

      DO J=1,INY
C Do interior points on each process
        DO K=2,INZ
          IF (B(J,K-1).NE.0.d0) THEN
            A(J,K)=-A(J,K)/B(J,K-1)
          ELSE
            A(J,K)=0.d0
          END IF
          B(J,K)=B(J,K)+A(J,K)*C(J,K-1)
          G(J,K)=G(J,K)+A(J,K)*G(J,K-1)
        END DO
      END DO
      
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE TDMA_B_COMP_DUCT_Z_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm BACKWARD sweep in Z-DIRECTION. Input data is in
C Fourier (COMPLEX) space and formatted for DUCT geometry in YZ.
C Requires no communication.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)

      DO J=1,INY
C Compute the solution at K=INZ      
        IF (B(J,INZ).NE.0.d0) THEN
          G(J,INZ)=G(J,INZ)/B(J,INZ)
        ELSE
          G(J,INZ)=0.d0
        END IF
C Do interior points on each process
        DO K=INZ-1,1,-1
          IF (B(J,K).NE.0.d0) THEN
            G(J,K)=(G(J,K)-C(J,K)*G(J,K+1))/B(J,K)
          ELSE
            G(J,K)=0.d0
          END IF
        END DO
      END DO

      RETURN
      END










C depreciated TDMA routines below
      
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_DUCT_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Tri-Diagonal Matrix Algorithm forward sweep in Y-direction. Input data is in
C physical (real) space and formatted for duct geometry as YZ slices.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
     &     , G(0:INY+1,0:INZ+1)
      REAL*8 OCPACK(4,0:INZ+1),ICPACK(4,0:INZ+1)

!       IF (RANKY.NE.0) THEN
! C If we aren't the lowest process, then wait for data
!         CALL MPI_RECV(OCPACK,4*(INZ+2),MPI_DOUBLE_PRECISION,RANKY-1,12
!      &               ,MPI_COMM_Y,status,ierror)
! C Unpack the data
!         DO K=0,INZ+1
!           A(1,K)=OCPACK(1,K)
!           B(1,K)=OCPACK(2,K)
!           C(1,K)=OCPACK(3,K)
!           G(1,K)=OCPACK(4,K)
!         END DO
! ! If we aren't the lowest process, start at J=2
!         DO J=2,INY
!           DO K=0,INZ+1
!             A(J,K)=-A(J,K)/B(J-1,K)
!             B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
!             G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
!           END DO
!         END DO
!       ELSE
! ! Here, we are the lowest process, start solving at J=1
!         DO J=1,INY
!           DO I=0,INX-1
!           A(I,J)=-A(I,J)/B(I,J-1)
!           B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
!           G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
!           END DO
!         END DO
!       END IF

!       IF (RANKY.NE.NPROCY-1) THEN
! C If we aren't the upper process, then wait for data
!         DO I=0,INX-1
!         ICPACK(I,1)=A(I,INY)
!         ICPACK(I,2)=B(I,INY)
!         ICPACK(I,3)=C(I,INY)
!         ICPACK(I,4)=G(I,INY)
!         END DO
!         CALL MPI_SEND(ICPACK,4*INX,MPI_DOUBLE_PRECISION,RANKY+1,12
!      &               ,MPI_COMM_Y,ierror)
!       ELSE
! ! Here, we are at the upper process, so solve one more row containing
! ! the boundary conditions
!         J=INY+1
!         DO I=0,INX-1
!         A(I,J)=-A(I,J)/B(I,J-1)
!         B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
!         G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
!         END DO
!       END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_COMPLEX_DUCT_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)
      COMPLEX*16 OCPACK(4),ICPACK(4)

      DO K=1,INZ
C If we aren't the lowest process, then wait for data from below
        IF (RANKY.NE.0) THEN
          CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_COMPLEX,RANKY-1,13
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
          A(1,K)=REAL(OCPACK(1))
          B(1,K)=REAL(OCPACK(2))
          C(1,K)=REAL(OCPACK(3))
          G(1,K)=OCPACK(4)
        END IF

C Do interior points on each process
        DO J=2,INY
          IF (B(J-1,K).NE.0.d0) THEN
            A(J,K)=-A(J,K)/B(J-1,K)
          ELSE
            A(J,K)=0.d0
          END IF
          B(J,K)=B(J,K)+A(J,K)*C(J-1,K)
          G(J,K)=G(J,K)+A(J,K)*G(J-1,K)
        END DO

C If we aren't the highest process, then send data upwards
        IF (RANKY.NE.NPROCY-1) THEN
          ICPACK(1)=A(INY,K)
          ICPACK(2)=B(INY,K)
          ICPACK(3)=C(INY,K)
          ICPACK(4)=G(INY,K)
          CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_COMPLEX,RANKY+1,13
     &               ,MPI_COMM_Y,ierror)
        END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_BACKWARD_REAL_DUCT_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_DUCT_MPI(A,B,C,G,INY,INZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,INY,INZ
      REAL*8 A(0:INY+1,0:INZ+1), B(0:INY+1,0:INZ+1), C(0:INY+1,0:INZ+1)
      COMPLEX*16 G(0:INY+1,0:INZ+1)

      DO K=1,INZ
C If we aren't the highest process, then wait for data from above
        IF (RANKY.NE.NPROCY-1) THEN
          CALL MPI_RECV(G(INY+1,K),1,MPI_DOUBLE_COMPLEX,RANKY+1,11
     &               ,MPI_COMM_Y,status,ierror)
          IF (B(INY,K).NE.0.d0) THEN
            G(INY,K)=(G(INY,K)-C(INY,K)*G(INY+1,K))/B(INY,K)
          END IF
        ELSE
C If we are the highest process, then compute the solution at J=INY
          IF (B(INY,K).NE.0.d0) THEN
            G(INY,K)=G(INY,K)/B(INY,K)
          ELSE
            G(INY,K)=0.d0
          END IF
        END IF

C Do interior points on each process
        DO J=INY-1,1,-1
          IF (B(J,K).NE.0.d0) THEN
            G(J,K)=(G(J,K)-C(J,K)*G(J+1,K))/B(J,K)
          ELSE
            G(J,K)=0.d0
          END IF
        END DO

C If we aren't the lowest process, then send data down
        IF (RANKY.NE.0) THEN
          CALL MPI_SEND(G(2,K),1,MPI_DOUBLE_COMPLEX,RANKY-1,11
     &               ,MPI_COMM_Y,ierror)
        END IF
      END DO

      RETURN
      END
      


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'

      INTEGER bl(2),disp(2),types(2)

      INTEGER IPROCS,TYPE1,COMM_CART
      INTEGER DIMS(2),PERDIM(2)

      INTEGER I,J,K,XI,ZI

      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      CHARACTER*55 FNAME

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,IPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      IF (NPROCS.NE.IPROCS) THEN
         IF (RANK.EQ.0) WRITE(*,*) ' NPROCS is not equal to',
     &        ' the number of processes which we run on. '
         CALL MPI_FINALIZE(ierror)
         stop
      END IF

      IF (MOD(NPROCS,NPROCY).NE.0) THEN
         IF (RANK.EQ.0) WRITE(*,*) ' Error. NPROCS is not a ',
     &        'multiple of NPROCY'
         CALL MPI_FINALIZE(ierror)
         stop
      END IF

      DIMS(2)=NPROCY
      DIMS(1)=NPROCZ

      NXPP=int(NX/(2*NPROCZ))

      call MPI_CART_CREATE(MPI_COMM_WORLD,2,DIMS,PERDIM,.FALSE.,
     &     COMM_CART,IERROR)
      ! In PERDIM I put the information for the remain_dims
      PERDIM=(/0,1/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Y,IERROR)
      PERDIM=(/1,0/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Z,IERROR)

      call MPI_COMM_RANK(MPI_COMM_Y,RANKY,IERROR)
      call MPI_COMM_RANK(MPI_COMM_Z,RANKZ,IERROR)

      if (RANK.eq.0) write(*,'(1A,1I8)') 'MPI Initialised. NPROCS: ',
     &     NPROCS
      call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      write(*,'(1A,4I8)') 'RANK,RANKY,RANKZ: ',RANK,RANKY,RANKZ



c     ------------------------------
c     Define datatypes
c     ------------------------------

      ! Box full x to z (1)
      call MPI_TYPE_VECTOR(NZP,NXP,NX/2,
     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror)
      call MPI_TYPE_COMMIT(TYPE1,ierror)
      bl(1:2)=(/1, 1/)
      disp(1:2)= (/0, NXP*16/)
      types=(/TYPE1, MPI_UB/)

      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_1,ierror)
      call MPI_TYPE_COMMIT(XY2ZY_1,ierror)
      call MPI_TYPE_FREE(TYPE1,ierror)

      ! Box full x to z (2)
C$$$      call MPI_TYPE_VECTOR(NZP,NXP,NXP+1,
C$$$     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror)
C$$$      call MPI_TYPE_COMMIT(TYPE1,ierror)

C$$$      bl(1:2)=(/1, 1/)
C$$$      disp(1:2)= (/0, NZP*(NXP+1)*16/)
C$$$      types=(/TYPE1, MPI_UB/)

C$$$      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_2,ierror)
C$$$      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)
C$$$      call MPI_TYPE_FREE(TYPE1,ierror)

      call MPI_TYPE_VECTOR(NZP,NXP,NXP,
     &     MPI_DOUBLE_COMPLEX,XY2ZY_2,ierror)
      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)
       
c     ///////////////////////////////////
c     OTHER POSSIBLE TRANSPOSES!!!!
c     ///////////////////////////////////
c$$$      ! Box full x to z (2)
c$$$      call MPI_TYPE_VECTOR(NXP,1,NZ,
c$$$     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror)
c$$$      call MPI_TYPE_COMMIT(TYPE1,ierror)
c$$$
c$$$      bl(1:2)=(/1, 1/)
c$$$      disp(1:2)= (/0, 16/)
c$$$      types=(/TYPE1, MPI_UB/)
c$$$
c$$$      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_2,ierror)
c$$$      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)

c$$$      call MPI_TYPE_VECTOR(NZP,NXP,NXP+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY_2,ierror)
c$$$      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)
c$$$      ! Box full x to z
c$$$      call MPI_TYPE_VECTOR(NZ,NXP,NX/2+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY_1,ierr)
c$$$      call MPI_TYPE_COMMIT(XY2ZY_1,ierr)
c$$$
c$$$      ! Box full z to x
c$$$      call MPI_TYPE_VECTOR(NX/2,NX/2,NX/2+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY,ierr)
c$$$      call MPI_TYPE_COMMIT(XY2ZY,ierr)


c     ------------------------------
c     Define datatypes - for DUCT
c     ------------------------------

      ! Box full x to z (1)
      call MPI_TYPE_VECTOR(NZP+2,NXP,NX/2,
     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror)
      call MPI_TYPE_COMMIT(TYPE1,ierror)
      bl(1:2)=(/1, 1/)
      disp(1:2)= (/0, NXP*16/)
      types=(/TYPE1, MPI_UB/)

      call MPI_TYPE_STRUCT(2,bl,disp,types,X2Z_1,ierror)
      call MPI_TYPE_COMMIT(X2Z_1,ierror)
      call MPI_TYPE_FREE(TYPE1,ierror)

    !   call MPI_TYPE_VECTOR(NZP+2,NXP,NX/2,
    !  &     MPI_DOUBLE_COMPLEX,X2Z_1,ierror)
    !   call MPI_TYPE_COMMIT(X2Z_1,ierror)

      ! Box full x to z (2)
      call MPI_TYPE_VECTOR(NZP+2,NXP,NXP,
     &     MPI_DOUBLE_COMPLEX,X2Z_2,ierror)
      call MPI_TYPE_COMMIT(X2Z_2,ierror)



C Set a string to determine which input/output files to use
C When MPI is used, each process will read/write to files with the
C number of their rank (+1) appended to the end of the file.
C The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCY.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,1000)/100+48)
     &             //CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,10000)/1000+48)
     &             //CHAR(MOD(RANKY+1,1000)/100+48)
     &             //CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here
      INCLUDE 'header'

      INTEGER J, N

      IF (RANK.EQ.0)
     &     write(*,*) '*******IN INIT_CHAN_MPI *********'

      PI=4.D0*ATAN(1.D0)

! Set the upper and lower bounds for timestepping
      IF (RANKY.eq.0) THEN
        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        JEND=NY
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=2
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
      ELSE IF (RANKY.eq.NPROCY-1) THEN
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        JSTART=2
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY-1
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      ELSE
! Here, we are on a middle process
        JSTART=2
        JEND=NY
        DO N=1,N_TH
           JSTART_TH(N)=2
           JEND_TH(N)=NY
        END DO
      END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'

      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
        END IF
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Y-sweep of the duct TH ADI step.
! Only modifies Y-boundary rows (J=0, J=1 on RANKY=0; J=NY, J=NY+1
! on RANKY=NPROCY-1). Z-direction BCs are not enforced here; they
! enter the Y-sweep only through the explicit TH(K-1)/TH(K+1) halo
! reads on the RHS (handled elsewhere).
! MATL_YZ/MATD_YZ/MATU_YZ/VEC_YZ are passed by reference; their
! shapes come from the declarations in 'header'.
      INCLUDE 'header'
      INTEGER N

      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_TH_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
      ELSE IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_TH_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_TH_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Z-sweep of the duct TH ADI step.
! Only modifies Z-boundary rows (K=0, K=1 on RANKZ=0; K=NZP, K=NZP+1
! on RANKZ=NPROCZ-1).
! MATL_YZ/MATD_YZ/MATU_YZ/VEC_YZ are passed by reference; their
! shapes come from the declarations in 'header'.
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
! We need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
      INCLUDE 'header'
      INTEGER N

      IF (RANKZ.EQ.0) THEN
        CALL APPLY_BC_TH_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
      ELSE IF (RANKZ.EQ.NPROCZ-1) THEN
        CALL APPLY_BC_TH_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ,N)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Y-sweep of the duct U1 ADI step.
! U1 is cell-centered in Y -- TH-style BC convention.
      INCLUDE 'header'

      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_U1_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_U1_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U1_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Z-sweep of the duct U1 ADI step.
! U1 is cell-centered in Z -- TH-style BC convention.
      INCLUDE 'header'

      IF (RANKZ.EQ.0) THEN
        CALL APPLY_BC_U1_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKZ.EQ.NPROCZ-1) THEN
        CALL APPLY_BC_U1_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Y-sweep of the duct U2 ADI step.
! U2 is face-centered in Y -- channel-style staggered BC convention.
      INCLUDE 'header'

      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_U2_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_U2_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U2_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Z-sweep of the duct U2 ADI step.
! U2 is cell-centered in Z -- TH-style BC convention.
      INCLUDE 'header'

      IF (RANKZ.EQ.0) THEN
        CALL APPLY_BC_U2_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKZ.EQ.NPROCZ-1) THEN
        CALL APPLY_BC_U2_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_DUCT_Y_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Y-sweep of the duct U3 ADI step.
! U3 is cell-centered in Y -- TH-style BC convention.
      INCLUDE 'header'

      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_U3_YLOWER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_U3_YUPPER_DUCT_Y(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_U3_DUCT_Z_MPI(MATL_YZ,MATD_YZ,
     &         MATU_YZ,VEC_YZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! BC dispatch for the Z-sweep of the duct U3 ADI step.
! U3 is face-centered in Z -- channel-style staggered BC convention.
      INCLUDE 'header'

      IF (RANKZ.EQ.0) THEN
        CALL APPLY_BC_U3_ZLOWER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      ELSE IF (RANKZ.EQ.NPROCZ-1) THEN
        CALL APPLY_BC_U3_ZUPPER_DUCT_Z(MATL_YZ,MATD_YZ,MATU_YZ,VEC_YZ)
      END IF

      RETURN
      END


      SUBROUTINE APPLY_BC_TH_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'

      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Apply the boundary conditions
        IF (RANKY.EQ.0) THEN
        DO I=0,NXP-1
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
          IF ((K.EQ.0).AND.(I.EQ.0).AND.(RANKZ.EQ.0)) THEN
C Otherwise the matrix will be singular
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Apply the boundary conditions
        IF (RANKY.EQ.NPROCY-1) THEN
          DO I=0,NXP-1
            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END DO
        END IF

        RETURN
        END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_REM_DIV_DUCT_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K

C First, apply BC in Z
      DO J=1,NY
C Use Neumann boundary conditions, dp/dz=0 at walls
        MATL_C_YZ(J,1)=0.
        MATD_C_YZ(J,1)=1.
        MATU_C_YZ(J,1)=-1.
        VEC_C_YZ(J,1)=(0.,0.)
        MATL_C_YZ(J,NZ)=1.
        MATD_C_YZ(J,NZ)=-1.
        MATU_C_YZ(J,NZ)=0.
        VEC_C_YZ(J,NZ)=(0.,0.)
      END DO


C Apply BC in Y
C Check which processor we are on, only apply BC on top and bottom process
      IF (RANKY.EQ.0) THEN
        DO K=1,NZ
          IF ((K.EQ.1).AND.(I.EQ.0).AND.(RANKZ.EQ.0)) THEN
C Use homogeneous Dirichlet BC for kx=0 component at bottom YZ corner
C Otherwise the matrix will be singular
            MATL_C_YZ(1,K)=0.
            MATD_C_YZ(1,K)=1.
            MATU_C_YZ(1,K)=0.
            VEC_C_YZ(1,K)=(0.,0.)
          ELSE
C Use Neumann boundary conditions, dp/dy=0 at walls
            MATL_C_YZ(1,K)=0.
            MATD_C_YZ(1,K)=1.
            MATU_C_YZ(1,K)=-1.
            VEC_C_YZ(1,K)=(0.,0.)
          END IF
        END DO
      ELSEIF (RANKY.EQ.NPROCY-1) THEN
        DO K=1,NZ
          MATL_C_YZ(NY,K)=1.
          MATD_C_YZ(NY,K)=-1.
          MATU_C_YZ(NY,K)=0.
          VEC_C_YZ(NY,K)=(0.,0.)
        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_REM_DIV_DUCT_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K

C First, apply BC in Z
      ! DO J=1,NY
      !   VEC_C_YZ(J,1)=(0.,0.)
      !   VEC_C_YZ(J,NZ)=(0.,0.)
      ! END DO

C Apply BC in Y
C Check which processor we are on, only apply BC on top and bottom process
      IF (RANKY.EQ.0) THEN
        DO K=1,NZ
          IF ((I.EQ.0).AND.(RANKZ.EQ.0)) THEN
C Use homogeneous Dirichlet BC for kx=0 component at bottom YZ corner
C Otherwise the matrix will be singular
            MATL_C_YZ(1,K)=0.
            MATD_C_YZ(1,K)=1.
            MATU_C_YZ(1,K)=0.
            VEC_C_YZ(1,K)=(0.,0.)
          ELSE
C Use Neumann boundary conditions, dp/dz=0 at walls
            MATL_C_YZ(1,K)=0.
            MATD_C_YZ(1,K)=1.
            MATU_C_YZ(1,K)=-1.
            VEC_C_YZ(1,K)=(0.,0.)
          END IF
        END DO
      ELSEIF (RANKY.EQ.NPROCY-1) THEN
        DO K=1,NZ
          MATL_C_YZ(NY,K)=1.
          MATD_C_YZ(NY,K)=-1.
          MATU_C_YZ(NY,K)=0.
          VEC_C_YZ(NY,K)=(0.,0.)
        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_REM_DIV_DUCT6_Y_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Apply modal Y-boundary conditions for REM_DIV_DUCT6.
C Z Neumann boundary conditions are handled by the DCT-I diagonalization.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,M

C Lower Y wall. Pin only the global null mode: kx=0 and z-constant.
      IF (RANKY.EQ.0) THEN
        DO M=1,NZ
          IF ((I.EQ.0).AND.(RANKZ.EQ.0).AND.(M.EQ.1)) THEN
            MATL_C_YZ(1,M)=0.d0
            MATD_C_YZ(1,M)=1.d0
            MATU_C_YZ(1,M)=0.d0
            VEC_C_YZ(1,M)=(0.d0,0.d0)
          ELSE
            MATL_C_YZ(1,M)=0.d0
            MATD_C_YZ(1,M)=1.d0
            MATU_C_YZ(1,M)=-1.d0
            VEC_C_YZ(1,M)=(0.d0,0.d0)
          END IF
        END DO
      END IF

C Upper Y wall.
      IF (RANKY.EQ.NPROCY-1) THEN
        DO M=1,NZ
          MATL_C_YZ(NY,M)=1.d0
          MATD_C_YZ(NY,M)=-1.d0
          MATU_C_YZ(NY,M)=0.d0
          VEC_C_YZ(NY,M)=(0.d0,0.d0)
        END DO
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_REM_DIV_DUCT4_MPI(I)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Apply boundary conditions to the /PENTA_VARS/ block-tridiagonal system
C for REM_DIV_DUCT4 at kx mode I.
C
C Z-walls (Neumann, every J on every rank):
C   Row K=1  of B_J: phi(J,1)  - phi(J,2)    = 0  -> diag=1, super=-1, sub=0
C   Row K=NZ of B_J: phi(J,NZ-1) - phi(J,NZ) = 0  -> sub=1, diag=-1, super=0
C   A_J(K=1)=A_J(K=NZ)=C_J(K=1)=C_J(K=NZ)=0; G_J(K=1)=G_J(K=NZ)=0.
C
C Y-floor (Neumann, J=1 on RANKY=0):
C   Block-row 1 entirely overwritten: B_1=I, A_1=0, C_1=-I, G_1=0.
C
C Y-ceiling (Neumann, J=NY on RANKY=NPROCY-1):
C   Block-row NY: B_NY=-I, A_NY=I, C_NY=0, G_NY=0.
C
C kx=0 anchor (Dirichlet phi=0 at the J=1, K=1 corner) on
C RANKY=0 .AND. RANKZ=0 .AND. I=0:
C   B_1(K=1,:)=identity row, A_1(1)=C_1(1)=0, G_1(1)=0.
C   This breaks the otherwise-singular pure-Neumann system for kx=0.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I, J, K

C Z-walls applied to every active block-row J=1..NY: rows K=1 and K=NZ
C of each B_J block. Halo block-rows (J=0, J=NY+1) are never consumed by
C BLOCK_TDMA_F/B, so we don't touch them here.
      DO J = 1, NY
C K=1 row: phi(J,1) - phi(J,2) = 0
        B_PENTA(J, 1, 1) = 0.d0
        B_PENTA(J, 1, 2) = 1.d0
        B_PENTA(J, 1, 3) = -1.d0
        A_PENTA(J, 1)    = 0.d0
        C_PENTA(J, 1)    = 0.d0
        G_PENTA(J, 1)    = (0.d0, 0.d0)
C K=NZ row: phi(J,NZ-1) - phi(J,NZ) = 0
        B_PENTA(J, NZ, 1) = 1.d0
        B_PENTA(J, NZ, 2) = -1.d0
        B_PENTA(J, NZ, 3) = 0.d0
        A_PENTA(J, NZ)    = 0.d0
        C_PENTA(J, NZ)    = 0.d0
        G_PENTA(J, NZ)    = (0.d0, 0.d0)
      END DO

C Y-floor (RANKY=0): replace block-row J=1 with B_1=I, A_1=0, C_1=-I, G_1=0
      IF (RANKY.EQ.0) THEN
        DO K = 1, NZ
          B_PENTA(1, K, 1) = 0.d0
          B_PENTA(1, K, 2) = 1.d0
          B_PENTA(1, K, 3) = 0.d0
          A_PENTA(1, K)    = 0.d0
          C_PENTA(1, K)    = -1.d0
          G_PENTA(1, K)    = (0.d0, 0.d0)
        END DO
      END IF

C Y-ceiling (RANKY=NPROCY-1): replace J=NY with B_NY=-I, A_NY=I, C_NY=0, G_NY=0
      IF (RANKY.EQ.NPROCY-1) THEN
        DO K = 1, NZ
          B_PENTA(NY, K, 1) = 0.d0
          B_PENTA(NY, K, 2) = -1.d0
          B_PENTA(NY, K, 3) = 0.d0
          A_PENTA(NY, K)    = 1.d0
          C_PENTA(NY, K)    = 0.d0
          G_PENTA(NY, K)    = (0.d0, 0.d0)
        END DO
      END IF

C kx=0 corner anchor (I=0, RANKY=0, RANKZ=0): Dirichlet phi=0 at (J=1, K=1).
C The Y-floor block already set B_1(K=1,:) to the identity row above; here
C we just need to zero the C_PENTA(1,1) entry (was -1) so the corner is
C truly decoupled from J=2.
      IF ((I.EQ.0).AND.(RANKY.EQ.0).AND.(RANKZ.EQ.0)) THEN
        B_PENTA(1, 1, 1) = 0.d0
        B_PENTA(1, 1, 2) = 1.d0
        B_PENTA(1, 1, 3) = 0.d0
        A_PENTA(1, 1)    = 0.d0
        C_PENTA(1, 1)    = 0.d0
        G_PENTA(1, 1)    = (0.d0, 0.d0)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_REM_DIV_DUCT_Z_MPI(MATL_C_YZ,MATD_C_YZ,
     &         MATU_C_YZ,VEC_C_YZ,I)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K

C First, apply BC in Y (including numerical processor boundary??)
      ! IF (RANKY.EQ.0) THEN
      !   DO K=1,NZ
      !     VEC_C_YZ(1,K)=(0.,0.)
      !   END DO
      ! ELSEIF (RANKY.EQ.NPROCY-1) THEN
      !   DO K=1,NZ
      !     VEC_C_YZ(NY,K)=(0.,0.)
      !   END DO
      ! END IF
      ! DO K=1,NZ
      !   VEC_C_YZ(1,K)=(0.,0.)
      !   VEC_C_YZ(NY,K)=(0.,0.)
      ! END DO

C Apply BC in Z
      DO J=1,NY
        IF ((I.EQ.0).AND.(RANKZ.EQ.0)) THEN
C Use homogeneous Dirichlet BC for kx=0 component at bottom YZ corner
C Otherwise the matrix will be singular
          MATL_C_YZ(J,1)=0.
          MATD_C_YZ(J,1)=1.
          MATU_C_YZ(J,1)=0.
          VEC_C_YZ(J,1)=(0.,0.)
        ELSE
C Use Neumann boundary conditions, dp/dz=0 at walls
          MATL_C_YZ(J,1)=0.
          MATD_C_YZ(J,1)=1.
          MATU_C_YZ(J,1)=-1.
          VEC_C_YZ(J,1)=(0.,0.)
        END IF
      END DO
      DO J=1,NY
        MATL_C_YZ(J,NZ)=1.
        MATD_C_YZ(J,NZ)=-1.
        MATU_C_YZ(J,NZ)=0.
        VEC_C_YZ(J,NZ)=(0.,0.)
      END DO

      RETURN
      END







      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANKY.eq.0) THEN
        DO I=0,NXP-1
          IF ((I.EQ.0).AND.(K.EQ.0).AND.(RANKZ.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANKY.eq.NPROCY-1) THEN
        DO I=0,NXP-1
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
        END DO
        END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_VEL_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      IF (NUM_PER_DIR.EQ.2) THEN
        IF (RANKY.EQ.0) THEN
          CALL APPLY_BC_VEL_LOWER
        ELSE IF (RANKY.EQ.NPROCY-1) THEN
          CALL APPLY_BC_VEL_UPPER
        END IF
      ELSE IF (NUM_PER_DIR.EQ.1) THEN
        ! in duct flow Fourier space, Z data is not split between procs
        CALL APPLY_BC_VEL_ZLOWER_DUCT
        CALL APPLY_BC_VEL_ZUPPER_DUCT
        ! apply Y BC second
        IF (RANKY.EQ.0) THEN
          CALL APPLY_BC_VEL_YLOWER_DUCT
        ELSE IF (RANKY.EQ.NPROCY-1) THEN
          CALL APPLY_BC_VEL_YUPPER_DUCT
        END IF
      END IF

      RETURN
      END


      SUBROUTINE APPLY_BC_SCALAR_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! Apply Boundary conditions to velocity field
      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_SCALAR_LOWER
      END IF
      IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_SCALAR_UPPER
      END IF

      RETURN
      END

!       SUBROUTINE APPLY_BC_TH_PHYS_MPI
! ! This subroutine applies the boundary conditions for the Poisson Eq.
! ! Note, MATL, MATD, etc. are dimensioned in header
!       INCLUDE 'header'

! ! Apply Boundary conditions to velocity field
!       IF (RANKY.EQ.0) THEN
!         CALL APPLY_BC_TH_PHYS_LOWER
!       END IF
!       IF (RANKY.EQ.NPROCY-1) THEN
!         CALL APPLY_BC_TH_PHYS_UPPER
!       END IF

!       RETURN
!       END

      SUBROUTINE APPLY_BC_VEL_PHYS_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! Apply Boundary conditions to velocity field
      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_VEL_PHYS_LOWER
      END IF
      IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_VEL_PHYS_UPPER
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_DUCT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize any constants here
      INCLUDE 'header'

      INTEGER I,J,K,N

      IF (RANK.EQ.0)
     &     write(*,*) '*******IN INIT_DUCT_MPI *********'

      PI=4.D0*ATAN(1.D0)

! Y - Set the lower bounds for timestepping
      IF (RANKY.eq.0) THEN
        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        JEND=NY
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=2
        ELSE
          JSTART=2
        END IF
! Y - Set the lower bounds for timestepping of the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
! Y - Set the upper bounds for timestepping
      ELSE IF (RANKY.eq.NPROCY-1) THEN
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        JSTART=2
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY-1
        ELSE
          JEND=NY-1
        END IF
! Y - Set the upper bounds for timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO
! Y - Here, we are on a middle process
      ELSE
        JSTART=2
        JEND=NY
        DO N=1,N_TH
           JSTART_TH(N)=2
           JEND_TH(N)=NY
        END DO
      END IF




! Z - Set the lower bounds for timestepping
      IF (RANKZ.eq.0) THEN
        write(*,*) 'U_BC_ZMIN: ',U_BC_ZMIN
        KEND=NZP
        IF (U_BC_ZMIN.EQ.0) THEN
          KSTART=2
        ELSE IF (U_BC_ZMIN.EQ.1) THEN
          KSTART=2
        ELSE
          KSTART=2
        END IF
! Z - Set the lower bounds for timestepping of the scalar equations
        DO N=1,N_TH
          KEND_TH(N)=NZP
          IF (TH_BC_ZMIN(N).EQ.0) THEN
            KSTART_TH(N)=2
          ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
            KSTART_TH(N)=1
          ELSE
            KSTART_TH(N)=2
          END IF
        END DO
! Z - Set the upper bounds for timestepping
      ELSE IF (RANKZ.eq.NPROCZ-1) THEN
        write(*,*) 'U_BC_ZMAX: ',U_BC_ZMAX
        KSTART=2
        IF (U_BC_ZMAX.EQ.0) THEN
          KEND=NZP-1
        ELSE IF (U_BC_ZMAX.EQ.1) THEN
          KEND=NZP-1
        ELSE
          KEND=NZP-1
        END IF
! Z - Set the upper bounds for timestepping of the scalar equations
        DO N=1,N_TH
        KSTART_TH(N)=2
        IF (TH_BC_ZMAX(N).EQ.0) THEN
          KEND_TH(N)=NZP-1
        ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
          KEND_TH(N)=NZP
        ELSE
          KEND_TH(N)=NZP-1
        END IF
        END DO
! Z - Here, we are on a middle process
      ELSE
        KSTART=2
        KEND=NZP
        DO N=1,N_TH
           KSTART_TH(N)=2
           KEND_TH(N)=NZP
        END DO
      END IF

C When variables are in Fourier space, data is not divided in Z.
C We need lower and upper bounds that are consistent across all
C processors - denoted by kstartC, kendC, etc.
      IF (U_BC_ZMIN.EQ.0) THEN
        KSTARTC=2
      ELSE IF (U_BC_ZMIN.EQ.1) THEN
        KSTARTC=2
      ELSE
        KSTARTC=2
      END IF
      DO N=1,N_TH
        IF (TH_BC_ZMIN(N).EQ.0) THEN
          KSTARTC_TH(N)=2
        ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
          KSTARTC_TH(N)=1
        ELSE
          KSTARTC_TH(N)=2
        END IF
      END DO

      IF (U_BC_ZMAX.EQ.0) THEN
        KENDC=NZ-1
      ELSE IF (U_BC_ZMAX.EQ.1) THEN
        KENDC=NZ-1
      ELSE
        KENDC=NZ-1
      END IF
      DO N=1,N_TH
        IF (TH_BC_ZMAX(N).EQ.0) THEN
          KENDC_TH(N)=NZ-1
        ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
          KENDC_TH(N)=NZ
        ELSE
          KENDC_TH(N)=NZ-1
        END IF
      END DO

      RETURN
      END


      SUBROUTINE TRANSPOSE_MPI_XZ_TO_XY()
! This subroutine starts with all arrays decomposed in x-z slabs
! and transposes the data so that it is decomposed in x-y slabs
! x-y slabs.
c$$$      include 'header'
c$$$
c$$$      integer i,j,k,N
c$$$
c$$$      real*8 A(0:NX+1,0:NZ+1,0:NY+1)
c$$$      real*8 buffer(0:NX+1,0:NY+1,0:NZ+1)
c$$$
c$$$      real*8 test1(1:NX,1:NY,1:NZ)
c$$$      real*8 test2(1:NX,1:NY,1:NZ)
c$$$
c$$$       do i=0,NXM
c$$$       do j=1,NY
c$$$         write(100+RANK,*) 'i,j,U1: ',i,j,A(I,0,J)
c$$$       end do
c$$$       end do
c$$$
c$$$      do k=0,NZ-1
c$$$        do j=1,NY
c$$$          do i=0,NX-1
c$$$            buffer(i,j,k)=A(i,k,j)
c$$$            test1(i+1,j,k+1)=A(i,k,j)
c$$$          end do
c$$$        end do
c$$$      end do
c$$$
c$$$      CALL MPI_AllToAll(test1,(NX)*(NY)*(NZ)/NPROCS
c$$$     &   ,MPI_DOUBLE_PRECISION,test2,(NX)*(NY)*(NZ)/NPROCS
c$$$     &       ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
c$$$      CALL MPI_Barrier(MPI_COMM_WORLD)
c$$$
c$$$! A should now be indexed by A(NX,NZ/nprocs,NY*nprocs)
c$$$
c$$$         n=1
c$$$         do k=0,NZ/NPROCS-1
c$$$         do i=0,NX-1
c$$$         do j=1,NY
c$$$           A(i,k,j)=test2(i+1,j,k+1)
c$$$         end do
c$$$         end do
c$$$         end do
c$$$
c$$$         do n=2,NPROCS
c$$$         do k=0,NZ/NPROCS-1
c$$$         do i=0,NX-1
c$$$         do j=2,NY
c$$$           A(i,k,NY+j-1)=test2(i+1,NY+j,k+1)
c$$$         end do
c$$$         end do
c$$$         end do
c$$$         end do
c$$$
c$$$         NX_t=NX
c$$$         NZ_t=NZ/NPROCS
c$$$         NY_t=NY*NPROCS-(NPROCS-1)
c$$$
c$$$! A should now be indexed from A(0:NX_t-1,0:NZ_t/NPROCS-1,1:NY*NPROCS-(NPROCS-1))
c$$$! (The new wall locations are 1 and NY-(NPROCS-1) )
c$$$
c$$$          do i=0,NXM
c$$$          do j=1,NY_T
c$$$            write(200+RANK,*) 'i,j,U1: ',i,j,A(i,0,J)
c$$$          end do
c$$$          end do
c$$$


      RETURN
      END




      SUBROUTINE get_minimum_mpi(val)

      include 'header'

      REAL*8 val,vmin

      call MPI_ALLREDUCE(val,vmin,1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN,MPI_COMM_WORLD,ierror)

      val=vmin

      RETURN

      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER(V,VV,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP1A(0:NX/2-1,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP2(0:NXP-1,0:(NZP+2)*NPROCZ-1)

      INTEGER I,J,K, N
      INTEGER JMIN,JMAX
      ! CHARACTER*55 FNAME

      
      ! FFT in X
      DO J=JMIN,JMAX
        CALL DFFTW_EXECUTE_DFT_R2C(FFTW_X_TO_F_PLAN,
     *     V(0,0,J), TMP(0,0,J))
        DO K=0,NZP+1
          DO I=0,NKX
            TMP(I,K,J)=TMP(I,K,J)/NX
            TMP1A(I,K,J)=TMP(I,K,J)
          END DO
          DO I=NKX+1,NX/2
            TMP(I,K,J)=cmplx(0.d0,0.d0)
          END DO
          DO I=NKX+1,NX/2-1
            TMP1A(I,K,J)=TMP(I,K,J)
          END DO
        END DO
      END DO

      ! WRITE(*,*) 'writing TMP1A on rank ', RANK  
      ! FNAME='TMP1A_'//char(RANK+48)//'.dat'
      ! OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      ! WRITE(15,151) real(TMP1A)
      ! CLOSE(15)

    !   WRITE(*,*) 'FFT CU2;   RANK',RANK,
    !  &  '; VV loc 1',LOC(VV(0,0,0)),'; VV loc 2',LOC(VV(NXP,NZ+1,NY+1))


    !   WRITE(*,*) 'writing TMP2 on rank ', RANK  
    !   FNAME='TMP2_'//char(RANK+48)//'.dat'
    !   OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
      
      

      ! reconfigure and send
      DO J=JMIN,JMAX
        CALL MPI_ALLTOALL(TMP1A(0,0,J),1,X2Z_1,TMP2,1,X2Z_2,
     &     MPI_COMM_Z,IERROR)

      ! WRITE(15,151) real(TMP2)

        ! first Z processor
        DO K=0,NZP
          DO I=0,NXP-1
            VV(I,K,J)=TMP2(I,K)
          END DO
          VV(NXP,K,J)=cmplx(0.d0,0.d0)
        END DO

        ! last Z processor
        DO K=0,NZP-1
          DO I=0,NXP-1
            VV(I,(NZP-1)*(NPROCZ-1)+2+K,J)
     &          =TMP2(I,(NZP+2)*(NPROCZ-1)+2+K)
          END DO
          VV(NXP,(NZP-1)*(NPROCZ-1)+2+K,J)=cmplx(0.d0,0.d0)
    !       IF (RANK.EQ.0) WRITE(*,*) 'RANK ',RANK,';  J ',J,
    !  &     ';  K ',K,';  ind ',(NZP-1)*(NPROCZ-1)+2+K
    !       IF (RANK.EQ.0) WRITE(*,*) '             '
    !  &     ,(NZP+2)*(NPROCZ-1)+2+K
        END DO
        
        ! intermediate Z processors (if NPROCZ>2)
        DO N=1,NPROCZ-2
          DO K=0,NZP-2
            DO I=0,NXP-1
              VV(I,(NZP-1)*N+2+K,J)=TMP2(I,(NZP+2)*N+2+K)
            END DO
            VV(NXP,(NZP-1)*N+2+K,J)=cmplx(0.d0,0.d0)
          END DO
        END DO

      END DO
!       CLOSE(15)

!       WRITE(*,*) 'writing VV on rank ', RANK  
!       FNAME='VV_'//char(RANK+48)//'.dat'
!       OPEN(UNIT=15,FILE=FNAME,STATUS="UNKNOWN",FORM="FORMATTED")
!       WRITE(15,151) real(VV)
!       CLOSE(15)

! 151   format(F30.20)

      END SUBROUTINE


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL(VV,V,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP1A(0:NX/2-1,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP2(0:NXP-1,0:(NZP+2)*NPROCZ-1)

      INTEGER I,J,K, N
      INTEGER JMIN,JMAX


      ! reconfigure and send
      DO J=JMIN,JMAX
        DO I=0,NXP-1

          ! first Z processor
          DO K=0,NZP
            TMP2(I,K)=VV(I,K,J)
          END DO
          TMP2(I,NZP+1)=VV(I,NZP+1,J)

          ! last Z processor
          DO K=0,NZP-1
            TMP2(I,(NZP+2)*(NPROCZ-1)+2+K)
     &         =VV(I,(NZP-1)*(NPROCZ-1)+2+K,J)
          END DO
          TMP2(I,(NZP+2)*(NPROCZ-1))=VV(I,(NZP-1)*(NPROCZ-1),J)
          TMP2(I,(NZP+2)*(NPROCZ-1)+1)=VV(I,(NZP-1)*(NPROCZ-1)+1,J)

          ! intermediate Z processors (if NPROCZ>2)
          DO N=1,NPROCZ-2
            DO K=0,NZP-2
              TMP2(I,(NZP+2)*N+2+K)=VV(I,(NZP-1)*N+2+K,J)
            END DO
            TMP2(I,(NZP+2)*N)=VV(I,(NZP-1)*N,J)
            TMP2(I,(NZP+2)*N+1)=VV(I,(NZP-1)*N+1,J)
            TMP2(I,(NZP+2)*N+NZP+1)=VV(I,(NZP-1)*N+NZP+1,J)
          END DO

        END DO
        CALL MPI_ALLTOALL(TMP2,1,X2Z_2,TMP1A(0,0,J),1,X2Z_1,
     &       MPI_COMM_Z,IERROR)
      END DO

      ! FFT in X
      DO J=JMIN,JMAX
        DO K=0,NZP+1
          DO I=0,NKX
            TMP(I,K,J)=TMP1A(I,K,J)
          END DO
          DO I=NKX+1,NX/2
            TMP(I,K,J)=cmplx(0.d0,0.d0)
          END DO
        END DO
        CALL DFFTW_EXECUTE_DFT_C2R(FFTW_X_TO_P_PLAN,
     *       TMP(0,0,J), V(0,0,J))
      END DO

      END SUBROUTINE


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_FOURIER(V,VV,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP2(0:NXP-1,0:NZ-1)
      COMPLEX*16 TMP_1A(0:NX/2-1,0:NZP-1,0:NY+1)

      INTEGER I,J,K
      INTEGER JMIN,JMAX

      ! FFT in X
      DO J=JMIN,JMAX
        CALL DFFTW_EXECUTE_DFT_R2C(FFTW_X_TO_F_PLAN,
     *       V(0,0,J), TMP(0,0,J))
        DO K=0,NZP-1
          DO I=0,NKX
            TMP(I,K,J)=TMP(I,K,J)/NX
            TMP_1A(I,K,J)=TMP(I,K,J)
          END DO
          DO I=NKX+1,NX/2
            TMP(I,K,J)=cmplx(0.d0,0.d0)
          END DO
          DO I=NKX+1,NX/2-1
            TMP_1A(I,K,J)=TMP(I,K,J)
          END DO
        END DO
      END DO

      DO J=JMIN,JMAX
        call mpi_alltoall(TMP_1A(0,0,J),1,XY2ZY_1,TMP2,1,XY2ZY_2,
     &       MPI_COMM_Z,IERROR)
        DO K=0,NZ-1
          DO I=0,NXP-1
            VV(I,K,J)=TMP2(I,K)
          END DO
          VV(NXP,K,J)=cmplx(0.d0,0.d0)
        END DO
        DO K=NZ,NZ+1
          DO I=0,NXP
            VV(I,K,J)=cmplx(0.d0,0.d0)
          END DO 
        END DO
      END DO

      ! FFT in Z
      DO J=JMIN,JMAX
        CALL DFFTW_EXECUTE_DFT_R2C(FFTW_Z_TO_F_PLAN,
     *       VV(0,0,J), VV(0,0,J))
        DO K=0,NKZ
          DO I=0,NXP
            VV(I,K,J)=VV(I,K,J)/NZ
          END DO
        END DO
        DO K=1,NKZ
          DO I=0,NXP
            VV(I,NKZ+K,J)=VV(I,NZ-1+K-NKZ,J)/NZ
          END DO
        END DO
      END DO

      END SUBROUTINE



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_PHYSICAL(VV,V,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 direction) planes JMIN-JMAX to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      COMPLEX*16 TMP2(0:NXP-1,0:NZ-1)
      COMPLEX*16 TMP_1A(0:NX/2-1,0:NZP-1,0:NY+1)

      INTEGER I,J,K
      INTEGER JMIN,JMAX

      ! FFT in Z
      DO J=JMIN,JMAX
        ! UNPACK
        DO K=NKZ,1,-1
          DO I=0,NXP
            VV(I,NZ-1+K-NKZ,J)=VV(I,NKZ+K,J)
          END DO
        END DO
        DO K=NKZ+1,NZ-NKZ-1
          DO I=0,NXP
            VV(I,K,J)=cmplx(0.d0,0.d0)
          END DO
        END DO
        CALL DFFTW_EXECUTE_DFT(FFTW_Z_TO_P_PLAN,
     *       VV(0,0,J), VV(0,0,J))
      END DO

      DO J=JMIN,JMAX
        DO K=0,NZ-1
          DO I=0,NXP-1
            TMP2(I,K)=VV(I,K,J)
          END DO
        END DO
        call mpi_alltoall(TMP2,1,XY2ZY_2,TMP_1A(0,0,J),1,XY2ZY_1,
     &       MPI_COMM_Z,IERROR)
        DO K=0,NZP-1
          DO I=0,NX/2-1
            TMP(I,K,J)=TMP_1A(I,K,J)
          END DO
          TMP(NX/2,K,J)=cmplx(0.d0,0.d0)
        END DO
        DO K=NZP,NZP+1
          DO I=0,NX/2
            TMP(I,K,J)=cmplx(0.d0,0.d0)
          END DO
        END DO
      END DO

      ! FFT in X
      DO J=JMIN,JMAX
        CALL DFFTW_EXECUTE_DFT_C2R(FFTW_X_TO_P_PLAN,
     *       TMP(0,0,J), V(0,0,J))
      END DO

      END SUBROUTINE


      
      SUBROUTINE INTEGRATE_Y_VAR(VAR,RES,COMM)

      INCLUDE 'header'

      INTEGER J,COMM
      REAL*8 VAR(0:NY+1),RES

! Integrat the instantaneous mean profile numerically at GY points
      RES=0.
      IF (USE_MPI) THEN
         do j=2,NY
            RES=RES+0.5*(VAR(j)+VAR(j-1))*DY(j)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,RES,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,COMM,ierror)
      ELSE
         do j=2,NY
            RES=RES+0.5*(VAR(j)+VAR(j-1))*DY(j)
         end do
      END IF
      RES=RES/LY

      END SUBROUTINE

      SUBROUTINE INTEGRATE_Z_VAR(VAR,RES,COMM)

      INCLUDE 'header'

      INTEGER i,k,j,COMM
      REAL*8 VAR(0:NX+1,0:NZP+1,0:NY+1),RES(0:NXM,1:NY)

! Integrat the instantaneous mean profile numerically at GY points
      IF (USE_MPI) THEN
        do i=0,NXM
        do j=1,NY
        RES(i,j)=0.
         do k=0,NZP
            RES(i,j)=RES(i,j)+VAR(i,k,j)*DZ(k)
         end do
        end do
        end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,RES,NX*NY,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,COMM,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      ELSE
      do i=0,NXM
      do j=1,NY
         do k=0,NZM
            RES(i,j)=RES(i,j)+VAR(i,k,j)*DZ(k)
         end do
      end do
      end do
      END IF
      RES=RES/LZ

      END SUBROUTINE


      SUBROUTINE END_RUN_MPI(FLAG)

      INCLUDE 'header'

      LOGICAL FLAG

      IF (RANK.EQ.0) THEN
         CALL END_RUN(FLAG)
      END IF
      CALL MPI_BCAST(FLAG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)

      END
