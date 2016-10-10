      PROGRAM MD
      IMPLICIT NONE
      CHARACTER*100 Buffer
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTCOORDS
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTFORCES
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTACCELS
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTVELS
      INTEGER, DIMENSION(:,:),  ALLOCATABLE :: CellSizes
      INTEGER  nCellRows,NCellCols,CellSize,NumIterations
      INTEGER  nParticles,Dimension,nPartPerCell,CellDim
      REAL*8 ran,centerx,centery
      REAL*8 rx1,ry1,rx,ry,r,rt
      REAL*8 TotPot,TotKin
      REAL*8 ARmass,AReps,ARsigma
      REAL*8 timestart,timeend
      INTEGER  nPartsPerCell,MAXPART,MINPART
      INTEGER TotalCells,TotSize
      INTEGER i,j,k,l,cx1,cy1,cx2,cy2,t
      INTEGER maxCellSize
      PARAMETER(CellDim=12.0)
      PARAMETER(NPartPerCell=2)
      PARAMETER(ARMass=39.948)


      IF(IARGC() < 2) THEN
         WRITE(6,*)"Incorrect syntax: should be two arguments"
         CALL EXIT(2)
      ENDIF

      CALL GETARG(1,Buffer)
      READ(Buffer,*)nParticles
      CALL GETARG(2,Buffer)
      READ(Buffer,*)NumIterations

      ! Lets Allocate the number of cells
      Dimension= int(Sqrt(NParticles/real(NPartPerCell)))
      TotalCells = Dimension**2
      NCellRows = Dimension
      NCellCols = Dimension


      ALLOCATE(CellSizes(NCellCols,NCellRows))
      CALL ComputeAtomsPerCell(CellSizes,NCellCols,NCellRows,
     $     NParticles,maxCellSize)

      WRITE(6,'(A,I6,A,I6,A,I6,A)')' The Total Number of Cells is ',
     $     TotalCells,' With ',maxCellSize,' particles per cell, and ',
     $     nParticles,' particles total in system'

!OK, now we can allocate the Particle Matrix
      ALLOCATE(PARTCOORDS(2,maxCellSize,NCellCols,NCellRows))
      ALLOCATE(PARTFORCES(2,maxCellSize,NCellCols,NCellRows))
      ALLOCATE(PARTACCELS(2,maxCellSize,NCellCols,NCellRows))
      ALLOCATE(PARTVELS(2,maxCellSize,NCellCols,NCellRows))
!     Zero all of the arrays
      DO i=1,NCellRows
         do j=1,NCellCols
            DO k=1,maxCellSize
               do l=1,2
                  PartCoords(l,k,j,i) = 0.0
                  PartForces(l,k,j,i) = 0.0
                  PartAccels(l,k,j,i) = 0.0
                  PartVels(l,k,j,i)   = 0.0
               enddo
            ENDDO
         enddo
      ENDDO

!     WE will be filling the cells, making sure than
!     No atom is less than 4 Angstroms from another

      DO i = 1,NCellRows
         DO j = 1,NCellCols
            Centerx = CellDim*real(j) + 0.5*CellDim
            Centery = CellDim*real(i) + 0.5*CellDim
            call random_number(ran)
            PartCoords(1,1,j,i) = centerx+((ran-0.5)*(CellDim-2))
            call random_number(ran)
            PartCoords(2,1,j,i) = centery+((ran-0.5)*(CellDim-2))
            DO k=2,CellSizes(j,i)
               R=0
               DO WHILE (R < 4.0)
                  R=4.0001
                  call random_number(ran)
                  RX1 = centerx+((ran-0.5)*(CellDim-2))
                  call random_number(ran)
                  RY1 = centery+((ran-0.5)*(CellDim-2))
!     check to make sure it is far enough away from the
!     others in the cell
                  DO l=1,k-1
                     RX=RX1-PartCoords(1,l,j,i)
                     RY=RY1-PartCoords(2,l,j,i)
                     RT = RX**2+RY**2
                     IF(RT < R) THEN
                        R=RT
                     ENDIF
                  ENDDO
               ENDDO
               PartCoords(1,k,j,i) = RX1
               PartCoords(2,k,j,i) = RY1
            ENDDO
         ENDDO
      ENDDO
      CALL WallTime(timestart)
      ! Lets Start Iterating
      DO t=1,NumIterations
         ! Zero Energy variables
         TotPot = 0.0
         TotKin = 0.0
         ! Loop over Cells
         DO CX1 = 1, NCellRows
            DO CY1 = 1,NCellCols
!     Cells must be neighbors
! Let's do the ones within a cell first
               DO I=1,CellSizes(CY1,CX1)
                  DO J=1,CellSizes(CY1,CX1)
                     if(i.ne.j)then
                        CALL INTERACT(
     $                       PartCoords(1,j,cy1,cx1),
     $                       PartCoords(2,j,cy1,cx1),
     $                       PartCoords(1,i,cy1,cx1),
     $                       PartCoords(2,i,cy1,cx1),
     $                       PartForces(1,j,cy1,cx1),
     $                       PartForces(2,j,cy1,cx1),
     $                       CellSize,TotPot)
                     ENDIF
                  ENDDO
               ENDDO
! Lets do the other cells
               DO CX2 = 1,NCellRows
                  DO CY2 = 1,NCellCols
                     if(((ABS(CX1-CX2) < 2).and.(ABS(CY1-CY2) < 2)).and.
     $                    ((CX1.ne.CX2).or.(CY1.ne.CY2))) THEN
                        DO i=1,CellSizes(CY2,CX2)
                           DO j=1,CellSizes(CY1,CX1)
                              CALL INTERACT(
     $                             PartCoords(1,j,cy1,cx1),
     $                             PartCoords(2,j,cy1,cx1),
     $                             PartCoords(1,i,cy2,cx2),
     $                             PartCoords(2,i,cy2,cx2),
     $                             PartForces(1,j,cy1,cx1),
     $                             PartForces(2,j,cy1,cx1),
     $                             CellDim,TotPot)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
! Now apply the forces to integrate to new positions
         DO CX1=1,NCellRows
            DO CY1=1,NCellCols
               DO I=1,CellSizes(CY1,CX1)
                  CALL UPDATE(
     $                 PARTFORCES(1,i,cy1,cx1),
     $                 PARTFORCES(2,i,cy1,cx1),
     $                 PARTACCELS(1,i,cy1,cx1),
     $                 PARTACCELS(2,i,cy1,cx1),
     $                 PARTVELS(1,i,cy1,cx1),
     $                 PARTVELS(2,i,cy1,cx1),
     $                 PARTCOORDS(1,i,cy1,cx1),
     $                 PARTCOORDS(2,i,cy1,cx1),
     $                 ARMass,TotKin)
               ENDDO
            ENDDO
         ENDDO
         WRITE(6,'(A,I10,A,E20.10,A)')'Iteration ',t,
     $        ' with Total Energy',(TotPot+TotKin)/NParticles
     $         ,' Per Particle'
      ENDDO
      CALL WallTime(timeend)
      WRITE(6,'(A,F20.10)')'The Iteration Time is ',timeend-timestart
      DEALLOCATE(PARTCOORDS)
      DEALLOCATE(PARTFORCES)
      DEALLOCATE(PARTACCELS)
      DEALLOCATE(PARTVELS)
      DEALLOCATE(CellSizes)


      END PROGRAM

      SUBROUTINE INTERACT(x1,y1,x2,y2,fx,fy,CS,TP)
      IMPLICIT NONE
      REAL*8 x1,y1,x2,y2,fx,fy
      REAL*8 rx,ry,r,f,Sig6,Sig12
      REAL*8 AReps,ARsigma,TP
      INTEGER CS
      PARAMETER(AReps=119.8)
      PARAMETER(ARsigma=3.405)

      Rx=x1-x2
      Ry=y1-y2

      R=Rx*Rx+Ry*Ry
      IF((R < 0.00000001))THEN
         R=0.0
         fx=0.0
         fy=0.0
         RETURN
      ENDIF
      R = SQRT(R)
!Derivative of the potential
      Sig6 = (ARsigma/R)**6
      Sig12 = Sig6**2
      F = ((Sig12 - 0.5*Sig6)*48.0*AReps)/R
      Fx = FX + F*(RX)
      Fy = FY + F*(RY)
      TP = TP + 4.0*AReps*(Sig12-Sig6)

      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATE(FX,FY,AX,AY,VX,VY,X,Y,MASS,TK)
      REAL*8 FX,FY,AX,AY,VX,VY,X,Y
      REAL*8 TP,TK,TE
      REAL*8 DT
      REAL*8 MASS
      ! We are using a 1.0 fs timestep, this is the conversion factor
      DT=0.000911633
      AX = FX/MASS
      AY = FY/MASS
      !Update Velocities
      VX = VX + 0.5*DT*AX
      VY = VY + 0.5*DT*AY
      !Update Positions
      X  = X + DT*VX
      Y  = Y + DT*VY
      !Update Energy
      TK = TK + 0.5*MASS*(VX**2+VY**2)
      FX = 0.0
      FY = 0.0

      RETURN
      END

      SUBROUTINE ComputeAtomsPerCell(CellSizes,NCols,NRows,
     $                               NParts,maxCellSize)
      IMPLICIT NONE
      INTEGER nCols,nRows,maxCellSize,nParts
      INTEGER CellSizes(NCols,NRows)
      INTEGER nPartsPerCell
      INTEGER molsum,i,j
      PARAMETER(nPartsPerCell=2)

      maxCellSize = nPartsPerCell
      DO i=1,NRows
         DO j=1,NCols
            CellSizes(j,i) =  nPartsPerCell
         ENDDO
      ENDDO
      molsum = NRows*NCols*NPartsPerCell
      !Diivy up the rest of the atoms amongst the cells
      if(molsum < NParts) THEN
 10      maxCellSize = maxCellSize + 1
         Do i=1,NRows
            DO j=1,NCols
               CellSizes(j,i) = CellSizes(j,i) + 1
               molsum = molsum + 1
               if(molsum >=nParts) Then
                  RETURN
               ENDIF
            ENDDO
         ENDDO
         if(molsum < nParts) then
            goto 10
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE

      subroutine WallTime(t)
      implicit none
      real*8 t
      integer t0, tr, tm
      call system_clock(t0, tr, tm)
      t = dble(t0)/dble(tr)
      return
      end
