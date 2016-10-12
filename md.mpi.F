      PROGRAM MD
      IMPLICIT NONE
      CHARACTER*100 Buffer
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTCOORDS
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTFORCES
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTACCELS
      REAL*8,  DIMENSION(:,:,:,:),ALLOCATABLE :: PARTVELS
      INTEGER, DIMENSION(:,:),  ALLOCATABLE :: CellSizes
      INTEGER  nCellRows,NCellCols,CellSize,NumIterations,REQUEST
      INTEGER  nParticles, nPartLoc,Dimen,nPartPerCell,CellDim
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
      INTEGER i_proc,j_proc, i, j, proc_col, proc_row,proc_dims 
      INTEGER my_rank, rem, sz, ierror, tag, stat(MPI_STATUS_SIZE), comm

!INITIALIZE MPI
      call MPI_INIT(ierror)

      comm = MPI_COMM_WORLD !Define World of proceses
      
!get my proc nr: my_rank
      call MPI_COMM_RANK(comm, my_rank, ierror)
      
!get total nr of proc.: sz
      call MPI_COMM_SIZE(comm, sz, ierror)

!SET PROC TOPOLOGY, DO IT WITH CART; THIS IS RUBBISH
!get procesor column and row 
      proc_dim = int( sqrt (sz) ) ! Number of proc in each row/col
      proc_row = int( my_rank / proc_dim )
      proc_col = mod( my_rank, proc_dim)   
  
!Read from console Input
      IF(my_rank .eq. 0 ) THEN            
          IF( (IARGC() < 2) ) THEN
             WRITE(6,*)"Incorrect syntax: should be two arguments"
             CALL EXIT(2)
          ENDIF

          CALL GETARG(1,Buffer)
          READ(Buffer,*)nParticles
          CALL GETARG(2,Buffer)
          READ(Buffer,*)NumIterations
      ENDIF

!Pass read variables to other procs
      call MPI_BCAST(nParticles, 1, MPI_INTEGER, 0, comm, ierror)
      call MPI_BCAST(NumIterations, 1, MPI_INTEGER, 0, comm, ierror)
      
!Calculate Number of particles in each proc
      nPartLoc = int( nParticles / sz )
      rem = mod(nParticles / sz)
!Give one particle more to each proc until rem = 0
      if( my_rank .ge. sz-rem ) then 
          nPartLoc = nPartLoc + 1
      end if

      ! Lets Allocate the number of cells in one proc 
      !NO GHOST CELLS IMPLEMENTED
      !Do it on 1 proc, so that grid is homogeneous over all space
      if (my_rank .eq. 0 ) then 
          Dimen = int(Sqrt(NPartLoc/real(NPartPerCell)))
      end if
!Pass Dimen to other procs
      call MPI_BCAST(Dimen, 1, MPI_INTEGER, 0, comm, ierror)
!Each Proc calculates TotalCells, NCellRows, NCellCols (all local)
      TotalCells = Dimen**2
      NCellRows = Dimen
      NCellCols = Dimen
          
      ALLOCATE(CellSizes(0:NCellCols+1,0:NCellRows+1))
      CALL ComputeAtomsPerCell(CellSizes,NCellCols,NCellRows,
     $     NPartLoc,maxCellSize)

      WRITE(6,'(A,I6,A,I6,A,I6,A)')' The Total Number of Cells is ',
     $     TotalCells,' With ',maxCellSize,' particles per cell, and ',
     $     nParticles,' particles total in system'

!OK, now we can allocate the Particle Matrix
      ALLOCATE(PARTCOORDS(2,maxCellSize,0:NCellCols+1,0:NCellRows+1))
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
! INITIAL POSITIONS DEPEND ON THE PROCESOR
      DO i = 1,NCellRows
         DO j = 1,NCellCols
            Centerx = proc_col * CellDim * real(NCellCols) +
     $                    CellDim*real(j) + 0.5*CellDim
            Centery = proc_row * CellDim * real(NCellRows) +
     $                    CellDim*real(i) + 0.5*CellDim
! Set first paricle in the box
            call random_number(ran)
            PartCoords(1,1,j,i) = centerx+((ran-0.5)*(CellDim-2))
            call random_number(ran)
            PartCoords(2,1,j,i) = centery+((ran-0.5)*(CellDim-2))
! Set initial position of the remaining particles in the cell
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

      !START MAIN TIME LOOP
      ! Lets Start Iterating
      DO t=1,NumIterations
         ! Zero Energy variables
         TotPot = 0.0
         TotKin = 0.0

!HERE I SHOULD COMUNICATE WITH NEIGHBORING PROCS
!Exchange cell_sizes to neighboring cells
      CALL CELL_SZ_EXCH

!Recieve all neighboring cells from 8 neighbouring procs      
      CALL PROC_DATA_EXCH
       
!$OMP PARALLEL  PRIVATE(CY1,I,J,CX2,CY2)  
! SHARED(PartCoords,CellSizes,NCellRows,NCellCols) 
!$OMP DO  REDUCTION (+: TotPot,PartForces)
         ! FORCE LOOP ~N*LOG(N) !Loop over Cells
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
         ENDDO !END OF FORCE LOOP

!$OMP END DO
!$OMP END PARALLEL

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
         WRITE(60,*) t, (TotPot+TotKin)/NParticles

!TO TEST CORRECTNESS OF THE CODE WRITE OUT TRAJECTORY
!         WRITE(77,*) t, PARTCOORDS(1,1,1,1),PARTCOORDS(2,1,1,1)


      ENDDO
      CALL WallTime(timeend)
      OPEN(UNIT=99,FILE='Time_to_Solution')
      WRITE(99,'(A,F20.10)')'The Iteration Time is ',timeend-timestart
      CLOSE(99)
      DEALLOCATE(PARTCOORDS)
      DEALLOCATE(PARTFORCES)
      DEALLOCATE(PARTACCELS)
      DEALLOCATE(PARTVELS)
      DEALLOCATE(CellSizes)
     
!END PARALLEL PROCESES
      call MPI_FINALIZE(ierror)

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

      SUBROUTINE ComputeAtomsPerCell(CellSizes, !Number of particles in this cell
     $                               NCols,NRows,
     $             NParts,maxCellSize) !Maximum number of particles in a single cell
!This routine gives out maxCellSize = maximum number of particles in one cell
! And CellSizes(j,i) = number of particles in cell i,j
      IMPLICIT NONE
      INTEGER nCols,nRows,maxCellSize,nParts
      INTEGER CellSizes(NCols,NRows)
      INTEGER nPartsPerCell !Parameter to initiate system
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

      SUBROUTINE PROC_DATA_EXCH(proc_row, proc_col, PartCoords)
!DEFINE ALL VARIABLES THAT IM USING
      INTEGER :: ne_rank,i_proc_NOPBC,j_proc_NOPBC, n_data
    
!MPI COMMUNICATION 
      do i_proc = proc_row-1,proc_row+1
          do j_proc = proc_col-1 ,proc_col+1
              i_proc_NOPBC = i_proc
              j_proc_NOPBC = j_proc

              !Take care of BOUNDARY CONDITIONS (PBC in this case)
               if(i_proc.lt.0) then
                   i_proc = proc_dim - 1
               else if(i_proc.ge.proc_dim) then
                   i_proc = 0
               endif 
               if(j_proc.lt.0) then
                   j_proc = proc_dim - 1  
               else if(j_proc.ge.proc_dim) then
                   j_proc = 0
               endif  
               
        !set rank of neighbouring cell: ne_rank
              ne_rank = i_proc * proc_dim + j_proc 

        !No self-communication, skip central cell
              if ( (i_proc.eq.proc_row).and.(j_proc.eq.proc_col)) cycle
            
! SAME ROW COMMUNICATION BETWEEN PROCS
              if ( i_proc.eq.proc_row ) then 
                  if (j_proc_NOPBC.gt.proc_col) then !Pass to right
                      !Set cells to pass
                      cel_col_send = NCellCols  
                      cel_col_recv = NCellCols + 1
                  else                               !Pass to left
                      cel_col_send = 1
                      cel_col_recv = 0
                  end if
                  
                  call MPI_ISEND(PartCoords(:,:,cel_col_send,:),
     $                sum(CellSizes(cel_col_send,:)), MPI_REAL, 
     $                ne_rank, tag,comm, REQUEST, ierror)
                  
                  call MPI_RECV(PartCoords(:,:,cel_col_recv,:), 
     $            sum(CellSizes(cel_col_recv,:)), MPI_REAL, ne_rank,
     $                  tag, comm, STAT, ierror)   
                  

! SAME COLUMN COMMUNICATION BETWEEN PROCS
              else if (j_proc.eq.proc_col) then
                   if (i_proc_NOPBC.gt.proc_row) then !Pass UP
                      !Set cells to pass
                      cel_row_send = NCellRows  
                      cel_row_recv = NCellRows + 1
                  else                               !Pass DOWN
                      cel_row_send = 1
                      cel_row_recv = 0
                  end if

                   call MPI_ISEND(PartCoords(:,:,:,cel_row_send),
     $                sum(CellSizes(:,cel_row_send)), MPI_REAL, 
     $                ne_rank, tag,comm, REQUEST, ierror)
                  
                  call MPI_RECV(PartCoords(:,:,:,cel_row_recv), 
     $            sum(CellSizes(:,cel_row_recv)), MPI_REAL, ne_rank,
     $                  tag, comm, STAT, ierror)   
                                  
! DIAGONAL COMMUNICATION BETWEEN PROCS
              else
        !Set cel_col and cel_row to send
                  if(i_proc_NOPBC.gt.proc_row) then
                      cel_row = Dimen
                  else
                      cel_row = 1
                  endif
                  if(j_proc_NOPBC.gt.proc_col) then
                      cel_col = Dimen
                  else
                      cel_row = 1
                  endif
 
                  call MPI_ISEND(PartCoords(:,:,cel_col,cel_row), 
     $                   CellSizes(cel_col,cel_row), MPI_REAL, ne_rank,
     $                                 tag,comm, REQUEST, ierror)  
                 
         !Set cel_col and cel_row to recieve
                  if(cel_row.eq.Dimen) then
                      cel_row = Dimen + 1
                  else
                      cel_row = 0
                  endif
                  if(proc_col.eq.Dimen) then
                      cel_col = Dimen + 1
                  else
                      cel_row = 0
                  endif
                   
                  call MPI_RECV(PartCoords(:,:,cel_col,cel_row),
     $           CellSizes(cel_col,cel_row) , MPI_REAL, ne_rank , 
     $                   tag, comm, STAT, ierror)
              endif
          end do
      end do 

      END SUBROUTINE

