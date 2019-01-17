!This program generates NE (set below) EoS candidates. The number of
!points (the grid) in the energy-density/pressure plane can be set
!below (Npr, Npe). The program outputs three files. The first one
! is called boundary_pts.dat. This file contains the coordinates of
!the boundaries dictated by causality and monotony for a given choice
!of boundary conditions. The number of points used for the boundary 
!is set by the parameter Nbd. The second output file is called
!EoS_grid.dat. It contains the e-density and pressure coordinates
!of the grid that is used to generate the EoS candidates. Finally,
!a file EoS_points.dat is created. It has NE+1 lines. The first line 
!contains the e-ensity coordinates. The second line corresponds to
!the coordinates of the first EoS candidates, the third line
!to the second EoS candidate, and so on. Each line has
!Ned+2 entries. The extra '2' comes from the start and end point
!of a candidate, that is, the boundary value, which is chosen at
!random (for now, only the pressure value is picked at random).

PROGRAM main

 USE aux

 IMPLICIT NONE


integer, parameter :: Nd=20, Np=20 !in low band, number of pts en densty (Nd),number of pts pressure (Np)                                                                         
integer, parameter :: Nep=200  !number of EoS                                                                                                                                    
integer, parameter :: Neh=12!number of points of density and pressure in high (pQCD) band 
double precision :: h1,E1,P0,Emax,Pmin,Pmax!sp1
double precision :: h,E(Nd),P(Nd), edensity (Nd)!Bandmax(Nd,2),Bandmin(Nd,2)
double precision :: Pin(Nep)!,sp2(Nep)
double precision :: Bandground (Nd),Bandroof (Nd),Bandfloor (Nd) 
double precision :: table(Neh,3), Bandhighground (Neh),Bandhighroof (Neh),Bandhighfloor(Neh), Ehigh(Neh),Phigh(Neh)
double precision :: factor


 INTEGER(KIND=4), PARAMETER :: NE = 1					!number of EoS candidates to be generated
 INTEGER(KIND=4), PARAMETER :: Npr = 20					!number of pts pressure 
 INTEGER(KIND=4), PARAMETER :: Ned = 20					!number of pts en densty
 INTEGER(KIND=4), PARAMETER :: Nbd = 1000				!number of pts for boundary plot
 REAL(KIND=8) :: P1x,P2x,P1y,P2y,P1intx, P1inty,P2intx,P2inty
 REAL(KIND=8), PARAMETER :: Q1x=15231.618285			!bndry pt (energy density at baryon chemical potential 2.6 GeV)
 REAL(KIND=8), PARAMETER :: Q1y=3465.932101				!bndry pt (lower pressure at baryon chemical ppotential 2.6 GeV)
 REAL(KIND=8), PARAMETER :: Q2x=15231.618285			!bndry pt (energy density at baryon chemical potential 2.6 GeV)
 REAL(KIND=8), PARAMETER :: Q2y=4642.5168				!bndry pt (higher pressure at baryon chemical potential 2.6 GeV) 
! REAL(KIND=8), PARAMETER :: Q1x=20305.13				!bndry pt (energy density at baryon chemical potential 2.8 GeV)
! REAL(KIND=8), PARAMETER :: Q1y= 5063.841				!bndry pt (lower pressure at baryon chemical ppotential 2.8 GeV)
! REAL(KIND=8), PARAMETER :: Q2x=20305.13				!bndry pt (energy density at baryon chemical potential 2.8 GeV)
! REAL(KIND=8), PARAMETER :: Q2y=6244.83 				!bndry pt (higher pressure at baryon chemical potential 2.8 GeV)
! REAL(KIND=8), PARAMETER :: P1intx=P1y-Q1y+Q1x			!bndry pt
! REAL(KIND=8), PARAMETER :: P1inty=P1y					!bndry pt
! REAL(KIND=8), PARAMETER :: P2intx=Q2y-P2y+P2x				!bndry pt
! REAL(KIND=8), PARAMETER :: P2inty=Q2y					!bndry pt
!  REAL(KIND=8) :: Q1x,Q1y,Q2x,Q2y
 REAL(KIND=8) :: P1(2),P2(2),Q1(2),Q2(2),P1int(2),P2int(2)		!bndry pts
 REAL(KIND=8) :: eqgrid(Ned,Npr,2)					!grid for EoS
 REAL(KIND=8) :: x, y							    !grid point variables
 REAL(KIND=8) :: sp(2),ep(2)   						!start and end points
 REAL(KIND=8) :: pt1(2),pt2(2)						!points for candidate generation
 REAL(KIND=8) :: expmin,expmax,expinc				!exponent for grid generation
 REAL(KIND=8) :: rn							        !aux random number variable
 INTEGER(KIND=4) :: ydim(Ned),xdim					!grid dimension
 INTEGER(KIND=4) :: candidate(NE,Ned)				!EoS candidate
 LOGICAL :: ACCEPT							        !accept points for EoS
 LOGICAL :: GRIDMASK(Ned,Npr)						!mask with upper and lower boundary encoded
 LOGICAL :: EOSMASK(Ned,Npr)						!mask used throughout generation of EoS candidate
 INTEGER(KIND=4) :: seed    						!seed for random number generator
 CHARACTER(LEN=100) :: cdate						!date used for seed
 CHARACTER(LEN=10) :: ctime 						!time used for seed
 INTEGER(KIND=4) :: i,j,l,m,n,lo,hi,idx,z			!aux ints
 character (len=90) :: EOSfilename
 integer :: jfile,counter,file

 !INITIALIZATIONS
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 !init random number generator with a date & time based seed
 cdate=""
 ctime=""
 CALL date_and_time(cdate,ctime)
 READ( ctime, '(I6)' ) seed
 CALL srand(seed)
 
 
 !Select upper band from: Sammarruca (NLO:case 1; N2LO: case 3; N3LO: case 5), Holt (NLO: case 7; N3LO: case 9;  Jinniu (NLO; case 13; N3LO; case 15; N4LO: case 17) and  Drischler (case 21)
 ! In all the files first column E in MeV/fm≥; second column P in MeV/fm≥
 write(*,*) 'Select up band:1(SammaNLO);3(SammaN2LO);5(SammaN3LO);7(HoltNLO);9(HoltN3LO);13(HuNLO);15(HuN3LO);17(HuN4LO);21(Dris)'
 read (*,*) file
 call select_case (file)
     
   do i=1,Nd
      read (40,*) edensity(i), Bandroof(i)
      E(i) = edensity(i)
   enddo
  close(40)
      E(1) = edensity(1)
      Pmax= Bandroof(Nd)
      Emax= edensity(Nd)
!Select lower band
 write(*,*) 'Select down band case: 2(SammaNLO);4(SammaN2LO);6(SammaN3LO);8(HoltNLO);&
             10(HoltN3LO);14(HuNLO);16(HuN3LO);18(HuN4LO);22(Dris)'
 read (*,*) file
 call select_case (file)
 
 
     do i=1,Nd
       read (40,*) edensity (i), Bandground(i)
     enddo
 close (40)
        P(1) = Bandground(1)
        Pmin= Bandground(Nd)
 

! INITIALIZATION FINISHED, START COMPUTING PRESSURES
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccEOSSammarrucaNLOcutoff1GeV

do jfile = 50+1,50+Nep  !the 30 is to avoid low reserved numbers (such as the screen...)                                                                      \
   write (EOSfilename, '( "EOSDATA.DIR/EOSSammarrucaNLOcutoff1GeV", I6.6, ".dat" )' )  jfile-50
  OPEN(unit=jfile,file=EOSfilename)
enddo  
do m=1, Nep   ! This loop generates each EoS one at a time
550 jfile=m+50
   Rewind (jfile)
   
!EOS points low band    
      factor=Rand()
      do i=2,Nd  
        bandfloor (i) = Min(bandground(i),P(i-1))
150     P(i)= factor*bandroof(i)+(1-factor)*bandfloor(i)
        IF (P(i) .GT. bandroof(i) .OR. (P(i) .LT. bandground(i))) THEN
         factor=Rand()
         goto 150
        ELSE
         IF(P(i) .GE. P(i-1))  THEN
           write(jfile,*) E(i), P(i)
         ELSE
           factor=Rand()
           goto 150
         endif 
        endif 
     
      enddo            

  close(40)
  close (44)

!EOS points intermediate  band   
 !store boundary points for easy use, don't change
 

 
  P1x=Emax
  P2x=Emax
  P1y=Pmin
  P2y=Pmax

  P1intx=P1y-Q1y+Q1x                         !bndry pt
  P1inty=P1y                                 !bndry pt
  P2intx=Q2y-P2y+P2x                         !bndry pt
  P2inty=Q2y                                 !bndry pt

  P1(1) = P1x; P1(2) = P1y
  P2(1) = P2x; P2(2) = P2y
  Q1(1) = Q1x; Q1(2) = Q1y
  Q2(1) = Q2x; Q2(2) = Q2y
  P1int(1) = P1intx; P1int(2) = P1inty 
  P2int(1) = P2intx; P2int(2) = P2inty 

 
 !print upper and lower bound as function of energy density (x)
 !and write it to a file.
  
  OPEN (UNIT=10,FILE='boundary_pts.dat',ACTION="WRITE",STATUS="REPLACE")
  DO i=1, Nbd
    x = P2x + REAL((i-1),8)/(REAL(Nbd-1,8))*(Q2x - P1x)
!   PRINT *,Real((i-1),8), x, up_bd(x,P2,Q2,P2int),  lo_bd(x,P1,Q1,P1int) 
   WRITE( 10,'(F18.5 X F18.5 X F18.5)')  x, up_bd(x,P2,Q2,P2int),  lo_bd(x,P1,Q1,P1int) 
  END DO
  CLOSE(10)


 !Construct the grid.
 !The variable ydim is an array that stores the number of available points in y-direction 
 !(pressure) for a given x-point (energy density). Available means that is lies within the
 !boundaries derived from causality and monotony. Initially, it is put to zero.
 !The data type is integer 4.
  ydim = 0  

 !The grid for the equation candidates that is to be constructed is also initialized to zero. 
 !Its dimensions are: Ned x Npr x 2. For every grid index i and j, there is an energy-density
 !value ( indices (i,j,1) ) and a pressure value ( indices (i,j,2) ).
 !The data type is real 8.

  eqgrid = 0D0
 
 !This logical mask keeps track of grid-points that are within the desired boundary.
 !It is used internally for the generation of the grid. It has dimensions
 ! Ned x Npr, and only those entries which satisfy the boundary constraints are
 !put to T later.
  GRIDMASK = .FALSE.

 !Next, fill the (logarithmic) grid with values that lie in the range
 !as given by the boundary points hardcoded above.
  DO i=1, Ned !loop over energy density points (x-axis)
   expmin = LOG10(P1x+50)
   expmax = LOG10(Q2x-1700)
   expinc = expmin + REAL((i-1),8)/(REAL(Ned-1,8))*(expmax-expmin)
   x = 10.0_8**(expinc)
   DO j=1,Npr !loop over pressure points (y-axis)
    expmin = LOG10(P1y+0.5)
    expmax = LOG10(Q2y-500)
    expinc = expmin + REAL((j-1),8)/(REAL(Ned-1,8))*(expmax-expmin)
    y = 10.0_8**(expinc)
    !PRINT *, x, y
    eqgrid(i,j,1) = x 
    eqgrid(i,j,2) = y
   END DO
  END DO

 !Now fill the logical mask for all values that lie within the boundaries. 
  DO i=1, Ned
   DO j=1, Npr
    x = eqgrid(i,j,1)
    y = eqgrid(i,j,2)
    IF( y .GT. lo_bd(x,P1,Q1,P1int)  .AND. y .LT. up_bd(x,P2,Q2,P2int)) THEN
     GRIDMASK(i,j) = .TRUE.
    END IF
   END DO
  END DO

 !Write the grid data to file.
  OPEN (UNIT=20,FILE='EoS_grid.dat',ACTION="WRITE",STATUS="REPLACE")
  DO i=1, Ned
   DO j=1, Npr
    IF( GRIDMASK(i,j) .EQV. .TRUE. ) THEN
     WRITE (20,'(F18.5 X F18.5)') eqgrid(i,j,1), eqgrid(i,j,2)
    END IF
   END DO
  END DO
  CLOSE(20)

 !Determine number of x-points (it general, this should match Ned).
  xdim = 0
  ydim = 0
  ydim = COUNT(GRIDMASK,2)
! PRINT *, ydim
  DO i=1, Ned !Count all e-density values for which a non-zero number of pressure points satisfy the b.c.
   IF( ydim(i) .GT. 0 ) xdim = xdim + 1
  END DO
!!  PRINT *, "Number of x-points: ", xdim 
  !DO i=1, Ned
 ! PRINT *, "Number of y-points at xpos = ", i, " : ", ydim(i)
 !END DO
 !PRINT *, ydim 
 !PRINT *, eqgrid(1,1,1), eqgrid(1,1,2)

!  PRINT *, "Done."
!  PRINT *, "================================================="
!  PRINT *, "Generating EoSs."

  !Now that the grid is set, let us create the EoS candidates
   sp(1)=P1x
   sp(2)=P(Nd)
!  l=1 !This parameter l numbers the candidates from l=1,...,NE
  l=1
   DO WHILE( l .LE. NE )  !loop over number of candidates
  
     EOSMASK = GRIDMASK    !initialize aux logical mask 
     ACCEPT = .FALSE.      !upper boundary condition flag
   !come up with a random starting point (boundary condition)
      sp(1) = P1x
      sp(2) = P(Nd)

  !come up with a random end point (boundary condition)
      ep(1) = Q1x
      ep(2) =RAND()*(Q2y-Q1y)+Q1y
 
  !initialze aux ints
      i = 1 
      lo = 0
      hi = 0
      idx = 0
 
   DO WHILE (i .LE. Ned+1) !loop over points in x-direction (energy density)
   !For a given position in energy density (x-direction), starting from
   !the randomly chosen starting point (boundary condition), determine
   !a possible paths forward by chosing points at random. If a chosen random 
   !point (a random step) satisfies the monoyony and causality constraints, 
   !accept the point and move on. If it does not, disregard the candidate and
   !start over. One could optimize this procedure by taking only one step back
   !and alter the last step choice, thereby systematically excluding candidates
   !that are to be excluded. However, even without this improvement the code
   !seems to run reasonably fast.    
!   print*, 'llegu√© al if '
     IF( i .EQ. 1 ) THEN !i=1 corresponds to the (randomly chosen) start point
      pt1 = sp
      do j=1,ydim(i)
       pt2 = eqgrid(i,j,:)
        IF( check_constraints(pt1, pt2) .EQV. .FALSE. ) THEN
        EOSMASK(i,j) = .FALSE.
        ELSE 
         IF( lo .EQ. 0 ) THEN 
          lo = j
         END IF
        hi = j 
        ENDIF
!     PRINT *, "Accepted!"
      END DO

    ELSE IF( i .LE. Ned) THEN !i <=Ned are the points that are within the two boundary points.
    CALL RANDOM_NUMBER(rn)
    idx = FLOOR((hi+1-lo)*rn) + lo 
    IF( MOD(l,2) .EQ. 0D0 ) THEN !This introduces the possibility to intorduce a bias for the random procedure.
     IF( idx .GT. lo+1 ) idx = FLOOR((hi+1-lo)*rn) + lo
    END IF
    lo = idx
    pt1 = eqgrid(i-1,idx,:)
    DO j=1, ydim(i)
     pt2 = eqgrid(i,j,:)
     IF( check_constraints(pt1, pt2) .EQV. .FALSE. ) THEN
      EOSMASK(i,j) = .FALSE.
     ELSE 
      hi = j 
     ENDIF
    END DO

   ELSE IF( i .EQ. Ned+1) THEN !Finally, check that also the last step (the choice of the end point) is OK.
    CALL RANDOM_NUMBER(rn)
    idx = FLOOR((hi+1-lo)*rn) + lo 
    pt1 = eqgrid(i-1,idx,:)
    pt2 = ep
    IF( check_constraints(pt1, pt2) .EQV. .TRUE. ) THEN 
     ACCEPT = .TRUE.
 !    PRINT *, "Accepted!"
    END IF
   END IF
  
   !Up to this point, we have determined wether there is a path forward for the chosen
   !random points.
   !PRINT *, "i   = ", i
   !PRINT *, "lo  = ", lo
   !PRINT *, "hi  = ", hi
   !PRINT *, "idx = ", idx
   !PRINT *, "-------------------------------"
 
  !If the step above has been found to satisfy the constraints, accept the point (for now). 

   IF (  ( (i .LE. Ned ) .AND. (COUNT(EOSMASK(i,:)) .GT. 0 ) ) .OR. (ACCEPT .EQV. .TRUE.) ) THEN
!   print*, 'Paso por la primera rama '
    IF ( i .GT. 1 ) THEN
     candidate(l,i-1) = idx
    END IF

   i = i+1
   
   ACCEPT = .FALSE.
   ELSE !Else kicks in if the above procedure failed to produce a candidate. Reset and restart.
!   print*, 'Paso por la segunda rama '
!   print*, 'comprobando bucle en i ',i, Ned+1
      EOSMASK = GRIDMASK    !initialize logical mask 
      ACCEPT = .FALSE.      !upper boundary condition flag
!      REWIND(jfile)
      counter=counter+1
     if (counter.LT.500) then
      sp(1) = P1x
      sp(2) = P(Nd)
      ep(1) = Q1x
      ep(2) = RAND()*(Q2y-Q1y)+Q1y
      i=1
      lo = 0
      hi = 0 
      idx = 0
     else
        counter=0
        Rewind (jfile)
       go to 550
     end if
   END IF  
  END DO
 
 !This point is only reached if a candidate has been produced successfully. 
  !The coordinates of the candidate are written to a file.
  IF ( l .EQ. 1 ) THEN
   OPEN (UNIT=30,FILE="EoS_pts.dat",ACTION="WRITE",STATUS="REPLACE")
   WRITE(30,*)
   WRITE (30,'(F18.5 )',ADVANCE='NO') sp(1)
   DO i=1, xdim
    WRITE (30,'(F18.5 )',ADVANCE='NO') eqgrid(i,1,1)
   END DO
   WRITE (30,'(F18.5 )',ADVANCE='NO') ep(1)
   WRITE(30,*)
  
   WRITE (30,'(F18.5 )',ADVANCE='NO') sp(2)
 
   DO i=1, xdim

    WRITE (30,'(F18.5 )',ADVANCE='NO') eqgrid(i,candidate(l,i),2)
   write(jfile,*) eqgrid(i,1,1), eqgrid(i,candidate(l,i),2)
   END DO
   WRITE (30,'(F18.5 )',ADVANCE='NO') ep(2)
   CLOSE(30)
     
  ELSE
   OPEN (UNIT=30,FILE="EoS_pts.dat",ACTION="WRITE",POSITION="APPEND",STATUS="OLD")
   WRITE (30,'(F18.5 )',ADVANCE='NO') sp(2)
      DO i=1, xdim
       WRITE (30,'(F18.5 )',ADVANCE='NO') eqgrid(i,candidate(l,i),2)
      END DO
   WRITE (30,'(F18.5 )',ADVANCE='NO') ep(2)
   WRITE(30,*)
   CLOSE(30)
  
  END IF
!  PRINT *, "l : ", l
  l = l + 1 !Increment the candidate id.
 
! print*, E(i),P(i)



! PRINT *, "Done."
! PRINT *, "================================================="
! PRINT *, "Output files:"
! PRINT *, "================================================="
! PRINT *, "boundary_pts.dat  (contains upper and lower bdry)"
! PRINT *, "-------------------------------------------------"
! PRINT *, "EoS_grid.dat (contains grid points)"
! PRINT *, "-------------------------------------------------"
! PRINT *, "EoS_pts.dat (contains coordinates of candidates: "
! PRINT *, "            1st line: x-coords, rest: y-coords,  "
! PRINT *, "            each line corresponds to 1 candidate)"
! PRINT *, "================================================="
! PRINT *, "Done."
! PRINT *, "================================================="



 !EOS points high band
 open(unit=43,file='tabla_High.dat')                 !pQCD table at muB=2,6 GeV in three columns (E, Pmin, Pmax) in MeV/fm≥
!open (unit=43,file='tabla_High1.dat')               !pQCD table at muB=2.8 GeV in three columns (E, Pmin, Pmax) in MeV/fm≥
do i=1,Neh
      read (43,*) table(i,1), table(i,2), table(i,3)
 !      if (stat /= 0) exit
      Ehigh(i) = table(i,1)
      Bandhighground(i) = table(i,2)
      Ehigh(1) = Q1x
      Bandhighroof(i) = table(i,3)
      Phigh(1) = ep(2)

   enddo
           
   write(jfile,*) Ehigh(1), Phigh(1)

   do i=2,Neh
        Bandhighfloor(i) = Min(bandhighground(i),Phigh(i-1))
350     Phigh(i)= factor*bandhighroof(i)+(1-factor)*Bandhighfloor(i)
        IF(Phigh(i) .GT. bandhighroof(i)  .OR.(Phigh(i) .LE. bandhighground(i))) THEN 
         factor=Rand()
         goto 350
        Else  
         IF(Phigh(i) .GE. Phigh(i-1))  THEN
          write(jfile,*) Ehigh(i), Phigh(i)
         Else
         factor=Rand()
          goto 350
!          print*, rn        
         endif                   
  
        endif
       
       enddo


  Close (43)
      
    enddo
 
 Close (jfile)

 Close (30)
 !Close (45)


! print*, 'aqui ', m
 enddo !finish the loop iterating over the equation of state


 
END PROGRAM main
