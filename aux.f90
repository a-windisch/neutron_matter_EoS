MODULE aux

 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION up_bd(x,P2,Q2,Pint)

  !This function returns the value of the upper boundary
  !(pressure value) for a given energy density value.

  IMPLICIT NONE

  REAL(KIND=8) :: x, P2(:), Q2(:), Pint(:)
  REAL(KIND=8) :: up_bd

  IF( x .LT. Pint(1) ) THEN
   up_bd = x + P2(2) - P2(1)
  ELSE
   up_bd = Q2(2)
  END IF

 END FUNCTION up_bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION lo_bd(x,P1,Q1,Pint)

  !This function returns the value of the lower boundary
  !(pressure value) for a given energy density value.

  IMPLICIT NONE

  REAL(KIND=8) :: x, P1(:), Q1(:), Pint(:)
  REAL(KIND=8) :: lo_bd

  IF( x .LT. Pint(1) ) THEN
   lo_bd = P1(2)
  ELSE
   lo_bd = x + Q1(2) - Q1(1)
  END IF

 END FUNCTION lo_bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION check_constraints(pt1, pt2)

  !This function returns T if, for two given points p1 and p2, monotony and
  !causality constraints are satisfied.

  IMPLICIT NONE

  REAL(KIND=8) :: pt1(:), pt2(:)
  LOGICAL ::check_constraints

  IF( (pt2(2) .GE. pt1(2))  .AND. (pt2(2)-pt1(2))/(pt2(1)-pt1(1)) .LE. 1D0 ) THEN
   check_constraints = .TRUE.
  ELSE
   check_constraints = .FALSE.
  END IF
  !PRINT *, (pt2(2)-pt1(2))/(pt2(1)-pt1(1))
  !PRINT *, pt1(:)
  !PRINT *, pt2(:)

 END FUNCTION check_constraints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE spline(x, y, n, yp1, ypn, y2)

  !This routine has been taken from Numerical Recipes in Fortran 90, p109
  !Given arrays x(1:n) and y(1:n) containing a tabulated function. i.e. yi=f(xi),
  !with x1<x2<...<xn, and given values yp1 and ypn for the first derivative of the
  !interpolating function at points 1 and n, respectively, this routine returns an
  !array y2(1:n) of length n which contains the second derivatives of the 
  !interpolating function at the tabulated points xi. If yp1 and/or ypn are equal to 1x10^30 
  !or larger, the routine is signaled to set the corresponding boundary condition for a
  !natural spline, with zero second derrivative on that boundary. Parameter NMAX is the 
  !largest anticipated value of n.

  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER :: NMAX = 4000
  INTEGER(KIND=4) :: n
  REAL(KIND=8) :: yp1, ypn, x(n), y(n), y2(n)
  INTEGER(KIND=4) :: i,k
  REAL(KIND=8):: p, qn, sig, un, u(NMAX)	
  
  IF (yp1 .GT. .99E30 ) THEN			!The lower boundary condition is set either to be "natural"
   y2(1) = 0.
   u(1) = 0.
  ELSE						!or else to have a specified first derivative
   y2(1)=-0.5
   u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  END IF
 
  DO i=2, n-1					!This is the decomposition loop of the tridiagonal algorithm.
   sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))		!y2 and u are used for temporary storage of the decomposed 
   p = sig*y2(i-1)+2.				!factors.
   y2(i) = (sig-1.)/p
   u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))	&
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  END DO

  IF(ypn .GT. .99E30) THEN			!The upper boundary condition is set either to be "natural"
   qn=0.
   un=0.
  ELSE						!or else to have a specified first derivative.
   qn=0.5	
   un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  END IF
  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
  DO k=n-1,1,-1					!This is the backsubstitution loop of the tridiagonal
   y2(k) = y2(k)*y2(k+1)+u(k)			!algorithm
  END DO
!  PRINT *, u(1), u(2), u(3)
 
 END SUBROUTINE spline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE splint(xa, ya, y2a, n, x, y)

  !This routine has been taken from Numerical Recipes in Fortran 90, p110
  !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
  !xai's in order), and given the array y2a(1:n), which is the output from spline above,
  !and given a value of x, this routine returns a cubic-spline interpolated value y.

  IMPLICIT NONE

  INTEGER(KIND=4) :: n
  INTEGER (KIND=4) :: k, khi, klo
  REAL(KIND=8) :: x,y,xa(n),y2a(n),ya(n)
  REAL(KIND=8) :: a,b,h
  
  klo = 1					!We will find the right place in the table by means of 
  khi = n					!bisection. This is optimal if sequential calls to this routine
  DO WHILE (khi-klo .GT. 1)  			!are at random values of x. If sequential calls are in order, and closely
   k=(khi+klo)/2				!spaced, one would do better to store previous values of klo and khi
   IF(xa(k) .GT. x) THEN			!and test if they remain appropriate on the next call.
    khi=k
   ELSE
    klo=k
   END IF
  END DO					!klo and khi now bracket the input value of x.

  h = xa(khi)-xa(klo)
  a = (xa(khi)-x)/h				!Cubic spline polynomial is now evaluated
  b = (x-xa(klo))/h
  y = a*ya(klo)+b*ya(khi) + &
      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  
 END SUBROUTINE splint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION linear_interpol(x,f,N,valx)

  IMPLICIT NONE

  REAL(KIND=8) :: linear_interpol
  REAL(KIND=8) :: x(:), f(:)
  REAL(KIND=8) :: valx, x1, x2
  INTEGER(KIND=4) :: N
  INTEGER(KIND=4) :: ct, x1i, x2i

  IF( valx .LT. x(1) ) THEN
   x1=x(1); x1i=1
   x2=x(2); x2i=2
  ELSE IF ( valx .GT. x(N) ) THEN
   x1=x(N-1); x1i=N-1
   x2=x(N); x2i=N
  ELSE
   DO ct=1, N-1
    IF( valx .GE. x(ct) .AND. valx .LE. x(ct+1) ) THEN
     x1=x(ct); x1i=ct
     x2=x(ct+1); x2i=ct+1
    END IF
   END DO
  END IF

  linear_interpol = f(x1i) + (f(x2i) - f(x1i))*(valx - x1)/(x2 - x1) 
 
 END FUNCTION linear_interpol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE create_candidate( xa, candidate, Ned, Npr )

  !This subroutine generates a random  EoS candidate.
   
  IMPLICIT NONE
 
  INTEGER(KIND=4) :: Npr                                                 !number of pts pressure 
  INTEGER(KIND=4) :: Ned                                                 !number of pts en densty
  REAL(KIND=8) :: P1x,P2x,P1y,P2y,P1intx, P1inty,P2intx,P2inty
 ! REAL(KIND=8), PARAMETER :: P1x=173.0_8                                 !bndry pt
 ! REAL(KIND=8), PARAMETER :: P1y=2.0_8                                   !bndry pt
 ! REAL(KIND=8), PARAMETER :: P2x=173.0_8                                 !bndry pt
 ! REAL(KIND=8), PARAMETER :: P2y=3.0_8                                   !bndry pt
  REAL(KIND=8), PARAMETER :: Q1x=14333.0_8                               !bndry pt
  REAL(KIND=8), PARAMETER :: Q1y=1900.0_8                                !bndry pt
  REAL(KIND=8), PARAMETER :: Q2x=14333.0_8                               !bndry pt
  REAL(KIND=8), PARAMETER :: Q2y=2900.0_8                                !bndry pt
!  REAL(KIND=8), PARAMETER :: P1intx=P1y-Q1y+Q1x                          !bndry pt
!  REAL(KIND=8), PARAMETER :: P1inty=P1y                                  !bndry pt
!  REAL(KIND=8), PARAMETER :: P2intx=Q2y-P2y+P2x                          !bndry pt
!  REAL(KIND=8), PARAMETER :: P2inty=Q2y                                  !bndry pt
  REAL(KIND=8) :: P1(2),P2(2),Q1(2),Q2(2),P1int(2),P2int(2)              !bndry pts
  REAL(KIND=8) :: eqgrid(Ned,Npr,2)                                      !grid for EoS
  REAL(KIND=8) :: x, y                                                   !grid point variables
  REAL(KIND=8) :: xa(:)                                                  !abscissa for EoS
  REAL(KIND=8) :: candidate(:)                                           !EoS candidate
  REAL(KIND=8) :: sp(2),ep(2)                                            !start and end points
  REAL(KIND=8) :: pt1(2),pt2(2)                                          !points for candidate generation
  REAL(KIND=8) :: expmin,expmax,expinc                                   !exponent for grid generation
  REAL(KIND=8) :: rn                                                     !aux random number variable
  REAL(KIND=8) :: xval                                                   !energy density value at which 
  INTEGER(KIND=4) :: ydim(Ned),xdim                                      !grid dimension
  LOGICAL :: ACCEPT                                                      !accept points for EoS
  LOGICAL :: GRIDMASK(Ned,Npr)                                           !mask with upper and lower boundary encoded
  LOGICAL :: EOSMASK(Ned,Npr)                                            !mask used throughout generation of EoS candidate
  LOGICAL :: INTERPOL                                                    !interpolation flag
  INTEGER(KIND=4) :: i,j,l,n,lo,hi,idx                                   !aux ints
                                                                                    
                                            
          

  !check dimesions
  IF ( (SIZE(xa) .NE. Npr+2) .AND. ( SIZE(candidate) .NE. Npr+2 ) ) THEN
   PRINT *, "The dimensions of the arrays 'x' and 'candidate' must be"
   PRINT *, "Npr + 2. Please confirm, compile and rerun."
   STOP
  END IF

  !store boundary points for easy use, don't change
  P1intx=P1y-Q1y+Q1x                         !bndry pt                                                  
  P1inty=P1y                                 !bndry pt                                                  
  P2intx=Q2y-P2y+P2x                         !bndry pt                                                  
  P2inty=Q2y 
  P1(1) = P1x; P1(2) = P1y
  P2(1) = P2x; P2(2) = P2y
  Q1(1) = Q1x; Q1(2) = Q1y
  Q2(1) = Q2x; Q2(2) = Q2y
  P1int(1) = P1intx; P1int(2) = P1inty 
  P2int(1) = P2intx; P2int(2) = P2inty 

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
   expmin = LOG10(P1x+5)
   expmax = LOG10(Q2x-50)
   expinc = expmin + REAL((i-1),8)/(REAL(Ned-1,8))*(expmax-expmin)
   x = 10.0_8**(expinc)
   DO j=1,Npr !loop over pressure points (y-axis)
    expmin = LOG10(P1y+0.1)
    expmax = LOG10(Q2y-50)
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
 
  !Determine number of x-points (it general, this should match Ned).
  xdim = 0
  ydim = 0
  ydim = COUNT(GRIDMASK,2)
  !PRINT *, ydim
  DO i=1, Ned !Count all e-density values for which a non-zero number of pressure points satisfy the b.c.
   IF( ydim(i) .GT. 0 ) xdim = xdim + 1
  END DO

  EOSMASK = GRIDMASK    !initialize aux logical mask 
  ACCEPT = .FALSE.      !upper boundary condition flag
  !come up with a random starting point (boundary condition)
  sp(1) = P1x
  sp(2) = RAND()*(P2y-P1y)+P1y
  !come up with a random end point (boundary condition)
  ep(1) = Q1x
  ep(2) = RAND()*(Q2y-Q1y)+Q1y
  !initialze aux ints
  i = 1 
  lo = 0
  hi = 0
  idx = 0 

  DO WHILE (i .LE. Ned+1  ) !loop over points in x-direction (energy density)
   !For a given position in energy density (x-direction), starting from
   !the randomly chosen starting point (boundary condition), determine
   !a possible paths forward by chosing points at random. If a chosen random 
   !point (a random step) satisfies the monoyony and causality constraints, 
   !accept the point and move on. If it does not, disregard the candidate and
   !start over. One could optimize this procedure by taking only one step back
   !and alter the last step choice, thereby systematically excluding candidates
   !that are to be excluded. However, even without this improvement the code
   !seems to run reasonably fast.    
   IF( i .EQ. 1 ) THEN !i=1 corresponds to the (randomly chosen) start point
    pt1 = sp
    DO j=1, ydim(i)
     pt2 = eqgrid(i,j,:)
     IF( check_constraints(pt1, pt2) .EQV. .FALSE. ) THEN
        EOSMASK(i,j) = .FALSE.
     ELSE 
      IF( lo .EQ. 0 ) THEN 
       lo = j
      END IF
      hi = j 
     ENDIF
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
     !PRINT *, "Accepted!"
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
    IF ( i .GT. 1 ) THEN
     candidate(i) =  eqgrid(i-1,idx,2) 
    END IF
    i = i+1
    ACCEPT = .FALSE.
   ELSE !Else kicks in if the above procedure failed to produce a candidate. Reset and restart.
    EOSMASK = GRIDMASK    !initialize logical mask 
    ACCEPT = .FALSE.      !upper boundary condition flag
    !come up with a new random start point (boundary condition)
    sp(1) = P1x
    sp(2) = RAND()*(P2y-P1y)+P1y
    !come up with a new random end point (boundary condition)
    ep(1) = Q1x
    ep(2) = RAND()*(Q2y-Q1y)+Q1y
    i = 1
    lo = 0
    hi = 0
    idx = 0 
    !PRINT *, "Reset"
   END IF  
  END DO  

  !This point is only reached if a candidate has been produced successfully. 
  !Now fill the arrays x and candidate such that also the starting point
  !and the end point are contained.
  xa(1)       = sp(1)         !store start point x-coordinate
  xa(2:Ned+1) = eqgrid(:,1,1) !store grid x-coordinates
  xa(Ned+2)   = ep(1)         !store end point x-coordinate
  candidate(1)       = sp(2)  !store start point y-coordinate
  candidate(Ned+2)   = ep(2)  !store end point y-coordinate 

 END SUBROUTINE create_candidate 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION interpol_EoS_splines( x, candidate, N, xpos )

  !This function uses the functions spline and splint
  !above to perform a cubic spline interpolation. 

  IMPLICIT NONE

  INTEGER(KIND=4) :: N                   !dimension of x and candidate
  REAL(KIND=8) :: x(:), candidate(:)     !abscissa and EoS candidate
  REAL(KIND=8) :: sd(N)                  !2nd derivative
  REAL(KIND=8) :: xpos                   !interpolate at this point
  REAL(KIND=8) :: res                    !interpolation result
  REAL(KIND=8) :: interpol_EoS_splines   !return value

  CALL spline( x, candidate, N , 0D0, 0D0, sd )  
  CALL splint( x, candidate, sd, N, xpos, res ) 

  interpol_EoS_splines = res

 END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION monotonic_interpolation( x, y, N, valx)

  !This function is an implementation of an algorithm presented in the paper of
  !M. Steffen, 'A simple method for monotonic interpolation in one dimension',
  !Astron. Astrophys. 239, 443-450 (1990)
  !The resulting interpolation is smooth (in the first order derivative, not
  !necessarily in second), and there are no extrema between any of the given points.
  !Local extrema are thus only to be found at the nodes of the original data.
  !Implemented by Andreas Windisch

  IMPLICIT NONE

  INTEGER(KIND=4) :: N                     !dimension of x and candidate
  REAL(KIND=8) :: x(:), y(:)               !abscissa and function
  REAL(KIND=8) :: valx                     !interpolate at this point
  REAL(KIND=8) :: res                      !interpolation result
  REAL(KIND=8) :: monotonic_interpolation  !return value
  REAL(KIND=8) :: xi,xip1,yi,yip1,si,hi    !aux values
  REAL(KIND=8) :: ypi,ypip1,ai,bi,ci,di    !aux values
  REAL(KIND=8) :: him1,sim1,xim1,yim1,pi   !aux values
  REAL(KIND=8) :: xip2,yip2,pip1,hip1,sip1 !aux values
  INTEGER(KIND=4) :: ct,i                  !aux ints

  !I do not extrapolate, but simply kepp the value constant beyond the boundaries
  IF( valx .LE. x(1) ) THEN
   res = y(1)
  ELSE IF ( valx .GE. x(N) ) THEN
   res = y(N)
  ELSE

  !This interpolation procedure is local, that is, it only depends on the
  !neighboring points. We thus have to find out where those neighboring points are.
   DO i=1, N-1
    IF( valx .GE. x(i) .AND. valx .LE. x(i+1) ) THEN
     xi=x(i); yi=y(i) 
     xip1=x(i+1); yip1=y(i+1)
     IF( i .GT. 1 ) THEN
      xim1=x(i-1); yim1=y(i-1)
     END IF
     IF( i .LT. N-2 ) THEN
      xip2 = x(i+2); yip2 = y(i+2)
     END IF
     EXIT
    END IF
   END DO

   IF( i .EQ. 1 ) THEN
    si    = (yip1 - yi  )/(xip1 - xi  )
    hi    = (xip1 - xi  )
    sip1  = (yip2 - yip1)/(xip2 - xip1)
    hip1  = (xip2 - xip1)
    pi    = si*(1.0_8+hi/(hi+hip1))-sip1*hi/(hi+hip1) 
    ypi   = ABS(SIGN(1.0_8,pi) + SIGN(1.0_8,si))*MIN(ABS(si),0.5_8*ABS(pi))
    pip1  = (si*hip1 + sip1*hi)/(hi + hip1)
    ypip1 = (SIGN(1.0_8,si) + SIGN(1.0_8,sip1))*MIN(ABS(si),ABS(sip1),0.5_8*ABS(pip1))
   ELSE IF ( i .EQ. N-1 ) THEN
    si    = (yip1 - yi  )/(xip1 - xi  )
    sim1  = (yi   - yim1)/(xi   - xim1)
    hi    = (xip1 - xi  )
    him1  = (xi   - xim1)
    pi    = (sim1*hi + si*him1)/(him1 + hi)
    ypi   = (SIGN(1.0_8,sim1) + SIGN(1.0_8,si))*MIN(ABS(sim1),ABS(si),0.5_8*ABS(pi))
    pip1  = si*(1.0_8+hi/(hi+him1))-sim1*hi/(hi+him1)
    ypip1 = ABS(SIGN(1.0_8,pip1) + SIGN(1.0_8,si))*MIN(ABS(si),0.5_8*ABS(pip1))
   ELSE
    si    = (yip1 - yi  )/(xip1 - xi  )
    sim1  = (yi   - yim1)/(xi   - xim1)
    hi    = (xip1 - xi  )
    him1  = (xi   - xim1)
    pi    = (sim1*hi + si*him1)/(him1 + hi)
    ypi   = (SIGN(1.0_8,sim1) + SIGN(1.0_8,si))*MIN(ABS(sim1),ABS(si),0.5_8*ABS(pi))
    sip1  = (yip2 - yip1)/(xip2 - xip1)
    hip1  = (xip2 - xip1)
    pip1  = (si*hip1 + sip1*hi)/(hi + hip1)
    ypip1 = (SIGN(1.0_8,si) + SIGN(1.0_8,sip1))*MIN(ABS(si),ABS(sip1),0.5_8*ABS(pip1))
   END IF
   ai    = (ypi  + ypip1 - 2*si )/hi**2
   bi    = (3*si - 2*ypi - ypip1)/hi
   ci    = ypi
   di    = yi
   res   = ai*(valx - xi)**3 + bi*(valx - xi)**2 + ci*(valx - xi) + di 
  END IF

  !return interpolated value 
  monotonic_interpolation = res 

END FUNCTION monotonic_interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE select_case (file)

  !This subroutine selects  upper/lower bands                                                                                                                                        

   IMPLICIT NONE
   integer :: file
   REAL(KIND=8) :: Q1x,Q1y,Q2x,Q2y

    SELECT CASE (File)
!! In all the files first column is E in MeV/fm³; second column  is P in MeV/fm³    
!Sammarruca bands    
     CASE (1)
     open(unit=40,file='EOS_Samma_NLO_600.dat')! up NLO Sammarruca band
     CASE (2)
      open(unit=40,file='EOS_Samma_NLO_450.dat')! down NLO Sammarruca  
     CASE (3)
       open(unit=40,file='EOS_Samma_N2LO_450.dat')!high N2LO Sammarruca band                                                                                                                                                                                                      
     CASE (4)
      open(unit=40,file='EOS_Samma_N2LO_600.dat')! down N2LO Sammarruca band
     CASE (5)
       open(unit=40,file='EOS_Samma_N3LO_450.dat')!  up N3LO Sammarruca band                                                                                                                                                                                                            
     CASE (6)
      open(unit=40,file='EOS_Samma_N3LO_600.dat')! down N3LO Sammarruca band
!Holt bands      
     CASE (7)
     open(unit=40,file='EOS_Holt_NLO_500.dat')! up NLO Holt band   
     CASE (8)
       open(unit=40,file='EOS_Holt_NLO_450.dat')!down NLO Holt band
     CASE (9)
      open(unit=40,file='EOS_Holt_N3LO_450.dat')! Holt band                                                                                                                                                                                                                    
     CASE (10)
      open(unit=40,file='EOS_Holt_N3LO_500.dat')! Holt                                                                                                                                                                                                                         
     CASE (11)
         open(unit=40,file='EOS_Holt_N2LO_450.dat')!higher Holt  band                                                                                                                                                                                                              
     CASE (12)
        open(unit=40,file='EOS_Holt_N2LO_500.dat')! very near N2LO450                                                                                                                                                                                                               
! Hu bands                                                                                                                                                                                                                                                          
    CASE (13)
      open(unit=40,file='EOS_Jinniu_10_NLO.dat')! up NLO Jinniu band
    CASE (14)
      open(unit=40,file='EOS_Jinniu_0_9_NLO.dat')! down NLO Jinniuband 
    CASE (15)
      open(unit=40,file='EOS_Jinniu_10_N3LO.dat')! up N3LO Jinniu band 
    CASE (16)
        open(unit=40,file='EOS_Jinniu_0_9_N3LO.dat')! down N3LO Jinniu
    CASE (17)
     open(unit=40,file='EOS_Jinniu_0_9_N4LO.dat')! up N4LO (higher Jinniu 0.9 fm  band)
    CASE (18)
       open(unit=40,file='EOS_Jinniu_10_N4LO.dat')! down N4LO Bans (higher Jinniu 1.0  fm  band)                                                                                                                                                                                                         
    CASE (19)
      open(unit=40,file='EOS_Jinniu_10_N2LO.dat')!lower Jinniu 0.9 fm band                                                                                                                                                                                                     
    CASE (20)
       open(unit=40,file='EOS_Jinniu_0_9_N2LO.dat')!lower Jinniu 0.9 fm band                                                                                                                                                                                                           
!Drischler band    
    CASE (21)
      open(unit=40,file='EOS_Drischler_high.dat')!!Drischler high 22 pts                                                                                                                                                                                                       
    CASE (22)
       open(unit=40,file='EOS_Drischler_down.dat') !Drischler down 22 pt)

     

    END SELECT
    
 END SUBROUTINE select_case 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE aux
