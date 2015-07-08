! ==============================================================================
! 		        MAH Example Calling Code version 1.1		       !
! ==============================================================================
!
! AUTHOR: David Kipping
!         Columbia University
!         Please report any problems to: d.kipping@columbia.edu
!
! CITATION: If using this code, please cite:
!           [1] Kipping, D. M., Spiegel, D. S. & Sasselov, D. D., 2013, 'A 
!               simple, quantitative method to infer the minimum atmospheric 
!               height of small exoplanets', MNRAS, 434, 1883-1888
!           [2] Zeng, L. & Sasselov, D. D., 2013, PASP, 125, 227
!
! DESCRIPTION: This is an example routine which calls the MAH.f90 module.
!              See the accompanying README file for details on its execution.
!
! CHANGES: v1.0 Initial version released
!          v1.1 Added new functions to simplify main code

PROGRAM example

USE MAHmod

 implicit none

 INTEGER :: i, k, m
 INTEGER, PARAMETER :: n = 1D5
 REAL(8), DIMENSION(n) :: RP, MP
 REAL(8), DIMENSION(n,2) :: RH2O, RMAH
 REAL(8), DIMENSION(2) :: PRMAH
 LOGICAL, DIMENSION(n,2) :: valid
 INTEGER, DIMENSION(2) :: nvalid
 REAL(8), DIMENSION(3) :: RH2O_temp, RMAH_temp
 REAL(8), DIMENSION(3,2) :: RH2O_sol, RMAH_sol
 REAL(8), DIMENSION(3) :: MP_sol, RP_sol

 ! === SECTION 1: BASIC EXECUTION ===
 ! Here, I show a simple example of calling the MAH module

 ! = Read in the data ===
 OPEN(UNIT=10,FILE='example_data.dat')
 DO i = 1,n
  read(10,*) MP(i),RP(i)
 END DO
 CLOSE(10)

 ! = Call the MAH subroutine =
 call MAH(n,RP,MP,RH2O,RMAH,PRMAH,nvalid,valid)

 ! = Output the results==
 OPEN(UNIT=11,FILE='100_H2O.dat')
 OPEN(UNIT=12,FILE='75_H2O_25_MgSiO3.dat')
 DO i = 1,n
  write(11,*) RH2O(i,1),RMAH(i,1),valid(i,1)
  write(12,*) RH2O(i,2),RMAH(i,2),valid(i,2)
 END DO
 CLOSE(11)
 CLOSE(12)

 ! === SECTION 2: SOME SIMPLE STATISTICS ===
 ! Here, a few simple statistics are calculated using the results. These may
 ! be useful when preparing a paper.

 ! = Calculate median, negian and posian of RH2O & RMAH =
 DO k=1,2
  call best_values(n,RH2O,RMAH,nvalid,valid,k,RH2O_sol(:,k),RMAH_sol(:,k))
 END DO

 ! = Calculate median, negian and posian of RP & MP =
 RP_sol = marginalize(RP,n)
 MP_sol = marginalize(MP,n)

 ! = Summarize the results =
 101 format (A)
 102 format (A,F6.3,A,F5.3,A)
 103 format (A,F4.2,A,F4.2,A,F4.2,A)
 write(*,101) '=== Observations ==='
 write(*,103) '$M_P =$ $',MP_sol(1),'_{-',MP_sol(2),'}^{+',MP_sol(3),'}$'
 write(*,103) '$R_P =$ $',RP_sol(1),'_{-',RP_sol(2),'}^{+',RP_sol(3),'}$'
 write(*,*) ' '
 DO k=1,2
  IF( k .EQ. 1 ) THEN
   write(*,101) '=== 100%-H2O model ==='
  ELSE IF( k .EQ. 2 ) THEN
   write(*,101) '=== 75%-H2O-25%-MgSiO3 model ==='
  END IF
  write(*,103) '$R_{H2O} =$ $',RH2O_sol(1,k),'_{-',RH2O_sol(2,k),'}^{+',&
               RH2O_sol(3,k),'}$'
  write(*,103) '$R_{MAH} =$ $',RMAH_sol(1,k),'_{-',RMAH_sol(2,k),'}^{+',&
               RMAH_sol(3,k),'}$'
  write(*,102) '$P(R_{MAH}>0) =$ ',PRMAH(k)*1D2,'\% (',&
               DSQRT(2.0D0)*inverf(PRMAH(k)),'-sigma)'
  write(*,*) ' ' 
 END DO

 CONTAINS

! ==============================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ==============================================================================

! ==============================================================================
FUNCTION marginalize(x,n)

 implicit none

 INTEGER, INTENT(IN) :: n
 REAL(8), DIMENSION(n), INTENT(IN) :: x
 REAL(8), DIMENSION(n) :: x_sort
 INTEGER :: median_I, negian_I, posian_I
 REAL(8), DIMENSION(3) :: marginalize

 ! Define useful constants
 negian_I = FLOOR( 0.1586552539314571D0*n )
 median_I = NINT( 0.5D0*n )
 posian_I = CEILING( 0.8413447460685429D0*n )

 ! Sort the list
 x_sort(:) = x(:)
 call sorter(x_sort,n)

 ! Extract negian, median and posian
 marginalize(1) = x_sort(median_I)
 marginalize(2) = x_sort(median_I) - x_sort(negian_I)
 marginalize(3) = x_sort(posian_I) - x_sort(median_I)

END FUNCTION
! ==============================================================================

! ==============================================================================
SUBROUTINE best_values(n,RH2O,RMAH,nvalid,valid,k,RH2O_sol,RMAH_sol)

 implicit none

 INTEGER, INTENT(IN) :: n, k
 INTEGER :: i, m
 REAL(8), DIMENSION(n,2), INTENT(IN) :: RH2O, RMAH
 LOGICAL, DIMENSION(n,2), INTENT(IN) :: valid
 INTEGER, DIMENSION(2), INTENT(IN) :: nvalid
 REAL(8), DIMENSION(nvalid(k)) :: RH2O_temp, RMAH_temp
 REAL(8), DIMENSION(3), INTENT(OUT) :: RH2O_sol, RMAH_sol

 m = 0
 DO i=1,n
   IF( valid(i,k) ) THEN
     m = m + 1
     RH2O_temp(m) = RH2O(i,k)
     RMAH_temp(m) = RMAH(i,k)
  END IF
 END DO
 RH2O_sol = marginalize(RH2O_temp,nvalid(k))
 RMAH_sol = marginalize(RMAH_temp,nvalid(k))

END SUBROUTINE best_values
! ==============================================================================

! ==============================================================================
SUBROUTINE sorter(asort,length)
! Bubble sorts vector asort(length) from smallest to largest

 implicit none

 INTEGER, INTENT(in) :: length
 REAL(8), DIMENSION(length), INTENT(inout) :: asort
 INTEGER :: swaps_made, counts
 REAL(8) :: temp

 DO ! Repeat this loop until we break out
   swaps_made = 0  ! Initially, we've made no swaps
   ! Make one pass of the bubble sort algorithm
   DO counts = 1,(length-1)
     ! If item is greater than the one after it,
     ! then we initiate a swap
     IF( asort(counts) > asort(counts+1) ) THEN
       temp = asort(counts)
       asort(counts) = asort(counts+1)
       asort(counts + 1) = temp
       swaps_made = swaps_made + 1
     END IF
   END DO
   ! If no swaps, break loop
   IF( swaps_made .EQ. 0 ) EXIT
 END DO

END SUBROUTINE sorter
! ==============================================================================

END PROGRAM example
