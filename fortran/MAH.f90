! ==============================================================================
! 				     MAH version 1.1			       !
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
! DESCRIPTION: MAH is designed to compute RMAH = the minimum atmospheric
!              height of an exoplanet. RMAH is calculated solely from an
!              inputted mass, MP, and radius, RP, of the planet (both of which
!              should be Earth units). Please see the original paper for details
!              on how this is done. Please see the accompanying README file
!              for instrucitons on using this code.
!
! CHANGES: v1.0 Initial version released
!
! ------------------------------------------------------------------------------
! INPUTS:
! ------------------------------------------------------------------------------
!
! n     = Integer length of the vectors MP and RP
! MP(i) = Vector containing planetary masses, in units of Earth masses
! RP(i) = Vector containing planetary radii, in units of Earth radii
!
! ------------------------------------------------------------------------------
! OUTPUTS:
! ------------------------------------------------------------------------------
!
! [Every output has k=2 different versions. k=1 indicates the calculation
!  assumed a water planet is composed of 100%-H20, k=2 indicates that that the
!  calculation instead used 75%-H2O-25%-MgSiO3 (more realistic)]
!
! RH2O(i,k)  = Radius of a planet with mass MP(i) and pure-water composition
! RMAH(i,k)  = Minimum atmospheric height of the i'th object, in Earth radii
! PRMAH(k)   = Fraction of trials yielding a positive RMAH
! nvalid(k)  = # of trials for which the MAH calculation is reliable
! valid(i,k) = Logical flag for the i'th object, specifying whether the 
!              calculation is reliable or not
!    
! ------------------------------------------------------------------------------
!
! ==============================================================================

MODULE MAHmod

 implicit none
 
 CONTAINS

 SUBROUTINE MAH(n,RP,MP,RH2O,RMAH,PRMAH,nvalid,valid)

 implicit none

 INTEGER, INTENT(IN) :: n
 INTEGER :: i, j, k
 REAL(8), DIMENSION(n), INTENT(IN) :: RP, MP
 REAL(8), DIMENSION(n,2), INTENT(OUT) :: RH2O, RMAH
 INTEGER, DIMENSION(2) :: nvalid, npositive
 LOGICAL, DIMENSION(n,2) :: valid
 REAL(8), DIMENSION(2) :: fvalid
 REAL(8), DIMENSION(2), INTENT(OUT) :: PRMAH

 ! Compute RMAH
 nvalid(:) = 0
 DO k=1,2
   DO i=1,n
     call R_water(MP(i),k,RH2O(i,k),valid(i,k))
     IF( valid(i,k) ) THEN
       nvalid(k) = nvalid(k) + 1 ! Increase the count of valid points
     END IF
     RMAH(i,k) = RP(i) - RH2O(i,k)
   END DO
 END DO

 ! What percentage of inputs were valid?
 DO k=1,2
   fvalid(k) = (REAL(nvalid(k))/REAL(n))
   !write(*,*) 'For model ',k,' ',fvalid(k)*1D2,'% of trials valid'
 END DO

 ! Compute P(RMAH>0)
 npositive(:) = 0
 DO k=1,2
   DO i=1,n
     IF( valid(i,k) .AND. RMAH(i,k) .GT. 0.0D0 ) THEN
       npositive(k) = npositive(k) + 1
     END IF
   END DO
   PRMAH(k) = REAL(npositive(k))/REAL(nvalid(k))
   !write(*,*) 'For model ',k,' P(RMAH>0) = ',PRMAH(k)
 END DO
   
 END SUBROUTINE MAH
! =======================================================

! ==============================================================================
SUBROUTINE R_water(M,flag,Rwater,valid)

 ! If flag = 1, then 100%-H2O model
 ! If flag = 2, then 75%-H2O-25%-MgSiO3 model

 implicit none

 INTEGER :: j
 REAL(8), INTENT(IN) :: M
 INTEGER, INTENT(IN) :: flag
 REAL(8) :: logM
 REAL(8), DIMENSION(8,2) :: a
 LOGICAL, INTENT(OUT) :: valid
 REAL(8), INTENT(OUT) :: Rwater
 ! === Polynomial coefficients for Zeng & Sasselov (2013) mass-radius models ===
 REAL(8), PARAMETER, DIMENSION(8) :: p = (/+1.409482429344076D+0,&
					   +3.942477439199774D-1,&
					   +5.014920562075038D-2,&
					   +2.513336445152832D-3,&
					   -4.557149516701664D-4,&
					   -9.717343773829619D-5,&
					   -3.900230062481005D-6,&
					   +1.776907877676882D-7/)
 REAL(8), PARAMETER, DIMENSION(8) :: q = (/+1.346431690428415D+0,&
					   +3.797476623436896D-1,&
					   +4.669407861902054D-2,&
					   +1.992173679766138D-3,&
					   -3.468932022678234D-4,&
					   -7.637885099727161D-5,&
					   -6.314645065750879D-6,&
					   -1.981417530222748D-7/)
 ! === Limiting range of the Zeng & Sasselov (2013) mass-radius models ===
 REAL(8), PARAMETER, DIMENSION(2) :: MPmin = (/4.854443348229839D-4,&
					       6.196844780385370D-5/)
 REAL(8), PARAMETER, DIMENSION(2) :: MPmax = (/4.864096638526573D+2,&
					       3.393232615628230D+2/)

 ! Useful variables
 a(:,1) = p(:)
 a(:,2) = q(:)
 logM = DLOG(M)

 ! Check if valid
 IF( M .GT. MPmin(flag) .AND. M .LT. MPmax(flag) ) THEN
   valid = .TRUE.       ! M is inside range of Zeng & Sasselov (2013) models
 ELSE
   valid = .FALSE. ! M is outside range of Zeng & Sasselov (2013) models
 END IF

 ! Compute R_water
 Rwater = a(1,flag)
 DO j=1,7
   Rwater = Rwater + a(j+1,flag)*(logM)**j
 END DO

END SUBROUTINE R_water
! ==============================================================================

END MODULE MAHmod
