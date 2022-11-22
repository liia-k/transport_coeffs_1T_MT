MODULE CONSTANT


IMPLICIT NONE

! Common constants definition: Boltzmann constant kb; 
! atomic mass unit amu; pi; Planck constant hp; Avogadro number
! navog; dimension factor ww for energy calculation cm^{-1} --> J
	
REAL :: kb, amu, pi, hp, navog, ww

PARAMETER (kb=1.3806504e-23, amu=1.6605402e-27, pi=3.14159265358979, hp=6.62606957e-34, &
		   navog=6.02214179e23, ww=1.60219e-19/8065.47)

! N2 spectroscopic data (we_O2, wexe_O2, J)

REAL :: we_N2,wexe_N2

PARAMETER (we_N2=ww*2358.57, wexe_N2=ww*14.32)

! O2 spectroscopic data (we_O2, wexe_O2, J) 

REAL :: we_O2,wexe_O2

PARAMETER (we_O2=ww*1580.19, wexe_O2=ww*11.98)

! NO spectroscopic data (we_NO, wexe_NO, J) 

REAL :: we_NO,wexe_NO

PARAMETER (we_NO=ww*1904.20, wexe_NO=ww*14.075)

! Species mass definition, mass, kg

REAL, DIMENSION(5) :: MASS=(/amu*28.0134, amu*31.998, amu*30.0061, amu*14.0067, amu*15.9994/) 
!REAL, DIMENSION(5) :: MASS=(/7.3080374202e-26, 5.31339653196e-26, 4.6511731002e-26, 2.656764687588e-26, 1.99447483422e-26/)

	! mass(1)= amu*28.0134   = mass N2, kg
	! mass(2)= amu*31.9988  = mass O2, kg
	! mass(3)= amu*30.0061   = mass NO, kg
	! mass(4)= amu*14.0067 = mass N, kg
	! mass(5)= amu*15.9994  = mass O, kg


! For SPECIFIC_HEAT module:

! Number of vibrational levels in CO2 modes (1-3), O2, CO

INTEGER :: L_N2, L_O2, L_NO
DATA L_N2, L_O2, L_NO / 67, 46, 52 /
	
! Arrays containing values of vibrational energy of CO2, O2, CO

REAL, DIMENSION(0:67) :: EN_N2
REAL, DIMENSION(0:46) :: EN_O2
REAL, DIMENSION(0:52) :: EN_NO


! For OMEGA_INTEGRALS and BRACKET_INTEGRALS modules:

! Species gaskinetic diameter definition, sigma, m (parameter for LJ potential)

REAL, DIMENSION(5) :: SIGMA=(/3.621e-10, 3.458e-10, 3.47e-10,&
						3.298e-10, 2.75e-10/) 

	! sigma(1)=3.621e-10 = diameter N2, m
	! sigma(2)=3.458e-10 = diameter O2, m
	! sigma(3)=3.47e-10  = diameter NO, m
	! sigma(4)=3.298e-10  = diameter N, m
	! sigma(5)=2.75e-10  = diameter O, m

! Species well depth definition, epsilon_k, K (parameter for LJ potential)

REAL, DIMENSION(5) :: EPS=(/97.5, 107.4, 119.0, 71.4, 80.0/) 

	! eps(1)=97.5	= epsilon/k N2, K
	! eps(2)=107.4	= epsilon/k O2, K
	! eps(3)=119.0 	= epsilon/k NO, K
	! eps(4)=71.4  	= epsilon/k N, K
	! eps(5)=80.0	= epsilon/k O, K

! Species formation enthalpy definition, h_form, J

REAL, DIMENSION(5) :: HFORM=(/0., 0., 9.029e4/navog,&
							 4.726e5/navog, 2.54e5/navog/) 
	
	! hform(1)= 0             = h_form N2, J
	! hform(2)= 0             = h_form O2, J
	! hform(3)=9.029e4/navog  = h_form NO, J
	! hform(4)=2.54e5/navog   = h_form N, J
	! hform(5)=2.54e5/navog   = h_form O, J
		  
!Omega-integrals and their ratios

REAL, DIMENSION(5,5) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
						AA, BB, CC

!Bracket integrals

REAL, DIMENSION(5,5) :: LAMBDA, LAMBDA00, LAMBDA01, LAMBDA11						

REAL, DIMENSION(5,5) :: ETA, H00, BETA11

REAL, DIMENSION(5,3) :: BETA01

REAL, DIMENSION(3) :: BETA0011

!internal heat conductivity coefficients (LAMBDA_INT)

REAL, DIMENSION(5) :: LAMBDA_INT=(/0., 0., 0., 0., 0./) 


! For TRANSPORT_AIR5_1T module:

!Matrices for the linear transport systems defining
!heat conductivity and thermal diffusion (LTH);
!bulk viscosity (BVISC);
!diffusion (LDIFF);
!shear viscisity (HVISC).

REAL, DIMENSION(10,10) :: LTH 

REAL, DIMENSION(8,8) :: BVISC 

REAL, DIMENSION(5,5) ::  LDIFF, HVISC, b1

!Vectors of right hand terms

REAL, DIMENSION(10,1) :: b

REAL, DIMENSION(5,1) :: b2

REAL, DIMENSION(8,1) :: b3

						
!Macroparameters required for transport coeffs. calculation:

!Vectors of species molar fractions (X), mass fractions (Y)

REAL, DIMENSION(5) :: X, Y

!Common variables: gas temperature (T); pressure (press); total number density (ntot); mixture density (rho)

REAL T, press, ntot, rho 

!REAL :: EPSILON=1e-10 


!Transport coefficients: 

!total heat conductivity; translational heat conductivity;
!N2, O2, NO rotational heat conductivity; ... vibrational heat conductivity;
!shear viscosity; bulk viscosity

REAL ltot, ltr, lint, visc, bulk_visc!, lrot_n2, lrot_o2, lrot_no, lvibr_n2, lvibr_o2, lvibr_no, 

!thermal diffusion coefficients (THDIFF);

REAL, DIMENSION(5) :: THDIFF 

!Diffusion coeffcients matrix

REAL, DIMENSION(5,5) :: DIFF


END MODULE 	CONSTANT

