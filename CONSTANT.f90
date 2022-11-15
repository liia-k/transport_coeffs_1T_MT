MODULE constant

IMPLICIT NONE

! Common constants definition: Boltzmann constant kb; 
! atomic mass unit amu; pi; Planck constant hp; Avogadro number
! navog; dimension factor ww for energy calculation cm^{-1} --> J
	
REAL :: kb, amu, pi, hp, navog, ww

PARAMETER (kb=1.3805e-23, amu=1.6605402e-27, pi=3.14156, hp=6.6254E-34, navog=6.0221e23, ww=1.60219e-19/8065.47)

! CO2 spectroscopic data (wi, J) and dissociation energy d_co2, K

REAL :: w1,w2,w3,wx11,wx12,wx22,wx13,wx23,wx33,wxll,d_CO2

PARAMETER (w1=1345.04*ww, w2=667.25*ww, &
	w3=2361.71*ww, wx11=-3.63*ww, wx12=3.44*ww, wx22=-0.635*ww, &
	wx13=-19.28*ww, wx23=-12.51*ww, wx33=-12.56*ww,&
	wxll=0.775*ww,d_co2=64017)

! O2 spectroscopic data (we_O2, wexe_O2, J) 

REAL :: we_O2,wexe_O2

PARAMETER (we_O2=ww*1580.9, wexe_O2=ww*11.98)

! CO spectroscopic data (we_CO, wexe_CO, J) 

REAL :: we_CO,wexe_CO

PARAMETER (we_CO=ww*2169.8, wexe_CO=ww*13.29)


! Species mass definition, mass, kg

REAL, DIMENSION(5) :: MASS=(/amu*44.01, amu*31.998, amu*28.01, amu*15.9994, amu*12.011/) 
!REAL, DIMENSION(5) :: MASS=(/7.3080374202e-26, 5.31339653196e-26, 4.6511731002e-26, 2.656764687588e-26, 1.99447483422e-26/)

	! mass(1)= amu*44.01   = mass CO2, kg
	! mass(2)= amu*31.998  = mass O2, kg
	! mass(3)= amu*28.01   = mass CO, kg
	! mass(4)= amu*15.9994 = mass O, kg
	! mass(5)= amu*12.011  = mass C, kg


! Species gaskinetic diameter definition, sigma, m

REAL, DIMENSION(5) :: SIGMA=(/3.763e-10, 3.458e-10, 3.65e-10,&
						2.75e-10, 2.75e-10/) 

	! sigma(1)=3.763e-10 = diameter CO2, m
	! sigma(2)=3.458e-10 = diameter O2, m
	! sigma(3)=3.65e-10  = diameter CO, m
	! sigma(4)=2.75e-10  = diameter O, m
	! sigma(5)=2.75e-10  = diameter C, m


! Species well depth definition, epsilon_k, K

REAL, DIMENSION(5) :: EPS=(/244.0, 107.4, 98.1, 80.0, 80.0/) 

	! eps(1)=244.0	= epsilon/k CO2, K
	! eps(2)=107.4	= epsilon/k O2, K
	! eps(3)=98.1 	= epsilon/k CO, K
	! eps(4)=80.0  	= epsilon/k O, K
	! eps(5)=80.0	= epsilon/k C, K

! Species formation enthalpy definition, h_form, J

REAL, DIMENSION(5) :: HFORM=(/-3.95e5/navog, 0., -1.15e5/navog,&
						2.54e5/navog, 7.17e5/navog/) 
	
	! hform(1)=-3.95e5/navog  = h_form CO2, J
	! hform(2)=0              = h_form O2, J
	! hform(3)=-1.15e5/navog  = h_form CO, J
	! hform(4)=2.54e5/navog   = h_form O, J
	! hform(5)=7.17e5/navog   = h_form C, J
		  

! Number of vibrational levels in CO2 modes (1-3), O2, CO

INTEGER :: L1, L2, L3, L_O2, L_CO
DATA L1, L2, L3, L_O2, L_CO / 40, 40, 40, 35, 45 /
	

! Arrays containing values of vibrational energy of CO2, O2, CO

REAL, DIMENSION(0:40,0:40,0:40) :: EN3
REAL, DIMENSION(0:35) :: EN_O2
REAL, DIMENSION(0:45) :: EN_CO


!Omega-integrals and their ratios

REAL, DIMENSION(5,5) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
						AA, BB, CC

!Bracket integrals

REAL, DIMENSION(5,5) :: LAMBDA, LAMBDA00, LAMBDA01, LAMBDA11						

REAL, DIMENSION(5,5) :: ETA, H00, BETA11

REAL, DIMENSION(5,3) :: BETA01

REAL, DIMENSION(3) :: BETA0011


!Diffusion coeffcients matrix

REAL, DIMENSION(5,5) :: DIFF						


!Vectors of species number densities (X); thermal diffusion coefficients (THDIFF);
!internal heat conductivity coefficients (LAMBDA_INT)

REAL, DIMENSION(5) :: X, THDIF, LAMBDA_INT=(/0., 0., 0., 0., 0./) 


!Common variables: gas temperature (T), CO2 combined mode temperature (T12);
!asymmetric mode temperature (T3); O2 vibrational temperature (TVO2); 
!CO vibrational temperature (TVCO); pressure (press); total number density (ntot);
!mixture density (rho)

REAL T, T12, T3, TVO2, TVCO, press, ntot, rho 

REAL :: EPSILON=1e-10 

!Transport coefficients: total heat conductivity; translational heat conductivity;
!CO2 rotational heat conductivity; O2, CO internal heat conductivity;
!heat conductivity of combined (1+2) and asymmetric (3) CO2 modes;
!shear viscosity; bulk viscosity

REAL ltot, ltr, lrot_co2, lrot_o2, lrot_co,  &
		  lvibr_o2,lvibr_co,lvibr_12,lvibr_3, visc, bulk_visc


END MODULE constant

