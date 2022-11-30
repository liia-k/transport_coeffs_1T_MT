MODULE CONSTANT


IMPLICIT NONE

! Common constants definition: Boltzmann constant kb; 
! atomic mass unit amu; pi; Planck constant hp; Avogadro number
! navog; dimension factor ww for energy calculation cm^{-1} --> J
	
REAL :: kb, amu, pi, hp, navog, ww, R

PARAMETER (kb=1.3806504e-23, amu=1.660538921e-27, pi=3.14159265358979, hp=6.62606957e-34, &
		   navog=6.02214179e23, ww=1.9864806390700107e-23, R = kb*navog)! 1.60219e-19/8065.47)

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

! Species molar mass definition, mass, kg

REAL, DIMENSION(5) :: MOLAR=(/28.0134, 31.998, 30.0061, 14.0067, 15.9994/) 

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

						
!Macroparameters required for transport coeffs. calculation:

type transport_in
real rho  
!real vel(1:3), not required
real temp
real mass_fractions(1:5)
end type

!Transport coefficients: 

type transport_out
real visc, bulk_visc, ltot ! shear viscosity,  bulk viscosity, thermal conductivity
real :: thdiff(1:5) ! thermal diffusion coeffs
real :: diff(1:5,1:5) ! diffusion coeffs
end type							 

END MODULE 	CONSTANT

