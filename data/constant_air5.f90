! In this module constants for air5 mixture are defined. 
! The mixture compoents are assumed to be in the following order: N2, O2, NO, N, O.

module constant_air5


implicit none

! Number of species, number of molecular mixture components

integer, parameter :: NUM_SP=5, NUM_MOL=3

! Vector of numbers of vibrational levels for N2, O2, NO

integer, dimension(NUM_MOL), parameter  :: L_vibr = (/67, 56, 52/) ! L_N2=67, L_O2=56, L_NO=52

! Constants for calculation of bracket int PHI:
real, dimension(NUM_MOL), parameter :: CONST_v = (/23.73, 20.72, 9.16/)

! Common constants definition: Boltzmann constant kb (J/K); 
! atomic mass unit amu (kg); pi; Planck constant hp (J*s); Avogadro number
! navog (mol); dimension factor ww for energy calculation m^{-1} --> J; gas constant R (J/K/mol)
	
real, parameter ::  Kb=1.3806504e-23, amu=1.660538921e-27, PI=3.14159265358979, hp=6.626068e-34, &
		   			NAvog=6.02214179e23, ww=1.9864806390700107e-21, R = 8.3144598 ! 1.60219e-19/8065.47)


! Spectroscopic data
					
real, dimension(NUM_MOL), parameter :: we = (/ww*2358.57, ww*1580.19, ww*1904.20/)
real, dimension(NUM_MOL), parameter :: wexe = (/ww*14.32, ww*11.98, ww*14.075/)

! Species mass definition, mass_spcs, kg

real, dimension(NUM_SP), parameter :: MASS_SPCS=(/4.6517344343135997e-26, 5.3135256633152e-26, 4.9826300488143997e-26, &
												  2.3258672171567998e-26, 2.6567628316576e-26/) 

! Species molar mass definition, molar, kg/mol

real, dimension(NUM_SP), parameter  :: MOLAR=(/28.0134e-3, 31.998e-3, 30.0061e-3, 14.0067e-3, 15.9994e-3/) 


! For OMEGA_INTEGRALS and BRACKET_INTEGRALS modules:

! Species gaskinetic diameter definition, sigma, m (parameter for LJ potential)

real, dimension(NUM_SP), parameter  :: SIGMA_LJ=(/3.621e-10, 3.458e-10, 3.47e-10,&
										 		  3.298e-10, 2.75e-10/) 

!Parameters sigma(i) and eps(i) of the Lenard-Jones potential are given by:
   !R.J.Kee, J.A.Miller, T.N. Jefferson, CHEMKIN: A General-Purpose, 
   !Problem-Independent, Transportable, Fortran Chemical Kinetics Code Package,
   !Sandia National Laboratories, SAND80-8003, 1980 

! Species well depth definition, epsilon_k, K (parameter for LJ potential)

real, dimension(NUM_SP), parameter  :: EPS_LJ=(/97.5, 107.4, 119.0, 71.4, 80.0/) 

! Species formation enthalpy definition, h_form, J

real, dimension(NUM_SP), parameter  :: H_FORM=(/0., 0., 1.507112761911427e-19,&
										  		7.818078240456827e-19, 4.098045681049634e-19/) 
	
end module 	constant_air5

