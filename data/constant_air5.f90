! In this module constants for air5 mixture are defined. 
! The mixture compoents are assumed to be in the following order: N2, O2, NO, N, O.

module constant_air5


implicit none

! Number of species, number of molecular mixture components

integer, parameter :: NUM_SP=5, NUM_MOL=3

! Vector of numbers of vibrational levels for N2, O2, NO

integer, dimension(NUM_MOL), parameter  :: L_vibr = (/67, 56, 52/) ! L_N2=67, L_O2=56, L_NO=52

! Common constants definition: Boltzmann constant kB (J/K), kB_eV (eV/K); 
! atomic mass unit amu (kg); pi; Planck constant hp (J*s); Avogadro number
! navog (mol); dimension factor ww for energy calculation m^{-1} --> J; gas constant R (J/K/mol)
	
real, parameter ::  kB=1.3806504e-23, amu=1.660538921e-27, PI=3.14159265358979, hp=6.626068e-34, &
		   			NAvog=6.02214179e23, ww=1.9864806390700107e-21, R = 8.3144598, kB_eV = 8.617333262e-5 ! 1.60219e-19/8065.47)


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
   !Problem-Independent, Transportable, Fortran Chemical Kinetics Codae Package,
   !Sandia National Laboratories, SAND80-8003, 1980 

! Species well depth definition, epsilon_k, K (parameter for LJ potential)

real, dimension(NUM_SP), parameter  :: EPS_LJ=(/97.5, 107.4, 119.0, 71.4, 80.0/) 

! VSS parameters: C and omega coefficients, which are taken from:
! Koura, K. and Matsumoto, H., 1992. Variable soft sphere molecular model for air species. 
! Physics of Fluids A: Fluid Dynamics, 4(5), pp.1083-1085.

real, dimension(NUM_SP,NUM_SP), parameter :: C_VSS = reshape((/ &
													216.1e-20, 192.0e-20, 203.8e-20, 240.3e-20, 167.1e-20, & 
													192.0e-20, 155.6e-20, 173.1e-20, 139.9e-20, 114.7e-20, &
													203.8e-20, 173.1e-20, 187.6e-20, 181.5e-20, 135.4e-20, &
													240.3e-20, 139.9e-20, 181.5e-20, 306.8e-20, 163.6e-20, &
													167.1e-20, 114.7e-20, 135.4e-20, 163.6e-20, 106.9e-20 /), &
													shape(C_VSS), order=(/2,1/)) ! m^2
real, dimension(NUM_SP,NUM_SP), parameter :: OMEGA_VSS = reshape((/ &
													0.235, 0.225, 0.230, 0.275, 0.239, & 
													0.225, 0.201, 0.213, 0.221, 0.199, &
													0.230, 0.213, 0.221, 0.244, 0.217, &
													0.275, 0.221, 0.244, 0.328, 0.268, &
													0.239, 0.199, 0.217, 0.268, 0.225 /), &
													shape(OMEGA_VSS), order=(/2,1/)) ! dimensionless

real, dimension(NUM_SP,NUM_SP), parameter :: C_VSS_D = reshape((/ &
													241.4e-20, 199.2e-20, 219.1e-20, 215.4e-20, 150.1e-20, &
													199.2e-20, 160.6e-20, 179.1e-20, 125.9e-20, 103.7e-20, &
													219.1e-20, 179.1e-20, 198.1e-20, 163.1e-20, 122.0e-20, &
													215.4e-20, 125.9e-20, 163.1e-20, 274.3e-20, 146.6e-20, &
													150.1e-20, 103.7e-20, 122.0e-20, 146.6e-20,  96.2e-20 /), &
													shape(C_VSS_D), order=(/2,1/)) ! m^2

real, dimension(NUM_SP,NUM_SP), parameter :: C_VSS_ETA = reshape((/ &
													142.5e-20, 128.6e-20, 135.5e-20, 165.4e-20, 114.6e-20, &
													128.6e-20, 103.9e-20, 115.7e-20,  95.8e-20,  78.5e-20, &
													135.5e-20, 115.7e-20, 125.0e-20, 124.6e-20,  92.7e-20, &
													165.4e-20,  95.8e-20, 124.6e-20, 212.9e-20, 112.4e-20, &
													114.6e-20,  78.5e-20,  92.7e-20, 112.4e-20,  73.3e-20 /), &
													shape(C_VSS_ETA), order=(/2,1/)) ! m^2

real, dimension(NUM_SP,NUM_SP), parameter :: OMEGA_VSS_D = reshape((/ &
													0.274, 0.252, 0.263, 0.295, 0.254, &
													0.252, 0.224, 0.238, 0.234, 0.211, &
													0.263, 0.238, 0.250, 0.260, 0.230, &
													0.295, 0.234, 0.260, 0.357, 0.286, &
													0.254, 0.211, 0.230, 0.286, 0.239 /), &
													shape(OMEGA_VSS_D), order=(/2,1/)) ! dimensionless

real, dimension(NUM_SP,NUM_SP), parameter :: OMEGA_VSS_ETA = reshape((/ &
													0.231, 0.222, 0.227, 0.276, 0.239, &
													0.222, 0.198, 0.210, 0.221, 0.199, &
													0.227, 0.210, 0.218, 0.244, 0.217, &
													0.276, 0.221, 0.244, 0.335, 0.269, &
													0.239, 0.199, 0.217, 0.269, 0.225 /), &
													shape(OMEGA_VSS_ETA), order=(/2,1/)) ! dimensionless

! Born-Mayer potential coefficients:

real, dimension(4,6), parameter :: A_CONST_BM = reshape((/ -267.0, 201.570, 174.672, 54.305, &  
															26700.0, -19226.5, -27693.8, -10860.9, &
															-8.90e5, 6.3201e5, 1.0227e6, 5.4304e5, &
															-33.0838, 20.0862, 72.1059, 68.5001, & 
 															101.571, -56.4472, -286.393, -315.4531, &
															-87.7036, 46.3130, 277.146, 363.1807 /), &
															shape(A_CONST_BM))

! Born-Mayer potential parameters: phi and beta coefficients, which are taken from:														

real, dimension(NUM_SP,NUM_SP), parameter :: PHI_BM = reshape((/ 415.7, 2316.0, 52.38, 184.9, 860.4, & 
																 2316.0, 1.485e4, 6373.0, 905.7, 4530.0, &
																 52.38, 6373.0, 2678.0, 428.6, 2142.0, &
																 184.9, 905.7, 428.6, 86.0, 348.2, &
																 860.4, 4530.0, 2142.0, 348.2, 1410.0 /), &
																 shape(PHI_BM)) ! eV

real, dimension(NUM_SP,NUM_SP), parameter :: BETA_BM = reshape((/ 2.573e+10, 3.267e+10, 1.761e+10, 2.614e+10, 3.331e+10, & 
																  3.267e+10, 3.964e+10, 3.644e+10, 3.320e+10, 4.039e+10, &
																  1.761e+10, 3.644e+10, 3.303e+10, 2.983e+10, 3.717e+10, &  
															  	  2.614e+10, 3.320e+10, 2.983e+10, 2.68e+10, 3.41e+10, &
																  3.331e+10, 4.039e+10, 3.717e+10, 3.41e+10, 4.14e+10 /), &
																 shape(PHI_BM))	! 1/m													



! ESA phenomenological model polynomials a_i (i=1,7) coefficients b_j (j=0,1,2) for different Omega-ints:

real, dimension(7,3), parameter :: OMEGA11_CONST_ESA = reshape((/ 7.884756E-01, -2.438494E-02, 0., &  
																  -2.952759E-01, -1.744149E-03, 0., &
																  5.020892E-01, 4.316985E-02, 0., &
																  -9.042460E-01, -4.017103E-02, 0., & 
																  -3.373058E+00, 2.458538E-01, -4.850047E-03, &
																  4.161981E+00, 2.202737E-01, -1.718010E-02, &
																  2.462523E+00, 3.231308E-01, -2.281072E-02/), &
																 shape(OMEGA11_CONST_ESA), order=(/2,1/))																 

real, dimension(7,3), parameter :: OMEGA12_CONST_ESA = reshape((/ 7.123565E-01, -2.688875E-02, 0.0, &  
																  -2.910530E-01, -2.065175E-03, 0.0, &
																  4.187065E-02, 4.060236E-02, 0.0, &
																  9.287685E-01, -2.342270E-02, 0.0, & 
																  -3.598542E+00, 2.545120E-01, -4.685966E-03, &
																  3.934824E+00, 2.699944E-01, -2.009886E-02, &
																  2.578084E+00, 3.449024E-01, -2.292710E-02/), &
																shape(OMEGA12_CONST_ESA), order=(/2,1/))

real, dimension(7,3), parameter :: OMEGA13_CONST_ESA = reshape((/ 6.606022E-01, -2.831448E-02, 0.0, &  
																 -2.870900E-01, -2.232827E-03, 0.0, &
																 -2.519690E-01, 3.778211E-02, 0.0, &
																 -9.173046E-01, -1.864476E-02, 0.0, & 
																 -3.776812E+00, 2.552528E-01, -4.237220E-03, &
																  3.768103E+00, 3.155025E-01, -2.218849E-02, &
																  2.695440E+00, 3.597998E-01, -2.267102E-02/), &
															  shape(OMEGA13_CONST_ESA), order=(/2,1/))

															  
real, dimension(7,3), parameter :: OMEGA22_CONST_ESA = reshape((/ 7.898524E-01, -2.114115E-02, 0.0, &  
																  -2.998325E-01, -1.243977E-03, 0.0, &
																  7.077103E-01, 3.583907E-02, 0.0, &
																  -8.946857E-01, -2.473947E-02, 0.0, & 
																  -2.958969E+00, 2.303358E-01, -5.226562E-03, &
																  4.348412E+00, 1.920321E-01, -1.496557E-02, &
																  2.205440E+00, 2.567027E-01, -1.861359E-02/), &
															  shape(OMEGA22_CONST_ESA), order=(/2,1/))

! ESA pptential parameters

real, dimension(NUM_SP,NUM_SP), parameter :: E_0_VSS = reshape((/ 11.443, 11.665, 11.635, 8.424, 7.677, & 
															  	  11.665, 11.972, 11.895, 8.633, 7.978, &
																  11.635, 11.895, 11.845, 8.578, 7.862, &
																  8.424, 8.633, 8.578, 6.432, 5.989, &
																  7.677, 7.978, 7.862, 5.989, 5.763 /), &
															  	shape(E_0_VSS), order=(/2,1/)) ! meV (milli eV)															   

real, dimension(NUM_SP,NUM_SP), parameter :: BETA_VSS = reshape((/ 8.07, 8.11, 8.08, 6.94, 7.25, & 
															  	   8.11, 8.14, 8.12, 6.94, 7.26, &
																   8.08, 8.12, 8.09, 6.94, 7.26, &
																   6.94, 6.94, 6.94, 6.61, 6.72, &
																   7.25, 7.26, 7.26, 6.72, 6.90 /), &
															  shape(BETA_VSS), order=(/2,1/)) 															  
! Species formation enthalpy definition, h_form, J

real, dimension(NUM_SP), parameter  :: H_FORM=(/0., 0., 1.507112761911427e-19,&
										  		7.818078240456827e-19, 4.098045681049634e-19/) 

! Molecular components Parker's formula z_inf coefficients

real, dimension(NUM_MOL), parameter :: Z_INF = (/20.38982, 20.84783, 24.0/)

end module 	constant_air5

