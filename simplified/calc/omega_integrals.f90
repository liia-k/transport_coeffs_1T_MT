! In this module, Omega-integrals and their ratios AA, BB, CC are calculated 
! using the Lennard-Jones potential for moderate temperatures 

module omega_integrals

use constant_air5

implicit none

! Omega-integrals and their ratios

type omega_int
	real, dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
									  AA, BB, CC
end type


contains

subroutine OmegaInt(T, omega_out)

	! Calculation of OMEGA-integrals at given temperature

	real,intent(in)   :: T ! temperature
	type(omega_int),intent(out) :: omega_out
  
	integer i, j
	real x11, eij, mij, sig_ij, sij, tx
  
	!Omega-integrals and their ratios
  
	real, dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
									  AA, BB, CC
  
	do i=1,NUM_SP
	  do j=1,NUM_SP
		sig_ij = (SIGMA_LJ(i)+SIGMA_LJ(j))/2
		mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E27)/(MASS_SPCS(i)*1E27 + MASS_SPCS(j)*1E27)
		eij = sqrt(EPS_LJ(i)*EPS_LJ(j)*(SIGMA_LJ(i)*1e10)**6 &
			  *(SIGMA_LJ(j)*1e10)**6)/(sig_ij*1e10)**6
		sij = PI*(sig_ij)**2*sqrt((Kb/mij)*t/2./pi) ! of 1e-7 and higher order
		tx = T/eij
			  ! Omega11, Omega12, Omega12, Omega13, B, C for 
			  ! moderate temperatures; Lennard-Jones potential

		! All the omega-ints are of 1e-6 and higher order

		x11 = log((tx)) + 1.4
		omega11(i,j) = sij/(-0.16845 - 2.25768e-2/x11/x11 &
		+ 0.19779/x11 + 0.64373*x11 - 9.26718e-2*x11*x11 &
		+ 7.1131e-3*x11**3)
		!WRITE (*,*) 'Omega11_ij = ', omega11(i,j)

		x11 = log((tx)) + 1.5
		omega22(i,j) = 2.*sij/(-0.40811 - 5.08552e-2/x11/x11 &
		+ 0.3401/x11 + 0.70375*x11 - 0.10699*x11*x11 &
		+ 7.62686e-3*x11**3)
		!WRITE (*,*) 'Omega22_ij = ', omega22(i,j)

		x11 = log((tx)) + 1.1
		omega12(i,j)=3.*sij/(0.40785 + 9.25303e-4/x11/x11 &
		+ 2.79680e-4/x11 + 0.44739*x11-6.27242e-2*x11*x11 &
		+ 5.98567e-3*x11**3)
		!WRITE (*,*) 'Omega12_ij = ', omega12(i,j)

		x11 = log((tx)) + 4.
		omega13(i,j) = 12.*sij/(25.04929 + 63.12444/x11/x11 &
		- 65.87398/x11 - 4.13758*x11 + 0.34999*x11*x11 &
		- 1.0096e-2*x11**3)
		!WRITE (*,*) 'Omega13_ij = ', omega13(i,j)
	
		bb(i,j) = (5./3.*omega12(i,j) - 4./12.*omega13(i,j))/omega11(i,j)

		cc(i,j) = omega12(i,j)/omega11(i,j)/3.0
		
		aa(i,j) = omega22(i,j)/omega11(i,j)/2.0
	  end do
	end do
  
	omega_out%OMEGA11 = omega11
	omega_out%OMEGA22 = omega22
	omega_out%OMEGA12 = omega12
	omega_out%OMEGA13 = omega13
	omega_out%AA = AA
	omega_out%BB = BB
	omega_out%CC = CC
  
	end subroutine OmegaInt

end module omega_integrals