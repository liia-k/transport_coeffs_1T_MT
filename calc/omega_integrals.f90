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

subroutine OmegaInt(T, omega_out, interactionType)

	! Calculation of OMEGA-integrals at given temperature

	real,intent(in)   :: T ! temperature
	type(omega_int),intent(out) :: omega_out
	character(len=*), intent(in) :: interactionType
  
	integer i, j, k
	real x11, eij, mij, sig_ij, sij, tx, term, term2, term3, term4, eij_eV, T_star, nu_star, r_star, x
	real, dimension(6) :: Ak_BM
	real, dimension(7) :: aj_ESA
  
	!Omega-integrals and their ratios
  
	real, dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
									  AA, BB, CC
	
	select case (interactionType)
	
	case ('VSS')

	do i=1,NUM_SP
		do j=1,NUM_SP
			mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E27)/(MASS_SPCS(i)*1E27 + MASS_SPCS(j)*1E27)
			term = sqrt(Kb/mij/8./pi) 
			term2 = C_VSS(i,j)*T**(0.5 - OMEGA_VSS(i,j))
			! All the omega-ints are of 1e-6 and higher order
	
			omega11(i,j) = term * term2 * gamma(3 - OMEGA_VSS(i,j))
	
			omega22(i,j) = term * term2 * gamma(4 - OMEGA_VSS(i,j))
	
			omega12(i,j) = term * term2 * (3 - OMEGA_VSS(i,j)) * gamma(3 - OMEGA_VSS(i,j))
	
			omega13(i,j) = term * term2 * (3 - OMEGA_VSS(i,j)) * (4 - OMEGA_VSS(i,j)) * gamma(3 - OMEGA_VSS(i,j))
		
			bb(i,j) = (5./3.*omega12(i,j) - 1./3.*omega13(i,j))/omega11(i,j)
	
			cc(i,j) = omega12(i,j) / omega11(i,j) / 3.0
			
			aa(i,j) = omega22(i,j) / omega11(i,j) / 2.0
		end do
	end do
	
	case ('Born-Mayer')

	do i=1,NUM_SP
		do j=1,NUM_SP
			sig_ij = (SIGMA_LJ(i)+SIGMA_LJ(j))/2
			mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E27)/(MASS_SPCS(i)*1E27 + MASS_SPCS(j)*1E27)
			eij = sqrt(EPS_LJ(i)*EPS_LJ(j)*(SIGMA_LJ(i)*1e10)**6 &
				 *(SIGMA_LJ(j)*1e10)**6)/(sig_ij*1e10)**6
			eij_eV = eij * kB_eV

			sij = PI*(sig_ij)**2*sqrt((Kb/mij)*T/2./pi) ! of 1e-7 and higher order

			T_star = T/eij
			nu_star = PHI_BM(i,j) / eij_eV
			r_star = 1 / (sig_ij * BETA_BM(i,j))

			term = 1 / (r_star * log(nu_star/10))**2
			term2 = 1 / log(nu_star/10)
			do k=1,6
				Ak_BM(k) = A_CONST_BM(1,k) + term * (A_CONST_BM(2,k) &
					    + A_CONST_BM(3,k) * term2 + A_CONST_BM(4,k) * term2**2 )
			end do
				
			! All the omega-ints are of 1e-6 and higher order
	
			term3 = (r_star * log(nu_star/T_star))**2
			term4 = log(T_star)

			omega11(i,j) = sij * term3 * (0.89 + Ak_BM(1)/T_star**2 +  Ak_BM(2)/T_star**4 +  Ak_BM(3)/T_star**6)
			! WRITE (*,*) 'Omega11_ij = ', omega11(i,j)
	
			omega22(i,j) = 2. * sij * term3 * (1.04 + Ak_BM(4)/term4**2 +  Ak_BM(5)/term4**3 +  Ak_BM(6)/term4**4)
			! WRITE (*,*) 'Omeg22_ij = ', omega22(i,j)
	
			omega12(i,j) = 2.5 * omega11(i,j) ! ToDo
	
			omega13(i,j) = 3.5 * omega12(i,j) ! ToDo
		
			bb(i,j) = (5./3.*omega12(i,j) - 1./3.*omega13(i,j))/omega11(i,j)
	
			cc(i,j) = omega12(i,j)/omega11(i,j)/3.0
			
			aa(i,j) = omega22(i,j)/omega11(i,j)/2.0
		end do
		end do

	case('ESA-Bruno')

	do i=1,NUM_SP
		do j=1,NUM_SP
			sig_ij = (SIGMA_LJ(i)+SIGMA_LJ(j))/2.
			mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E27)/(MASS_SPCS(i)*1E27 + MASS_SPCS(j)*1E27)
			sij = PI * (sig_ij)**2 * sqrt((Kb/mij)*T/2./pi) ! of 1e-7 and higher order

			x = log(kB_eV * T / (E_0_VSS(i,j) / 1.0e3)) ! the latter is done to convert meV to eV
			! WRITE (*,*) 'x in ESa model ', x
	
			! All the omega-ints are of 1e-6 and higher order
			
			! omega11
			do k=1,7
				aj_ESA(k) = OMEGA11_CONST_ESA(k,1) + OMEGA11_CONST_ESA(k,2)*BETA_VSS(i,j) &
							+ OMEGA11_CONST_ESA(k,3)*BETA_VSS(i,j)**2
			end do
			! WRITE (*,*) 'a_j polynomials ', (aj_ESA(k), k=1,7)

			term = (x - aj_ESA(3))/aj_ESA(4)
			term2 = (x - aj_ESA(6))/aj_ESA(7)
			term3 = aj_ESA(1) + aj_ESA(2)*x
			
			omega11(i,j) = sij * exp( term3 * exp(term)/ (exp(term) + exp(-1.*term)) &
						   + aj_ESA(5) * exp(term2)/ (exp(term2) + exp(-1.*term2)))
			! WRITE (*,*) 'Omega11_ij = ', omega11(i,j)
			! omega12
			do k=1,7
				aj_ESA(k) = OMEGA12_CONST_ESA(k,1) + OMEGA12_CONST_ESA(k,2)*BETA_VSS(i,j) &
							+ OMEGA12_CONST_ESA(k,3)*BETA_VSS(i,j)**2
			end do

			term = (x - aj_ESA(3))/aj_ESA(4)
			term2 = (x - aj_ESA(6))/aj_ESA(7)
			term3 = aj_ESA(1) + aj_ESA(2)*x
			omega12(i,j) = 3 * sij * exp( term3 * exp(term)/ (exp(term) + exp(-term)) &
						   + aj_ESA(5) * exp(term2)/ (exp(term2) + exp(-term2)))
			
			! omega13
			do k=1,7
				aj_ESA(k) = OMEGA13_CONST_ESA(k,1) + OMEGA13_CONST_ESA(k,2)*BETA_VSS(i,j) &
							+ OMEGA13_CONST_ESA(k,3)*BETA_VSS(i,j)**2
			end do

			term = (x - aj_ESA(3))/aj_ESA(4)
			term2 = (x - aj_ESA(6))/aj_ESA(7)
			term3 = aj_ESA(1) + aj_ESA(2)*x
			omega13(i,j) = 12 * sij * exp( term3 * exp(term)/ (exp(term) + exp(-term)) &
						   + aj_ESA(5) * exp(term2)/ (exp(term2) + exp(-term2)))

			! omega22
			do k=1,7
				aj_ESA(k) = OMEGA22_CONST_ESA(k,1) + OMEGA22_CONST_ESA(k,2)*BETA_VSS(i,j) &
							+ OMEGA22_CONST_ESA(k,3)*BETA_VSS(i,j)**2
			end do

			term = (x - aj_ESA(3))/aj_ESA(4)
			term2 = (x - aj_ESA(6))/aj_ESA(7)
			term3 = aj_ESA(1) + aj_ESA(2)*x
			omega22(i,j) = 2 * sij * exp( term3 * exp(term)/ (exp(term) + exp(-term)) &
						   + aj_ESA(5) * exp(term2)/ (exp(term2) + exp(-term2)))
		
			bb(i,j) = (5./3.*omega12(i,j) - 1./3.*omega13(i,j))/omega11(i,j)
	
			cc(i,j) = omega12(i,j)/omega11(i,j)/3.0
			
			aa(i,j) = omega22(i,j)/omega11(i,j)/2.0
		end do
		end do

	case ('Lennard-Jones')

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
	
		bb(i,j) = (5./3.*omega12(i,j) - 1./3.*omega13(i,j))/omega11(i,j)

		cc(i,j) = omega12(i,j)/omega11(i,j)/3.0
		
		aa(i,j) = omega22(i,j)/omega11(i,j)/2.0
	  end do
	end do

	end select
  
	omega_out%OMEGA11 = omega11
	omega_out%OMEGA22 = omega22
	omega_out%OMEGA12 = omega12
	omega_out%OMEGA13 = omega13
	omega_out%AA = AA
	omega_out%BB = BB
	omega_out%CC = CC
  
	end subroutine OmegaInt

end module omega_integrals