!In this module, bracket integrals in the transport linear systems are calculated.

module bracket_integrals

use constant_air5
use specific_heat_sp

implicit none

! Bracket integrals 

type bracket_int
	real(8), dimension(NUM_SP, NUM_SP) :: LAMBDA0000, LAMBDA0100, LAMBDA1100						
	real(8), dimension(NUM_SP, NUM_SP) :: H00, BETA1100
	real(8), dimension(NUM_SP, NUM_MOL) :: BETA0110
	real(8), dimension(NUM_MOL) :: beta0011
	real(8), dimension(NUM_SP) :: lambda_int
end type

contains

	subroutine BracketInt(T, ntot, x, omega_in, cv_in, bracket_out)

	use omega_integrals

	real, intent(in) :: T ! temperature
	real, intent(in) :: ntot ! number density
	real, dimension(NUM_SP), intent(in) :: x ! molar fractions
	type(omega_int), intent(in) :: omega_in
	type(SpHeatVOut), intent(in) :: cv_in
	type(bracket_int), intent(out)   :: bracket_out
	

	integer i, j, k
	real(8) z_rot, F

	! Supporting bracket integrals
	real(8), dimension(NUM_SP,NUM_SP) :: diff, eta, red_mass_sp ! diffusion coefficients, viscosity coefficients, reduced mass for species
	real(8), dimension(NUM_SP) :: zeta, phi ! bulk viscosity for species (fictional), term in beta integrals
	real(8), dimension(NUM_MOL) :: int_p_rel ! relaxation time * pressure, s*Pa

	real(8), dimension(NUM_SP,NUM_SP) :: LAMBDA0000, LAMBDA0100, LAMBDA1100						
	real(8), dimension(NUM_SP,NUM_SP) :: H00, BETA1100
	real(8), dimension(NUM_SP,NUM_MOL) :: BETA0110
	real(8), dimension(NUM_MOL) :: beta0011

	!internal heat conductivity coefficients (lambda_int)

	real(8), dimension(NUM_SP) :: lambda_int

	real, dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
										AA, BB, CC


	OMEGA11 = omega_in%omega11
	OMEGA22 = omega_in%omega22
	OMEGA12 = omega_in%omega12
	OMEGA13 = omega_in%omega13
	AA = omega_in%aa
	BB = omega_in%bb
	CC = omega_in%cc
	
	LAMBDA0000 = 0.0
	LAMBDA0100 = 0.0
	LAMBDA1100 = 0.0
	H00 = 0.0
	BETA1100 = 0.0
	BETA0110 = 0.0
	BETA0011 = 0.0
	beta0011 = 0.0
	lambda_int = 0.0
	zeta = 0.0
	phi = 0.0

	do i = 1, NUM_SP
		do j = 1, NUM_SP
			red_mass_sp(i,j) = MASS_SPCS(j) * (MASS_SPCS(i) * 1E26) / (MASS_SPCS(i) * 1E26 + MASS_SPCS(j) * 1E26)
			diff(i,j) = (3. / 16.) * T * (Kb / (red_mass_sp(i,j) * ntot)) / OMEGA11(i,j) 
			eta(i,j) = (5. / 8.) * Kb * T / OMEGA22(i,j)
		end do
	end do

	! Parker's formula for relaxation time:
	do i = 1, NUM_MOL
		F = 1 + PI**(3. / 2.) / 2. / (T / EPS_LJ(i))**(1. / 2.) + (PI**2. / 4. + 2.) / (T / EPS_LJ(i)) &
		    + PI**(3. / 2.) / (T / EPS_LJ(i))**(3. / 2.)
		z_rot = Z_INF(i) / F
		int_p_rel(i) = (PI / 4.) * eta(i,i) * z_rot
		! Bulk viscocity for species
		zeta(i) = (4. / PI) * int_p_rel(i) / eta(i,i)
		phi(i) = MASS_SPCS(i) * cv_in%cv_int_sp(i) / zeta(i)
	end do
	
	! Calculation of bracket integrals
	do i = 1, NUM_SP
		do j = 1, NUM_SP
		if (i == j) then
			lambda0000(i,j) = 0
			lambda0100(i,j) = 0
			lambda1100(i,j) = 4. * x(i)**2 *MASS_SPCS(i) * omega22(i,i)
			h00(i,j) = x(i)**2 / eta(i,i)
			beta1100(i,j) = (4. * T / PI) * phi(i) * x(i)**2 / eta(i,i)
			do k = 1, NUM_SP 
				if (k .NE. i) then 
				lambda0000(i,j) = lambda0000(i,j) + 8. * x(i) * x(k) * red_mass_sp(i,k) * omega11(i,k)
				lambda0100(i,j) = lambda0100(i,j) -  4. * x(i) * x(k) * red_mass_sp(i,k) &
								* (MASS_SPCS(i) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(k) * 1e27) &
								* (2 * omega12(i,k) - 5 * omega11(i,k))

				lambda1100(i,j) = lambda1100(i,j) + 8. * x(i) * x(k) * red_mass_sp(i,k) &
								  * ((red_mass_sp(i,k) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(k) * 1e27)) &
								  * ((15. / 2.) * (MASS_SPCS(i) * 1e27) / (MASS_SPCS(k) * 1e27) * omega11(i,k)  &
								  + (25. / 4.) * (MASS_SPCS(k) * 1e27) / (MASS_SPCS(i) * 1e27) * omega11(i,k) &
								  + (MASS_SPCS(k) * 1e27) / (MASS_SPCS(i) * 1e27) * (omega13(i,k) - 5. * omega12(i,k)) &
								  + 2. * omega22(i,k)) 
				
				h00(i,j) = h00(i,j) + 2. * (x(i) * x(k) / eta(i,k)) * (MASS_SPCS(i) * 1e27) * (MASS_SPCS(k) * 1e27) &
						   / (MASS_SPCS(i) * 1e27 + MASS_SPCS(k) * 1e27)**2 & 
						   * (5. / 3. / aa(i,k) + (MASS_SPCS(k) * 1e27) / (MASS_SPCS(i) * 1e27))

				beta1100(i,j) = beta1100(i,j) + (x(i) * x(k) / eta(i,k)) * (MASS_SPCS(k) * 1e27) &
								* (MASS_SPCS(i) * 1e27) / (MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2 &
								* (5. * Kb * T / aa(i,k) + 4. * T / PI * (MASS_SPCS(k) * 1e27)  &
								/ (MASS_SPCS(i) * 1e27) * (phi(i) + phi(k)))
				
				end if
			end do
			
		else
			lambda0000(i,j) = - 8. * x(i) * x(j) * red_mass_sp(i,j) * omega11(i,j)
			lambda0100(i,j) = 4. * x(i) * x(j) * red_mass_sp(i,j) &
							  * (MASS_SPCS(i) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(j) * 1e27) &
							  * (2 * omega12(i,j) - 5 * omega11(i,j))

			lambda1100(i,j) = - 8. * x(i) * x(j) * red_mass_sp(i,j) &
								  * ((red_mass_sp(i,j) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(j) * 1e27)) &
								  * ((55. / 4.) * omega11(i,j)  - 5. * omega12(i,j) + omega13(i,j) - 2. * omega22(i,j))

			h00(i,j) = 2. * (x(i) * x(j) / eta(i,j)) * (MASS_SPCS(i) * 1e27) * (MASS_SPCS(j) * 1e27) &
						   / (MASS_SPCS(i) * 1e27 + MASS_SPCS(j) * 1e27)**2 & 
						   * (- 5. / 3. / aa(i,j) + 1.)
			
			beta1100(i,j) = (x(i) * x(j) / eta(i,j)) * (MASS_SPCS(j) * 1e27) &
							* (MASS_SPCS(i) * 1e27) / (MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)**2 &
							* (- 5. * Kb * T / aa(i,j) + 4. * T / PI * (phi(i) + phi(j)))

		end if
		lambda_int(i) = lambda_int(i) + x(j) * red_mass_sp(i,j) * omega11(i,j) 
		end do
	end do

	do i = 1, NUM_SP
		do j = 1, NUM_MOL
		if(i == j) then
			beta0110(i,j) = -(4. * T / PI) * x(i)**2 * phi(i) / eta(i,i) 
			do k = 1, NUM_SP 
				if (k .NE. i) then 
				beta0110(i,j) = beta0110(i,j) - (4. * T / PI) * x(i) * x(k) * phi(i) / eta(i,k) &
								* (MASS_SPCS(k) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(k) * 1e27)
				end if
			end do
		else
			beta0110(i,j) = -(4. * T / PI) * x(i) * x(j) * phi(j) / eta(i,j) &
							*(MASS_SPCS(j) * 1e27) / (MASS_SPCS(i) * 1e27 + MASS_SPCS(j) * 1e27)
		
		end if
		beta0011(j) = beta0011(j) + (4. * T / PI) * x(i) * x(j) * phi(j) / eta(i,j)
		end do
	end do

	bracket_out%LAMBDA0000 = LAMBDA0000
	bracket_out%LAMBDA0100 = LAMBDA0100
	bracket_out%LAMBDA1100 = LAMBDA1100
	bracket_out%H00 = H00
	bracket_out%BETA1100 = BETA1100
	bracket_out%BETA0110 = BETA0110
	bracket_out%BETA0011 = BETA0011
	bracket_out%LAMBDA_INT = LAMBDA_INT

	end subroutine BracketInt


end module bracket_integrals