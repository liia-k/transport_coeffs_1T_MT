!In this module, bracket integrals in the transport linear systems are calculated.

module bracket_integrals

use constant_air5

implicit none

! Bracket integrals 

type bracket_int
	real(8), dimension(NUM_SP,NUM_SP) :: LAMBDA0000, LAMBDA0100, LAMBDA1100						
	real(8), dimension(NUM_SP,NUM_SP) :: H00, BETA1100
	real(8), dimension(NUM_SP,NUM_MOL) :: BETA0110
	real(8), dimension(NUM_MOL) :: beta0011
	real(8), dimension(NUM_SP) :: lambda_int
end type

contains

	subroutine BracketInt(T, x, omega_in, bracket_out)

	use omega_integrals

	real(8),intent(in) :: T ! temperature
	real(8), dimension(NUM_SP), intent(in) :: x ! molar fractions
	type(omega_int), intent(in)      :: omega_in
	type(bracket_int), intent(out)   :: bracket_out
	

	integer i, j, k
	real(8) mij
	real(8), dimension(NUM_SP) :: PHI

	real(8), dimension(NUM_SP,NUM_SP) :: LAMBDA, LAMBDA0000, LAMBDA0100, LAMBDA1100						
	real(8), dimension(NUM_SP,NUM_SP) :: ETA, H00, BETA1100
	real(8), dimension(NUM_SP,NUM_MOL) :: BETA0110
	real(8), dimension(NUM_MOL) :: beta0011

	!internal heat conductivity coefficients (lambda_int)

	real(8), dimension(NUM_SP) :: lambda_int

	real(8), dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
										AA, BB, CC


	OMEGA11 = omega_in%omega11
	OMEGA22 = omega_in%omega22
	OMEGA12 = omega_in%omega12
	OMEGA13 = omega_in%omega13
	AA = omega_in%aa
	BB = omega_in%bb
	CC = omega_in%cc
	
	phi = 0
	beta0011 = 0
	lambda_int = 0

	do i=1,NUM_SP
		do j=1,NUM_SP
			mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E26)/(MASS_SPCS(i)*1E26 + MASS_SPCS(j)*1E26)
			lambda(i,j) = (75./64.)*Kb*T*(Kb/mij)/OMEGA22(i,j) ! of 1e-11 order
			eta(i,j) = 5./8.*Kb*T/OMEGA22(i,j) ! of 1e-14 order
		end do
	end do
	
	phi(1:NUM_MOL) = Kb/Z_INF*(1. + (PI**(3./2.)/2.)*sqrt(EPS_LJ(1:NUM_MOL)/T) &
						+ (PI*PI/4.+2.)*EPS_LJ(1:NUM_MOL)/T + (PI*EPS_LJ(1:NUM_MOL)/T)**(3./2.)) ! of 1e-22 order

	! PHI(1) = Kb/23.73*(1. + (PI**(3./2.)/2.)*SQRT(EPS_LJ(1)/T) + (PI*PI/4.+2.)*EPS_LJ(1)/T &
	! 		 + (PI*EPS_LJ(1)/T)**(3./2.))
	! PHI(2) = Kb/20.72*(1. + (PI**(3./2.)/2.)*SQRT(EPS_LJ(2)/T) + (PI*PI/4.+2.)*EPS_LJ(2)/T &
	! 		 + (PI*EPS_LJ(2)/T)**(3./2.))
	! PHI(3) = Kb/9.16*(1. + (PI**(3./2.)/2.)*SQRT(EPS_LJ(3)/T) + (PI*PI/4.+2.)*EPS_LJ(3)/T &
	! 		 + (PI*EPS_LJ(3)/T)**(3./2.))


	do i=1,NUM_SP
		do j=1,NUM_SP
		mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E26)/(MASS_SPCS(i)*1E26 + MASS_SPCS(j)*1E26)
		if(i==j) then
			lambda0000(i,j) = 0
			lambda0100(i,j) = 0
			lambda1100(i,j) = (x(i))**2/lambda(i,i)
			h00(i,j) = (x(i))**2/eta(i,i)
			beta1100(i,j) = 4.*T*phi(i)*x(i)**2/eta(i,i)/PI
			do k=1,NUM_SP 
				if(k.NE.i) then 
				lambda0000(i,j) = lambda0000(i,j) + x(i)*x(k)/lambda(i,k)/2./aa(i,k)
				lambda0100(i,j) = lambda0100(i,j) - x(i)*x(k) &
								/lambda(i,k)/4./aa(i,k)*(MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27) &
								*(6.*cc(i,k)-5.)

				lambda1100(i,j) = lambda1100(i,j) + x(i)*x(k)/lambda(i,k)/2./aa(i,k) &
								*((15./2.)*(MASS_SPCS(i)*1e27)**2 + (25./4.)*(MASS_SPCS(k)*1e27)**2 &
								- 3.*(MASS_SPCS(k)*1e27)**2*bb(i,k) + 4.*(MASS_SPCS(i)*1e27)*(MASS_SPCS(k)*1e27)*aa(i,k)) &
								/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2
				
				h00(i,j) = h00(i,j) + 2.*(x(i)*x(k)/eta(i,k))*(MASS_SPCS(i)*1e27)*(MASS_SPCS(k)*1e27) &
							/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2 & 
							*(5./3./aa(i,k) + (MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27))

				beta1100(i,j) = beta1100(i,j) + (x(i)*x(k)/eta(i,k))*(MASS_SPCS(k)*1e27) &
								*(5.*Kb*T*(MASS_SPCS(i)*1e27)/aa(i,k) + 4.*(phi(i)+phi(k))*T*(MASS_SPCS(k)*1e27)/PI) &
								/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2
				
				end if
			end do
			
		else
			lambda0000(i,j) = -x(i)*x(j)/lambda(i,j)/2./aa(i,j)
			lambda0100(i,j) = x(i)*x(j)*(6.*cc(i,j) - 5.)/4./aa(i,j)/lambda(i,j) &
							*(MASS_SPCS(i)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)

			lambda1100(i,j) = -x(i)*x(j)*(55./4. - 3.*bb(i,j) - 4.*aa(i,j))/lambda(i,j)/2./aa(i,j) &
							*mij/(MASS_SPCS(i) + MASS_SPCS(j))
			h00(i,j) = -2.*x(i)*x(j)*(5./3./aa(i,j) - 1.)/eta(i,j) &
						*(mij*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)
			
			beta1100(i,j) = x(i)*x(j)*(-5.*Kb*T/aa(i,j) + 4.*T*(phi(i)+phi(j))/PI)/eta(i,j) &
							*(MASS_SPCS(i)*1e27)*(MASS_SPCS(j)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)**2

		end if

		lambda_int(i) = lambda_int(i) + x(j)*mij*omega11(i,j) 

		end do
	end do

	do i=1,NUM_SP
		do j=1,NUM_MOL
		if(i==j) then
			beta0110(i,j) = -4.*T*x(i)**2*phi(i)/eta(i,i)/PI
			do k=1,NUM_SP 
				if(k.NE.i) then 
				beta0110(i,j) = beta0110(i,j) - 4.*T*x(i)*x(k)*phi(i)/eta(i,k)/PI &
								*(MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)
				end if
			end do

		else
			beta0110(i,j) = -4.*T*x(i)*x(j)*phi(j)/eta(i,j)/PI &
							*(MASS_SPCS(j)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)
		
		end if

		beta0011(j) = beta0011(j) + 4.*T*x(i)*phi(j)*x(j)/eta(i,j)/PI

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