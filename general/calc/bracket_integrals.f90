!In this module, bracket integrals in the transport linear
!systems are calculated.
!It uses modules constant.f90, Specific_heat.f90, Omega_integrals.f90
!Input variables: X(i),T


module bracket_integrals

	use constant_air5
	
	implicit none
	
	! Bracket integrals 
	
	type bracket_int
		real, dimension(NUM_SP,NUM_SP) :: LAMBDA, LAMBDA00, LAMBDA01, LAMBDA11						
		
		real, dimension(NUM_SP,NUM_SP) :: ETA, H00, BETA11
		
		real, dimension(NUM_SP,NUM_MOL) :: BETA01
		
		real, dimension(NUM_MOL) :: BETA0011
				
		real, dimension(NUM_SP) :: LAMBDA_INT=0
	end type
	
	contains
	
	  subroutine BracketInt(T, x, omega_in, bracket_out)

		use omega_integrals
	
		real,intent(in) :: T ! temperature
		real, dimension(NUM_SP), intent(in) :: x ! molar fractions
		type(omega_int), intent(in)      :: omega_in
		type(bracket_int), intent(out)   :: bracket_out
	
		integer i, j, k
		real mij
		real, dimension(NUM_SP) :: PHI
	
		real, dimension(NUM_SP,NUM_SP) :: LAMBDA, LAMBDA00, LAMBDA01, LAMBDA11						
		real, dimension(NUM_SP,NUM_SP) :: ETA, H00, BETA11
		real, dimension(NUM_SP,NUM_MOL) :: BETA01
		real, dimension(NUM_MOL) :: BETA0011
	
		!internal heat conductivity coefficients (LAMBDA_INT)
	
		real, dimension(NUM_SP) :: LAMBDA_INT
	
		real, dimension(NUM_SP,NUM_SP) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
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
				!WRITE (*,*) 'i = ', i 
				!WRITE (*,*) 'j = ', j
				!WRITE (*,*) 'm_ij = ', mij
				!WRITE (*,*) 'lambda_ij = ', lambda(i,j)
				!WRITE (*,*) 'eta_ij = ', eta(i,j)
			end do
		end do
		
		phi(1:NUM_MOL) = Kb/CONST_v*(1. + (PI**(3./2.)/2.)*sqrt(EPS_LJ(1:NUM_MOL)/T) &
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
				lambda00(i,j) = 0
				lambda01(i,j) = 0
				lambda11(i,j) = (x(i))**2/lambda(i,i)
				h00(i,j) = (x(i))**2/eta(i,i)
				beta11(i,j) = 4.*T*phi(i)*x(i)**2/eta(i,i)/PI
				do k=1,NUM_SP 
					if(k.NE.i) then 
					lambda00(i,j) = lambda00(i,j) + x(i)*x(k)/lambda(i,k)/2./aa(i,k)
					lambda01(i,j) = lambda01(i,j) - x(i)*x(k) &
									/lambda(i,k)/4./aa(i,k)*(MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27) &
									*(6.*cc(i,k)-5.)
	
					lambda11(i,j) = lambda11(i,j) + x(i)*x(k)/lambda(i,k)/2./aa(i,k) &
									*((15./2.)*(MASS_SPCS(i)*1e27)**2 + (25./4.)*(MASS_SPCS(k)*1e27)**2 &
									- 3.*(MASS_SPCS(k)*1e27)**2*bb(i,k) + 4.*(MASS_SPCS(i)*1e27)*(MASS_SPCS(k)*1e27)*aa(i,k)) &
									/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2
					
					h00(i,j) = h00(i,j) + 2.*(x(i)*x(k)/eta(i,k))*(MASS_SPCS(i)*1e27)*(MASS_SPCS(k)*1e27) &
							   /(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2 & 
							   *(5./3./aa(i,k) + (MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27))
	
					beta11(i,j) = beta11(i,j) + (x(i)*x(k)/eta(i,k))*(MASS_SPCS(k)*1e27) &
								  *(5.*Kb*T*(MASS_SPCS(i)*1e27)/aa(i,k) + 4.*(phi(i)+phi(k))*T*(MASS_SPCS(k)*1e27)/PI) &
								  /(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)**2
					
					end if
				end do
				!WRITE (*,*) 'i = ', i 
				!WRITE (*,*) 'j = ', j
				!WRITE (*,*) 'lambda00(i,j) = ', lambda00(i,j)
				!WRITE (*,*) 'lambda01(i,j) = ', lambda01(i,j)
				!WRITE (*,*) 'lambda11(i,j) = ', lambda11(i,j)
				!WRITE (*,*) 'h00(i,j) =      ', h00(i,j)
				!WRITE (*,*) 'beta11(i,j) =   ', beta11(i,j)
			else
				lambda00(i,j) = -x(i)*x(j)/lambda(i,j)/2./aa(i,j)
				lambda01(i,j) = x(i)*x(j)*(6.*cc(i,j) - 5.)/4./aa(i,j)/lambda(i,j) &
								*(MASS_SPCS(i)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)
	
				lambda11(i,j) = -x(i)*x(j)*(55./4. - 3.*bb(i,j) - 4.*aa(i,j))/lambda(i,j)/2./aa(i,j) &
								*mij/(MASS_SPCS(i) + MASS_SPCS(j))
				h00(i,j) = -2.*x(i)*x(j)*(5./3./aa(i,j) - 1.)/eta(i,j) &
							*(mij*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)
				
				beta11(i,j) = x(i)*x(j)*(-5.*Kb*T/aa(i,j) + 4.*T*(phi(i)+phi(j))/PI)/eta(i,j) &
							  *(MASS_SPCS(i)*1e27)*(MASS_SPCS(j)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)**2
	
				!WRITE (*,*) 'i = ', i 
				!WRITE (*,*) 'j = ', j
				!WRITE (*,*) 'lambda00(i,j) = ', lambda00(i,j)
				!WRITE (*,*) 'lambda01(i,j) = ', lambda01(i,j)
				!WRITE (*,*) 'lambda11(i,j) = ', lambda11(i,j)
				!WRITE (*,*) 'h00(i,j) =      ', h00(i,j)
				!WRITE (*,*) 'beta11(i,j) =   ', beta11(i,j)
	
			end if
	
			lambda_int(i) = lambda_int(i) + x(j)*mij*omega11(i,j) 
	
			!WRITE (*,*) 'i = ', i 
			!WRITE (*,*) 'lambda_int(i) = ', lambda_int(i)
	
			end do
		end do
	
		  do i=1,NUM_SP
			do j=1,NUM_MOL
			if(i==j) then
				beta01(i,j) = -4.*T*x(i)**2*phi(i)/eta(i,i)/PI
				do k=1,NUM_SP 
					if(k.NE.i) then 
					beta01(i,j) = beta01(i,j) - 4.*T*x(i)*x(k)*phi(i)/eta(i,k)/PI &
								  *(MASS_SPCS(k)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(k)*1e27)
					end if
				end do
				!WRITE (*,*) 'i = ', i 
				!WRITE (*,*) 'j = ', j
				!WRITE (*,*) 'beta01(i,j) =   ', beta01(i,j)
			else
				beta01(i,j) = -4.*T*x(i)*x(j)*phi(j)/eta(i,j)/PI &
							  *(MASS_SPCS(j)*1e27)/(MASS_SPCS(i)*1e27 + MASS_SPCS(j)*1e27)
				!WRITE (*,*) 'i = ', i 
				!WRITE (*,*) 'j = ', j
				!WRITE (*,*) 'beta01(i,j) =   ', beta01(i,j)
			end if
	
			beta0011(j) = beta0011(j) + 4.*T*x(i)*phi(j)*x(j)/eta(i,j)/PI
	
			end do
		end do
	
	
		bracket_out%LAMBDA = LAMBDA
		bracket_out%LAMBDA00 = LAMBDA00
		bracket_out%LAMBDA01 = LAMBDA01
		bracket_out%LAMBDA11 = LAMBDA11
		bracket_out%ETA = ETA
		bracket_out%H00 = H00
		bracket_out%BETA11 = BETA11
		bracket_out%BETA01 = BETA01
		bracket_out%BETA0011 = BETA0011
		bracket_out%LAMBDA_INT = LAMBDA_INT
	
	  end subroutine BracketInt
	
	
	end module bracket_integrals