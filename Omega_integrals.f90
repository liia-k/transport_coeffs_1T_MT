!In this module, Omega-integrals and their ratios AA, BB, CC
!are calculated using the Lennard-Jones potential for 
!moderate temperatures 
!It uses the module constant.f90 containing main constants and 
!variables definition 

!Input variable: T


MODULE OMEGA_INTEGRALS

USE CONSTANT

IMPLICIT NONE

!Omega-integrals and their ratios

type omega_int
	REAL, DIMENSION(5,5) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
						AA, BB, CC
end type



CONTAINS

  SUBROUTINE OMEGA(T, omega_out)

  real,intent(in)   :: T
  type(omega_int),intent(out) :: omega_out

  ! Calculation of OMEGA-integrals for given T

  INTEGER i, j
  REAL x11, eij, mij, sig_ij, sij, tx

  !Omega-integrals and their ratios

  REAL, DIMENSION(5,5) :: OMEGA11, OMEGA22, OMEGA12, OMEGA13, &
  AA, BB, CC

 !Parameters sigma(i) and eps(i) of the Lenard-Jones potential are
 !defined in module constant.f90. Data given by:
 !R.J.Kee, J.A.Miller, T.N. Jefferson, CHEMKIN: A General-Purpose, 
 !Problem-Independent, Transportable, Fortran Chemical Kinetics Code Package,
 !Sandia National Laboratories, SAND80-8003, 1980 

  DO i=1,5
	DO j=1,5
  		sig_ij=(sigma(i)+sigma(j))/2
		mij=mass(j)*(mass(i)*1E27)/(mass(i)*1E27 + mass(j)*1E27)
		eij=sqrt(eps(i)*eps(j)*(sigma(i)*1e10)**6 &
			*(sigma(j)*1e10)**6)/(sig_ij*1e10)**6
		sij=pi*(sig_ij)**2*sqrt(kb*t/2.0/pi/mij)
		!WRITE (*,*) 'i = ', i 
		!WRITE (*,*) 'j = ', j
		!WRITE (*,*) 'e_ij = ', eij
		!WRITE (6,'(1x, A15, E50.40)') 'm_ij = ', mij
		!WRITE (*,*) 'Omega_rs_ij = ', sij
	    
		tx=t/eij
			! Omega11, Omega12, Omega12, Omega13, B, C for 
			! moderate temperatures; Lennard-Jones potential
			x11=log((tx))+1.4
			omega11(i,j)=1/(-0.16845-2.25768e-2/x11/x11+ &
			0.19779/x11+0.64373*x11-9.26718e-2*x11*x11+ &
			7.1131e-3*x11**3)*sij
			!WRITE (*,*) 'Omega11_ij = ', omega11(i,j)

			x11=log((tx))+1.5
			omega22(i,j)=1/(-0.40811-5.08552e-2/x11/x11+ &
			0.3401/x11+0.70375*x11-0.10699*x11*x11+ &
			7.62686e-3*x11**3)*sij*2.0
			!WRITE (*,*) 'Omega22_ij = ', omega22(i,j)

			x11=log((tx))+1.1
			omega12(i,j)=1/(0.40785+9.25303e-4/x11/x11+ &
			2.79680e-4/x11+0.44739*x11-6.27242e-2*x11*x11+ &
			5.98567e-3*x11**3)*sij*3.0
			!WRITE (*,*) 'Omega12_ij = ', omega12(i,j)

			x11=log((tx))+4.0
			omega13(i,j)=1/(25.04929+63.12444/x11/x11- &
			65.87398/x11-4.13758*x11+0.34999*x11*x11- &
			1.0096e-2*x11**3)*sij*12.0
			!WRITE (*,*) 'Omega13_ij = ', omega13(i,j)
		
			bb(i,j)=(5./3.*omega12(i,j)-4./12.*omega13(i,j)) &
			/omega11(i,j)

			cc(i,j)=1/omega11(i,j)*omega12(i,j)/3.0
		
		aa(i,j)=1/omega11(i,j)*omega22(i,j)/2.0
			
	    	 
	END DO
  END DO

  omega_out%OMEGA11 = omega11
  omega_out%OMEGA22 = omega22
  omega_out%OMEGA12 = omega12
  omega_out%OMEGA13 = omega13
  omega_out%AA = AA
  omega_out%BB = BB
  omega_out%CC = CC

  END SUBROUTINE OMEGA

END MODULE OMEGA_INTEGRALS