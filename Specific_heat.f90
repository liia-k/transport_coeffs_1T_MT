!In this module, vibrational energy, non-equilibrium vibrational
!partition functions, vibrational specific heat capacities,
! and internal heat capacities c_int of molecular species are calculated.
!It uses the module constant.f90 containing main constants and 
!variables definition 

!Input variables: T - temperature, y - molar fractions, ntot - total number density, rho - density

MODULE Specific_heat

USE CONSTANT

IMPLICIT NONE

! For SPECIFIC_HEAT module:

! Number of vibrational levels in CO2 modes (1-3), O2, CO

INTEGER :: L_N2, L_O2, L_NO
DATA L_N2, L_O2, L_NO / 67, 46, 52 /
	
! Arrays containing values of vibrational energy of CO2, O2, CO

REAL, DIMENSION(0:67) :: EN_N2
REAL, DIMENSION(0:46) :: EN_O2
REAL, DIMENSION(0:52) :: EN_NO

Double precision zvibr_n2, zvibr_o2, zvibr_no, c_vibr_n2, c_vibr_o2, c_vibr_no &
                  , cv_int, cv_tot

real, dimension(5) :: c_int 

CONTAINS

  SUBROUTINE ENERGY

    INTEGER I1

    do i1=1,5
      c_int(i1)=0
    end do

  ! Calculation of vibrational energy levels for N2, O2, NO
    
      DO i1=0,L_N2
        en_N2(i1)=we_N2*(i1) - wexe_N2*i1**2
      END DO

      DO i1=0,L_O2
        en_O2(i1)=we_O2*(i1) - wexe_O2*i1**2
      END DO

      DO i1=0,L_NO
        en_NO(i1)=we_NO*(i1) - wexe_NO*i1**2
      END DO

  END SUBROUTINE ENERGY


  SUBROUTINE PART_FUNC_N2(T)
	
    ! Calculation of non-equilibrium N2 partition function Z_N2(T)
  
      INTEGER I1
      REAL EE,T
    
      zvibr_N2=0 
  
      DO i1=0,L_N2
        ee=en_N2(i1)
        zvibr_N2=zvibr_N2+exp(-ee/kb/T)
      end do
  
  END SUBROUTINE PART_FUNC_N2

  SUBROUTINE PART_FUNC_O2(T)
	
	! Calculation of non-equilibrium O2 partition function Z_O2(T)

    INTEGER I1
	  REAL EE,T
	
	  zvibr_O2=0 

    DO i1=0,L_O2
      ee=en_O2(i1)
      zvibr_O2=zvibr_O2+exp(-ee/kb/T)
    end do

  END SUBROUTINE PART_FUNC_O2

  SUBROUTINE PART_FUNC_NO(T)
	
	! Calculation of non-equilibrium NO partition function Z_NO(T)

    INTEGER I1
    REAL EE, T

    zvibr_NO=0 
      
    DO i1=0,L_NO
      ee=en_NO(i1)
      zvibr_NO=zvibr_NO+exp(-ee/kb/T)
    end do

  END SUBROUTINE PART_FUNC_NO

  SUBROUTINE s_heat_N2

    ! Calculation of non-equilibrium O2 vibrational specific heat
    ! Cvibr_N2
  
    integer I1
    real ee
    double precision ppp,s,s0
  
    s=0; s0=0;
  
      DO i1=0,L_N2
         ee=en_N2(i1)
         ppp=exp(-ee/kb/T)/zvibr_N2
         s=s+ee/kb/T*ppp;
         s0=s0+ee*ee/kb/T/kb/T*ppp
      END DO
      c_vibr_n2=(s0-s*s)
      c_int(1)=c_vibr_n2 + kb/mass(1)
    
  END SUBROUTINE s_heat_N2

  SUBROUTINE s_heat_O2

	! Calculation of non-equilibrium O2 vibrational specific heat
	! Cvibr_O2

	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_O2
       ee=en_O2(i1)
       ppp=exp(-ee/kb/T)/zvibr_O2
       s=s+ee/kb/T*ppp;
       s0=s0+ee*ee/kb/T/kb/T*ppp;
    END DO
    c_vibr_o2=(s0-s*s)
    c_int(2)=c_vibr_o2 + kb/mass(2)
	
  END SUBROUTINE s_heat_O2

  SUBROUTINE s_heat_NO
	! Calculation of non-equilibrium CO vibrational specific heat
	! Cvibr_NO
	
	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_NO
       ee=en_NO(i1)
       ppp=exp(-ee/kb/T)/zvibr_NO
       s=s+ee/kb/T*ppp
       s0=s0+ee*ee/kb/T/kb/T*ppp
    END DO
    c_vibr_no=(s0-s*s)
    c_int(3)=c_vibr_no + kb/mass(3)
	
  END SUBROUTINE s_heat_NO

  SUBROUTINE s_heat
    ! Calculation of internal and total specific heats
    ! CV_int, CV_tot
    
    cv_int = y(1)*c_int(1) + y(2)*c_int(2) + y(3)*c_int(3)
    cv_tot = 3/2*ntot*kb/rho + cv_int
    
  END SUBROUTINE s_heat


END	MODULE Specific_heat


