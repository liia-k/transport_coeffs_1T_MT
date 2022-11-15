!In this module, vibrational energy, non-equilibrium vibrational
!partition functions, vibrational specific heat capacities,
! and mean vibrational energy en_int of molecular species are calculated.
!It uses the module constant.f90 containing main constants and 
!variables definition 
!Input variables: T

MODULE Specific_heat

USE CONSTANT

IMPLICIT NONE

Double precision zv_n2, zv_o2, zv_no, c_v_n2, c_v_o2, c_v_no

real, dimension(5) :: en_int 

INTEGER I1,I2,I3,IL,G

CONTAINS

  SUBROUTINE ENERGY

    do i1=1,5
      en_int(i1)=0
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
    
      zv_N2=0 
  
      DO i1=0,L_N2
        ee=en_N2(i1)
        zv_N2=zv_N2+exp(-ee/kb/T)
      end do
  
  END SUBROUTINE PART_FUNC_N2

  SUBROUTINE PART_FUNC_O2(T)
	
	! Calculation of non-equilibrium O2 partition function Z_O2(T)

    INTEGER I1
	  REAL EE,T
	
	  zv_O2=0 

    DO i1=0,L_O2
      ee=en_O2(i1)
      zv_O2=zv_O2+exp(-ee/kb/T)
    end do

  END SUBROUTINE PART_FUNC_O2

  SUBROUTINE PART_FUNC_NO(T)
	
	! Calculation of non-equilibrium NO partition function Z_NO(T)

    INTEGER I1
    REAL EE, T

    zv_NO=0 
      
    DO i1=0,L_NO
      ee=en_NO(i1)
      zv_NO=zv_NO+exp(-ee/kb/T)
    end do

  END SUBROUTINE PART_FUNC_NO

  SUBROUTINE s_heat_N2

    ! Calculation of non-equilibrium O2 vibrational specific heat
    ! CV_N2
  
    integer I1
    real ee
    double precision ppp,s,s0
  
    s=0;s0=0;
  
      DO i1=0,L_N2
         ee=en_N2(i1)
         ppp=exp(-ee/kb/T)/zv_N2
         s=s+ee/kb/T*ppp;
         s0=s0+ee*ee/kb/T/kb/T*ppp;
       en_int(1)=en_int(1)+ee*ppp
      END DO
      c_v_n2=(s0-s*s)
    
  END SUBROUTINE s_heat_N2

  SUBROUTINE s_heat_O2

	! Calculation of non-equilibrium O2 vibrational specific heat
	! CV_O2

	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_O2
       ee=en_O2(i1)
       ppp=exp(-ee/kb/T)/zv_O2
       s=s+ee/kb/T*ppp;
       s0=s0+ee*ee/kb/T/kb/T*ppp;
	   en_int(2)=en_int(2)+ee*ppp
    END DO
    c_v_o2=(s0-s*s)
	
  END SUBROUTINE s_heat_O2

  SUBROUTINE s_heat_NO
	! Calculation of non-equilibrium CO vibrational specific heat
	! CV_NO
	
	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_NO
       ee=en_NO(i1)
       ppp=exp(-ee/kb/T)/zv_NO
       s=s+ee/kb/T*ppp
       s0=s0+ee*ee/kb/T/kb/T*ppp
       en_int(3)=en_int(3)+ee*ppp
    END DO
    c_v_no=(s0-s*s)
	
  END SUBROUTINE s_heat_NO


END	MODULE Specific_heat


