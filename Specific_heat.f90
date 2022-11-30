!In this module, vibrational energy, non-equilibrium vibrational
!partition functions, vibrational specific heat capacities,
! and internal heat capacities c_int of molecular species are calculated.
!It uses the module constant.f90 containing main constants and 
!variables definition 

!Input variables -- primitive variables: T - temperature, y - mass fractions, rho - density

MODULE Specific_heat

USE CONSTANT

IMPLICIT NONE

! For SPECIFIC_HEAT module:

! Number of vibrational levels for N2, O2, NC

INTEGER :: L_N2, L_O2, L_NO
PARAMETER (L_N2 = 67, L_O2 = 56, L_NO = 52)

type cv_out
  Double precision cv_int, cv_tot
  real, dimension(5) :: c_int
end type


CONTAINS
SUBROUTINE S_Heat(data_in, c_out)
 
type(transport_in),intent(in)   :: data_in
type(cv_out),intent(out) :: c_out

! Arrays containing values of vibrational energy of N2, O2, NO

REAL, DIMENSION(0:67) :: EN_N2
REAL, DIMENSION(0:56) :: EN_O2
REAL, DIMENSION(0:52) :: EN_NO

Double precision zvibr_n2, zvibr_o2, zvibr_no, c_vibr_n2, c_vibr_o2, c_vibr_no &
                  , cv_int, cv_tot, ppp,s,s0

real, dimension(5) :: c_int, y

INTEGER I1

real T, EE, rho, M, M1

T = data_in%temp
rho = data_in%rho
y = data_in%mass_fractions

M1 = 0 ! 1/total molar mass
do i1=1,5
  M1 = M1 + y(i1)/molar(i1)
end do
M = 1/M1


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

	
  ! Calculation of non-equilibrium N2 partition function Z_N2(T)

    zvibr_N2=0 

    DO i1=0,L_N2
      ee=en_N2(i1)
      zvibr_N2=zvibr_N2+exp(-ee/kb/T)
    end do

! Calculation of non-equilibrium O2 partition function Z_O2(T)


  zvibr_O2=0 

  DO i1=0,L_O2
    ee=en_O2(i1)
    zvibr_O2=zvibr_O2+exp(-ee/kb/T)
  end do


! Calculation of non-equilibrium NO partition function Z_NO(T)

  zvibr_NO=0 
    
  DO i1=0,L_NO
    ee=en_NO(i1)
    zvibr_NO=zvibr_NO+exp(-ee/kb/T)
  end do

  ! Calculation of non-equilibrium O2 vibrational specific heat
  ! Cvibr_N2

  s=0; s0=0;

    DO i1=0,L_N2
        ee=en_N2(i1)
        ppp=exp(-ee/kb/T)/zvibr_N2
        s=s+ee/kb/T*ppp;
        s0=s0+ee*ee/kb/T/kb/T*ppp
    END DO
    c_vibr_n2=(s0-s*s)
    c_int(1)=c_vibr_n2 + kb/mass(1)

! Calculation of non-equilibrium O2 vibrational specific heat
! Cvibr_O2

s=0;s0=0;

  DO i1=0,L_O2
      ee=en_O2(i1)
      ppp=exp(-ee/kb/T)/zvibr_O2
      s=s+ee/kb/T*ppp;
      s0=s0+ee*ee/kb/T/kb/T*ppp;
  END DO
  c_vibr_o2=(s0-s*s)
  c_int(2)=c_vibr_o2 + kb/mass(2)

! Calculation of non-equilibrium CO vibrational specific heat
! Cvibr_NO

s=0;s0=0;

  DO i1=0,L_NO
      ee=en_NO(i1)
      ppp=exp(-ee/kb/T)/zvibr_NO
      s=s+ee/kb/T*ppp
      s0=s0+ee*ee/kb/T/kb/T*ppp
  END DO
  c_vibr_no=(s0-s*s)
  c_int(3)=c_vibr_no + kb/mass(3)


  ! Calculation of internal and total specific heats
  ! CV_int, CV_tot
  
  cv_int = y(1)*c_int(1) + y(2)*c_int(2) + y(3)*c_int(3)
  cv_tot = 3/2*R/M + cv_int!3/2*ntot*kb/rho + cv_int

  c_out%cv_int = cv_int
  c_out%cv_tot = cv_tot
  c_out%c_int = c_int
  
END SUBROUTINE S_Heat


END	MODULE Specific_heat


