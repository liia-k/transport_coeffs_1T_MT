!In this module, vibrational energy, non-equilibrium vibrational
!partition functions, vibrational specific heat capacities,
! and internal heat capacities cv_int_sp of molecular species are calculated.
!It uses the module constant.f90 containing main constants and 
!variables definition 

!Input variables -- primitive variables: T - temperature, y - mass fractions

module specific_heat_sp

  use constant_air5
    
  implicit none

! For SPECIFIC_HEAT module:

type cv_out
  double precision cv_int_tot, cv_tot
  real, dimension(NUM_SP) :: cv_int_sp, cv_vibr_sp
end type

type dim
  real, allocatable :: vibr_en_sp(:)
end type

contains

subroutine VibrEn(temp, en_out)
  
  real, intent(in)   :: temp
  type(dim), dimension(NUM_MOL), intent(out) :: en_out

  integer i, j

  real T
  
  do i=1,NUM_MOL
      allocate(en_out(i)%vibr_en_sp(1 + L_vibr(i)))
  end do

  T = temp
  
  ! Calculation of vibrational energy levels for N2, O2, NO
  
  do i=1,NUM_MOL
      do j=0,L_vibr(i)
          en_out(i)%vibr_en_sp(j) = we(i)*j - wexe(i)*j**2
      end do
  end do

end subroutine VibrEn

subroutine SpHeat(temp, mass_fr, c_out)

  real,intent(in)   :: temp
  real, dimension(num_sp), intent(in) :: mass_fr
  type(cv_out),intent(out) :: c_out

  integer i
  
  ! Arrays containing values of vibrational energy of molecules

  type(dim), dimension(NUM_MOL) :: en_vibr

  double precision, dimension(NUM_MOL) :: zvibr
  
  double precision cv_int_tot, cv_tot, s,s0
  
  real, dimension(NUM_SP) :: cv_int_sp, cv_vibr_sp, y
  
  real T, M
  
  T = temp
  y = mass_fr
  
  M = 1/dot_product(y,1/MOLAR) ! total molar mass

  cv_vibr_sp = 0
  cv_int_sp = 0 ! initial zero internal specific heats
  
  ! Calculation of vibrational energy levels for N2, O2, NO
  
  call VibrEn(T, en_vibr)

  
  ! Calculation of non-equilibrium partition functions Z_c(T)

  do i=1,NUM_MOL
    zvibr = SUM(exp(-en_vibr(i)%vibr_en_sp/(Kb*T)))
  end do
    
  ! Calculation of non-equilibrium vibrational specific heats
  
  do i=1,NUM_MOL
    s = SUM(en_vibr(i)%vibr_en_sp*exp(-en_vibr(i)%vibr_en_sp/(Kb*T))/zvibr(i)/(Kb*T))
    s0 = SUM((en_vibr(i)%vibr_en_sp/(Kb*T))**2*exp(-en_vibr(i)%vibr_en_sp/(Kb*T))/zvibr(i))
    cv_vibr_sp(i) = s0 - s*s
  end do
  
  ! Calculation of internal and total specific heats
  ! cv_int_tot, CV_tot
  
  cv_int_sp(1:NUM_MOL) = cv_vibr_sp(1:NUM_MOL) + Kb/MASS_SPCS(1:NUM_MOL)

  cv_int_tot = dot_product(y, cv_int_sp)
  cv_tot = (3./2.)*R/M + cv_int_tot!3/2*ntot*kb/rho + cv_int_tot

  c_out%cv_int_tot = cv_int_tot
  c_out%cv_tot = cv_tot
  c_out%cv_int_sp = cv_int_sp
  c_out%cv_vibr_sp = cv_vibr_sp
    
end subroutine SpHeat
    
end module specific_heat_sp
