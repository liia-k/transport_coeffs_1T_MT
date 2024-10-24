
module defs_models

use constant_air5

type transport_in

! structure of macroparameters required for transport coeffs. calculation:

    real(8) rho  ! density
    !real(8) vel(1:3), not required
    real(8) temp ! temperature
    real(8) mass_fractions(1:NUM_SP) ! vector of mass fractions
end type

type transport_in_additional

! structure of macroparameters required for transport coeffs. calculation:

    real(8) rho  ! density
    real(8) ntot ! total number density
    real(8) temp ! temperature
    real(8) M ! gas molar mass
    real(8) mass_fractions(1:NUM_SP) ! vector of mass fractions
    real(8) num_fractions(1:NUM_SP) ! vector of number fractions
end type

!Transport coefficients: 

type transport_out

! structure of transport coeffs.:

    real(8) visc, bulk_visc, ltot ! shear viscosity,  bulk viscosity, thermal conductivity
    real(8) :: thdiff(1:NUM_SP) ! vector of thermal diffusion coeffs
    real(8) :: effDiff(1:NUM_SP) ! vector of effective diffusion coeffs
    real(8) :: DIFF(1:NUM_SP,1:NUM_SP) ! matrix of diffusion coeffs
end type

end module defs_models