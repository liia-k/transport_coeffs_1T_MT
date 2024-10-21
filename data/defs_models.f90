
module defs_models

use constant_air5

type transport_in

! structure of macroparameters required for transport coeffs. calculation:

    real rho  ! density
    !real vel(1:3), not required
    real temp ! temperature
    real mass_fractions(1:NUM_SP) ! vector of mass fractions
end type

type transport_in_additional

! structure of macroparameters required for transport coeffs. calculation:

    real rho  ! density
    real ntot ! total number density
    real temp ! temperature
    real M ! gas molar mass
    real mass_fractions(1:NUM_SP) ! vector of mass fractions
    real num_fractions(1:NUM_SP) ! vector of number fractions
end type

!Transport coefficients: 

type transport_out

! structure of transport coeffs.:

    real visc, bulk_visc, ltot ! shear viscosity,  bulk viscosity, thermal conductivity
    real :: thdiff(1:NUM_SP) ! vector of thermal diffusion coeffs
    real :: effDiff(1:NUM_SP) ! vector of effective diffusion coeffs
    real :: DIFF(1:NUM_SP,1:NUM_SP) ! matrix of diffusion coeffs
end type

end module defs_models