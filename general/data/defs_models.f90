!Macroparameters required for transport coeffs. calculation:
module defs_models

use constant_air5

type transport_in
    real rho  
    !real vel(1:3), not required
    real temp
    real mass_fractions(1:NUM_SP)
end type

!Transport coefficients: 

type transport_out
    real visc, bulk_visc, ltot ! shear viscosity,  bulk viscosity, thermal conductivity
    real :: thdiff(1:NUM_SP) ! thermal diffusion coeffs
    real :: diff(1:NUM_SP,1:NUM_SP) ! diffusion coeffs
end type

end module defs_models