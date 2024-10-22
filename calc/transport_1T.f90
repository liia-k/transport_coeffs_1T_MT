!Module for calculation of transport coefficients. Uses the module constant_...(name of mixture).f90 containing main constants and 
!variables definition, modules for calculation of specific heats, omega and bracket integrals.

module transport_1t

    ! Use statements for required modules
    use defs_models
    use constant_air5
    use specific_heat_sp
    use omega_integrals
    use bracket_integrals
    use qr_decomposition
    use, intrinsic :: ieee_arithmetic 
    
    implicit none
    ! Define constants
    real, parameter :: SMALL_VALUE = 1e-5

    contains

    subroutine Transport1TGeneral(data_in, data_out, interactionType)
        implicit none
        type(transport_in), intent(in) :: data_in
        character(len=*), intent(in), optional :: interactionType
        type(transport_out), intent(out) :: data_out

        ! Local variables
        type(transport_in_additional) :: data_add
        real :: ltot, visc, bulk_visc
        real, dimension(NUM_SP) :: effDiff, thdiff
        real, dimension(NUM_SP, NUM_SP) :: DIFF
        type(SpHeatVOut) :: cv
        type(omega_int) :: omega_out
        type(bracket_int) :: bracket_out

        ! Initial calculations
        call initialize(data_in, data_add)

        ! Calculation of specific heats and omega + bracket integrals
        call SpHeat(data_in%temp, data_in%mass_fractions, cv)
        if (present(interactionType)) then
            call OmegaInt(data_in%temp, omega_out, interactionType)
        else
            call OmegaInt(data_in%temp, omega_out)
        end if
        call BracketInt(data_in%temp, data_add%num_fractions, omega_out, bracket_out)

        ! Calculations of transport coefficients
        call calculateThermalCondAndDiff(data_add, cv, bracket_out, ltot, thdiff)
        call calculateDiffCoeffs(data_add, bracket_out, DIFF)
        call calculateShearVisc(data_add, bracket_out, visc)
        call calculateBulkVisc(data_add, cv, bracket_out, bulk_visc)

        effDiff = 0
        ! Output results
        call outputResults(data_in, data_out, ltot, DIFF, visc, bulk_visc, effDiff)

    end subroutine Transport1TGeneral

    subroutine Transport1TSimpl(data_in, data_out, interactionType)
        implicit none
        type(transport_in), intent(in) :: data_in
        character(len=*), intent(in), optional :: interactionType
        type(transport_out), intent(out) :: data_out

        ! Local variables
        type(transport_in_additional) :: data_add
        real :: ltot, visc, bulk_visc
        real, dimension(NUM_SP) :: effDiff, thdiff
        real, dimension(NUM_SP, NUM_SP) :: DIFF
        type(SpHeatVOut) :: cv
        type(omega_int) :: omega_out
        type(bracket_int) :: bracket_out

        ! Initial calculations
        call initialize(data_in, data_add)

        ! Calculation of specific heats and omega + bracket integrals
        call SpHeat(data_in%temp, data_in%mass_fractions, cv)
        if (present(interactionType)) then
            call OmegaInt(data_in%temp, omega_out, interactionType)
        else
            call OmegaInt(data_in%temp, omega_out)
        end if
        call BracketInt(data_in%temp, data_add%num_fractions, omega_out, bracket_out)

        ! Calculations of transport coefficients
        call calculateThermalCond(data_add, cv, bracket_out, ltot)
        call calculateEffDiffCoeffs(data_add, omega_out, effDiff)
        ! call calculateDiffCoeffs(data_add, bracket_out, DIFF)
        call calculateShearVisc(data_add, bracket_out, visc)

        thdiff = 0
        DIFF = 0
        bulk_visc = 0

        ! Output results
        call outputResults(data_in, data_out, ltot, DIFF, visc, bulk_visc, effDiff)

    end subroutine Transport1TSimpl

    subroutine initialize(data_in, data_add)
        implicit none
        type(transport_in), intent(in) :: data_in
        type(transport_in_additional), intent(out) :: data_add
        real, dimension(NUM_SP) :: y, x
        real :: ntot, M
        integer :: i, max_index

        y = data_in%mass_fractions

        ! Adjust small mass fractions
        max_index = max(maxloc(y, dim=1), 0)
        do i = 1, NUM_SP
            if (abs(y(i)) < SMALL_VALUE) then
                y(max_index) = y(max_index) + y(i) - SMALL_VALUE
                y(i) = SMALL_VALUE
            end if
        end do

        M = 1 / dot_product(y, 1 / MOLAR)
        ntot = sum(data_in%rho * y / MASS_SPCS)
        x = (data_in%rho / ntot) * y / MASS_SPCS

        ! Set additional macroparameters
        data_add%mass_fractions = y
        data_add%rho = data_in%rho
        data_add%temp = data_in%temp
        data_add%ntot = ntot
        data_add%M = M
        data_add%num_fractions = x
    end subroutine initialize

    subroutine outputResults(data_in, data_out, ltot, DIFF, visc, bulk_visc, effDiff)
        implicit none
        type(transport_in), intent(in) :: data_in
        type(transport_out), intent(out) :: data_out
        real, intent(in) :: ltot, visc, bulk_visc
        real, dimension(NUM_SP, NUM_SP), intent(in) :: DIFF
        real, dimension(NUM_SP), intent(in) :: effDiff
        integer :: i, j

        ! Output thermal conductivity
        if (exception(ltot)) then
            data_out%ltot = ltot
        else
            print *, "thermal conductivity is not calculated for the set:"
            call macro_output(data_in)
            data_out%ltot = 0
        end if

        ! Output diffusion coefficients
        do i = 1, NUM_SP
            do j = 1, NUM_SP
                if (exception(DIFF(i, j))) then
                    data_out%DIFF(i, j) = DIFF(i, j)
                else
                    print *, "diffusion coeffs are not calculated for the set:"
                    call macro_output(data_in)
                    data_out%DIFF(i, j) = 0
                end if
            end do
        end do

        ! Output shear viscosity
        if (exception(visc)) then
            data_out%visc = visc
        else
            print *, "shear viscosity is not calculated for the set:"
            call macro_output(data_in)
            data_out%visc = 0
        end if

        ! Output bulk viscosity
        if (exception(bulk_visc)) then
            data_out%bulk_visc = bulk_visc
        else
            print *, "bulk viscosity is not calculated for the set:"
            call macro_output(data_in)
            data_out%bulk_visc = 0
        end if

        ! Output effective diffusion coefficients
        do i = 1, NUM_SP
            if (exception(effDiff(i))) then
                data_out%effDiff(i) = effDiff(i)
            else
                print *, "effective diffusion coeffs are not calculated for the set:"
                call macro_output(data_in)
                data_out%effDiff(i) = 0
            end if
        end do

    end subroutine outputResults

    logical function exception(var)
    ! has a true value if var in not NaN and Infinity
        real, intent(in) :: var
        exception = .not.(ieee_is_nan(var)) .and. (ieee_is_finite(var))
    end function

    subroutine macro_output(data_in)
    ! macroparameters given values output
        type(transport_in),intent(in)   :: data_in
        print *, "temperature = ", data_in%temp
        print *, "density = ", data_in%rho
        print *, "mass fractions: ", data_in%mass_fractions
    end subroutine


    subroutine calculateThermalCondAndDiff(data_in, cv, bracket_ints, ltot, thdiff)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(SpHeatVOut), intent(in) :: cv
        type(bracket_int), intent(in) :: bracket_ints

        real, intent(out) :: ltot
        real, dimension(NUM_SP), intent(out) :: thdiff

        ! Translational and internal heat conductivity
        real :: ltr, lint

       ! Matrices for the linear transport system defining heat conductivity and thermal diffusion (LTH)
        real, dimension(2*NUM_SP, 2*NUM_SP) :: LTH, Q_LTH, R_LTH, Inverse_LTH

        ! Vector of the RHS of the system
        real, dimension(2*NUM_SP, 1) :: b
        
        ! Define matrix LTH for calculation of thermal conductivity and thermal diffusion coefficients      
        LTH(1:NUM_SP, 1:NUM_SP) = bracket_ints%Lambda00  
        LTH(1:NUM_SP, NUM_SP+1 : 2*NUM_SP) = bracket_ints%Lambda01
        LTH(NUM_SP+1 : 2*NUM_SP, 1:NUM_SP) = LTH(1 : NUM_SP, NUM_SP+1 : 2*NUM_SP)
        LTH(NUM_SP+1 : 2*NUM_SP, NUM_SP+1 : 2*NUM_SP) = bracket_ints%Lambda11    

        LTH(1, 1 : NUM_SP) = data_in%mass_fractions    
        LTH(1, NUM_SP+1 : 2*NUM_SP) = 0

       ! Define vector b (right-hand side of the system for calculation of thermal conductivity and thermal diffusion coefficients)
        b(: NUM_SP, 1) = 0
        b(NUM_SP+1 : 2*NUM_SP, 1) = (4./5./Kb) * data_in%num_fractions
        
        ! Linear system solution using QR decomposition 
        call QRDecomposition(LTH, Q_LTH, R_LTH, 2*NUM_SP)
        call InvertQR(Q_LTH, R_LTH, Inverse_LTH, 2*NUM_SP)
        b = matmul(Inverse_LTH, b)

        ! Calculate thermal diffusion coefficients (thdiff)
        thdiff = -(1./2./data_in%ntot) * b(1:NUM_SP, 1)

        ! Calculate thermal conductivity coefficient associated with translational energy (ltr)
        ltr = (5./4.) * Kb * sum(data_in%num_fractions * b(NUM_SP+1:2*NUM_SP, 1))

        ! Calculate thermal conductivity coefficients associated with internal energies (lint)
        lint = (3./16.) * data_in%temp * sum(data_in%num_fractions(:NUM_MOL) * cv%cv_int_sp(:NUM_MOL) &
                * ((Kb) * (MASS_SPCS(:NUM_MOL)*1e30) / (bracket_ints%lambda_int(:NUM_MOL)*1e30)))

        ! Total thermal conductivity coefficient at the translational temperature gradient
        ltot = ltr + lint

    end subroutine calculateThermalCondAndDiff

    subroutine calculateThermalCond(data_in, cv, bracket_ints, ltot)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(SpHeatVOut), intent(in) :: cv
        type(bracket_int), intent(in) :: bracket_ints

        real, intent(out) :: ltot

        ! Translational and internal heat conductivity
        real :: ltr, lint

       ! Matrices for the linear transport system defining heat conductivity and thermal diffusion (LTH)
        real, dimension(NUM_SP, NUM_SP) :: LTH, Q_LTH, R_LTH, Inverse_LTH

        ! Vector of the RHS of the system
        real, dimension(NUM_SP, 1) :: b
        
        ! Define matrix LTH for calculation of thermal conductivity coefficients      
        LTH(:NUM_SP, :NUM_SP) = bracket_ints%Lambda11  

       ! Define vector b (right-hand side of the system for calculation of thermal conductivity coefficients)
        b(: NUM_SP, 1) = (4./5./Kb) * data_in%num_fractions
        
        ! Linear system solution using QR decomposition 
        call QRDecomposition(LTH, Q_LTH, R_LTH, NUM_SP)
        call InvertQR(Q_LTH, R_LTH, Inverse_LTH, NUM_SP)
        b = matmul(Inverse_LTH, b)

        ! Calculate thermal conductivity coefficient associated with translational energy (ltr)
        ltr = (5./4.) * Kb * sum(data_in%num_fractions * b(: NUM_SP, 1))


        ! Calculate thermal conductivity coefficients associated with internal energies (lint)
        lint = (3./16.) * data_in%temp * sum(data_in%num_fractions(:NUM_MOL) * cv%cv_int_sp(:NUM_MOL) &
                * ((Kb) * (MASS_SPCS(:NUM_MOL)*1e30) / (bracket_ints%lambda_int(:NUM_MOL)*1e30)))

        ! Total thermal conductivity coefficient at the translational temperature gradient
        ltot = ltr + lint

    end subroutine calculateThermalCond

    subroutine calculateDiffCoeffs(data_in, bracket_ints, DIFF)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(bracket_int), intent(in) :: bracket_ints
        real, dimension(NUM_SP, NUM_SP), intent(out) :: DIFF

        integer i, j, delta

        ! Matrices for the linear transport systems defining diffusion (LDIFF)
        real, dimension(NUM_SP, NUM_SP) :: LDIFF, Q_LDIFF, R_LDIFF, Inverse_LDIFF

        ! Vectors of the RHS of systems
        real, dimension(NUM_SP, NUM_SP) :: B1

        ! Define matrix LDIFF for calculation of diffusion coefficients
        LDIFF = bracket_ints%Lambda11

        ! Define matrix B1 (right-hand side of the system for calculation of diffusion coefficients)
        do i = 1, NUM_SP
            do j = 1, NUM_SP
                if (i == j) then
                    delta = 1
                else
                    delta = 0
                end if
                B1(i, j) = 8./25./Kb * (delta - data_in%mass_fractions(i))
            end do
        end do

        B1(1, 1:NUM_SP) = 0

        ! Linear system solution using QR decomposition
        call QRDecomposition(LDIFF, Q_LDIFF, R_LDIFF, NUM_SP)
        call InvertQR(Q_LDIFF, R_LDIFF, Inverse_LDIFF, NUM_SP)
        B1 = matmul(Inverse_LDIFF, B1)

        ! Calculate diffusion coefficients, DIFF(i, j)
        DIFF = (1./2./data_in%ntot) * B1

    end subroutine calculateDiffCoeffs

    subroutine calculateEffDiffCoeffs(data_in, omega_ints, effDiff)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(omega_int), intent(in) :: omega_ints
        real, dimension(NUM_SP), intent(out) :: effDiff

        integer :: i, j
        real :: mij
        real, dimension(NUM_SP, NUM_SP) :: binDiff

        ! Calculate binary diffusion coefficients
        do i = 1, NUM_SP
            do j = 1, NUM_SP
                mij = MASS_SPCS(j) * (MASS_SPCS(i) * 1E27) / (MASS_SPCS(i) * 1E27 + MASS_SPCS(j) * 1E27)
                binDiff(i, j) = 3.0 * Kb * data_in%temp / 16.0 / (data_in%ntot * mij) &
                                / omega_ints%OMEGA11(i, j)
            end do
        end do

        ! Calculate effective diffusion coefficients
        do i = 1, NUM_SP
            effDiff(i) = (1.0 - data_in%num_fractions(i)) &
                         / sum(data_in%num_fractions(1:NUM_MOL) / binDiff(i, 1:NUM_MOL))
        end do

    end subroutine calculateEffDiffCoeffs

    subroutine calculateShearVisc(data_in, bracket_ints, visc)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(bracket_int), intent(in) :: bracket_ints
        real, intent(out) :: visc

        ! Matrices for the linear transport systems defining shear viscosity (HVISC)
        real, dimension(NUM_SP, NUM_SP) :: HVISC, Q_HVISC, R_HVISC, Inverse_HVISC

        ! Vectors of the RHS of systems
        real, dimension(NUM_SP, 1) :: b2

        ! Define matrix HVISC for calculation of shear viscosity coefficient
        HVISC = bracket_ints%H00

        ! Define vector b2 (right-hand side of the system for calculation of shear viscosity coefficient)
        b2(1:NUM_SP, 1) = (2.0 / Kb /  data_in%temp) * data_in%num_fractions

        ! Linear system solution using QR decomposition
        call QRDecomposition(HVISC, Q_HVISC, R_HVISC, NUM_SP)
        call InvertQR(Q_HVISC, R_HVISC, Inverse_HVISC, NUM_SP)
        b2 = matmul(Inverse_HVISC, b2)

        ! Calculate shear viscosity coefficient, visc
        visc = (Kb *  data_in%temp / 2.0) * sum(data_in%num_fractions * b2(1:NUM_SP, 1))

    end subroutine calculateShearVisc

    subroutine calculateBulkVisc(data_in, cv, bracket_ints, bulk_visc)
        implicit none

        type(transport_in_additional), intent(in) :: data_in
        type(SpHeatVOut), intent(in) :: cv
        type(bracket_int), intent(in) :: bracket_ints
        real, intent(out) :: bulk_visc

        integer :: i

        ! Matrices for the linear transport systems defining bulk viscosity (BVISC)
        real, dimension(NUM_SP+NUM_MOL, NUM_SP+NUM_MOL) :: BVISC, Q_BVISC, R_BVISC, Inverse_BVISC

        ! Vectors of the RHS of systems
        real, dimension(NUM_SP+NUM_MOL, 1) :: b3

        ! Define matrix BVISC for calculation of bulk viscosity coefficients
        do i = 1, NUM_SP
            BVISC(i, 1:NUM_SP) = data_in%num_fractions * bracket_ints%beta11(i, 1:NUM_SP)
        end do
        
        do i = 1, NUM_SP
            BVISC(i, NUM_SP+1:NUM_SP+NUM_MOL) = data_in%mass_fractions(:NUM_MOL) * bracket_ints%beta01(i, :NUM_MOL)
        end do
    
        BVISC(NUM_SP+1 : NUM_SP+NUM_MOL, 1:NUM_SP) = transpose(BVISC(1:NUM_SP, NUM_SP+1:NUM_SP+NUM_MOL))  
        BVISC(NUM_SP+1 : NUM_SP+NUM_MOL, NUM_SP+1 : NUM_SP+NUM_MOL) = 0
        
        do i = NUM_SP+1, NUM_SP+NUM_MOL
            BVISC(i, i) = data_in%mass_fractions(i-5) * bracket_ints%beta0011(i-5)
        end do
        
        BVISC(1, 1:NUM_SP) = ((3.0 / 2.0) * R / data_in%M) * data_in%num_fractions ! Kb*ntot/rho    
        BVISC(1, NUM_SP+1:NUM_SP+NUM_MOL) = data_in%mass_fractions(1:NUM_MOL) * (Kb / MASS_SPCS(1:NUM_MOL))
        
        ! Define vector b3 (right-hand side of the system for calculation of bulk viscosity coefficients)
        b3(1:NUM_SP, 1) = -data_in%num_fractions * cv%cv_int_sp / cv%cv_tot
        b3(NUM_SP+1:NUM_SP+NUM_MOL, 1) = data_in%mass_fractions(:NUM_MOL) * cv%cv_int_sp(:NUM_MOL) &
                                        / cv%cv_tot

        b3(1, 1) = 0.0

    
        ! Linear system solution using QR decomposition
        call QRDecomposition(BVISC, Q_BVISC, R_BVISC, NUM_SP+NUM_MOL)
        call InvertQR(Q_BVISC, R_BVISC, Inverse_BVISC, NUM_SP+NUM_MOL)
        b3 = matmul(Inverse_BVISC, b3)
        
        ! call gaussj(BVISC,8,8,b3,1,1)

        ! Calculate bulk viscosity coefficient, bulk_visc
        bulk_visc = -Kb * data_in%temp * sum(data_in%num_fractions * b3(1:NUM_SP, 1))

    end subroutine calculateBulkVisc


end module transport_1t
