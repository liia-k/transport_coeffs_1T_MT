program test_potential_models

    ! This program calculates transport coefficients for a given mixture of air with varied temperature.
    ! Simplified model is applied: only shear viscosity, thermal conductivity
    ! and effective diffusion coefficients are calculated.
    ! Different interaction potentials are used: Lennard-Jones, Born-Mayer, VSS, ESA-Bruno.
    use constant_air5
    use defs_models
    use specific_heat_sp
    use omega_integrals
    use bracket_integrals
    use transport_1t
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k
    integer, parameter :: start = 500, end = 9000, step = 100
    integer, parameter :: size = (end - start) / step + 1
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    type(bracket_int) :: bracket_out
    type(omega_int) :: omega_in
    type(SpHeatVOut) :: cv_in
    real, dimension(4, size) :: shearVisc, thermCond
    real, dimension(size) :: tempArray
    character(len=25), dimension(4) :: interaction 

    ! Initialize interaction models
    interaction(1) = 'Lennard-Jones'
    interaction(2) = 'Born-Mayer'
    interaction(3) = 'VSS'
    interaction(4) = 'ESA-Bruno'

    ! Initialize molar fractions
    x(1) = 0.77999 
    x(2) = 0.19999 
    x(3) = 0.01999
    x(4) = 0.00086999 
    x(5) = 0.00099 

    press = 100000

    j = 0
    do k = start, end, step
        j = j + 1
        T = k
        ntot = press / T / Kb
        rho = sum(MASS_SPCS * ntot * x)
        y = (ntot / rho) * x * MASS_SPCS
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        tempArray(j) = T
        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        do i=1,4
            call Transport1TSimpl(transport, transport_coeff, interaction(i))
            shearVisc(i,j) = transport_coeff%visc
            thermCond(i,j) = transport_coeff%ltot
        end do

    end do

    ! Write shear viscosity coefficients to file
    open(6, file='../res/air5-1T-shear-visc.txt', status='unknown')
    write (6, *)
    write (6, '(1x, A40)') 'Shear viscosity coefficients:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)
    do k = 1, size
        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), shearVisc(:, k)
    end do
    close(6)

    ! Write thermal conductivity coefficients to file
    open(6, file='../res/air5-1T-therm-cond.txt', status='unknown')
    write (6, *)
    write (6, '(1x, A40)') 'Thermal conductivity coefficient:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)
    do k = 1, size
        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), thermCond(:, k)
    end do
    close(6)


    T = 1000 ! K
    ntot = press / T / Kb
    rho = sum(MASS_SPCS * ntot * x)
    y = (ntot / rho) * x * MASS_SPCS
    
    tempArray(j) = T
    transport%temp = T
    transport%mass_fractions = y
    transport%rho = rho

    call OmegaInt(T, omega_in, interaction(1))
    call SpHeat(T, y, cv_in)
    call BracketInt(T, ntot, x, omega_in, cv_in, bracket_out)

    open(6, file='../res/air5-1T-bracket-integrals.txt', status='unknown')
    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Lambda1100 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%lambda1100(i, j), j = 1, NUM_SP)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Lambda0000 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%lambda0000(i, j), j = 1, NUM_SP)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Lambda0100 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%lambda0100(i, j), j = 1, NUM_SP)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals H00 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%h00(i, j), j = 1, NUM_SP)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Beta1100 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%beta1100(i, j), j = 1, NUM_SP)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Beta0110 at T = 1000 K:'
    write (6, *)

    do i = 1, NUM_SP
        write (6, '(1x, 5E15.6)') (bracket_out%beta0110(i, j), j = 1, NUM_MOL)
    end do

    write (6, *)
    write (6, '(1x, A40)') 'Bracket integrals Beta0011 at T = 1000 K:'
    write (6, *)

    write (6, '(1x, 5E15.6)') (bracket_out%beta0011(i), i = 1, NUM_MOL)

    close(6)
end program test_potential_models