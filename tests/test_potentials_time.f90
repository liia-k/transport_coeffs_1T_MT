program test_potentials_time

    ! This program calculates transport coefficients for a given mixture of air with varied temperature.
    ! Simplified model is applied: only shear viscosity, thermal conductivity
    ! and effective diffusion coefficients are calculated.
    ! Different interaction potentials are used: Lennard-Jones, Born-Mayer, VSS, ESA-Bruno.
    ! The time taken to calculate the transport coefficients is recorded.
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
    integer, parameter :: begin = 300, end = 9000, step = 1
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
     type(bracket_int) :: bracket_out
    type(omega_int) :: omega_in
    type(SpHeatVOut) :: cv_in
    real(8), dimension(4) :: timeArray
    real(8) :: start_time, end_time
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

    ! Loop over interaction models
    do i = 1, 4
        call cpu_time(start_time)
        do k = begin, end, step
            T = k
            ntot = press / T / Kb
            rho = sum(MASS_SPCS * ntot * x)
            y = (ntot / rho) * x * MASS_SPCS

            transport%temp = T
            transport%mass_fractions = y
            transport%rho = rho

            call Transport1TSimpl(transport, transport_coeff, interaction(i))
        end do
        call cpu_time(end_time)
        timeArray(i) = end_time - start_time
    end do

    write (6, '(4x, A15, A15, A15, A15)') (interaction(k), k=1,4)
    write (6, '(8x, E13.6, E13.6, E13.6, E13.6)') (timeArray(k), k=1,4)

    T = 1000 ! K
    ntot = press / T / Kb
    rho = sum(MASS_SPCS * ntot * x)
    y = (ntot / rho) * x * MASS_SPCS
    
    transport%temp = T
    transport%mass_fractions = y
    transport%rho = rho

    call OmegaInt(T, omega_in, interaction(1))
    call SpHeat(T, y, cv_in)
    call BracketInt(T, x, omega_in, bracket_out)

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

end program test_potentials_time