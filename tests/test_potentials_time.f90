program test_potentials_time

    ! This program calculates transport coefficients for a given mixture of air with varied temperature.
    ! Simplified model is applied: only shear viscosity, thermal conductivity
    ! and effective diffusion coefficients are calculated.
    ! Different interaction potentials are used: Lennard-Jones, Born-Mayer, VSS, ESA-Bruno.
    ! The time taken to calculate the transport coefficients is recorded.
    use constant_air5
    use defs_models
    use transport_1t
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k
    integer, parameter :: begin = 300, end = 9000, step = 1
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
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

end program test_potentials_time