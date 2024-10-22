program test_simplified_random

    ! This program calculates transport coefficients for a random mixture of air.
    ! Simplified model is applied: only shear viscosity, thermal conductivity 
    ! and effective diffusion coefficients are calculated
    use defs_models
    use constant_air5
    use transport_1t
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y
    integer :: i, k, n
    real :: t1, t2, total_time
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    character(len=*), parameter :: interaction = 'ESA-Bruno' ! 'VSS', 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'

! Number of iterations
    n = 10000

    ! Start timing
    call cpu_time(t1)
    print *, "Start time: ", t1

    ! Open file for results
    open(10, file='../res/air5-1T-random.txt', status='unknown')

    do k = 1, n
        ! Generate random mass fractions
        call random_number(y(2:))
        y(3:4) = y(3:4) / 1e2
        y(1) = 1 - sum(y(2:5))
        if (y(1) < 0) cycle

        ! Generate random pressure and temperature
        call random_number(press)
        press = press * 1e5
        call random_number(T)
        T = T * 1e4

        ! Calculate molar mass and density
        M = 1 / dot_product(y, 1 / MOLAR)
        rho = M * press / T / R

        ! Set macroparameters for transport coefficients calculation
        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        ! Write input data to file
        write (10, *) 'INPUT DATA:'
        write (10, *)
        write (10, '(A25,E13.6)') 'Pressure, Pa:      ', press
        write (10, '(A25,E13.6)') 'Temperature, K:    ', T
        write (10, '(A25,E13.6)') 'N2 mass fraction:  ', y(1)
        write (10, '(A25,E13.6)') 'O2 mass fraction:  ', y(2)
        write (10, '(A25,E13.6)') 'NO mass fraction:  ', y(3)
        write (10, '(A25,E13.6)') 'N mass fraction:   ', y(4)
        write (10, '(A25,E13.6)') 'O mass fraction:   ', y(5)

        ! Calculate transport coefficients
        call Transport1TSimpl(transport, transport_coeff, interaction)

        ! Write transport coefficients to file
        write (10, *) 'TRANSPORT COEFFICIENTS:'
        write (10, *)
        write (10, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa*S             ', transport_coeff%visc
        write (10, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        write (10, *)
        write (10, *) 'EFFECTIVE DIFFUSION COEFFICIENTS D_ij, m^2/s'
        write (10, *)
        write (10, '(1x, 5E15.6)') (transport_coeff%effDiff(i), i = 1, NUM_SP)
        write (10, *)
    end do

    close(10)
    
    ! End timing
    call cpu_time(t2)
    write(6, *) "End time: ", t2

    total_time = t2 - t1
    write(6, '(a,f10.3,a)') "Total execution time: ", total_time, " seconds"
    ! print '(a,f10.3,a)', "Average time per iteration: ", total_time / n, " seconds"

end program test_simplified_random