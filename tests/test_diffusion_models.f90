program test_diffusion_models

    ! This program calculates transport coefficients for a specific mixture of air.
    ! Simplified model is applied: only shear viscosity, thermal conductivity 
    ! and effective diffusion coefficients are calculated
    use constant_air5
    use defs_models
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k, n
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    character(len=*), parameter :: interaction = 'ESA-Bruno' ! 'VSS', 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'

    ! Number of iterations
    n = 1000

    ! Open file for results
    open(10, file='../res/air5-1T-diffusion-models.txt', status='unknown')

    write (6, *)
    write (6, *) 'Potential model: ', interaction
    write (6, *)

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
        write (10, *) 'Input data:'
        write (10, *)
        write (10, '(A25,E13.6)') 'Pressure, Pa:      ', press
        write (10, '(A25,E13.6)') 'Temperature, K:    ', T
        write (10, '(A25,E13.6)') 'N2 mass fraction:  ', y(1)
        write (10, '(A25,E13.6)') 'O2 mass fraction:  ', y(2)
        write (10, '(A25,E13.6)') 'NO mass fraction:  ', y(3)
        write (10, '(A25,E13.6)') 'N mass fraction:   ', y(4)
        write (10, '(A25,E13.6)') 'O mass fraction:   ', y(5)

        ! Calculate transport coefficients
        call Transport1TGeneral(transport, transport_coeff, interaction)

        write (10, *) 'Diffusion coefficients D_ij, m^2/s'
        write (10, *)
        do i = 1, NUM_SP
            write (10, '(1x, 5E15.6)') (transport_coeff%DIFF(i, j), j = 1, NUM_SP)
        end do
        write (10, *)

        call Transport1TSimpl(transport, transport_coeff, interaction)
        
        write (10, *) 'Effective diffusion coefficients, m^2/s'
        write (10, *)
        write (10, '(1x, 5E15.6)') (transport_coeff%effDiff(i), i = 1, NUM_SP)   
        write (10, *)
    end do

    close(10)

end program test_diffusion_models