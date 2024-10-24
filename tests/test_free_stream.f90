program test_free_stream

! Test program for transport properties calculation bsed on 1T model (general procedure)
    use constant_air5
    use defs_models
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t_new
    
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k, n
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    character(len=*), parameter :: interaction = 'Lennard-Jones' ! 'VSS', 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'
    ! type(SpHeatVOut) :: cv
    ! type(omega_int) :: omega_test
    ! type(bracket_int) :: bracket_test


    ! Critical values for test, free-stream conditions
    y(1) = 0.768
    y(2) = 0.232
    y(3) = 1.e-5
    y(4) = 1.e-5
    y(5) = 1.e-5

    rho = 3.09674 * 1.e-3
    T = 247
    n = 5

    ! Open file for results
    open(6, file='../res/air5-1T-free-stream.txt', status='unknown')

    do k = 1, n
        M = 1 / dot_product(y, 1 / Molar)
        ntot = rho * R / Kb / M

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        call Transport1TGeneral(transport, transport_coeff, interaction)

        ! Write input data to file
        write (6, *) 'INPUT DATA:'
        write (6, *)
        write (6, *) 'Temperature, K         ', transport%temp
        write (6, *) 'Molar mass, kg         ', M
        write (6, *) 'Density, kg/m^3        ', transport%rho
        write (6, *) 'N2 mass fraction       ', transport%mass_fractions(1)
        write (6, *) 'O2 mass fraction       ', transport%mass_fractions(2)
        write (6, *) 'NO mass fraction       ', transport%mass_fractions(3)
        write (6, *) 'N mass fraction        ', transport%mass_fractions(4)
        write (6, *) 'O mass fraction        ', transport%mass_fractions(5)
        write (6, *)
        write (6, *) 'TRANSPORT COEFFICIENTS:'
        write (6, *)
        write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa*S             ', transport_coeff%visc
        write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        write (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa*S              ', transport_coeff%bulk_visc
        write (6, *)
        write (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        write (6, *)

        do i = 1, NUM_SP
            write (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i, j), j = 1, NUM_SP)
        end do

        write (6, *)
        write (6, *) 'THERMAL DIFFUSION COEFFICIENTS, m^2/s'
        write (6, *)
        write (6, '(1x, 5E15.6)') (transport_coeff%thdiff(i), i = 1, NUM_SP)
        write (6, *)

        ! Update mass fractions for next iteration
        y(3) = y(3) / 10
        y(4) = y(3)
        y(5) = y(3)
    end do

    close(6)

end program test_free_stream