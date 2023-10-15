program test_free_stream

    use constant_air5
    use specific_heat_sp
    use omega_integrals
    use bracket_integrals
    use transport_1t
    use defs_models
    
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k, n

    type(transport_in) :: transport
    type(cv_out) :: cv
    type(omega_int) :: omega_test
    type(bracket_int) :: bracket_test
    type(transport_out) :: transport_coeff

    ! critical values for test, free-stream conditions
    y(1) = 0.768
    y(2) = 0.232
    y(3) = 1.e-5
    y(4) = 1.e-5
    y(5) = 1.e-5

    rho = 3.09674*1.e-3
    T = 247


    n = 20

    do k = 1, n

        M = 1/dot_product(y, 1/Molar)
        !ntot = rho*R/kb/M

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        ntot = rho*R/Kb/M
        
        !CALL S_Heat(transport[k], cv)
        !CALL Omega(transport%temp, omega_test)
        !CALL Bracket(transport%temp, x, omega_test, bracket_test)
        call Transport1T(transport, transport_coeff)

        open(6, file='tests/air5_1T_test_free-stream.txt', status='unknown')

        write (6, *) 'INPUT DATA:'

        write (6, *)


        write (6, *) 'Temperature, K         ',transport%temp
        write (6, *) 'Molar mass, kg         ',M
        write (6, *) 'Density, kg/m^3        ',transport%rho
        write (6, *) 'Number density, 1/m^3  ',ntot
        write (6, *) 'N2 mass fraction       ',transport%mass_fractions(1)
        write (6, *) 'O2 massr fraction      ',transport%mass_fractions(2)
        write (6, *) 'NO mass fraction       ',transport%mass_fractions(3)
        write (6, *) 'N mass fraction        ',transport%mass_fractions(4)
        write (6, *) 'O mass fraction        ',transport%mass_fractions(5)

        ! write (6, *)


        ! write (6, *) 'Total internal Specific heat  ',cv%cv_int
        ! write (6, *) 'Total  Specific heat          ',cv%cv_tot

        ! do i1=1,5
        !     write (6, *) cv%c_int(i1)
        ! end do

        ! write (6, *)

        ! write (6, *) 'Omega Integrals Omega11_ij'

        ! write (6, *)

        ! do i=1,5
        !     write (6, '(1x, 5e15.6)') (omega_test%omega11(i,j), j=1,5)
        !     write (6, '(1x, 5E15.6)') (omega_test%omega12(i,j), j=1,5)
        ! end do

        ! write (6, *) 'Bracket Integrals Omega11_ij'

        ! write (6, *)

        ! do i=1,5
        !     write (6, '(1x, 5e15.6)') (bracket_test%lambda(i,j), j=1,5)
        !     write (6, '(1x, 5E15.6)') (bracket_test%beta01(i,j), j=1,5)
        ! end do

        !write (6, '(1x, 5E15.6)') (y(j), j=1,5)


        write (6, *) 'TRANSPORT COEFFICIENTS:'
        write (6, *)

        write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff%visc
        write (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa.s              ', transport_coeff%bulk_visc
        write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N2, m^2/s          ', transport_coeff%THDIFF(1)
        write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O2, m^2/s          ', transport_coeff%THDIFF(2)
        write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of NO, m^2/s          ', transport_coeff%THDIFF(3)
        write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N, m^2/s           ', transport_coeff%THDIFF(4)
        write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O, m^2/s           ', transport_coeff%THDIFF(5)

        write (6, *)
        write (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        write (6, *)


        do i=1,5
            write (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i,j), j=1,5)
        end do

        write (6, *)

        y(3) = y(3)/10
        y(4) = y(4)/10
        y(5) = y(5)/10

    end do

    close(6)

end program