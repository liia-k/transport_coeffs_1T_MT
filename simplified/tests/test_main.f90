program test_main
    
    use defs_models
    use constant_air5
    use transport_1t_simpl
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y

    integer :: i, j, k, n

    type(transport_in) :: transport
    type(transport_out) :: transport_coeff


    n = 10000

    open(6, file='../res/air5_1T_test_main.txt', status='unknown')

    do k = 1, n
        
        call random_number(y(2:))
        y(2:3) = y(2:3)/1e2
        y(1) = 1 - y(2) - y(3) - y(4) - y(5)
        if (y(1) < 0) then
            cycle
        end if
        call random_number(press)
        press = press*1e6
        call random_number(T)
        T = T*1e4

        M = 1/dot_product(y,1/MOLAR)
        rho = M*press/T/R

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        write (6, *) 'INPUT DATA:'
        write (6, *)

        write (6, '(A25,E13.6)') 'Pressure, Pa:      ',press
        write (6, '(A25,E13.6)') 'Temperature, K:    ',T
        write (6, '(A25,E13.6)') 'N2 mass fraction:  ',y(1)
        write (6, '(A25,E13.6)') 'O2 mass fraction:  ',y(2)
        write (6, '(A25,E13.6)') 'NO mass fraction:  ',y(3)
        write (6, '(A25,E13.6)') 'N mass fraction:   ',y(4)
        write (6, '(A25,E13.6)') 'O mass fraction:   ',y(5)

        call Transport1TSimpl(transport, transport_coeff)

        write (6, *) 'TRANSPORT COEFFICIENTS:'
        write (6, *)

        write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff%visc
        write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        
        write (6, *)
        write (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        write (6, *)


        do i=1,NUM_SP
            write (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i,j), j=1,NUM_SP)
        end do

        write (6, *)

    end do

    close(6)

end program