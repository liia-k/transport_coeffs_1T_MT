program test_n2n

    use constant_n2n
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t
    ! use transport_1t_simpl
    use defs_models
    
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y

    integer :: i, j, k, n

    type(transport_in) :: transport
    ! type(cv_out) :: cv
    ! type(omega_int) :: omega_test
    ! type(bracket_int) :: bracket_test
    type(transport_out) :: transport_coeff

    y(1) = 0.72
    y(2) = 0.28
    transport%mass_fractions = y

    press = 100000


    n = 5000

    open(6, file='../res/air5_1T_test_n2n.txt', status='unknown')

    write (6, *) 'INPUT DATA, N2, N:'
    write (6, *)

    write (6, '(A25,E13.6)') 'Pressure, Pa:       ',press


    write (6, *)
    write (6, *) 'TRANSPORT COEFFICIENTS:'
    write (6, *)

    write (6, '(1x, A12, A12, A12, A12)') 'Temp. ', 'Sh. visc. ', 'B. visc. ', 'Th. cond. '
    write (6, *)

    do k = 400, n, 100

        T = k
        M = 1/dot_product(y,1/MOLAR)
        rho = M*press/R/T
        
        transport%temp = T
        transport%rho = rho

        call Transport1T(transport, transport_coeff)

        write (6, '(1x, F12.0, E13.6, E13.6, E13.6)') transport%temp, transport_coeff%visc, &
                                                      transport_coeff%bulk_visc, transport_coeff%ltot

 

    end do

    close(6)

end program