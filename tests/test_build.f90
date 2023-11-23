program test_build

    use constant_air5
    use defs_models
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t_simpl
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k, n

    type(transport_in) :: transport
    ! type(SpHeatVOut) :: cv
    ! type(omega_int) :: omega_test
    ! type(bracket_int) :: bracket_test
    type(transport_out) :: transport_coeff

    x(1) = 0.77999 
    x(5) = 0.19999 
    x(3) = 0.01999
    x(4) = 0.00086999 
    x(2) = 0.00099 
   

    press = 100000


    n = 7000

    open(6, file='../res/air5_1T_test_new.txt', status='unknown')

    write (6, *) 'INPUT DATA:'
    write (6, *)

    write (6, '(A25,E13.6)') 'Pressure, Pa:       ',press
    write (6, '(A25,E13.6)') 'N2 molar fraction:  ',x(1)
    write (6, '(A25,E13.6)') 'O2 molar fraction:  ',x(2)
    write (6, '(A25,E13.6)') 'NO molar fraction:  ',x(3)
    write (6, '(A25,E13.6)') 'N molar fraction:   ',x(4)
    write (6, '(A25,E13.6)') 'O molar fraction:   ',x(5)

    write (6, *)
    write (6, *) 'TRANSPORT COEFFICIENTS:'
    write (6, *)

    write (6, '(1x, A12, A12, A12, A12)') 'Temp. ', 'Sh. visc. ', 'Th. cond. '
    write (6, *)

    do k = 500, n, 100

        T = k
        ntot = press/T/Kb
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        call Transport1TSimpl(transport, transport_coeff)

        write (6, '(1x, F12.0, E13.6, E13.6, E13.6)') transport%temp, transport_coeff%visc, &
                                                      transport_coeff%ltot


    end do

    close(6)


    open(6, file='../res/air5-1Ttest-effDiff.txt', status='unknown')

    write (6, *) 'INPUT DATA:'
    write (6, *)

    write (6, '(A25,E13.6)') 'Pressure, Pa:       ',press
    write (6, '(A25,E13.6)') 'N2 molar fraction:  ',x(1)
    write (6, '(A25,E13.6)') 'O2 molar fraction:  ',x(2)
    write (6, '(A25,E13.6)') 'NO molar fraction:  ',x(3)
    write (6, '(A25,E13.6)') 'N molar fraction:   ',x(4)
    write (6, '(A25,E13.6)') 'O molar fraction:   ',x(5)
    
    do k = 500, n, 100

        T = k
        ntot = press/T/Kb
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        call Transport1TSimpl(transport, transport_coeff)

        write (6, *)
        write (6, '(A25,E13.6)') 'Temperature:       ',T
        write (6, *)

        write (6, *) 'EFFECTIVE DIFFUSION COEFFICIENTS D_i, m^2/s'
        write (6, *)

        write (6, '(1x, E13.6)') (transport_coeff%effDiff(i), i=1,NUM_SP)

        write (6, *)
    end do

    close(6)

end program