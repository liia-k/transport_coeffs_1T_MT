program test_current

    use constant_air5
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t
    ! use transport_1t_simpl
    use defs_models
    
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k, n

    type(transport_in) :: transport
    ! type(cv_out) :: cv
    ! type(omega_int) :: omega_test
    ! type(bracket_int) :: bracket_test
    type(transport_out) :: transport_coeff

    x(1) = 0.77999 
    x(2) = 0.19999 
    x(3) = 0.01999
    x(4) = 0.00086999 
    x(5) = 0.00099 
   
    ! y(1) = 0.756656E+00   
    ! y(2) = 0.221602E+00   
    ! y(3) = 0.207714E-01   
    ! y(4) = 0.421982E-03   
    ! y(5) = 0.548507E-03

    press = 100000


    n = 5000

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

    write (6, '(1x, A12, A12, A12, A12)') 'Temp. ', 'Sh. visc. ', 'B. visc. ', 'Th. cond. '
    write (6, *)

    do k = 400, n, 100

        T = k
        ntot = press/T/Kb
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        call Transport1T(transport, transport_coeff)

        write (6, '(1x, F12.0, E13.6, E13.6, E13.6)') transport%temp, transport_coeff%visc, &
                                                      transport_coeff%bulk_visc, transport_coeff%ltot

        ! write (6, *) 'Temperature, K         ',transport%temp
        ! write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff%visc
        ! write (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa.s              ', transport_coeff%bulk_visc
        ! write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        ! write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, tr , W/m/K        ', ltr
        ! write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, int , W/m/K       ', lint
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_N2, W/m/K     ', lvibr_n2
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_O2, W/m/K     ', lvibr_o2
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_NO, W/m/K     ', lvibr_no
        ! write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N2, m^2/s          ', transport_coeff%THDIFF(1)
        ! write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O2, m^2/s          ', transport_coeff%THDIFF(2)
        ! write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of NO, m^2/s          ', transport_coeff%THDIFF(3)
        ! write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N, m^2/s           ', transport_coeff%THDIFF(4)
        ! write (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O, m^2/s           ', transport_coeff%THDIFF(5)

        ! write (6, *)
        ! write (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        ! write (6, *)


        ! do i=1,5
        !     write (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i,j), j=1,5)
        ! end do

        ! write (6, *)

    end do

    close(6)

end program