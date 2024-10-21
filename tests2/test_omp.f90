program test_omp

! 6 flows, 1e6 iterations at about 18 min...

    use constant_air5
    use defs_models
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t

    use omp_lib ! export OMP_NUM_THREADS=5
    
    implicit none

    real :: M, ntot, press, T, rho

    real :: t1, t2, total_time ! time

    real, dimension(5) :: y, x

    integer :: i, j, k, n

    integer :: num_threads

    type(transport_in), dimension(:), allocatable :: transport
    ! type(cv_out) :: cv
    ! type(omega_int) :: omega_test
    ! type(bracket_int) :: bracket_test
    type(transport_out), dimension(:), allocatable :: transport_coeff

    character(len=*), parameter :: interaction = 'ESA-Bruno' ! 'VSS' 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'


    x = [0.77999, 0.19999, 0.01999, 0.00086999, 0.00099]

    n = 10000 ! iterations
    press = 100000.0

    allocate(transport(n))
    allocate(transport_coeff(n))

    ! y(1) = 0.756656E+00   
    ! y(2) = 0.221602E+00   
    ! y(3) = 0.207714E-01   
    ! y(4) = 0.421982E-03   
    ! y(5) = 0.548507E-03

    ! get the number of threads
    !$omp parallel
    num_threads = omp_get_num_threads()
    !$omp end parallel

    print *, "Number of threads: ", num_threads

    call cpu_time(t1)

    !$omp parallel do private(t, rho, y) schedule(dynamic, 100)

    do k = 500, n

        T = k * 1.
        ntot = press/kb/T
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

        transport(k)%temp = T
        transport(k)%mass_fractions = y
        transport(k)%rho = rho

        
        !call SpHeat(transport[k], cv)
        !call OmegaInt(transport%temp, omega_test)
        !call BracketInt(transport%temp, x, omega_test, bracket_test)
        call Transport1TSimpl(transport(k), transport_coeff(k), interaction)

        ! write (6, *) 'Process num. ', OMP_GET_THREAD_NUM()

    end do

    !$omp end parallel do


    call cpu_time(t2)

    total_time = t2 - t1

    print '(a,f10.3,a)', "Total execution time: ", total_time, " seconds"
    print '(a,f10.3,a)', "Average time per iteration: ", total_time / (n-499), " seconds"

    open(6, file='../res2/air5_1T_test_omp.txt', status='unknown')

    do k = 500, n

        T = k * 1.

        write (6, *) 'INPUT DATA:'

        write (6, *)


        write (6, *) 'Temperature, K         ',T
        write (6, *) 'Pressure, Pa           ',press
        write (6, *) 'N2 molar fraction      ',x(1)
        write (6, *) 'O2 molar fraction      ',x(2)
        write (6, *) 'NO molar fraction      ',x(3)
        write (6, *) 'N molar fraction       ',x(4)
        write (6, *) 'O molar fraction       ',x(5)

        write (6, *)

        write (6, *) 'TRANSPORT COEFFICIENTS:'
        write (6, *)

        write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff(k)%visc
        write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff(k)%ltot
        ! write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, tr , W/m/K        ', ltr
        ! write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, int , W/m/K       ', lint
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_N2, W/m/K     ', lvibr_n2
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_O2, W/m/K     ', lvibr_o2
        ! write (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_NO, W/m/K     ', lvibr_no
        
        ! write (6, *)
        ! write (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        ! write (6, *)


        ! do i=1,5
        !     write (6, '(1x, 5E15.6)') (transport_coeff(k)%DIFF(i,j), j=1,5)
        ! end do

        ! write (6, *)

        write (6, *)
        write (6, *) 'EFFECTIVE DIFFUSION COEFFICIENTS D_i, m^2/s'
        write (6, *)


        write (6, '(1x, 5E15.6)') (transport_coeff(k)%effDiff(i), i=1,NUM_SP)

        write (6, *)


    end do

    close(6)

end program