program test_omp_simplified

    ! This program calculates transport coefficients for a given mixture of air with varied temperature.
    ! Simplified model is applied: only shear viscosity, thermal conductivity 
    ! and effective diffusion coefficients are calculated.
    ! The program is parallelized with OpenMP.
    
    use constant_air5
    use defs_models
    ! use specific_heat_sp
    ! use omega_integrals
    ! use bracket_integrals
    use transport_1t
    use omp_lib ! export OMP_NUM_THREADS=5
    
    implicit none

    ! Variables
    real(8) :: M, ntot, press, T, rho
    real :: t1, t2, total_time
    real(8), dimension(5) :: y, x
    integer :: i, k, n, num_threads
    type(transport_in), dimension(:), allocatable :: transport
    type(transport_out), dimension(:), allocatable :: transport_coeff
    character(len=*), parameter :: interaction = 'ESA-Bruno' ! 'VSS', 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'

    ! Initialize molar fractions
    x = [0.77999, 0.19999, 0.01999, 0.00086999, 0.00099]

    ! Number of iterations and pressure
    n = 10000
    press = 100000.0

    ! Allocate arrays
    allocate(transport(n))
    allocate(transport_coeff(n))

    ! Get the number of threads
    !$omp parallel
    num_threads = omp_get_num_threads()
    !$omp end parallel

    print *, "Number of threads: ", num_threads

    ! Start timing
    call cpu_time(t1)

    ! Parallel loop for transport coefficient calculations
    !$omp parallel do private(T, rho, y) schedule(dynamic, 100)
    do k = 500, n
        T = k * 1.0
        ntot = press / kb / T
        rho = sum(MASS_SPCS * ntot * x)
        y = (ntot / rho) * x * MASS_SPCS

        transport(k)%temp = T
        transport(k)%mass_fractions = y
        transport(k)%rho = rho

        call Transport1TSimpl(transport(k), transport_coeff(k), interaction)
    end do
    !$omp end parallel do
    
    ! End timing
    call cpu_time(t2)
    total_time = t2 - t1

    print '(a,f10.3,a)', "Total execution time: ", total_time, " seconds"
    print '(a,f10.3,a)', "Average time per iteration: ", total_time / (n - 499), " seconds"

    ! Open file for results
    open(6, file='../res/air5-1T-res-omp.txt', status='unknown')


    ! Write results to file
    do k = 500, n
        T = k * 1.0

        write (6, *) 'INPUT DATA:'
        write (6, *)
        write (6, *) 'Temperature, K         ', T
        write (6, *) 'Pressure, Pa           ', press
        write (6, *) 'N2 molar fraction      ', x(1)
        write (6, *) 'O2 molar fraction      ', x(2)
        write (6, *) 'NO molar fraction      ', x(3)
        write (6, *) 'N molar fraction       ', x(4)
        write (6, *) 'O molar fraction       ', x(5)
        write (6, *)
        write (6, *) 'TRANSPORT COEFFICIENTS:'
        write (6, *)
        write (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff(k)%visc
        write (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff(k)%ltot
        write (6, *)
        write (6, *) 'EFFECTIVE DIFFUSION COEFFICIENTS D_i, m^2/s'
        write (6, *)
        write (6, '(1x, 5E15.6)') (transport_coeff(k)%effDiff(i), i = 1, NUM_SP)
        write (6, *)
    end do

    close(6)

end program test_omp_simplified