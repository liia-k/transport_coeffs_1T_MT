program test_potential_models

    ! This program calculates transport coefficients for a given mixture of air with varied temperature.
    ! Simplified model is applied: only shear viscosity, thermal conductivity
    ! and effective diffusion coefficients are calculated.
    ! Different interaction potentials are used: Lennard-Jones, Born-Mayer, VSS, ESA-Bruno.
    use constant_air5
    use defs_models
    use transport_1t
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k
    integer, parameter :: start = 500, end = 9000, step = 100
    integer, parameter :: size = (end - start) / step + 1
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    real, dimension(4, size) :: shearVisc, thermCond
    real, dimension(size) :: tempArray
    character(len=25), dimension(4) :: interaction 

    ! Initialize interaction models
    interaction(1) = 'Lennard-Jones'
    interaction(2) = 'Born-Mayer'
    interaction(3) = 'VSS'
    interaction(4) = 'ESA-Bruno'

    ! Initialize molar fractions
    x(1) = 0.77999 
    x(2) = 0.19999 
    x(3) = 0.01999
    x(4) = 0.00086999 
    x(5) = 0.00099 

    press = 100000

    j = 0
    do k = start, end, step
        j = j + 1
        T = k
        ntot = press / T / Kb
        rho = sum(MASS_SPCS * ntot * x)
        y = (ntot / rho) * x * MASS_SPCS
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        tempArray(j) = T
        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        do i=1,4
            call Transport1TSimpl(transport, transport_coeff, interaction(i))
            shearVisc(i,j) = transport_coeff%visc
            thermCond(i,j) = transport_coeff%ltot
        end do

    end do

    ! Write shear viscosity coefficients to file
    open(6, file='../res/air5-1T-shear-visc.txt', status='unknown')
    write (6, *)
    write (6, '(1x, A40)') 'Shear viscosity coefficients:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)
    do k = 1, size
        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), shearVisc(:, k)
    end do
    close(6)

    ! Write thermal conductivity coefficients to file
    open(6, file='../res/air5-1T-therm-cond.txt', status='unknown')
    write (6, *)
    write (6, '(1x, A40)') 'Thermal conductivity coefficient:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)
    do k = 1, size
        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), thermCond(:, k)
    end do
    close(6)

end program test_potential_models