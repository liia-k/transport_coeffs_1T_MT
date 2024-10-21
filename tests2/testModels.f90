program testModels

    use constant_air5
    use defs_models
    use transport_1t
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k
    integer, parameter :: start = 500, end = 9000, step = 100
    integer, parameter :: size = (end - start)/step + 1

    type(transport_in) :: transport
    type(transport_out) :: transport_coeff

    real, dimension(4, size) :: shearVisc, thermCond
    real, dimension(size) :: tempArray

    character(len=25), dimension(4) :: interaction 
                      
    interaction(1) = 'Lennard-Jones'
    interaction(2) = 'Born-Mayer'
    interaction(3) = 'VSS'
    interaction(4) = 'ESA-Bruno'

    x(1) = 0.77999 
    x(2) = 0.19999 
    x(3) = 0.01999
    x(4) = 0.00086999 
    x(5) = 0.00099 
   

    press = 100000

    ! do i=1,NUM_SP
	! 	do j=1,NUM_SP
    !         write (*, *) gamma(3 - OMEGA_VSS(i,j))
    !     end do
    ! end do

    j = 0
    do k = start, end, step

        j = j + 1
        T = k
        ntot = press/T/Kb
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

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

    open(6, file='../res2/shearVisc.txt', status='unknown')

    write (6, *)
    write (6, '(1x, A40)') 'Shear viscosity coefficients:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)

    do k = 1, size

        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), shearVisc(:,k)

    end do

    close(6)

    open(6, file='../res2/thermCond.txt', status='unknown')

    write (6, *)
    write (6, '(1x, A40)') 'Thermal conductivity coefficient:'
    write (6, *)
    write (6, '(1x, A6, 4(A14, 1X))') 'Temp.', interaction(:)

    do k = 1, size

        write (6, '(1x, F5.0, 4(E13.6, 1X))') tempArray(k), thermCond(:,k)

    end do

    close(6)

    

end program