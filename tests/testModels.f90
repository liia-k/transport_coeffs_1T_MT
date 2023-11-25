program testModels

    use constant_air5
    use defs_models
    use transport_1t_simpl
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k
    integer, parameter :: start = 500, end = 9000, step = 100
    integer, parameter :: size = (end - start)/step + 1

    type(transport_in) :: transport
    type(transport_out) :: transport_coeff

    real, dimension(4, size) :: shearVisc, thermCond

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

    j = 0
    do k = start, end, step

        j = j + 1
        T = k
        ntot = press/T/Kb
        rho = sum(MASS_SPCS*ntot*x)
        ! M = dot_product(X,MOLAR)
        ! rho = M*press/R/T
        
        y = (ntot/rho)*x*MASS_SPCS

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        do i=1,4
            call Transport1TSimpl(transport, transport_coeff, interaction(i))
            shearVisc(i,j) = transport_coeff%visc
            thermCond(i,j) = transport_coeff%ltot
        end do

    end do

    open(6, file='../res/shearVisc.txt', status='unknown')

    write (6, *)
    write (6, '(1x, A40)') 'Shear viscosity coefficient:'
    write (6, *)
    write (6, '(4x, A15, A15, A15, A15)') interaction(:)
    write (6, *)

    do k = 1, size

        write (6, '(4x, E13.6, E13.6, E13.6, E13.6)') shearVisc(:,k)

    end do

    close(6)

    open(6, file='../res/thermCond.txt', status='unknown')

    write (6, *)
    write (6, '(1x, A40)') 'Thermal conductivity coefficient:'
    write (6, *)
    write (6, '(4x, A15, A15, A15, A15)') interaction(:)
    write (6, *)

    do k = 1, size

        write (6, '(4x, E13.6, E13.6, E13.6, E13.6)') thermCond(:,k)

    end do

    close(6)


end program