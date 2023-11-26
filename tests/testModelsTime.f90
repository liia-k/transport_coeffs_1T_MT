program testModelsTime

    use constant_air5
    use defs_models
    use transport_1t_simpl
    
    implicit none

    real :: M, ntot, press, T, rho

    real, dimension(NUM_SP) :: y, x

    integer :: i, j, k
    integer, parameter :: begin = 300, end = 9000, step = 1

    type(transport_in) :: transport
    type(transport_out) :: transport_coeff

    real(8), dimension(4) :: timeArray

    real(8) :: start_time, end_time

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

    do i=1,4
        call cpu_time(start_time)
        do k = begin, end, step

            T = k
            ntot = press/T/Kb
            rho = sum(MASS_SPCS*ntot*x)
            ! M = dot_product(X,MOLAR)
            ! rho = M*press/R/T
            
            y = (ntot/rho)*x*MASS_SPCS

            transport%temp = T
            transport%mass_fractions = y
            transport%rho = rho

            
            call Transport1TSimpl(transport, transport_coeff, interaction(i))
    
        end do
        call cpu_time(end_time)
        timeArray(i) = end_time - start_time
    end do

    write (6, '(4x, A15, A15, A15, A15)') (interaction(k), k=1,4)
    write (6, '(8x, E13.6, E13.6, E13.6, E13.6)') (timeArray(k), k=1,4)

end program