program test_n2_bulk_visc

! Test program for transport properties calculation based on 1T model (general procedure)
    use constant_air5
    use defs_models
    use transport_1t
    
    
    implicit none

    ! Variables
    real :: M, ntot, press, T, rho
    real, dimension(NUM_SP) :: y, x
    integer :: i, j, k, n
    type(transport_in) :: transport
    type(transport_out) :: transport_coeff
    character(len=*), parameter :: interaction = 'Lennard-Jones' ! 'VSS', 'Lennard-Jones', 'Born-Mayer', 'ESA-Bruno'

    ! Mostly N2 mixture
    y(1) = 0.9999
    y(2) = 0.00002
    y(3) = 0.00002
    y(4) = 0.00002
    y(5) = 0.00002

    press = 100000 ! Pa
    n = 4000

    ! Open file for results
    open(6, file='../res/n2-1T-bulk-visc.txt', status='unknown')

    ! Write input data to file
    write (6, *) 'Input data'
    write (6, *)
    write (6, '(A25,E13.6)') 'Pressure, Pa:       ', press
    ! write (6, '(A25,E13.6)') 'N2 molar fraction:  ', x(1)
    ! write (6, '(A25,E13.6)') 'O2 molar fraction:  ', x(2)
    ! write (6, '(A25,E13.6)') 'NO molar fraction:  ', x(3)
    ! write (6, '(A25,E13.6)') 'N molar fraction:   ', x(4)
    ! write (6, '(A25,E13.6)') 'O molar fraction:   ', x(5)
    write (6, *)
    write (6, *) 'Potential model: ', interaction
    write (6, *)
    write (6, '(1x, A13, A13, A13)') 'Temp. ', 'Bulk visc. ', 'Shear visc.'
    write (6, *)
    
    ! Loop over temperature range
    do k = 400, n, 100
        T = k
        M = 1 / sum(y / Molar)
        rho = M * press / T / R 

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        call Transport1TGeneral(transport, transport_coeff, interaction)

        write (6, '(1x, F13.6, 1x, E13.6, 1x, E13.6)') T, transport_coeff%bulk_visc, transport_coeff%visc
    end do

    close(6)

end program test_n2_bulk_visc