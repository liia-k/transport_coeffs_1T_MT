PROGRAM test_int

    USE CONSTANT
    USE SPECIFIC_HEAT
    USE OMEGA_INTEGRALS
    USE BRACKET_INTEGRALS
    USE TRANSPORT_AIR5_1T

    USE OMP_LIB !OMP_NUM_THREADS=5
    
    
    IMPLICIT NONE

    REAL :: M,M1,ntot,press,T,rho

    REAL, DIMENSION(5) :: y, x

    INTEGER :: I1, I, J, K, N

    type(transport_in) :: transport
    type(cv_out) :: cv
    type(omega_int) :: omega_test
    type(bracket_int) :: bracket_test
    type(transport_out) :: transport_coeff


    x(1)=0.77999
    x(2)=0.19999
    x(3)=0.01999
    x(4)=0.00086999
    x(5)=0.00099

    N = 10000

    ! y(1) = 0.756656E+00   
    ! y(2) = 0.221602E+00   
    ! y(3) = 0.207714E-01   
    ! y(4) = 0.421982E-03   
    ! y(5) = 0.548507E-03

    press = 100000.

    !ntot = 0.724296E+25
    !rho = 0.347314E+00

    open(6,file='air5_1Ttest_data.txt',status='unknown')

    WRITE (6, *) 'INPUT DATA:'

        WRITE (6, *)

        WRITE (6, *) 'Pressure, Pa           ',press
        WRITE (6, *) 'N2 molar fraction      ',x(1)
        WRITE (6, *) 'O2 molar fraction      ',x(2)
        WRITE (6, *) 'NO molar fraction      ',x(3)
        WRITE (6, *) 'N molar fraction       ',x(4)
        WRITE (6, *) 'O molar fraction       ',x(5)

        WRITE (6, *)
        WRITE (6, *) 'TRANSPORT COEFFICIENTS:'
        WRITE (6, *)

        WRITE (6, *) 'Temperature, K ',  'Shear viscosity, Pa.S ', 'Bulk viscosity, Pa.s ', &
         'Thermal cond., W/m/K ', 'Thermal diffusion N2, m^2/s ', 'Thermal diffusion O2, m^2/s ', &
         'Thermal diffusion NO, m^2/s ', 'Thermal diffusion N, m^2/s ', 'Thermal diffusion O, m^2/s '

    DO k = 500, N, 500

        T = k * 1.
        ntot = press/kb/T
        rho = 0
        do i1 = 1,5
            rho = rho + x(i1)*mass(i1)*ntot
        end do

        do i1 = 1,5
            y(i1) = x(i1)*mass(i1)*ntot/rho
        end do

        transport%temp = T
        transport%mass_fractions = y
        transport%rho = rho

        CALL TRANSPORT_1T(transport, transport_coeff)


        WRITE (6, '(1x, 9E15.6)') t, transport_coeff%visc, transport_coeff%bulk_visc, transport_coeff%ltot, &
        transport_coeff%THDIFF(1), transport_coeff%THDIFF(2), transport_coeff%THDIFF(3), transport_coeff%THDIFF(4), &
        transport_coeff%THDIFF(5)

        ! WRITE (6, *)
        ! WRITE (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        ! WRITE (6, *)


        ! do i=1,5
        !     WRITE (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i,j), j=1,5)
        ! end do

        ! WRITE (6, *)

    END DO

    close(6)


END PROGRAM