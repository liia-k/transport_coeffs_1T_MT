PROGRAM test_int

    USE CONSTANT
    USE SPECIFIC_HEAT
    USE OMEGA_INTEGRALS
    USE BRACKET_INTEGRALS
    USE TRANSPORT_AIR5_1T

    USE OMP_LIB
    
    
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

    ! y(1) = 0.756656E+00   
    ! y(2) = 0.221602E+00   
    ! y(3) = 0.207714E-01   
    ! y(4) = 0.421982E-03   
    ! y(5) = 0.548507E-03

    press = 100000

    !ntot = 0.724296E+25
    !rho = 0.347314E+00

    N = 10000

    DO k = 500, N

        T = 10000 * k*1./N
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

        
        !CALL S_Heat(transport[k], cv)
        !CALL Omega(transport%temp, omega_test)
        !CALL Bracket(transport%temp, x, omega_test, bracket_test)
        CALL TRANSPORT_1T(transport, transport_coeff)

        ! WRITE (6, *) 'INPUT DATA:'

        ! WRITE (6, *)


        ! WRITE (6, *) 'Temperature, K         ',transport%temp
        ! WRITE (6, *) 'Molar mass, kg         ',M
        ! WRITE (6, *) 'Density, kg/m^3        ',transport%rho
        ! WRITE (6, *) 'Number density, 1/m^3  ',ntot
        ! WRITE (6, *) 'N2 mass fraction       ',transport%mass_fractions(1)
        ! WRITE (6, *) 'O2 massr fraction      ',transport%mass_fractions(2)
        ! WRITE (6, *) 'NO mass fraction       ',transport%mass_fractions(3)
        ! WRITE (6, *) 'N mass fraction        ',transport%mass_fractions(4)
        ! WRITE (6, *) 'O mass fraction        ',transport%mass_fractions(5)

        ! WRITE (6, *)


        ! WRITE (6, *) 'Total internal Specific heat  ',cv%cv_int
        ! WRITE (6, *) 'Total  Specific heat          ',cv%cv_tot

        ! do i1=1,5
        !     WRITE (6, *) cv%c_int(i1)
        ! end do

        ! WRITE (6, *)

        ! WRITE (6, *) 'Omega Integrals Omega11_ij'

        ! WRITE (6, *)

        ! do i=1,5
        !     WRITE (6, '(1x, 5e15.6)') (omega_test%omega11(i,j), j=1,5)
        !     WRITE (6, '(1x, 5E15.6)') (omega_test%omega12(i,j), j=1,5)
        ! end do

        ! WRITE (6, *) 'Bracket Integrals Omega11_ij'

        ! WRITE (6, *)

        ! do i=1,5
        !     WRITE (6, '(1x, 5e15.6)') (bracket_test%lambda(i,j), j=1,5)
        !     WRITE (6, '(1x, 5E15.6)') (bracket_test%beta01(i,j), j=1,5)
        ! end do

        !WRITE (6, '(1x, 5E15.6)') (y(j), j=1,5)

        open(6,file='air5_1Ttest_old.txt',status='unknown')

        WRITE (6, *) 'INPUT DATA:'

        WRITE (6, *)


        WRITE (6, *) 'Temperature, K         ',t
        WRITE (6, *) 'Pressure, Pa           ',press
        WRITE (6, *) 'N2 molar fraction      ',x(1)
        WRITE (6, *) 'O2 molar fraction      ',x(2)
        WRITE (6, *) 'NO molar fraction      ',x(3)
        WRITE (6, *) 'N molar fraction       ',x(4)
        WRITE (6, *) 'O molar fraction       ',x(5)

        WRITE (6, *)

        WRITE (6, *) 'TRANSPORT COEFFICIENTS:'
        WRITE (6, *)

        WRITE (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', transport_coeff%visc
        WRITE (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa.s              ', transport_coeff%bulk_visc
        WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', transport_coeff%ltot
        !WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, tr , W/m/K        ', ltr
        !WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, int , W/m/K       ', lint
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_N2, W/m/K     ', lvibr_n2
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_O2, W/m/K     ', lvibr_o2
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_NO, W/m/K     ', lvibr_no
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N2, m^2/s          ', transport_coeff%THDIFF(1)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O2, m^2/s          ', transport_coeff%THDIFF(2)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of NO, m^2/s          ', transport_coeff%THDIFF(3)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N, m^2/s           ', transport_coeff%THDIFF(4)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O, m^2/s           ', transport_coeff%THDIFF(5)

        WRITE (6, *)
        WRITE (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        WRITE (6, *)


        do i=1,5
            WRITE (6, '(1x, 5E15.6)') (transport_coeff%DIFF(i,j), j=1,5)
        end do

        WRITE (6, *)

    END DO

    close(6)

END PROGRAM