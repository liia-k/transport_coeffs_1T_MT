!Test for calculation of all pre-required data: specific heats, omega and bracket integrals

PROGRAM test1

    USE CONSTANT
    USE SPECIFIC_HEAT
    USE OMEGA_INTEGRALS
    USE BRACKET_INTEGRALS
    USE TRANSPORT_AIR5_1T
    
    
    IMPLICIT NONE

    INTEGER I, J, k
    
    CALL ENERGY

! Input parameters: species molar fractions, temperatures,
! pressure
 

    x(1)=0.77999
    x(2)=0.19999
    x(3)=0.01999
    x(4)=0.00086999
    x(5)=0.00099

    press=100000!101325

    DO k = 1, 6

        T=1000*k

        ntot=press/kb/T

        rho=0
        do i=1,5
            rho = rho+x(i)*mass(i)*ntot
        end do
        do i=1,5
            y(i) = x(i)*mass(i)*ntot/rho
        end do


    ! Calculation of vibrational energy,
    ! partition functions and specific heats


        CALL PART_FUNC_N2(T)
        CALL PART_FUNC_O2(T)
        CALL PART_FUNC_NO(T)

        CALL S_HEAT_N2
        CALL S_HEAT_O2
        CALL S_HEAT_NO
        CALL S_HEAT


    ! Calculation of bracket integrals

        CALL OMEGA
        CALL BRACKET
        
    ! Calculation of transport coeffs.

        CALL TRANSPORT_1T

    ! WRITE (6, *) 'Matrices of LS:'   

    !     WRITE (6, *) 'LTH:'
    !     do i=1,10
    !         !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
    !         WRITE (6, '(1x, 10E15.6)') (LTH(i,j), j=1,10)
    !     end do
    !     WRITE (6, *) 'BVICS:'
    !     do i=1,8
    !         !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
    !         WRITE (6, '(1x, 8E15.6)') (BVISC(i,j), j=1,8)
    !     end do
    !     WRITE (6, *) 'LDIFF:'
    !     do i=1,5
    !         !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
    !         WRITE (6, '(1x, 5E15.6)') (LDIFF(i,j), j=1,5)
    !     end do
    !     WRITE (6, *) 'HVISC:'
    !     do i=1,5
    !         !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
    !         WRITE (6, '(1x, 5E15.6)') (HVISC(i,j), j=1,5)
    !     end do

    !     WRITE (6, *) 'Vectors of right hand terms:' 

    !     WRITE (6, *) 'b1:'
    !     do i=1,5
    !         !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
    !         WRITE (6, '(1x, 5E15.6)') (b1(i,j), j=1,5)
    !     end do

    !     WRITE (6, *) 'b:'
    !     WRITE (6, '(1x, 10E15.6)') (b(j,1), j=1,10)

    !     WRITE (6, *) 'b2:'
    !     WRITE (6, '(1x, 5E15.6)') (b2(j,1), j=1,5)

    !     WRITE (6, *) 'b3:'
    !     WRITE (6, '(1x, 8E15.6)') (b3(j,1), j=1,8)
        
    ! Output

        open(6,file='air5_1Ttransport.txt',status='unknown')

        WRITE (6, *) 'INPUT DATA:'

        WRITE (6, *)


        WRITE (6, *) 'Temperature, K         ',t
        WRITE (6, *) 'Pressure, Pa           ',press
        WRITE (6, *) 'N2 molar fraction      ',x(1)
        WRITE (6, *) 'O2 molar fraction      ',x(2)
        WRITE (6, *) 'NO molar fraction      ',x(3)
        WRITE (6, *) 'N molar fraction       ',x(4)
        WRITE (6, *) 'O molar fraction       ',x(5)

        ! WRITE (6, *)

        ! WRITE (6, *) 'Calculation parameters required:'
        ! WRITE (6, *)

        ! !WRITE (6, '(1x, A45, E13.5)') 'Internal energy N2, J              ', EN_N2
        ! !WRITE (6, '(1x, A45, E13.5)') 'Internal energy O2, J              ', EN_O2
        ! !WRITE (6, '(1x, A45, E13.5)') 'Internal energy NO, J              ', EN_NO
        ! WRITE (6, '(1x, A45, E13.5)') 'Omega13_11, J                      ', Omega13(1,1)


        ! WRITE (6, *)

        ! WRITE (6, *) 'Omega Integrals Omega11_ij'

        ! WRITE (6, *)

        ! do i=1,5
        !     !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
        !     WRITE (6, '(1x, 5E15.6)') (omega12(i,j), j=1,5)
        ! end do

        ! WRITE (6, *)

        ! WRITE (6, *) 'Backet integrals b01_ij'
        
        ! WRITE (6, *)

        ! do i=1,5
        !     !WRITE (6, '(0x, 5e15.6)') (DIFF(i,j), j=1,5)
        !     WRITE (6, '(1x, 5E15.6)') (beta01(i,j), j=1,5)
        ! end do

        ! WRITE (6, *)

        ! WRITE (6, *) 'Backet integrals b11_ij'
        
        ! WRITE (6, *)

        ! do i=1,5
        !     WRITE (6, '(1x, 5E15.6)') (beta11(i,j), j=1,5)
        ! end do
        
        ! WRITE (6, *)

        ! WRITE (6, *) 'Specific internal heats:'
        ! WRITE (6, *)

        ! WRITE (6, '(1x, A45, E13.5)') 'C_vibr,N2              ', C_vibr_N2
        ! WRITE (6, '(1x, A45, E13.5)') 'C_v_int,N2             ', c_int(1)
        ! WRITE (6, '(1x, A45, E13.5)') 'C_v_int,O2             ', c_int(2)
        ! WRITE (6, '(1x, A45, E13.5)') 'C_v_int,NO             ', c_int(3)
        ! WRITE (6, '(1x, A45, E13.5)') 'C_v_int,               ', cv_int
        ! WRITE (6, '(1x, A45, E13.5)') 'C_v,                   ', CV_tot

        WRITE (6, *)

        WRITE (6, *) 'TRANSPORT COEFFICIENTS:'
        WRITE (6, *)

        WRITE (6, '(1x, A45, E13.5)') 'Shear viscosity coefficient, Pa.S             ', visc
        WRITE (6, '(1x, A45, E13.5)') 'Bulk viscosity coefficient, Pa.s              ', bulk_visc
        WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, W/m/K             ', ltot
        WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, tr , W/m/K        ', ltr
        WRITE (6, '(1x, A45, E13.5)') 'Thermal cond. coef. lambda, int , W/m/K       ', lint
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_N2, W/m/K     ', lvibr_n2
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_O2, W/m/K     ', lvibr_o2
        !WRITE (6, '(1x, A45, E13.5)') 'Vibr. therm. cond. coef. lambda_NO, W/m/K     ', lvibr_no
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N2, m^2/s          ', THDIFF(1)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O2, m^2/s          ', THDIFF(2)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of NO, m^2/s          ', THDIFF(3)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of N, m^2/s           ', THDIFF(4)
        WRITE (6, '(1x, A45, E13.5)') 'Thermal diffusion coef. of O, m^2/s           ', THDIFF(5)

        WRITE (6, *)
        WRITE (6, *) 'DIFFUSION COEFFICIENTS D_ij, m^2/s'
        WRITE (6, *)


        do i=1,5
            WRITE (6, '(1x, 5E15.6)') (DIFF(i,j), j=1,5)
        end do

    END DO

END