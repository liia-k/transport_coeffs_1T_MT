! Module for simplified (no thermal diffusion and bulk viscocity) calculation of transport coefficients. 

module transport_1t_simpl

    use defs_models

    use, intrinsic :: ieee_arithmetic 
    
    implicit none

    contains

    logical function exception(var)

    ! has a true value if var in not NaN and Infinity

        real, intent(in) :: var

        exception = .not.(ieee_is_nan(var)) .and. (ieee_is_finite(var))

    end function

    subroutine macro_output(data_in)
    
    ! macroparameters given values output

        type(transport_in),intent(in)   :: data_in

        print *, "temperature = ", data_in%temp
        print *, "density = ", data_in%rho
        print *, "mass fractions: ", data_in%mass_fractions

    end subroutine



    subroutine Transport1TSimpl(data_in, data_out)

    ! Simplified calcultaion of transport coeffs

        use constant_air5
        use specific_heat_sp
        use omega_integrals
        use bracket_integrals
        use qr_decomposition


        type(transport_in),intent(in)   :: data_in ! macroparametrs
        type(transport_out),intent(out) :: data_out ! set of transport coeffs
        
        integer i, j, delta

        ! total heat conductivity; translational heat conductivity; internal heat conductivity;
        
        real ltot, ltr, lint
    
        ! shear viscosity

        real visc

        real, dimension(NUM_SP) :: effDiff

        ! diffusion coeffcients matrix

        real, dimension(NUM_SP,NUM_SP) :: DIFF, binDiff

        ! matrices for the linear transport systems defining:
        ! heat conductivity and thermal diffusion (LTH);
        ! bulk viscosity (BVISC);
        ! diffusion (LDIFF);
        ! shear viscisity (HVISC).

        real, dimension(NUM_SP,NUM_SP) :: LTH, Q_LTH, R_LTH, Inverse_LTH

        real, dimension(NUM_SP,NUM_SP) ::  LDIFF, Q_LDIFF, R_LDIFF, Inverse_LDIFF, &
                                           HVISC, Q_HVISC, R_HVISC, Inverse_HVISC, &
                                           b1

        ! vectors of the RHS of systems 

        real, dimension(NUM_SP,1) :: b

        real, dimension(NUM_SP,1) :: b2

        ! molar and mass fractions

        real, dimension(NUM_SP) :: x, Y

        ! macroparameters

        real T, ntot, rho, M, mij

        type(cv_out) :: cv
        type(omega_int) :: omega_out
        type(bracket_int) :: bracket_out


        T = data_in%temp
        rho = data_in%rho
        y = data_in%mass_fractions

        ! Since the transport coefficient are barely affected if one of the mass fractions is of order 1e-5 and smaller,
        ! the values of such mass fractions are set to be equal to the mentioned order:
        do i=1,NUM_SP
            if (abs(y(i)) < 1e-5) then
                y(i) = 1e-5
            end if
        end do

        M = 1/dot_product(Y,1/MOLAR)
        ntot = sum(rho*Y/MASS_SPCS)
        x = (rho/ntot)*Y/MASS_SPCS

        call SpHeat(T, y, cv)
        call OmegaInt(T, omega_out)
        call BracketInt(T, x, omega_out, bracket_out)

        ! Calculation of binary diffusion coefficients, for which only
        ! omega integrals and molar fractions are required:

        do i=1,NUM_SP
            do j=1,NUM_SP
              mij = MASS_SPCS(j)*(MASS_SPCS(i)*1E27)/(MASS_SPCS(i)*1E27 + MASS_SPCS(j)*1E27)
              binDiff(i,j) = 3.*kb*T / 16./ ntot / mij / omega_out%OMEGA11(i,j)
            end do
        end do

        ! Definition of matrix LTH for calculation of 
        ! thermal conductivity and thermal diffuaion coefficients
        ! The system has a form:
        ! LTH times a = b, a is the vector of unknowns
        
        LTH(: NUM_SP, : NUM_SP) = bracket_out%Lambda11
        
        ! End of matrix LTH definition


        ! Definition of matrix LDIFF for calculation of 
        ! diffuaion coefficients
        ! The system has a form:
        ! LDIFF times D = B1, D a is the matrix of unknowns
        
        LDIFF = LTH(1:NUM_SP, 1:NUM_SP)
        
        ! End of matrix LDIFF definition


        ! Definition of matrix HVISC for calculation of 
        ! shear viscocity coefficient
        ! The system has a form:
        ! HVISC times h = b2, h a is the vector of unknowns
        
        HVISC = bracket_out%H00
        
        ! End of matrix HVISC definition
            
            
        ! Definition of vector b (right hand side of system for calculation of
        ! thermal conductivity and thermal diffuaion coefficients)

        b(: NUM_SP, 1) = (4./5./Kb)*x
      
        ! End of vector b definition
            
            
        ! Definition of matrix B1 (right hand side of system for calculation of
        ! diffuaion coefficients)

        do i=1,NUM_SP
            do j=1,NUM_SP
            if (i==j) then 
                delta = 1
            else
                delta = 0
            end if
            B1(i,j) = 8./25./Kb*(delta - y(i))! 3*kb*T*(delta-y(i));
            end do
        end do
        
        B1(1, 1:NUM_SP) = 0
       
        ! End of matrix b1 definition
            
        ! Definition of vector b2 (right hand side of system for calculation of
        ! shear viscocity coefficient)

        b2(1:NUM_SP, 1) = (2./Kb/T)*x
       
        ! End of vector b2 definition
            
            
        ! Linear system solution using QR decomposition
        
        call QRDecomposition(LTH,Q_LTH,R_LTH,NUM_SP)
        call QRDecomposition(LDIFF,Q_LDIFF,R_LDIFF,NUM_SP)
        call QRDecomposition(HVISC,Q_HVISC,R_HVISC,NUM_SP)

        call InvertQR(Q_LTH,R_LTH,Inverse_LTH,NUM_SP)
        call InvertQR(Q_LDIFF,R_LDIFF,Inverse_LDIFF,NUM_SP)
        call InvertQR(Q_HVISC,R_HVISC,Inverse_HVISC,NUM_SP)

        ! Solutions:
        
        b = matmul(Inverse_LTH,b)
        B1 = matmul(Inverse_LDIFF,B1)
        b2 = matmul(Inverse_HVISC,b2)
        
       
        ! Thermal conductivity coefficient associated to translational energy, ltr

        ltr = (5./4.)*kb*sum(x*b(: NUM_SP,1))

        ! Thermal conductivity coefficients associated to internal energies, lint 
        
        lint = (3./16.)*T*sum(x(:NUM_MOL)*cv%cv_int_sp(:NUM_MOL)*((Kb)*(MASS_SPCS(:NUM_MOL)*1e30) &
                /(bracket_out%lambda_int(:NUM_MOL)*1e30)))

        ! Total thermal conductivity coefficient at the translational temperature gradient

        ltot = ltr + lint

        if (exception(ltot)) then
            data_out%ltot = ltot
        else
            print *, "thermal conductivity is not calculated for the set:"
            call macro_output(data_in)
        end if

        ! Diffusion coefficients DIFF(i,j)

        DIFF = (1./2./ntot)*b1

        do i=1,NUM_SP
            do j=1,NUM_SP
                if (exception(DIFF(i,j))) then
                    data_out%DIFF(i,j) = DIFF(i,j)
                else
                    print *, "diffusion coeffs are not calculated for the set:"
                    call macro_output(data_in)
                end if
            end do
        end do

        ! Shear viscosity coefficient VISC

        visc = (Kb*T/2.)*sum(x*b2(:NUM_SP,1))

        if (exception(visc)) then
            data_out%visc = visc
        else
            print *, "shear viscosity is not calculated for the set:"
            call macro_output(data_in)
        end if

        ! Calculation of effective diffusion coefficients

        do i=1,NUM_SP
            effDiff(i) = (1 - x(i)) / sum(x(:NUM_MOL)/binDiff(i,:NUM_MOL))
        end do

        do i=1,NUM_SP
            if (exception(effDiff(i))) then
                data_out%effDiff(i) = effDiff(i)
            else
                print *, "effective diffusion coeffs are not calculated for the set:"
                call macro_output(data_in)
            end if
        end do

    
    end subroutine Transport1TSimpl
    

end module transport_1t_simpl
