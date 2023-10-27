! Module for simplified (no thermal diffusion and bulk viscocity) calculation of transport coefficients. 
! Uses the module constant_...(name of mixture).f90 containing main constants and 
! variables definition, modules for calculation of specific heats, omega and bracket integrals.

module transport_1t_simpl

    use defs_models
    
    implicit none

    contains

    subroutine Transport1TSimpl(data_in, data_out)

        use constant_air5
        use specific_heat
        use bracket_integrals
        use qr_decomposition

        type(transport_in),intent(in)   :: data_in
        type(transport_out),intent(out) :: data_out 
        
        integer i, j, delta
        ! REAL CU, CUT

        !total heat conductivity; translational heat conductivity;
        !internal heat conductivity;
        !shear viscosity

        real ltot, ltr, lint, visc

        !Diffusion coeffcients matrix

        real, dimension(NUM_SP,NUM_SP) :: DIFF

        !Matrices for the linear transport systems defining
        !heat conductivity and thermal diffusion (LTH);
        !bulk viscosity (BVISC);
        !diffusion (LDIFF);
        !shear viscisity (HVISC).

        real, dimension(NUM_SP,NUM_SP) :: LTH, Q_LTH, R_LTH, Inverse_LTH

        real, dimension(NUM_SP,NUM_SP) ::  LDIFF, Q_LDIFF, R_LDIFF, Inverse_LDIFF, &
                                           HVISC, Q_HVISC, R_HVISC, Inverse_HVISC, &
                                           b1

        !Vectors of right hand terms

        real, dimension(NUM_SP,1) :: b

        real, dimension(NUM_SP,1) :: b2

        real, dimension(NUM_SP) :: X, Y

        real T, ntot, rho, M

        type(SpHeatVOut) :: cv
        type(omega_int) :: omega_out
        type(bracket_int) :: bracket_out

        T = data_in%temp
        rho = data_in%rho
        y = data_in%mass_fractions

        ! Since the transport coefficient are barely affected if one of the mass fractions is of order 1e-6 and smaller,
        ! the values of such mass fractions are set to be equal to the mentioned order:
        do i=1,NUM_SP
            if (abs(y(i)) < 1e-6) then
                y(i) = 1e-6
            end if
        end do

        M = 1/dot_product(Y,1/MOLAR)
        ntot = sum(rho*Y/MASS_SPCS)
        x = (rho/ntot)*Y/MASS_SPCS

        call SpHeat(T, y, cv)
        call OmegaInt(T, omega_out)
        call BracketInt(T, x, omega_out, bracket_out)

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
            
            
        ! Definition of matrix b1 (right hand side of system for calculation of
        ! diffuaion coefficients)

        do i=1,NUM_SP
            do j=1,NUM_SP
            if (i==j) then 
                delta = 1
            else
                delta = 0
            end if
            B1(i,j) = 8./25./Kb*(delta - y(i))!3*kb*T*(delta-y(i));
            end do
        end do
        
        B1(1, 1:NUM_SP) = 0
       
        ! End of matrix b1 definition
            
        ! Definition of vector b2 (right hand side of system for calculation of
        ! shear viscocity coefficient)

        b2(1:NUM_SP, 1) = (2./Kb/T)*x
       
        ! End of vector b2 definition
            
            
        ! Linear system solution using the Gauss method
        ! The solutions a, d, h, f are written to b, b1, b2, b3, respectively 
        
        call QRDecomposition(LTH,Q_LTH,R_LTH,NUM_SP)
        call QRDecomposition(LDIFF,Q_LDIFF,R_LDIFF,NUM_SP)
        call QRDecomposition(HVISC,Q_HVISC,R_HVISC,NUM_SP)

        call InvertQR(Q_LTH,R_LTH,Inverse_LTH,NUM_SP)
        call InvertQR(Q_LDIFF,R_LDIFF,Inverse_LDIFF,NUM_SP)
        call InvertQR(Q_HVISC,R_HVISC,Inverse_HVISC,NUM_SP)

        b = matmul(Inverse_LTH,b)
        B1 = matmul(Inverse_LDIFF,B1)
        b2 = matmul(Inverse_HVISC,b2)
        
        ! call gaussj(LTH,10,10,b,1,1)
        ! call gaussj(Ldiff,5,5,b1,5,5)
        ! call gaussj(HVISC,5,5,b2,1,1)

        ! Thermal diffusion coefficients THDIF(i)

        ! thdiff = -(1./2./ntot)*b(:NUM_SP,1)


        ! Thermal conductivity coefficient associated to translational
        ! energy, LTR 

        Ltr = (5./4.)*kb*sum(x*b(: NUM_SP,1))

        ! Thermal conductivity coefficients associated to internal
        ! energies 
        
        Lint = (3./16.)*T*sum(x(:NUM_MOL)*cv%cv_int_sp(:NUM_MOL)*((Kb)*(MASS_SPCS(:NUM_MOL)*1e30) &
                /(bracket_out%lambda_int(:NUM_MOL)*1e30)))

        ! Total thermal conductivity coefficient at the translational
        ! temperature gradient

        ltot = ltr + lint


        ! Diffusion coefficients	DIFF(i,j)

        diff = (1./2./ntot)*b1

        ! Shear viscosity coefficient VISC

        visc = (Kb*T/2.)*sum(x*b2(:NUM_SP,1))

       

        data_out%visc = visc
        data_out%ltot = ltot
        data_out%diff = diff
    
    end subroutine Transport1TSimpl
    

end module transport_1t_simpl
