!Module for calculation of transport coefficients. Uses the module constant_...(name of mixture).f90 containing main constants and 
!variables definition, modules for calculation of specific heats, omega and bracket integrals.

module transport_1t

    use defs_models
    
    implicit none

    contains

    subroutine Transport1T(data_in, data_out)

        use constant_air5
        use specific_heat_sp
        use omega_integrals
        use bracket_integrals
        use qr_decomposition

        type(transport_in),intent(in)   :: data_in
        type(transport_out),intent(out) :: data_out 
        
        integer i, j, delta
        ! REAL CU, CUT

        !total heat conductivity; translational heat conductivity;
        !internal heat conductivity;
        !shear viscosity; bulk viscosity

        real ltot, ltr, lint, visc, bulk_visc!, lrot_n2, lrot_o2, lrot_no, lvibr_n2, lvibr_o2, lvibr_no, 

        !thermal diffusion coefficients (THDIFF);

        real, dimension(NUM_SP) :: THDIFF 

        !Diffusion coeffcients matrix

        real, dimension(NUM_SP,NUM_SP) :: DIFF

        !Matrices for the linear transport systems defining
        !heat conductivity and thermal diffusion (LTH);
        !bulk viscosity (BVISC);
        !diffusion (LDIFF);
        !shear viscisity (HVISC).

        real, dimension(2*NUM_SP,2*NUM_SP) :: LTH, Q_LTH, R_LTH, Inverse_LTH

        real, dimension(NUM_SP+NUM_MOL,NUM_SP+NUM_MOL) :: BVISC, Q_BVISC, R_BVISC, Inverse_BVISC

        real, dimension(NUM_SP,NUM_SP) ::  LDIFF, Q_LDIFF, R_LDIFF, Inverse_LDIFF, &
                                           HVISC, Q_HVISC, R_HVISC, Inverse_HVISC, &
                                           b1

        !Vectors of right hand terms

        real, dimension(2*NUM_SP,1) :: b

        real, dimension(NUM_SP,1) :: b2

        real, dimension(NUM_SP+NUM_MOL,1) :: b3


        real, dimension(NUM_SP) :: x, y

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

        M = 1/dot_product(y,1/MOLAR)
        ntot = sum(rho*y/MASS_SPCS)
        x = (rho/ntot)*y/MASS_SPCS

        call SpHeat(T, y, cv)
        call OmegaInt(T, omega_out)
        call BracketInt(T, x, omega_out, bracket_out)

        ! Definition of matrix LTH for calculation of 
        ! thermal conductivity and thermal diffuaion coefficients
        ! The system has a form:
        ! LTH times a = b, a is the vector of unknowns
        
        LTH(1:NUM_SP, 1:NUM_SP) = bracket_out%Lambda00
        ! DO i=1,5
        !     DO j=1,5
        !         LTH(i,j)=bracket_out%Lambda00(i,j)
        !     END DO
        ! END DO
        
        LTH(1:NUM_SP, NUM_SP+1 : 2*NUM_SP) = bracket_out%Lambda01
        ! DO i=1,5
        !     DO j=6,10
        !         LTH(i,j)=bracket_out%Lambda01(i,j-5)
        !     END DO
        ! END DO
        
        LTH(NUM_SP+1 : 2*NUM_SP, 1:NUM_SP) = LTH(1:NUM_SP, NUM_SP+1:2*NUM_SP)
        ! DO i=6,10
        !     DO j=1,5
        !         LTH(i,j)=LTH(j,i)
        !     END DO
        ! END DO
        
        LTH(NUM_SP+1 : 2*NUM_SP, NUM_SP+1 : 2*NUM_SP) = bracket_out%Lambda11
        ! DO i=6,10
        !     DO j=6,10
        !         LTH(i,j)=bracket_out%Lambda11(i-5,j-5)
        !     END DO
        ! END DO
        
        LTH(1, 1:NUM_SP) = y
        ! DO j=1,5
        !     LTH(1,j)=y(j)!x(j)*MASS_SPCS(j)*ntot/rho
        ! END DO
        
        LTH(1, NUM_SP+1 : 2*NUM_SP) = 0
        ! DO j=6,10
        !     LTH(1,j)=0.
        ! END DO
        ! End of matrix LTH definition


        ! Definition of matrix LDIFF for calculation of 
        ! diffuaion coefficients
        ! The system has a form:
        ! LDIFF times D = B1, D a is the matrix of unknowns
        
        LDIFF = LTH(1:NUM_SP, 1:NUM_SP)
        ! DO i=1,5
        !     DO j=1,5
        !     LDIFF(i,j)=LTH(i,j)
        !     END DO
        ! END DO
        ! End of matrix LDIFF definition


        ! Definition of matrix HVISC for calculation of 
        ! shear viscocity coefficient
        ! The system has a form:
        ! HVISC times h = b2, h a is the vector of unknowns
        
        HVISC = bracket_out%H00
        ! DO i=1,5
        !     DO j=1,5
        !     HVISC(i,j)=bracket_out%H00(i,j)
        !     END DO
        ! END DO
        ! End of matrix HVISC definition
            
            
        ! Definition of matrix BVISC for calculation of 
        ! bulk viscosity coefficients
        ! The system has a form:
        ! BVISC times f = b3, f is the vector of unknowns

        do i=1,NUM_SP
            BVISC(i, 1:NUM_SP) = x*bracket_out%beta11(i, 1:NUM_SP)
            ! DO j=1,5
            !     BVISC(i,j) = x(j)*bracket_out%beta11(i,j)
            ! END DO
        end do
        
        do i=1,NUM_SP
            BVISC(i, NUM_SP+1 : NUM_SP+NUM_MOL) = y(:NUM_MOL)*bracket_out%beta01(i, :NUM_MOL)
            ! DO j=6,8
            !     BVISC(i,j) = y(j-5)*bracket_out%beta01(i,j-5)
            ! END DO
        end do
        
        BVISC(NUM_SP+1 : NUM_SP+NUM_MOL, 1:NUM_SP) = transpose(BVISC(:NUM_SP, NUM_SP+1 : NUM_SP+NUM_MOL))
        ! DO i=6,8
        !     DO j=1,5
        !         BVISC(i,j)=BVISC(j,i)!x(j)*beta11(j,i)
        !     END DO
        ! END DO
        
        BVISC(NUM_SP+1 : NUM_SP+NUM_MOL, NUM_SP+1 : NUM_SP+NUM_MOL) = 0
        ! DO i=6,8
        !     DO j=6,8
        !         BVISC(i,j)=0
        !     END DO
        ! END DO
        
        do i=NUM_SP+1,NUM_SP+NUM_MOL
            BVISC(i,i) = y(i-5)*bracket_out%beta0011(i-5)
        end do
        ! BVISC(6,6)=y(1)*bracket_out%beta0011(1)
        ! BVISC(7,7)=y(2)*bracket_out%beta0011(2)
        ! BVISC(8,8)=y(3)*bracket_out%beta0011(3)
        
        BVISC(1, 1:NUM_SP) = ((3./2.)*R/M)*x ! Kb*ntot/rho
        ! DO j=1,5
        !     BVISC(1,j)=x(j)*3./2.*Kb*ntot/rho ! R/M
        ! END DO
        
        BVISC(1, NUM_SP+1 : NUM_SP+NUM_MOL) = y(1:NUM_MOL)*(kb/MASS_SPCS(1:NUM_MOL))

        ! DO j=6,6
        !     BVISC(1,j)=y(1)*kb/MASS_SPCS(1)
        ! END DO
        
        ! DO j=7,7
        !     BVISC(1,j)=y(2)*kb/MASS_SPCS(2)
        ! END DO
        
        ! DO j=8,8
        !     BVISC(1,j)=y(3)*kb/MASS_SPCS(3)
        ! END DO
        
        ! End of matrix BVISC definition
            
            
        ! Definition of vector b (right hand side of system for calculation of
        ! thermal conductivity and thermal diffuaion coefficients)

        b(1:NUM_SP,1) = 0
        b(NUM_SP+1 : 2*NUM_SP, 1) = (4./5./Kb)*x
        ! DO i=1,5
        !     b(i,1)=0.
        ! END DO
        ! DO i=6,10
        !     b(i,1)=4./5./kb*x(i-5)!15./2.*T/kb*x(i-5)
        ! END DO
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
        ! DO j=1,NUM_SP
        !     B1(1,j)=0
        ! END DO
        ! End of matrix b1 definition
            
        ! Definition of vector b2 (right hand side of system for calculation of
        ! shear viscocity coefficient)

        b2(1:NUM_SP, 1) = (2./Kb/T)*x
        ! DO i=1,5
        !     b2(i,1)=2./kb/T*x(i)
        ! END DO
        ! End of vector b2 definition
            
            
        ! Definition of vector b3 (right hand side of system for calculation of
        ! bulk viscosity coefficients)
        
        ! cu=kb*ntot/rho*(3./2.+x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))
        ! cut=kb*ntot/rho*(x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))
        
        !cu = R/M*(3./2. + SUM(x(:NUM_MOL))) ! kb*ntot/rho
        !cut = R/M*SUM(x(:NUM_MOL)) ! kb*ntot/rho

        b3(1:NUM_SP, 1) = -x*cv%cv_int_sp/cv%cv_tot
        ! DO i=1,5
        !     b3(i,1)=-x(i)*cv%cv_int_sp(i)/cv%cv_tot
        !     !b3(i,1)=-x(i)*cut/cu
        ! END DO

        b3(NUM_SP+1 : NUM_SP+NUM_MOL, 1) = y(:NUM_MOL)*cv%cv_int_sp(:NUM_MOL)/cv%cv_tot
        ! DO i=6,8
        !     b3(i,1)=y(i-5)*cv%cv_int_sp(i-5)/cv%cv_tot
        ! END DO
        
        !b3(6,1)=x(1)*ntot/rho*kb/cu
        
        !b3(7,1)=x(2)*ntot/rho*kb/cu !*(1+c_v_O2)
        
        !b3(8,1)=x(3)*ntot/rho*kb/cu !*(1+c_v_co)
        
        b3(1,1)=0

        ! write (*, *) 'RHS of system for bulk visc coeffs'
        ! write (*, *) (b3(i,1), i=1,NUM_SP + NUM_MOL)
        
        ! End of vector b3 definition
        

        ! Linear system solution using the Gauss method
        ! The solutions a, d, h, f are written to b, b1, b2, b3, respectively 
        
        call QRDecomposition(LTH,Q_LTH,R_LTH,2*NUM_SP)
        call QRDecomposition(LDIFF,Q_LDIFF,R_LDIFF,NUM_SP)
        call QRDecomposition(HVISC,Q_HVISC,R_HVISC,NUM_SP)
        ! call QRDecomposition(BVISC,Q_BVISC,R_BVISC,NUM_SP+NUM_MOL)

        call InvertQR(Q_LTH,R_LTH,Inverse_LTH,2*NUM_SP)
        call InvertQR(Q_LDIFF,R_LDIFF,Inverse_LDIFF,NUM_SP)
        call InvertQR(Q_HVISC,R_HVISC,Inverse_HVISC,NUM_SP)
        ! call InvertQR(Q_BVISC,R_BVISC,Inverse_BVISC,NUM_SP+NUM_MOL)

        b = matmul(Inverse_LTH,b)
        B1 = matmul(Inverse_LDIFF,B1)
        b2 = matmul(Inverse_HVISC,b2)
        ! b3 = matmul(Inverse_BVISC,b3)
        ! print *, 'bulk visc coeffs ', LTH
        ! call gaussj(LTH,10,10,b,1,1)
        ! call gaussj(Ldiff,5,5,b1,5,5)
        ! call gaussj(HVISC,5,5,b2,1,1)
        call gaussj(BVISC,8,8,b3,1,1)


        ! Thermal diffusion coefficients THDIF(i)

        thdiff = -(1./2./ntot)*b(:NUM_SP,1)
        ! DO i=1,5
        !     thdiff(i)=-1./2./ntot*b(i,1)
        ! END DO


        ! Thermal conductivity coefficient associated to translational
        ! energy, LTR 

        ltr = (5./4.)*kb*sum(x*b(NUM_SP+1 : 2*NUM_SP,1))
        ! DO i=6,10
        !     LTR=LTR+5./4.*kb*x(i-5)*b(i,1)
        ! END DO


        ! Thermal conductivity coefficients associated to internal
        ! energies 
        
        lint = (3./16.)*T*sum(x(:NUM_MOL)*cv%cv_int_sp(:NUM_MOL)*((Kb)*(MASS_SPCS(:NUM_MOL)*1e30) &
                /(bracket_out%lambda_int(:NUM_MOL)*1e30)))

        ! DO i = 1,3
        !     lint = lint + 3./16.*kb*T*x(i)/bracket_out%lambda_int(i)*cv%cv_int_sp(i)*MASS_SPCS(i) ! division by zero, bracket_out%lambda_int(i)
        !     WRITE (*,*) 'bracket_out%lambda_int(i) = ', bracket_out%lambda_int(i)
        ! END DO

        ! Total thermal conductivity coefficient at the translational
        ! temperature gradient

        ltot = ltr + lint


        ! Diffusion coefficients	DIFF(i,j)

        diff = (1./2./ntot)*b1
        ! DO i=1,5
        !     DO j=1,5
        !         diff(i,j)=1./2./ntot*b1(i,j)
        !     END DO
        ! END DO


        ! Shear viscosity coefficient VISC

        visc = (Kb*T/2.)*sum(x*b2(:NUM_SP,1))
        ! DO i=1,5
        !     visc=visc+kb*t/2.*b2(i,1)*x(i)
        ! END DO

        ! Bulk viscosity coefficient BULK_VISC

        bulk_visc = -Kb*T*sum(x*b3(:NUM_SP,1))
        ! DO i=1,5
        !     bulk_visc=bulk_visc-kb*t*b3(i,1)*x(i)
        ! END DO

        data_out%visc = visc
        data_out%bulk_visc = bulk_visc
        data_out%ltot = ltot
        data_out%thdiff = thdiff
        data_out%diff = diff
    
    end subroutine Transport1T
    

end module transport_1t
