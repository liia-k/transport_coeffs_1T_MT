!Module for calculation of transport coefficients. Uses the module constant.f90 containing main constants and 
!variables definition, modules for calculation of specific heats, omega and bracket integrals along with 

MODULE TRANSPORT_AIR5_1T

    USE CONSTANT
    USE SPECIFIC_HEAT
    USE BRACKET_INTEGRALS
    
    IMPLICIT NONE

    
    CONTAINS

    SUBROUTINE TRANSPORT_1T
    
        INTEGER I, J, K, DELTA

        REAL CU, CUT

        !Matrices for the linear transport systems defining
        !heat conductivity and thermal diffusion (LTH);
        !bulk viscosity (BVISC);
        !diffusion (LDIFF);
        !shear viscisity (HVISC).

        REAL, DIMENSION(10,10) :: LTH 

        REAL, DIMENSION(8,8) :: BVISC 

        REAL, DIMENSION(5,5) ::  LDIFF, HVISC, b1

        !Vectors of right hand terms

        REAL, DIMENSION(10,1) :: b

        REAL, DIMENSION(5,1) :: b2

        REAL, DIMENSION(8,1) :: b3

    ! Definition of matrix LTH for calculation of 
    ! thermal conductivity and thermal diffuaion coefficients
    ! The system has a form	(1)
    ! LTH times a = b, a is the vector of unknowns
        DO i=1,5
            DO j=1,5
                LTH(i,j)=Lambda00(i,j)
            END DO
        END DO
        
        DO i=1,5
            DO j=6,10
                LTH(i,j)=Lambda01(i,j-5)
            END DO
        END DO
        
        DO i=6,10
            DO j=1,5
                LTH(i,j)=LTH(j,i)
            END DO
        END DO
        
        DO i=6,10
            DO j=6,10
                LTH(i,j)=Lambda11(i-5,j-5)
            END DO
        END DO
        
        DO j=1,5
            LTH(1,j)=y(j)!x(j)*mass(j)*ntot/rho
        END DO
        
        DO j=6,10
            LTH(1,j)=0.
        END DO
    ! End of matrix LTH definition


    ! Definition of matrix LDIFF for calculation of 
    ! diffuaion coefficients
    ! The system has a form	(2)
    ! LDIFF times D = B1, D a is the matrix of unknowns
        
        DO i=1,5
            DO j=1,5
            LDIFF(i,j)=LTH(i,j)
            END DO
        END DO
    ! End of matrix LDIFF definition


    ! Definition of matrix HVISC for calculation of 
    ! shear viscocity coefficient
    ! The system has a form	(3)
    ! HVISC times h = b2, h a is the vector of unknowns
        
        DO i=1,5
            DO j=1,5
            HVISC(i,j)=H00(i,j)
            END DO
        END DO
    ! End of matrix HVISC definition
        
        
    ! Definition of matrix BVISC for calculation of 
    ! bulk viscosity coefficients
    ! The system has a form	(4)
    ! BVISC times f = b3, f is the vector of unknowns
        DO i=1,5
            DO j=1,5
                BVISC(i,j)=x(j)*beta11(i,j)
            END DO
        END DO
        
        DO i=1,5
            DO j=6,8
                BVISC(i,j)=y(j-5)*beta01(i,j-5)
            END DO
        END DO
        
        DO i=6,8
            DO j=1,5
                BVISC(i,j)=BVISC(j,i)!x(j)*beta11(j,i)
            END DO
        END DO
        
        DO i=6,8
            DO j=6,8
                BVISC(i,j)=0
            END DO
        END DO
        
        BVISC(6,6)=y(1)*beta0011(1)
        BVISC(7,7)=y(2)*beta0011(2)
        BVISC(8,8)=y(3)*beta0011(3)
        
        
        DO j=1,5
            BVISC(1,j)=x(j)*3./2.*kb*ntot/rho
        END DO
        
        DO j=6,6
            BVISC(1,j)=y(1)*kb/mass(1)
        END DO
        
        DO j=7,7
            BVISC(1,j)=y(2)*kb/mass(2)
        END DO
        
        DO j=8,8
            BVISC(1,j)=y(3)*kb/mass(3)
        END DO
        
    ! End of matrix BVISC definition
        
        
    ! Definition of vector b (right hand side of system (1))
        DO i=1,5
            b(i,1)=0.
        END DO
        DO i=6,10
            b(i,1)=4./5./kb*x(i-5)!15./2.*T/kb*x(i-5)
        END DO
    ! End of vector b definition
        
        
    ! Definition of matrix b1 (right hand side of system (2))
        DO i=1,5
            DO j=1,5
            if(i==j) then 
                delta=1
            else
                delta=0
            end if
            B1(i,j)=8./25./kb*(delta - y(i))!3*kb*T*(delta-y(i));
            END DO
        END DO
        
        DO j=1,5
            B1(1,j)=0
        END DO
    ! End of matrix b1 definition
        
    ! Definition of vector b2 (right hand side of system (3))
        DO i=1,5
            b2(i,1)=2./kb/T*x(i)
        END DO
    ! End of vector b2 definition
        
        
    ! Definition of vector b3 (right hand side of system (4))
        
        !cu=kb*ntot/rho*(3./2.+x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))
        !cut=kb*ntot/rho*(x(1)+x(2)*(1+c_v_o2)+x(3)*(1+c_v_co))
        
        cu=kb*ntot/rho*(3./2.+x(1)+x(2)+x(3))
        cut=kb*ntot/rho*(x(1)+x(2)+x(3))
        
        DO i=1,5
            b3(i,1)=-x(i)*c_int(i)/cv_tot
            !b3(i,1)=-x(i)*cut/cu
        END DO

        DO i=6,8
            b3(i,1)=y(i-5)*c_int(i-5)/cv_tot
        END DO
        
        !b3(6,1)=x(1)*ntot/rho*kb/cu
        
        !b3(7,1)=x(2)*ntot/rho*kb/cu !*(1+c_v_O2)
        
        !b3(8,1)=x(3)*ntot/rho*kb/cu !*(1+c_v_co)
        
        b3(1,1)=0
        
    ! End of vector b3 definition
        

        ! Linear system solution using the Gauss method
    ! The solutions a, d, h, f are written to b, b1, b2, b3, respectively 
        
    CALL gaussj(LTH,10,10,b,1,1)
    CALL gaussj(Ldiff,5,5,b1,5,5)
    CALL gaussj(HVISC,5,5,b2,1,1)
    CALL gaussj(BVISC,8,8,b3,1,1)


! Thermal diffusion coefficients THDIF(i)

    DO i=1,5
        thdiff(i)=-1./2./ntot*b(i,1)
    END DO


! Thermal conductivity coefficient associated to translational
! energy, LTR 

    Ltr=0
    DO i=6,10
        LTR=LTR+5./4.*kb*x(i-5)*b(i,1)
    END DO


! Thermal conductivity coefficients associated to internal
! energies 
    Lint = 0
    DO i = 1,3
        lint = lint + 3./16.*kb*T*x(i)/lambda_int(i)*c_int(i)*mass(i)
    END DO
    !lrot_n2=3.*kb*t*x(1)/16./lambda_int(1)*kb
    !lrot_o2=3.*kb*t*x(2)/16./lambda_int(2)*kb
    !lrot_no=3.*kb*t*x(3)/16./lambda_int(3)*kb
    !lvibr_n2=3.*kb*t*x(1)/16./lambda_int(2)*kb*(c_vibr_n2)
    !lvibr_o2=3.*kb*t*x(2)/16./lambda_int(2)*kb*(c_vibr_o2)
    !lvibr_no=3.*kb*t*x(3)/16./lambda_int(3)*kb*(c_vibr_no)

! Total thermal conductivity coefficient at the translational
! temperature gradient

    ltot = ltr + lint


!Diffusion coefficients	DIFF(i,j)

    DO i=1,5
        DO j=1,5
            diff(i,j)=1./2./ntot*b1(i,j)
        END DO
    END DO


!Shear viscosity coefficient VISC

    visc=0
    DO i=1,5
        visc=visc+kb*t/2.*b2(i,1)*x(i)
    END DO

!Bulk viscosity coefficient BULK_VISC

    bulk_visc=0
    DO i=1,5
        bulk_visc=bulk_visc-kb*t*b3(i,1)*x(i)
    END DO

    END SUBROUTINE TRANSPORT_1T
    

END MODULE TRANSPORT_AIR5_1T
