!In this module, vibrational energy, non-equilibrium vibrational
!partition functions, vibrational specific heat capacities,
! and mean vibrational energy en_int of molecular species are calculated.
!It uses the module constant.f90 containing main constants and 
!variables definition 
!Input variables: T, T12, T3

MODULE Specific_heat

USE CONSTANT

IMPLICIT NONE

Double precision zv_12,zv_3,zv_o2,zv_co,c_v_t12,c_v_t3,c_v_o2,c_v_co

real, dimension(5) :: en_int 

REAL ee1,ee2,ee3,ee4

INTEGER I1,I2,I3,IL,G

CONTAINS

  SUBROUTINE ENERGY

  do i1=1,5
     en_int(i1)=0
  end do

! Calculation of vibrational energy levels for CO2, O2, CO
  
    DO i1=0,L1
      DO i2=0,L2
	    DO i3=0,L3
          en3(i1,i2,i3)=w1*i1+w2*i2+w3*i3
		  !+wx11*i1*i1+wx12*i1*i2+wx13*i1*i3+&
          !wx22*i2*i2+wx23*i2*i3+wx33*i3*i3+wxll*il*il
        END DO
      END DO
    END DO

    DO i1=0,L_O2
       en_O2(i1)=we_O2*(i1)
	   !-wexe_O2*i1**2
    END DO

    DO i1=0,L_CO
       en_CO(i1)=we_CO*(i1)
	   !-wexe_CO*i1**2
    END DO

  END SUBROUTINE ENERGY


  SUBROUTINE PART_FUNC3(T12,T3)

	! Calculation of non-equilibrium CO2 partition function
	! Z_CO2(T12,T3)

	REAL EE, T12, T3

    zv_12=0 
    zv_3=0 

    DO i1=0,L1
      DO i2=0,L2
			    g=i2+1
                ee=en3(i1,i2,0)
				if(ee<d_co2*kb) then
				zv_12=zv_12+&
					g*exp(-(i1*ee1+i2*ee2)/kb/t12)
				end if            
      end do
    end do

        DO i3=0,L3
                ee=en3(0,0,i3)
				if(ee<d_co2*kb) then
				zv_3=zv_3+&
					exp(-i3*ee3/kb/t3)
				end if                
        end do

  END SUBROUTINE PART_FUNC3

  SUBROUTINE PART_FUNC_O2(TVO2)
	
	! Calculation of non-equilibrium O2 partition function Z_O2(TVO2)

    INTEGER I1
	  REAL EE,TVO2
	
	  zv_O2=0 

    DO i1=0,L_O2
		ee=en_O2(i1)
		zv_O2=zv_O2+exp(-ee/kb/TVO2)
    end do

  END SUBROUTINE PART_FUNC_O2

  SUBROUTINE PART_FUNC_CO(TVCO)
	
	! Calculation of non-equilibrium CO partition function Z_CO(TVCO)

    INTEGER I1
	REAL EE, TVCO

	zv_CO=0 
    
    DO i1=0,L_CO
		ee=en_CO(i1)
		zv_CO=zv_CO+exp(-ee/kb/TVCO)
    end do

  END SUBROUTINE PART_FUNC_CO

  SUBROUTINE s_heat3

	! Calculation of non-equilibrium CO2 vibrational specific heats
	! CVT12, CVT3

	integer i1,i2,i3
	real ee
	double precision ppp,s12,s3,ss12,ss3

	s12=0;s3=0;
    ss12=0;ss3=0
		
	DO i1=0,L1
      DO i2=0,L2
			g=i2+1
            ee=en3(i1,i2,0)
			  if(ee<d_co2*kb) then
		    	ppp=g*exp(-(i1*ee1+i2*ee2)/kb/t12)/zv_12

                s12=s12+(i1*ee1+i2*ee2)/kb/t12*ppp
                ss12=ss12+((i1*ee1+i2*ee2)/kb/t12)**2*ppp
				en_int(1)=en_int(1)+ee*ppp
	    	  end if
      end do
    end do

        DO i3=0,L3
            ee=en3(0,0,i3)
			  if(ee<d_co2*kb) then
		    	ppp=exp(-i3*ee3/kb/t3)/zv_3

                s3=s3+i3*ee3/kb/t3*ppp
                ss3=ss3+i3*ee3/kb/t3*i3*ee3/kb/t3*ppp
				en_int(1)=en_int(1)+ee*ppp
	    	  end if
        end do


    c_v_t12=(ss12-s12*s12)
    c_v_t3=(ss3-s3*s3)


  END SUBROUTINE s_heat3

  SUBROUTINE s_heat_O2

	! Calculation of non-equilibrium O2 vibrational specific heat
	! CV_O2

	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_O2
       ee=en_O2(i1)
       ppp=exp(-ee/kb/TVO2)/zv_O2
       s=s+ee/kb/TVO2*ppp;
       s0=s0+ee*ee/kb/TVO2/kb/TVO2*ppp;
	   en_int(2)=en_int(2)+ee*ppp
    END DO
    c_v_o2=(s0-s*s)
	
  END SUBROUTINE s_heat_O2

  SUBROUTINE s_heat_CO
	! Calculation of non-equilibrium CO vibrational specific heat
	! CV_CO
	
	integer I1
	real ee
	double precision ppp,s,s0

	s=0;s0=0;

    DO i1=0,L_CO
       ee=en_CO(i1)
       ppp=exp(-ee/kb/TVCO)/zv_CO
       s=s+ee/kb/TVCO*ppp
       s0=s0+ee*ee/kb/TVCO/kb/TVCO*ppp
       en_int(3)=en_int(3)+ee*ppp
    END DO
    c_v_co=(s0-s*s)
	
  END SUBROUTINE s_heat_CO


END	MODULE Specific_heat


