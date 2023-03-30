# build obj files
cd data/
gfortran -c constant_air5.f90
gfortran -c defs_models.f90
cd ../solvers/
gfortran -c qr.f90
cd ../calc/
gfortran -c specific_heat_sp.f90 -I ../data
gfortran -c omega_integrals.f90 -I ../data
gfortran -c bracket_integrals.f90 -I ../data
gfortran -c transport_1t_simpl.f90 -I ../data -I ../solvers
# build particular tests
cd ../tests/
# 1
gfortran -c test_build.f90 -I ../data -I ../solvers -I ../calc
gfortran -o test test_build.o ../solvers/qr.o ../data/constant_air5.o ../data/defs_models.o ../calc/specific_heat_sp.o ../calc/omega_integrals.o ../calc/bracket_integrals.o ../calc/transport_1t_simpl.o -ffree-form -ffree-line-length-0 -fdefault-double-8  -cpp
# 2
gfortran -c test_f-s.f90 -I ../data -I ../solvers -I ../calc
gfortran -o test_f-s test_f-s.o ../solvers/qr.o ../data/constant_air5.o ../data/defs_models.o ../calc/specific_heat_sp.o ../calc/omega_integrals.o ../calc/bracket_integrals.o ../calc/transport_1t_simpl.o -ffree-form -ffree-line-length-0 -fdefault-double-8  -cpp
# 3
gfortran -c test_main.f90 -I ../data -I ../solvers -I ../calc
gfortran -o test_main test_main.o ../solvers/qr.o ../data/constant_air5.o ../data/defs_models.o ../calc/specific_heat_sp.o ../calc/omega_integrals.o ../calc/bracket_integrals.o ../calc/transport_1t_simpl.o -ffree-form -ffree-line-length-0 -fdefault-double-8  -cpp
# 4
gfortran -c test_omp.f90 -fopenmp -I ../data -I ../solvers -I ../calc
gfortran -o test_omp -fopenmp test_omp.o ../solvers/qr.o ../data/constant_air5.o ../data/defs_models.o ../calc/specific_heat_sp.o ../calc/omega_integrals.o ../calc/bracket_integrals.o ../calc/transport_1t_simpl.o -ffree-form -ffree-line-length-0 -fdefault-double-8  -cpp

