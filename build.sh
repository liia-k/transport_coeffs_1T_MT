# build obj files
cd data/
gfortran -c constant_air5.f90
gfortran -c defs_models.f90
cd ../solvers/
gfortran -c invers.f90
gfortran -c qr.f90
cd ../calc/
gfortran -c specific_heat_sp.f90 -I ../data
gfortran -c omega_integrals.f90 -I ../data
gfortran -c bracket_integrals.f90 -I ../data
gfortran -c transport_1T.f90 -I ../data -I ../solvers
# build particular test
cd ../tests/
gfortran -c test_current.f90 -I ../data -I ../solvers -I ../calc
gfortran -o test test_current.o ../solvers/invers.o ../solvers/qr.o ../data/constant_air5.o ../data/defs_models.o ../calc/specific_heat_sp.o ../calc/omega_integrals.o ../calc/bracket_integrals.o ../calc/transport_1T.o

rm **/*.mod **/*.o

