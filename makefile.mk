CFLAGS = -ffree-form -ffree-line-length-0 -fdefault-double-8 -cpp -fbounds-check -finit-real=snan -finit-integer=-ftrapv -O -Wuninitialized
OMPFLAGS = -fopenmp
COMPILER = gfortran

data_objects = data/constant_air5.o data/defs_models.o
solvers_objects = solvers/qr.o
calc_objects = calc/specific_heat_sp.o calc/omega_integrals.o calc/bracket_integrals.o calc/transport_1t_simpl.o
calc_objects_new = calc/specific_heat_sp.o calc/omega_integrals.o calc/bracket_integrals.o calc/transport_1t.o
test_objects = tests/test_build.o tests/test_f-s.o tests/test_main.o tests/test_omp.o tests/testModels.o tests/testModelsTime.o
test_objects_new = tests2/test_build.o tests2/test_f-s.o tests2/test_main.o tests2/test_omp.o tests2/testModels.o tests2/testModelsTime.o

all: test test_f-s test_main test_omp testModels testModelsTime test2 test_f-s2 test_main2 test_omp2 testModels2 testModelsTime2

test: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_build.o
	$(COMPILER) $(CFLAGS) -o tests/test $(solvers_objects) $(data_objects) $(calc_objects) tests/test_build.o

test2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_build.o
	$(COMPILER) $(CFLAGS) -o tests2/test $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_build.o

test_f-s: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_f-s.o
	$(COMPILER) $(CFLAGS) -o tests/test_f-s $(solvers_objects) $(data_objects) $(calc_objects) tests/test_f-s.o

test_f-s2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_f-s.o
	$(COMPILER) $(CFLAGS) -o tests2/test_f-s $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_f-s.o

test_main: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_main.o
	$(COMPILER) $(CFLAGS) -o tests/test_main $(solvers_objects) $(data_objects) $(calc_objects) tests/test_main.o

test_main2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_main.o
	$(COMPILER) $(CFLAGS) -o tests2/test_main $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_main.o

test_omp: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_omp.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -c tests/test_omp.f90 -o tests/test_omp.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -o tests/test_omp $(solvers_objects) $(data_objects) $(calc_objects) tests/test_omp.o

test_omp2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_omp.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -c tests2/test_omp.f90 -o tests2/test_omp.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -o tests2/test_omp $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/test_omp.o

testModels: $(solvers_objects) $(data_objects) $(calc_objects) tests/testModels.o
	$(COMPILER) $(CFLAGS) -o tests/testModels $(solvers_objects) $(data_objects) $(calc_objects) tests/testModels.o

testModels2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/testModels.o
	$(COMPILER) $(CFLAGS) -o tests2/testModels $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/testModels.o

testModelsTime: $(solvers_objects) $(data_objects) $(calc_objects) tests/testModelsTime.o
	$(COMPILER) $(CFLAGS) -o tests/testModelsTime $(solvers_objects) $(data_objects) $(calc_objects) tests/testModelsTime.o

testModelsTime2: $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/testModelsTime.o
	$(COMPILER) $(CFLAGS) -o tests2/testModelsTime $(solvers_objects) $(data_objects) $(calc_objects_new) tests2/testModelsTime.o

data/%.o: data/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@

solvers/%.o: solvers/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@

calc/%.o: calc/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@ -I data -I solvers

tests/%.o: tests/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@ -I data -I solvers -I calc

tests2/%.o: tests2/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@ -I data -I solvers -I calc

clean: 
	rm -f data/*.o data/*.mod solvers/*.o solvers/*.mod calc/*.o calc/*.mod tests/*.o tests/*.mod 
	rm *.o 
	rm *.mod