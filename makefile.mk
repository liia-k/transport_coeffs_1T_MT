CFLAGS = -ffree-form -ffree-line-length-0 -fdefault-double-8 -cpp -fbounds-check -finit-real=snan -finit-integer=-ftrapv -O -Wuninitialized
OMPFLAGS = -fopenmp
COMPILER = gfortran

data_objects = data/constant_air5.o data/defs_models.o
solvers_objects = solvers/qr.o
calc_objects = calc/specific_heat_sp.o calc/omega_integrals.o calc/bracket_integrals.o calc/transport_1t.o
calc_objects_new = calc/specific_heat_sp.o calc/omega_integrals.o calc/bracket_integrals_new.o calc/transport_1t_new.o
test_objects = tests/test_simplified_procedure.o tests/test_free_stream.o tests/test_simplified_random.o tests/test_omp_simplified.o tests/test_potential_models.o tests/test_potentials_time.o

all: test_simplified_procedure test_free_stream test_simplified_random test_omp_simplified test_potential_models test_potentials_time

test_simplified_procedure: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_simplified_procedure.o
	$(COMPILER) $(CFLAGS) -o tests/test_simplified_procedure $(solvers_objects) $(data_objects) $(calc_objects) tests/test_simplified_procedure.o

test_free_stream: $(solvers_objects) $(data_objects) $(calc_objects_new) tests/test_free_stream.o
	$(COMPILER) $(CFLAGS) -o tests/test_free_stream $(solvers_objects) $(data_objects) $(calc_objects_new) tests/test_free_stream.o

test_simplified_random: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_simplified_random.o
	$(COMPILER) $(CFLAGS) -o tests/test_simplified_random $(solvers_objects) $(data_objects) $(calc_objects) tests/test_simplified_random.o

test_omp_simplified: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_omp_simplified.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -c tests/test_omp_simplified.f90 -o tests/test_omp_simplified.o
	$(COMPILER) $(CFLAGS) $(OMPFLAGS) -o tests/test_omp_simplified $(solvers_objects) $(data_objects) $(calc_objects) tests/test_omp_simplified.o

test_potential_models: $(solvers_objects) $(data_objects) $(calc_objects_new) tests/test_potential_models.o
	$(COMPILER) $(CFLAGS) -o tests/test_potential_models $(solvers_objects) $(data_objects) $(calc_objects_new) tests/test_potential_models.o

test_potentials_time: $(solvers_objects) $(data_objects) $(calc_objects) tests/test_potentials_time.o
	$(COMPILER) $(CFLAGS) -o tests/test_potentials_time $(solvers_objects) $(data_objects) $(calc_objects) tests/test_potentials_time.o

data/%.o: data/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@

solvers/%.o: solvers/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@

calc/%.o: calc/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@ -I data -I solvers

tests/%.o: tests/%.f90
	$(COMPILER) $(CFLAGS) -c $< -o $@ -I data -I solvers -I calc

clean: 
	rm -f data/*.o data/*.mod solvers/*.o solvers/*.mod calc/*.o calc/*.mod tests/*.o tests/*.mod 
	rm *.o *.mod