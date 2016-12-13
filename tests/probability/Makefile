install:
	mkdir -p $(LIB)
	mv libprobability.a $(LIB)

lib: mvndstpack.o ccpdata.o unbiased_correlation.o
	ar cr libprobability.a $^

hslib: mvndstpack.o hsccpdata.o hsfuncs.o
	ar cr libprobability.a $^

derivatives_cost: mvndstpack.f ccpdata.f90 derivatives_cost.f90
	gfortran  $^ -o $@

mvndstpack.o: mvndstpack.f
	gfortran -c $< -o $@

hsfuncs.o: hsfuncs.f
	gfortran -c $< -o $@

%.o: %.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *.o *.mod *.a *~ *.out derivatives_cost

.PHONY: lib clean
