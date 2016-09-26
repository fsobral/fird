install:
	mkdir -p $(LIB)
	mv libprobability.a $(LIB)

lib: mvndstpack.o ccpdata.o unbiased_correlation.o
	ar cr libprobability.a $^

mvndstpack.o: mvndstpack.f
	gfortran -c $< -o $@

%.o: %.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *.o *.mod *.a *~ *.out

.PHONY: lib clean
