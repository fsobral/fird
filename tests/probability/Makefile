install:
	mkdir -p $(LIB)
	mv libprobability.a $(LIB)

lib: mvndstpack.o ccpdata.o
	ar cr libprobability.a $^

%.o: %.f
	gfortran -c $< -o $@

ccpdata.o: ccpdata.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *.o *.mod *.a *~

.PHONY: lib clean
