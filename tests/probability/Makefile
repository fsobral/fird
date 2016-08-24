install:
	mkdir -p $(LIB)
	mv libprobability.a $(LIB)

lib: mvndstpack.o
	ar cr libprobability.a $^

%.o: %.f
	gfortran -c $< -o $@

clean:
	rm -vf *.o *.mod *.a *~

.PHONY: lib clean
