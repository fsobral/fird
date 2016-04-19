lib: dfoirfilter.o dfoirfilter.mod restoration.mod restoration.o

%.o: %.f90
	$(FC) -c $^ -o $@

clean:
	rm -vf *~ *.o *.mod

.PHONY: clean
