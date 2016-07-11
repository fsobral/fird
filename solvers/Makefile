install: $(RESTORATION_INTERFACE).o $(OPTIMIZATION_INTERFACE).o

restoration.o: restoration.f90 algencan_restoration.o
	$(FC) -c $(FCC) -J$(LIB) $< -o $@

trdf_solver.o: trdf_solver.f90 trdf_algencan_solver.o
	$(FC) -c $(FCC) -J$(LIB) $< -o $@

optimization.o: optimization.f90 trdf_solver.o
	$(FC) -c $(FCC) -J$(LIB) $< -o $@

%.o: %.f90
	$(FC) -c $(FCC) -J$(LIB) $< -o $@

clean:
	rm -vf *.o *.mod *~
