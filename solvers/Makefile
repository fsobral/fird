install: $(RESTORATION_INTERFACE).o $(OPTIMIZATION_INTERFACE).o

restoration.o: restoration.f90 algencan_restoration.o
	$(FC) -c $(FCC) -I$(MODULES) -J$(LIB) $< -o $@

algencan_restoration.o: algencan_restoration.f90
	$(FC) -c $(FCC) -I$(MODULES) -J$(LIB) $< -o $@

trdf_solver.o: trdf_solver.f90
	$(FC) -c $(FCC) -I$(MODULES) -J$(LIB) $< -o $@

trdf_solver_interface.o: trdf_solver_interface.f90 trdf_solver.o
	$(FC) -c $(FCC) -I$(MODULES) -J$(LIB) $< -o $@

clean:
	rm -vf *.o *.mod *~
