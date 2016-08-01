base: userinterface.mod rinterface.mod ointerface.mod filters.mod

install:
	mkdir -p $(LIB)
	mv libdfoirfilter.a $(LIB)

all: dfoirfilter.o fird.o filters.o
	ar crv libdfoirfilter.a $^

fird.o: fird.f90
	$(FC) -c -J$(LIB) $< -o $@

dfoirfilter.o: dfoirfilter.f90
	$(FC) -c -J$(LIB) $< -o $@

userinterface.mod: userinterface.o

rinterface.mod: rinterface.o

ointerface.mod: filters.mod ointerface.o 

filters.mod: filters.o

%.o: %.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *~ *.o *.mod

.PHONY: all base clean install
