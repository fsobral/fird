base: userinterface.mod

install:
	mkdir -p $(LIB)
	mv libdfoirfilter.a $(LIB)

all: dfoirfilter.o restoration.o fkss.o userinterface.o
	ar crv libdfoirfilter.a $^

fkss.o: fkss.f90 dfoirfilter.o
	$(FC) -c $< -o $@

dfoirfilter.o: dfoirfilter.f90 restoration.o userinterface.o
	$(FC) -c $< -o $@

restoration.o: restoration.f90
	$(FC) -c $^ -o $@

userinterface.mod: userinterface.o

%.o: %.f90
	$(FC) -c -J$(LIB) $^ -o $@

clean:
	rm -vf *~ *.o *.mod

.PHONY: all clean install
