install:
	mkdir -p $(LIB)
	mv libdfoirfilter.a $(LIB)

all: dfoirfilter.o restoration.o dfoirfilter.mod restoration.mod
	ar crv libdfoirfilter.a $^

dfoirfilter.o: dfoirfilter.f90 restoration.o
	$(FC) -c $< -o $@

restoration.o: restoration.f90
	$(FC) -c $^ -o $@

%.mod: %.o

clean:
	rm -vf *~ *.o *.mod

.PHONY: all clean install
