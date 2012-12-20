
itebd : module.o itebd.o
	gfortran -o itebd module.o itebd.o -llapack

module.o : module.f95
	gfortran -c module.f95

itebd.o : itebd.f95
	gfortran -c itebd.f95


.PHONY: clean

clean:
	rm *.o
	[ -f itebd ] && rm itebd
	[ -f mymodule.mod ] && rm mymodule.mod
