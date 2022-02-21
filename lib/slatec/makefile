.PHONY: all clean

all: libslatec.a

libslatec.a:
	cd slatec-bessel-cpp ; env CXX=g++ make libslatec-f2c.a ; cp libslatec-f2c.a ../libslatec.a

clean:
	rm -rf libslatec.a
	cd slatec-bessel-cpp ; make clean
