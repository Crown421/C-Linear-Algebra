#all: gaussj #gmres1 #
all: test_matrix
#all: LUsolveO3

#For debugging
OPT=-g -Wall --std=c++11
#For optimistaion
#OPT=-O3

#All objects (except main) come from cpp and hpp
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<

#matrix.o:   .FORCE

#.FORCE:

LUsolve:	LUsolve.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -c -o matrix.o matrix.cpp
		g++ ${OPT} -c -o LUd.o LUd.cpp
		g++ ${OPT} -o LUsolve LUsolve.cpp matrix.o Exception.o LUd.o

LUsolveSAFE:	LUsolve.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -D SAFETYOFF -c -o LUd.o LUd.cpp
		g++ ${OPT} -D SAFETYOFF -o LUsolveSAFE LUsolve.cpp matrix.o Exception.o LUd.o

LUsolveO3:	LUsolve.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -O3 -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -O3 -D SAFETYOFF -c -o LUd.o LUd.cpp
		g++ ${OPT} -O3 -D SAFETYOFF -o LUsolveO3 LUsolve.cpp matrix.o Exception.o LUd.o

gaussj:	gaussj.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -c -o matrix.o matrix.cpp
		g++ ${OPT} -o gaussj gaussj.cpp matrix.o Exception.o LUd.o

gaussjSAFE:	gaussj.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -D SAFETYOFF -o gaussjSAFE gaussj.cpp matrix.o Exception.o LUd.o

gaussjO3:	gaussj.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -O3 -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -O3 -D SAFETYOFF -o gaussjO3 gaussj.cpp matrix.o Exception.o LUd.o

gmres1:	gmres1.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -c -o matrix.o matrix.cpp
		g++ ${OPT} -o gmres1 gmres1.cpp matrix.o Exception.o LUd.o

gmres1SAFE:	gmres1.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -D SAFETYOFF -o gmres1SAFE gmres1.cpp matrix.o Exception.o LUd.o

gmres1O3:	gmres1.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -O3 -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -O3 -D SAFETYOFF -o gmres1O3 gmres1.cpp matrix.o Exception.o LUd.o

gmres2:	gmres2.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -c -o matrix.o matrix.cpp
		g++ ${OPT} -o gmres2 gmres2.cpp matrix.o Exception.o LUd.o

gmres2SAFE:	gmres2.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -D SAFETYOFF -o gmres2SAFE gmres2.cpp matrix.o Exception.o LUd.o

gmres2O3:	gmres2.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -O3 -D SAFETYOFF -c -o matrix.o matrix.cpp
		g++ ${OPT} -O3 -D SAFETYOFF -o gmres2O3 gmres2.cpp matrix.o Exception.o LUd.o


test_matrix:	test_matrix.cpp matrix.o Exception.o LUd.o
		g++ ${OPT} -c -o matrix.o matrix.cpp
		g++ ${OPT} -o test_matrix test_matrix.cpp matrix.o Exception.o LUd.o

run:
		gmres1 random
		#test_matrix

clean:
	rm -f *.o *~ test_matrix.exe gmres1.exe
