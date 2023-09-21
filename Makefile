CFLAGS = -fdefault-real-8

CubicSolveForK: CubicSolveForK.f90 MatrixUtilities.o Solvers.o
	 gfortran  $(CFLAGS) -o CubicSolveForK CubicSolveForK.f90 Functions.o MatrixUtilities.o Solvers.o

Solvers.o: Solvers.f90
	 gfortran  $(CFLAGS) -c Solvers.f90

MatrixUtilities.o: MatrixUtilities.f90 Functions.o
	 gfortran  $(CFLAGS) -c MatrixUtilities.f90 Functions.o

Functions.o: Functions.f90
	 gfortran  $(CFLAGS) -c Functions.f90

.PHONY: test clean

test: CubicSolveForK
	./CubicSolveForK

clean:
	rm -f CubicSolveForK Functions.o functions.mod Solvers.o solvers.mod MatrixUtilities.o matrixutilities.mod
