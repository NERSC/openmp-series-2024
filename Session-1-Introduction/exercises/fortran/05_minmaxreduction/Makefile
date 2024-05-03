.PHONY: clean run go debug release

PROG=MinMaxReduction
FC = ftn
FLAGS_OPENMP = -fopenmp

release:
	$(FC) -g $(FLAGS_OPENMP) -O3 ${PROG}.f90 -o ${PROG}.exe

debug:
	$(FC) -g $(FLAGS_OPENMP) -O0 ${PROG}.f90 -o ${PROG}.exe

run go:
	./${PROG}.exe

clean:
	rm -f ${PROG}.exe ${PROG}.mod *~
