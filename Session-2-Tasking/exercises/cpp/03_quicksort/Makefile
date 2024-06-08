.PHONY: clean run go debug release

PROG=quicksort
CXX = CC
FLAGS_OPENMP = -fopenmp

release:
	$(CXX) -Wall -g $(FLAGS_OPENMP) -O3 ${PROG}.cpp -o ${PROG}.exe

debug:
	$(CXX) -Wall -g $(FLAGS_OPENMP) -O0 ${PROG}.cpp -o ${PROG}.exe

run go:
	./${PROG}.exe

clean:
	rm -f ${PROG}.exe ${PROG}.o *~
