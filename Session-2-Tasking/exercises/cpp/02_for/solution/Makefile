.PHONY: clean run go debug release

PROG=for-taskloop
CC = cc
FLAGS_OPENMP = -fopenmp

release:
	$(CC) -Wall -g $(FLAGS_OPENMP) -O3 ${PROG}.c -o ${PROG}.exe -lm

debug:
	$(CC) -Wall -g $(FLAGS_OPENMP) -O0 ${PROG}.c -o ${PROG}.exe -lm

run go:
	./${PROG}.exe

clean:
	rm -f ${PROG}.exe ${PROG}.o *~ *.er gmon.out gprof.txt
