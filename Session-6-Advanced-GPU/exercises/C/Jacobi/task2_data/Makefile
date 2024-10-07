all: jacobi.gpu jacobi.sol.gpu

jacobi.gpu: jacobi.c
	cc -fast -mp=gpu -Minfo=mp,accel -o jacobi.gpu jacobi.c

jacobi.sol.gpu: Makefile jacobi.sol.c
	cc -fast -mp=gpu -Minfo=mp,accel -o jacobi.sol.gpu jacobi.sol.c

clean:
	rm -f *.o
	rm -f jacobi.gpu
	rm -f jacobi.sol.gpu

realclean: clean
	rm -f *.bak *~

.PHONY: all clean realclean

