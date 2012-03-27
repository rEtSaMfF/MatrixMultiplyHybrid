all: main.c print.c
	mpicc -Wall -lpthread main.c
	mpicc -Wall -o print print.c
