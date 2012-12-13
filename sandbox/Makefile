all: findf

a.out: test.c
	gcc -Iinclude -Llib -Wall -g test.c -lgsl -lm -lgslcblas

findf: findf.c
	gcc -Wall -g findf.c -lgsl -lgslcblas -o findf
