# Gera o teste de execucao da funcao corcorr.c

all: debug

debug:
	gcc teste_corcorr.c -g -ggdb -lm
	gdb ./a.out

clean:
	-rm -rf a.out corcorr.o corcorr.so

corcorr.o: corcorr.c
	gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -O3 -pipe  -g  -c corcorr.c -o corcorr.o

corcorr.so: corcorr.o
	gcc -std=gnu99 -shared -o corcorr.so corcorr.o -L/usr/lib/R/lib -lR
