#include "corcorr.c"

int main () {
	// declaracoes
	double *vars;
	double *cor;
	int N=25, M=4, l=3, FLAGSTOP=0;
	int i,j;
	// inicializacoes
	srand(42);
	vars = (double *) calloc(N*M, sizeof(double));
	cor = (double *) calloc(M*M, sizeof(double));
	for (i =0; i<M;i++)
		for(j=0;j<N;j++)
			vars[i+N*j] = rand();

	cor[0+0*M] = 1.0; cor [1+0*M]=0.1; cor[2+0*M] = 0.2; cor[3+0*M] = 0.3;
	cor[0+1*M] = 0.1; cor [1+1*M]=1.0; cor[2+1*M] = 0.4; cor[3+1*M] = 0.5;
	cor[0+2*M] = 0.2; cor [1+2*M]=0.4; cor[2+2*M] = 1.0; cor[3+2*M] = 0.6;
	cor[0+3*M] = 0.3; cor [1+3*M]=0.5; cor[2+3*M] = 0.6; cor[3+3*M] = 1.0;

	corcorr(vars,cor,&N,&M,&l,&FLAGSTOP);

	return 0;
}
