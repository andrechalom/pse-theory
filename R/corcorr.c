#include<stdlib.h>
#include<float.h>
#include<math.h>
#define thisV vars[i+(*N)*j]
#define thisR R[i+(*N)*j]

void sw (double *obj, int i, int j, int k) {
	double tmp;
	tmp = obj[i+ k];
	obj[i + k] = obj[j+ k];
	obj[j + k] = tmp;
}
/* Funcao procura a menor correlacao entre as variaveis dadas apos um switch */
void corcorr (double *vars, double *cor, int * N, int * M, int * l) {
	/* R: Matriz temporaria para guardar os diferentes cenarios. Tj: vetor temporario */
	double *R;
	R = (double *) calloc((*N)*(*M), sizeof (double));
	/* Normaliza as variaveis */
	int i, j, k, m;
	double mean, sd, sum=0, sq_sum=0;
	for (j=0;j< *l;j++) {
		sum=0; sq_sum=0;
		for(i = 0; i < *N; ++i) {
			sum += thisV;
			sq_sum += thisV * thisV;
		}
		mean = sum / (*N);
		sd = sqrt(sq_sum / *N - mean * mean);
		for (i = 0; i< *N; i++)
			thisR = (thisV-mean)/sd;
	}
	/* min* usados para guardar o minimo ateh agora */
	double minE = DBL_MAX, E, Tj, tmp;
	int mini=0, minj=3;
	for (i =0; i < *N-1; i++) 
		for (j = i+1; j < *N; j++) {
			E =0;
			// Troca o valor de R[i] e R[j]
			sw(R, i, j, (*l-1)*(*N));
			for (m=0; m < *l-1;m++) {
				Tj = 0;
				for (k=0;k<*N;k++)
					Tj += R[k + (*N) * (*l-1)] * R[k+(*N)*m];
			E += (Tj - cor[*l + m*(*M)])*(Tj - cor[*l + m*(*M)])/(*N)/(*N);
			}
			// trabalha com E aqui
			if (E < minE) {
				mini = i;
				minj = j;
				minE = E;
			}
			// Finalmente, destroca i e j
			sw(R, i, j, (*l-1)*(*N));

	}
	/* Troca o valor da variavel na posicao i e j e finaliza */
	sw(vars, mini, minj, (*l-1)*(*N));
	free(R);
	return;
}
