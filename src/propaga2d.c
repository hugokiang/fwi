#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "string.h"
#include "utils.h"


int main(int argc, char **argv)
{

	Config config;
	char arqparam[2000];
	strcpy(arqparam, argv[1]);
	FILE *fp;
	fp = fopen(arqparam, "r");
	if(!fp) {
		perror(arqparam);
		exit(1);
	}	

	leParametros(&config, fp);

	int *ir2;
	int *ir1;
	ir2 = alocaMatrizInt(10000,1);
	ir1 = alocaMatrizInt(10000,1);
	leAquisicao(&config, ir2, ir1);
	if(config.nrec>10000) {
		printf("Problema com o numero de receptores. nrec=%d o maximo eh 10000\n", config.nrec);
	}

	fprintf(stderr,"nx = %d nz %d ns %d\n", config.nx, config.nz,  config.length_conv_operator); 

	double **modelo;
	modelo = double_matrix(config.nx, config.nz);
	read_model2(config.name_vp, config, modelo);
	



}