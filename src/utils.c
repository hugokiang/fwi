#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "libfwi.h"
#include "math.h"


double **double_matrix(int nlinhas, int ncolunas) {
	double **matriz;
	int i;
	// aloca as linhas
	matriz = malloc(sizeof(double)*nlinhas);
	// sempre teste se a alocacao foi realizada com sucesso
	if( !matriz) {
		fprintf(stderr,"Erro ao alocar linhas da matriz\n");
		fprintf(stderr,"nl = %d nc = %d\n", nlinhas, ncolunas);
		exit(1);
	}
	for(i=0; i<nlinhas;i++) {
		matriz[i] =  malloc(sizeof(double)*ncolunas);
		// sempre teste se a alocacao foi realizada com sucesso
		if( !matriz[i]) {
			fprintf(stderr,"Erro ao alocar colunas da matriz\n");
			exit(1);
		}
	}
	return matriz;
}

void free_matrix(double **m, int nlinhas, int ncolunas)
{
	int i;
	for(i=0; i<nlinhas; i++) {
		free(m[i]);
	}
	free(m);
}

float *float_vector(int nlinhas) {
	float *matriz;
	int i;
	// aloca as linhas
	matriz = malloc(sizeof(float)*nlinhas);
	// sempre teste se a alocacao foi realizada com sucesso
	if( !matriz) {
		fprintf(stderr,"Erro ao alocar linhas da matriz\n");
		exit(1);
	}
	return matriz;
}

double *double_vector(int nlinhas) {
	double *matriz;
	int i;
	// aloca as linhas
	matriz = malloc(sizeof(double)*nlinhas);
	// sempre teste se a alocacao foi realizada com sucesso
	if( !matriz) {
		fprintf(stderr,"Erro ao alocar double_vector\n");
		fprintf(stderr,"memoria necessaria = %.1lf MB\n", sizeof(double)*nlinhas/1000./1000.);
		exit(1);
	}
	return matriz;
}


size_t remove_comments(char *buffer)
{
char *p = strchr(buffer, '#');
if (p != NULL) *p = '\0';
return strlen(buffer);
}


void leParametros(Config *c, FILE *fp) {
		//Model names
	int i,j;
	char *linha;
	size_t comp;
	linha = malloc(2000);
	char name_vp[200], name_rho[200], name_acqui[200];
	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	int r = sscanf(linha,"%s",c->name_vp);
	printf("arquivo de velocidades: %s\n", c->name_vp);


	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	//Model dimensions
	printf("Tamanho do grid e espacamento (Nx Nz H)\n");
	r = sscanf(linha, "%d %d %lf", &c->nx, &c->nz, &c->h);
	printf("Nx Nz H: %d %d %lf\n\n", c->nx, c->nz, c->h);

	c->dx = c->h;
	c->dz = c->h;
	//Time-dependant paremeters

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%d %lf", &c->nt, &c->dt);
	printf("NT DT: %i %f\n\n", c->nt, c->dt);
	c->T= c->nt*c->dt;



	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%d", &c->iss);
	

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%lf",&c->fm);

	//Acquisition

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%s",c->name_acqui);
	printf("arquivo de aquisicao = %s \n", c->name_acqui);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%i", &c->pmlx);
	c->nsp = c->pmlx;
	printf("PMLX = %i\n\n", c->nsp);
	c->pmlz = c->pmlx;




	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%i", &c->oop);
	printf("OOP = %i\n\n", c->oop);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha, "%i", &c->istep);
	printf("istep = %i\n\n", c->istep);

	c->ndt = 1;
	c->n2e=c->nz+2*c->nsp;
	c->n1e=c->nx+2*c->nsp;
	int w = ceil((double)c->nt/c->ndt);
	c->namostras = w;
	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha,"%s %s",c->name_vp0, c->name_rho0);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha,"%d",&c->nshots);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha,"%d",&c->niter);

	printf("Numero de tiros %d\n", c->nshots);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha,"%lf",&c->vpmin);

	i = getline(&linha, &comp, fp);
	remove_comments(linha);
	r = sscanf(linha,"%lf",&c->vpmax);

	printf("vp min  = %lf   vpmax = %lf\n", c->vpmin, c->vpmax);

	if(c->iss==1){
		c->nts=5;
	} else{
		double td=sqrt(6)/(M_PI*c->fm);
		c->nts=(int)(5*td/c->dt)+1;
	}
	if(c->nts>c->nt){
		c->nts=c->nt;
	}

	c->alpha = 0.54;
	c->csi = 6.0;
	c->length_conv_operator=12;

}
int *alocaMatrizInt(int n1e, int n2e) {
	int *v = (int *)malloc(n1e*n2e*sizeof(int));
	if(!v) {
		printf("erro alocando memoria int\n");
		exit(1);
	}
	return v;
}

void leAquisicao(Config *c, int *ir2, int *ir1) {
	FILE *facqui;
	double zs, xs;
	int i;
	facqui = fopen(c->name_acqui, "r");
	if(facqui == NULL){
		perror(c->name_acqui);
		exit(1);
    }
	int r = fscanf(facqui, "%lf %lf", &zs, &xs);
	int nrec=-1;
	char ch;
	while((ch = fgetc(facqui)) != EOF ){
		if(ch == '\n'){
				nrec++;
		}
	}
	c->nrec = nrec;

	double t3, t4, t5;
	rewind(facqui);
	
	r = fscanf(facqui, "%lf %lf %lf %lf %lf", &zs, &xs, &t3, &t4, &t5);
	printf("zs = %lf  e xs = %lf nrec = %d\n", zs,xs, nrec);
	double zr, xr;
	double h = c->h;
	int is1, is2;
	for(i=0; i<nrec; i++) {
		r = fscanf(facqui, "%lf %lf %lf %lf %lf", &zr, &xr, &t3, &t4, &t5);
		is2=((int)(xr/h + 0.5));
		is1=((int)(zr/h + 0.5));
		ir2[i]=is2;
		ir1[i]=is1;
	//	printf("%d %d %d %lf %d\n", i, ir1[i], ir2[i], h, c->nsp);
	}
	// SOURCE POSITION
	is2=((int)(xs/h + 0.5))+c->nsp;
	is1=((int)(zs/h + 0.5))+c->nsp;
	c->is1 = is1;
	c->is2=is2;
	c->sx = xs;
	c->sz = zs;
	printf("posicao da fonte : %i %i\n", is1, is2);
	fclose(facqui);
}

void read_model2(char filename[], Config c, double **model)
{
	FILE *fp;
	int i, j;
	fp = fopen(filename,"rb");
	printf("lendo arquivo %s\n", filename);
	if(!fp) {
		perror(filename);
		exit(1);
	}
	float v;


	for(j=0; j<c.nx; j++) {	  
		for(i=0; i<c.nz; i++) {	
			int dumm = fread(&v, sizeof(float), 1, fp);
			//	printf("%d %d %lf\n", i, j, v);
			model[j][i] = v;
		}
	}
	fclose(fp); 
}

void save_model(char filename[], double **model, int nx, int nz)
{
	FILE *fp;
	int i;
	fp = fopen(filename,"wb");
	if(!fp) {
		perror(filename);
		exit(1);
	}
	float *m;
	int j;
	m = float_vector(nz);

	for(i=0; i<nx; i++) {
		for(j=0; j<nz; j++) {
			m[j] = (float)model[i][j];
			//  m[i][j] =  1.0;
		}
		fwrite(m, sizeof(float), nz, fp);
	}
	fclose(fp); 
}

void save_model2(char filename[], double *model, int nx, int nz)
{
	FILE *fp;
	int i;
	fp = fopen(filename,"wb");
	if(!fp) {
		perror(filename);
		exit(1);
	}
	float *m;
	int j;
	m = float_vector(nz);

	for(i=0; i<nx; i++) {
		for(j=0; j<nz; j++) {
			m[j] = (float)model[i*nz+j];
			//  m[i][j] =  1.0;
		}
		fwrite(m, sizeof(float), nz, fp);
	}
	fclose(fp); 
}


