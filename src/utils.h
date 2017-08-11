
#include "libfwi.h"
#include "math.h"

double **double_matrix(int nlinhas, int ncolunas);

void free_matrix(double **m, int nlinhas, int ncolunas);

float *float_vector(int nlinhas);

double *double_vector(int nlinhas); 
void leParametros(Config *c, FILE *fp);
size_t remove_comments(char *buffer);
void leAquisicao(Config *c, int *ir2, int *ir1);
void read_model2(char filename[], Config c, double **model);
void save_model(char filename[], double **model, int nx, int nz);
void save_model2(char filename[], double *model, int nx, int nz);
int *alocaMatrizInt(int n1e, int n2e);
