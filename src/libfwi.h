


typedef struct {
	double h; //espacamento da rede
	double dx, dz;
	int oop;  // ordem do operador derivada
	int nt;   // numero de passos no tempo
	int ndt;  // salva campos de onda a cada ndt
	int n1e, n2e;  // dimensoes expandidas, n1e = nx +28pml
	double dt;   // intervalo de tempo
	int istep;  
	int nts;
	int is1, is2;
	int nsp;
	int pmlx, pmlz;
	int nx, nz;
	char name_vp[2000], name_rho[2000], name_acqui[2000];
	char name_vp0[2000], name_rho0[2000];
	int iss;    // pico de frequencia da Ricker
	double fm; //frequencia maxima da Ricker
	int nrec;
	int namostras; // numero de amostras nos sismogramas
	int nshots;	   // numero de shots a processar
	int niter;    // numero de iteracoes maximo
	double T;
	int length_conv_operator;
	double alpha;
	double csi;
	double sx, sz;
	int rz;
	double dr;
	int ntraces;
	double vpmin, vpmax;
} Config;
