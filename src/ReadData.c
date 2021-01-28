#include "def.h"

void read_instance (const char *name, int trans_fact, double coll, double trans, double distr, int APset)
{
  int i,j,p,k,e, ee, value,s, incl, q, fake;
  int terminate, seed_initial, pp;
  FILE *in;
  double value1, value2, coeff_var;
  //ORDER *ord_comm;
  int   *close_edge;
  int all_count_F = 0;
  double  red_lev, rou;
  int  min_i;
  double min_f;
  int min_b, max_b;
  double center_x, center_y, max_dc, min_dc, b_max, f_0;
  double *OD;
  double *dc;
  char path [100];

  Q = 1;       // Number of capacity levels
  //Capacitated_instances = 0;  //Capcitated == 1,  uncapacitated == 0
  sprintf(path,"Data/");
  strcat(path,name);
  in=open_file (path,"r");
  fscanf(in,"%d",&NN);
  initialize_memory();
  OD = (double *)calloc(NN, sizeof(double));
  dc = (double *)calloc(NN, sizeof(double));
  for(p=0;p<NN;p++){
    if(fscanf(in,"%lf %lf",&pts[p].x,&pts[p].y) != 2){
      fprintf(stderr,"ERROR: Can't read coordinates\n");
      exit(1);
    }
    pts[p].i = p;
	OD[p] = 0;
  }
  //comm.dim = 0;
  AD = 0;
  for(i=0;i<NN;i++){
    for(j=0;j<NN;j++){
      if(fscanf(in,"%lf",&W[i][j]) != 1){
	fprintf(stderr,"ERROR: Can't read flow matrix\n");
        exit(1);
      }
	  AD += W[i][j];
	  OD[i] += W[i][j];
	  OD[j] += W[i][j];
    }
  }
  fscanf(in,"%d %lf %lf %lf",&fake,&collect,&transfer,&distribute);
  /*collect = 3.0;
  transfer = 0.75;
  distribute = 2.0;*/
  for(j=0;j<NN;j++) {
    for(i=0;i<NN;i++){
      c[i][j] = trans_fact*(sqrt(pow((pts[i].x-pts[j].x),2)+pow((pts[i].y-pts[j].y),2)))/1000;
      c_c[i][j] = collect*c[i][j];
      c_t[i][j] = transfer*c[i][j];
      c_d[i][j] = distribute*c[i][j];
      //W[i][j] *= 0.1;
    }
  }
  AggregatedDemand = 0;
  menor_O = MAX_DOUBLE;
  for (i = 0; i<NN; i++) {
	  O[i] = 0;
	  D[i] = 0;
	  for (j = 0; j<NN; j++){
		  O[i] += W[i][j];
		  D[i] += W[j][i];
	  }
	  if (O[i] < menor_O)
		  menor_O = O[i];
	  AggregatedDemand += O[i];
	 // printf("O[%d]= %.2f D[%d}= %.2f \n", i+1, O[i], i+1, D[i]);
  }
  Ordenar_costos();
  if (APset == 1){
	  for (i = 0; i<NN; i++){
		  if (fscanf(in, "%lf", &f[i][0]) != 1){
			  fprintf(stderr, "ERROR: Can't read fixed costs\n");
			  exit(1);
		  }
	  }
	  for (i = 0; i<NN; i++){
		  if (fscanf(in, "%lf", &b[i][0]) != 1){
			  fprintf(stderr, "ERROR: Can't read capacities\n");
			  exit(1);
		  }
	  }
  }
  else{
	  for (i = 0; i < NN; i++){
		  O[i] = 0;
		  D[i] = 0;
		  for (j = 0; j < NN; j++){
			  O[i] += W[i][j];
			  D[i] += W[j][i];
		  }
		  ord_O[i].i = i;
		  ord_O[i].W = O[i];
		  ord_D[i].i = i;
		  ord_D[i].W = D[i];
		  //printf("O[%d]: %.1f  b[%d]: %1f \n", i+1,O[i],i+1,b[i][0]);
	  }
	  qsort((ORD *)ord_O, NN, sizeof(ord_O[0]), Compareval);
	  qsort((ORD *)ord_D, NN, sizeof(ord_D[0]), Compareval);
	  center_x = 0;
	  center_y = 0;
	  for (i = 0; i < NN; i++){
		  center_x += OD[i] * pts[i].x;
		  center_y += OD[i] * pts[i].y;
	  }
	  center_x /= AD;
	  center_y /= AD;
	  max_dc = 0;
	  min_dc = MAX_DOUBLE;
	  min_i = -1;
	  for (i = 0; i<NN; i++){
		  dc[i] = (sqrt(pow((pts[i].x - center_x), 2) + pow((pts[i].y - center_y), 2))) / 1000;
		  if (dc[i] > max_dc)
			  max_dc = dc[i];
		  if (dc[i] < min_dc){
			  min_dc = dc[i];
			  min_i = i;
		  }
	  }
	  //New generated fixed costs and capacities
	  pp = 5;
	  b_max = 0;
	  for (i = 0; i<NN; i++){
		  b[i][0] = (NN / pp + (3 * dc[i] * O[i]) / (5 * max_dc*ord_O[NN - 1].W))*O[i];
		  if (b[i][0] > b_max)
			  b_max = b[i][0];
	  }
	  f_0 = 0;
	  for (i = 0; i < NN; i++){
		  for (j = 0; j < NN; j++){
			  f_0 += (c_c[i][min_i] + c_d[min_i][j])*W[i][j];
			  f_0 -= (c_t[i][j])*W[i][j];
		  }
	  }
	  f_0 /= pp;
	  //printf("f_0= %.2f \n",f_0);
	  for (i = 0; i < NN; i++){
		  f[i][0] = f_0*(1 - (3 * dc[i]) / (5 * max_dc));
		  //printf("f[%d]= %.2f \n", i+1, f[i][0]);
		  //f[i][0] = f_0*(5*(b[i][0] + O[i])/(b_max + ord_O[NN-1].W) + 0.5);
	  }
	  for (i = 0; i<NN; i++){
		  if (fscanf(in, "%lf", &f[i][0]) != 1){
			  fprintf(stderr, "ERROR: Can't read fixed costs \n");
			  exit(1);
		  }
	  }
	  if (Q > 1){
		  red_lev = 0.7;
		  rou = 1.3;
		  for (i = 0; i<NN; i++){
			  b[i][Q - 1] = b[i][0];
			  f[i][Q - 1] = f[i][0];
			  for (q = Q - 2; q >= 0; q--){
				  b[i][q] = red_lev*b[i][q + 1];
				  //f[i][q] = rou*f[i][q+1];
				  //f[i][q] = rou*red_lev*f[i][q+1];
				  f[i][q] = rou*b[i][q] * (f[i][q + 1]) / (b[i][q + 1]);
			  }
		  }
	  }
  }
  if(hybrid==0){  //Are we only solving the one without fixed costs
	  for(i=0;i<NN;i++){
		  f[i][0]=0;
	  }
  }
  if (Capacitated_instances == 0){
	  for (i = 0; i < NN; i++)
		  b[i][0] = AggregatedDemand;
  }
  free(OD);
  free(dc);
 // for(i=0;i<NN;i++) {
	//O[i]=0;
	//D[i]=0;
	//for(j=0;j<NN;j++){
	//  O[i]+=W[i][j];
	//  D[i]+=W[j][i];
	//  //AD += W[i][j];
	//}
 // }
 // initialize_memory2();
 fclose(in);
}

void initialize_memory(void) {
  int i;
  old_objval = 0;
  count_same_node = 0;
  cpuFacLocIni = 0;
  cpuGenAss = 0;
  /****Carlos's modifications for initializing the previous solutions*****/
  countsols=0;
  prevsols=(Solpool *) calloc(NN, sizeof(Solpool));
  prevsols[countsols].indhub=create_int_vector(NN);
  FlagHubLPsupport = create_int_vector(NN);
  not_eligible_hub = create_int_matrix(NN, NN);
  index_hub_oi = create_int_vector(NN);
  index_hub_dj = create_int_vector(NN);
  eli_per_com = create_int_vector(NN);
  for (i = 0; i < NN; i++) {
	  eli_per_com[i] = NN;
  }
  /*********************************************************/
  best_sol_facilities = create_int_vector(NN);
  best_sol_assignments = create_int_vector(NN);
  ord_O = (ORD *)calloc(NN, sizeof(ORD));
  ord_D = (ORD *)calloc(NN, sizeof(ORD));
  c = create_double_matrix(NN,NN);
  c_c = create_double_matrix(NN,NN);
  c_t = create_double_matrix(NN,NN);
  c_d = create_double_matrix(NN,NN);
  f = create_double_matrix(NN, Q);
  W = create_double_matrix(NN,NN);
  pts = (coordinate *)malloc((NN)*sizeof(coordinate));
  pos_z = create_int_matrix(NN,NN);
  //pos_eta = create_int_vector(NN);
  core = create_double_matrix(NN,NN);
  O = create_double_vector(NN);
  D = create_double_vector(NN);
  b = create_double_matrix(NN, Q);
  sepcut.cutind = create_int_vector(NN + 1 + (NN-1)*NN);
  sepcut.cutval = create_double_vector(NN + 1 + (NN - 1)*NN);
  alpha = (double ***) calloc(NN, sizeof(double **));
  beta = (double ***) calloc(NN, sizeof(double **));
  for(i=0;i<NN;i++){
	    alpha[i] = create_double_matrix(NN,NN);
        beta[i] = create_double_matrix(NN,NN);
  }
  initial_x = create_double_vector(NN*NN);
  best_assigmnent = create_int_vector(NN);
  open_plants = create_int_vector(NN);
  allocation = create_int_vector(NN);
  best_allocation = create_int_vector(NN);
  capacity = create_double_vector(NN);
  avail_capacity = create_double_vector(NN);
  best_capacity = create_double_vector(NN);
  costoso = (PORD *)calloc(NN, sizeof(PORD));
  orden_O = (SELEC *)calloc(NN, sizeof(SELEC));
  for (i = 0; i<NN; i++)
	  costoso[i].A = (CORD *)calloc(NN, sizeof(CORD));
  cover.z = create_int_vector(NN);
  ord_nodes = (NORD *)calloc(NN, sizeof(NORD));
}

void free_memory(void)
{
  int i,j;
  free(index_hub_oi);
  free(index_hub_dj);
  free(FlagHubLPsupport);
  free(eli_per_com);
  for (i = 0; i < NN; i++){
	  free(f[i]);
	  free(b[i]);
	  free(W[i]);
	  free(c[i]);
	  free(c_c[i]);
	  free(c_t[i]);
	  free(c_d[i]);
	  free(pos_z[i]);
	  free(core[i]);
	  for (j = 0; j < NN; j++){
		  free(alpha[i][j]);
		  free(beta[i][j]);
	  }
	  free(alpha[i]);
	  free(beta[i]);
	  free(not_eligible_hub[i]);
  }
  free(not_eligible_hub);
  free(W);
  free(pts);
  free(c);
  free(c_c);
  free(c_t);
  free(c_d);
  free(f);
  free(alpha);
  free(beta);
  free(pos_z);
  //free(pos_eta);
  free(core);
  free(sepcut.cutind);
  free(sepcut.cutval);
  free(O);
  free(D);
  free(b);
  free(ord_O);
  free(ord_D);
  free(initial_x);
  free(best_assigmnent);
  free(allocation);
  free(best_allocation);
  free(capacity);
  free(best_capacity);
  free(avail_capacity);
  for (i = 0; i<NN; i++)
	  free(costoso[i].A);
  free(cover.z);
  free(prevsols);
  free(best_sol_facilities);
  free(best_sol_assignments);
}

void Ordenar_costos(void)
{
	int i, j;
	for (i = 0; i<NN; i++) {
		costoso[i].i = i;
		for (j = 0; j<NN; j++) {
			costoso[i].A[j].j = j;
			costoso[i].A[j].delta = c[i][j];
		}
		qsort((CORD *)costoso[i].A, NN, sizeof(costoso[i].A[0]), Compararcostos);
	}
	for (i = 0; i < NN; i++){
		orden_O[i].hub = i;
		orden_O[i].cociente = O[i];
	}
	qsort((SELEC *)orden_O, NN, sizeof(orden_O[0]), Compararvalor_fb);
}

int Compareval(const void *a, const void *b)
{
	if (((ORD *)a)->W<((ORD *)b)->W)
		return 1;
	if (((ORD *)a)->W>((ORD *)b)->W)
		return -1;
	return 0;
}

int Compararvalor_fb(const void *a, const void *b)
{
	if (((SELEC *)a)->cociente>((SELEC *)b)->cociente)
		return 1;
	if (((SELEC *)a)->cociente<((SELEC *)b)->cociente)
		return -1;
	return 0;
}

int Compararcostos(const void *a, const void *b)
{
	if (((CORD *)a)->delta>((CORD *)b)->delta)
		return 1;
	if (((CORD *)a)->delta<((CORD *)b)->delta)
		return -1;
	return 0;
}

FILE *open_file (const char *name, const char *mode)
{
 FILE *file;
 if((file=fopen(name,mode))==NULL) {
    printf("\nError: Failed to open file\n");
    exit(8);
 }
 return file;
}

int **create_int_matrix (int rows, int Columns)
{
 int i;
 int **ptr;
 if((ptr=(int **) calloc (rows, sizeof(int *)) ) == NULL) {
   printf ("\nError: Memoria insuficiente\n");
   exit (8);
 }
 for(i=0;i<rows;i++)
   ptr[i] = create_int_vector(Columns);
 return ptr;
}

double **create_double_matrix (int rows, int Columns)
{
 int i;
 double **ptr;
 if((ptr=(double **) calloc (rows, sizeof(double *)))==NULL) {
    printf("\nError: Memoria insuficiente\n");
    exit(8);
  }
  for(i=0;i<rows;i++) {
    ptr[i]=create_double_vector(Columns);
  }
  return ptr;
}

int *create_int_vector (int dim)
{
 int *ptr;
 if((ptr=(int *) calloc (dim, sizeof(int)))==NULL) {
    printf("\nError: Memoria insuficiente\n");
    exit(8);
 }
 return ptr;
}

double *create_double_vector (int dim)
{
 double *ptr;
 if((ptr=(double *) calloc (dim, sizeof(double)))==NULL) {
    printf("\nError: Memoria insuficiente\n");
    exit(8);
 }
 return ptr;
}

void i_vector(int **vector,int n,char *s)
{
if((*vector=(int *)calloc(n,sizeof(int)))==NULL)
  //error(s);
  printf("Error \n");
return;
}

void d_vector(double **vector,int n,char *s)
{
if((*vector=(double *)calloc(n,sizeof(double)))==NULL)
 // error(s);
 printf("Error \n");
return;
}

void c_vector(char **vector,int n,char *s)
{
if((*vector=(char *)calloc(n,sizeof(char)))==NULL)
  //error(s);
  printf("Error \n");
return;
}

void free_and_null(char** ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */