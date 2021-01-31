#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include <cplex.h>
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define rand_double()(rand()/(double) RAND_MAX)
#define MAX_DOUBLE  10000000000
// Tolerances that we keep for the LP bounds
#define epsilon_LP_Pre  0.0005;
#define epsilon_LP_BB  0.5;

struct cutinfo {
   CPXLPptr lp;
   int      numcols;
   int      num;
   double   *x;
   int      *beg;
   int      *ind; 
   double   *val;
   double   *rhs;
   int      nodeid;
   double   nodeobjval;
   int      objsen;
};

typedef struct cutinfo CUTINFO, *CUTINFOptr;

typedef struct INCUTS {
	char    *sense;
	double  *rhs;
	int     *beg;
	int     *ind;
	double  *val;
	double  *origval;
	int     index;
} INCUTS;

typedef struct SCUTS {
  double value;
  int cutnz;
  int *cutind;
  double *cutval;
 } SCUTS;

typedef struct COMMOD {
  int  *i;
  int  *j;
  int  dim;
} COMMOD;

typedef struct EDGES {
  int  e1;
  int  e2;
  int  e;
} EDGES;

typedef struct COV {
	int rhs;
	int *z;
} COV;

typedef struct NORD {
	int i;
	int removed;
	double coeff;
} NORD;

typedef struct COMM {
  double diff;
  int    comm;
} COMM;

typedef struct ZVAL {
	int k;
	double value;
} ZVAL;

typedef struct CCUTS {
  int node;
  double val;
} CCUTS;

typedef struct{
  int i;                        /* index in original */
  int new_i;                    /* new index */
  double x,y;
} coordinate;

typedef struct CORD {
	int j;
	double delta;
} CORD;

typedef struct PORD {
	int i;
	CORD *A;
} PORD;

typedef struct SELEC {
	int hub;
	double cociente;
} SELEC;

typedef struct ORD {
	int i;
	double W;
} ORD;

/********Carlos's modifications*******/
typedef struct Solpool{
	double preeta;
	int * indhub;
	int num_comb;
}Solpool;
char** create_stringarray(int n, int m);
void	free_stringarray(char*** ptr, int n);
void	free_char_vector(char** ptr);
int		add_assignment_constr(CPXENVptr env, CPXLPptr lp);
int		reduced_cost_var_elimination(CPXENVptr env, CPXLPptr lp, double value, double* dj, int* count_fixed, int* total_assign_fixed, int* flag_fixed, int* assign_fixed);
int		pe_fix_to_zero(CPXENVptr env, CPXLPptr lp, double* x, double value, double* dj, int* count_fixed);
int		pe_fix_to_one(CPXENVptr env, CPXLPptr lp, double* x, double value, double* dj, int* count_fixed);
int		partial_enumeration(CPXENVptr env, CPXLPptr lp, double* x, double value, double* dj, int* count_fixed, int flag_fix_to_zero, int flag_fix_to_one);
int		solve_Benders_subproblem(CPXENVptr env, CPXLPptr lp, double* x, double value);
double  solve_as_LP(CPXENVptr env, CPXLPptr lp);
int		set_mip_start(CPXENVptr env, CPXLPptr lp, int cur_numcols);

time_t		t; //Time stamps
struct tm	*tm;//time pointer	
int Compareval(const void *, const void *);

double Det_Iterated_Local_Search(void);

int Comparevalue_zo(const void *, const void *);
int Comparevalue_zc(const void *, const void *);
int Reassign_nodes(int *);
void Ordenar_costos(void);
int Compararcostos(const void *, const void *);
int Compararvalor_fb(const void *, const void *);
double Evaluate_Quadratic_Objective(int *, int *, double *);
void SSCFLP_model(void);
int clients_swap1_f(int *, double *);
int clients_shift1_f(int *, double *);
double sum2_F_ijkm(int , int , int , int , int *);
double sum_F_ijkm(int , int , int *);
int open_hub(int *, double *);
int close_hub(int *, double *);
int open_close_hub(int *, double *);
int NetFlow_TP(double *, int , int );
static int buildNetwork (CPXENVptr , CPXNETptr , double *, int , int);
void Define_Core_Point(void);
void Update_Core_Point(double *);
void PopulateLPSupport(double* x);
int CompareLPSupport(double* x);
int SetBranchandCutParam(CPXENVptr env, CPXLPptr lp);
int CPXPUBLIC mycutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                  void *cbhandle, int *useraction_p);
int CPXPUBLIC mycheckcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                  void *cbhandle, int *useraction_p);
/****Carlos modification****/
int CPXPUBLIC Heur(CPXCENVptr env, void *cbdata, int wherefrom,  void *cbhandle, double *objval_p, double *x, int *checkfeas_p, int *useraction_p);
void Update_CP_MW(double *z_sol,int i, int j);
void  free_and_null (char **ptr);
void read_heur_cl_param(const char*);
 void i_vector(int **vector,int n,char *s);
 void d_vector(double **vector,int n,char *s);
 void c_vector(char **vector,int n,char *s);
 void free_memory(void);
 void initialize_memory(void);
 void read_instance(const char *, int, double, double, double, int);
 FILE *open_file(const char *, const char *);
 int **create_int_matrix(int, int);
 double **create_double_matrix (int, int);
 int *create_int_vector(int);
 double *create_double_vector(int);
 void i_vector(int **vector,int n,char *s);
 void d_vector(double **vector,int n,char *s);
 void c_vector(char **vector,int n,char *s);
 void Check_CP_MW(double *z_sol,int i, int j);
 int clients_shift_red(int *, int *, double *, double *);
 int clients_swap_red(int *, int *, double *, double *);
 void Benders_root_node_heur(void);
 int Construct_Feasible_Solution(double *, double *);
 int CFLP_reduced_model(int, ZVAL *, int *, int *);
 int Reassign_nodes_red(int *, int *);
 void Facility_Change_Phase(int *, int *, double *, ZVAL *, int, double *, double *);
 int Assignment_Change_Phase(int *, int *, double *, ZVAL *, double *);
 int open_close_hub_red(int *, int *, double *, double *);
 int open_hub_red(int *, int *, double *, double *);
 int close_hub_red(int *, int *, double *, double *);
 /****Carlos Modfications****/
 int  *best_open_plants;
 int	use_firstsolution;
 double Prev_incumbent;
 int        *best_assigmnent;
 char output_text[40];
 Solpool *prevsols; //pool of previous solutions
 int countsols; //coutner of how many solutions we have so far.
 double combtol; //combinatorial cut tolerance.
 /**Variable definitions***/
 int APset;
 int Capacitated_instances; //Indicator of whether the instance is capacitated or not.
 int p_hubs;//How many hubs
 int hybrid;//Are we doing both p_hubs and fixed costs? if hybrid==1 then we have only fixedcosts if 0 then have no fixed costs if hybrid==3 then we have both
 int missed;//Number of times we couldn't solve subproblem
 int vers; //To control which version is being tested vers=-1 is use everything.
 double UFLPrelval;//Value of solving the UFLP.
 double UBsearchUFLP;//Upper bound obtained from the search (Check ANY integer solution found)
 double UBOptUFLP;//Upper bound obtained from optimal of UFLP
 int FlagsingInd; //Flag that will let us know that we are solving Single Index.
 int  *FlagHubLPsupport; // array of flags stating whether the potential hub is being considered or not for the reduced problem.
 int* globvarind;  //Array containing all the Types for the complete model
 char* globvarctype;  //Array containing all the Types for the complete model
 int  glob_numcols;
 int* best_sol_facilities;
 int* best_sol_assignments;
 int** not_eligible_hub; //Matrix that keeps track which of the hubs is a feasible assignment or not.
 int* nom_com_assign; //NUmber of commodities assigned per hub
 int* nom_com_assign_few; //NUmber of commodities assigned per hub which there are few candidates.
 int* eli_per_com; //Number of eligible per com
 //For the reduced size of the subproblem
 int Breakpoint_O_D;
 int* index_hub_oi;
 int* index_hub_dj;
 //Additional times
 double cpuFacLocIni; //CPU time of the initial linear facility location problem solved to get an initial feasible solution.
 double cpuGenAss; //CPU time of the Generalized assignment problem.
 //Flags for the heuristic parameters
 int FlagHeuristic;
 int FlagLocalSearch;
 int FlagIteratedLocalSearch;

 /*Variables initially declared in the main*/
 double** c, ** c_c, ** c_t, ** c_d, ** f, ** W, * O, * D, ** b;
 double     collect, transfer, distribute, AD;
 int        NN, Q, p_hubs;
 double     UpperBound;
 coordinate* pts;
 int** pos_z;
 int *priority;
 int* indices;
 int        pos_eta;
 double     old_lower_bound;
 double*** alpha;
 double*** beta;
 double** core;
 double     LP_lower_bound;
 SCUTS      sepcut;
 INCUTS* initial_cuts;
 double     sum_core;
 int        MG;
 double* initial_x;
 int* allocation;
 int* best_allocation;
 int* open_plants;
 double* capacity;
 double* best_capacity;
 double* avail_capacity;
 double     menor_O;
 PORD* costoso;
 SELEC* orden_O;
 double     AggregatedDemand;
 ZVAL* z_open;
 ZVAL* z_closed;
 int* cand_hubs;
 int        count_cand_hubs;
 int* fixed_zero;
 int* fixed_one;
 int        count_added;
 COV        cover;
 NORD* ord_nodes;
 double* multiple;
 double* coeff_ES;
 ORD* ord_O;
 ORD* ord_D;
 double     sum_supply_i, sum_supply_j;
 double     old_objval;
 int        count_same_node;
