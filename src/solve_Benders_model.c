#include "def.h"

extern double     **c, **c_c, **c_t, **c_d, **f, **W, *O, *D, **b;
extern double     collect,transfer,distribute;
extern int        NN, Q, p_hubs;
extern double     MAX_DOUBLE;
extern double     UpperBound;
extern coordinate *pts;
extern int        **pos_z;
extern int        pos_eta;
extern double     old_lower_bound;
extern double     ***alpha;
extern double     ***beta;
extern double     **core;
extern double     LP_lower_bound;
extern SCUTS      sepcut;
extern INCUTS     *initial_cuts;
extern double     sum_core;
extern double     old_objval;
extern int        MG;
extern double     *initial_x;
extern    int     count_same_node;
extern int        *cand_hubs;
extern int        count_cand_hubs;
extern int        *fixed_zero;
extern int        *fixed_one;
extern int        count_added;
extern  double     sum_supply_i, sum_supply_j;
//extern int Capacitated_instances;
//extern int        *best_assigmnent;

