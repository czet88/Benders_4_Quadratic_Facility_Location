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

int SetBranchandCutParam(CPXENVptr env, CPXLPptr lp) {	//Here we will set the branch and cut parameters.
	int acumstatus = 0;
	//Here we are changing the variable types so that we can solve as a MIP
	/************************************************/
	acumstatus += CPXchgctype(env, lp, glob_numcols, globvarind, globvarctype);

	//CPLEX Branch and cut parameters (Not Fine-Tuned yet)
	/************************************************/
	CPXchgobjsen(env, lp, CPX_MIN);
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); // Number of threads to use
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetdblparam(env, CPX_PARAM_TILIM, 86400); // time limit
	//CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.0000001); // e-optimal solution (%gap)
	CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	CPXsetintparam(env, CPX_PARAM_PREIND, 0);
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UpperBound + 0.001); // provide an initial upper bound
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetdblparam(env, CPX_PARAM_TILIM, 86400); // time limit
	CPXsetdblparam(env, CPX_PARAM_TILIM, 2000); // time limit
	//CPXsetdblparam(env, CPX_PARAM_TRELIM, 14000); // B&B memory limit
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.0000001); // e-optimal solution (%gap)
											  //CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0); // add cuts or not
	//CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, 0.05);
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); //output display
	acumstatus += CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);											/* Do not use presolve */
	acumstatus += CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 					/* Turn on traditional search for use with control callbacks */
	acumstatus += CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);									 /* Let MIP callbacks work on the original model */



	return 0;
}



int CPXPUBLIC 
mycutcallback (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               int        *useraction_p)
{
   int status = 0;


   CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;

   int      numcols  = cutinfo->numcols;
   int      numcuts  = cutinfo->num;
   double   *x       = cutinfo->x;
   int      *beg     = cutinfo->beg;
   int      *ind     = cutinfo->ind;
   double   *val     = cutinfo->val;
   double   *rhs     = cutinfo->rhs;
   int      *cutind  = NULL;
   double   *cutval  = NULL;
   int      *feas    = NULL;
   double   cutvio;
   int      addcuts = 0;
   int      i, j, k, m, cutnz;
   int      non_int_feasible;
   double   objval;
   double   new_lower_bound;
   int      stop_cutgen, count_b;
   int      count_violcuts = 0;
   int      flag_solve_SP = 0;
   double   tolerance_sep;
   double   lhs, coeff_z;
   int      count_added, depth;
   double   EPSOBJ = 0.1;
   double   epsilon_LP = 0.5;
   int      oldnodeid = cutinfo->nodeid;
   double   oldnodeobjval = cutinfo->nodeobjval;
   int		Flag_Integer=0;

   *useraction_p = CPX_CALLBACK_DEFAULT; 



   // Get current depth in BB
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
   if (status) {
	   fprintf(stdout, "Can't get depth for node.");
	   goto TERMINATE;
   }

   //if (wherefrom == CPX_CALLBACK_MIP_CUT_LOOP ||
	  // wherefrom == CPX_CALLBACK_MIP_CUT_LAST) {
	  // int    oldnodeid = cutinfo->nodeid;
	  // double oldnodeobjval = cutinfo->nodeobjval;

	   /* Retrieve nodeid and node objval of the current node */

	   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
		   CPX_CALLBACK_INFO_NODE_SEQNUM,
		   &cutinfo->nodeid);
	   if (status) {
		   fprintf(stderr, "Failed to get node id.\n");
		   goto TERMINATE;
	   }

	   if (oldnodeid == cutinfo->nodeid){
		   count_same_node++;
	   }
	   else{
		   count_same_node = 0;
	   }

	  // status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
		 //  CPX_CALLBACK_INFO_NODE_OBJVAL,
		 //  &cutinfo->nodeobjval);
	  // if (status) {
		 //  fprintf(stderr, "Failed to get node objval.\n");
		 //  goto TERMINATE;
	  // }

	  // /* Abort the cut loop if we are stuck at the same node
	  // as before and there is no progress in the node objval */

	  // if (oldnodeid == cutinfo->nodeid && depth % 5 == 0) {
		 //  double objchg = (cutinfo->nodeobjval - oldnodeobjval);
		 //  /* Multiply objchg by objsen to normalize
		 //  the change in the objective function to
		 //  the case of a minimization problem */
		 //  objchg *= cutinfo->objsen;
		 //  if (objchg <= EPSOBJ) {
			//   *useraction_p = CPX_CALLBACK_ABORT_CUT_LOOP;
			//   goto TERMINATE;
		 //  }
	  // }
   //}

   /* If we reached this point, we are
   .. in a lazyconstraint callback, or
   .. in a user cut callback, and cuts seem to help
   improving the node objval.
   In both cases, we retrieve the x solution and
   look for violated cuts. */



   //Get lower bound at current node
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_OBJVAL, &new_lower_bound);
    if ( status ) {
      fprintf (stdout, "Can't get depth for node.");
      goto TERMINATE;
   }

   status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, &objval);
 // status = CPXgetcallbackgloballb (env, cbdata, wherefrom, &objval, 0, 0); 
   if ( status ) {
      fprintf (stdout, "Can't get objective value for node.");
      goto TERMINATE;
   }

   // We must do something in every integer solution:
   tolerance_sep = -100.0;
   //tolerance_sep = -0.01;
   flag_solve_SP = 0;

   if(wherefrom == CPX_CALLBACK_MIP_CUT_FEAS){
	 flag_solve_SP = 1;
	 tolerance_sep = -0.000000000001;
	 Flag_Integer = 1;
	 //tolerance_sep = -0.001;
   }
   else{
	  // printf("imp rel: %.6f objval: %.2f old_objval: %.2f \n", (double) ABS(objval - old_objval)/objval*100, objval, old_objval);
	   /*if (depth == 0 && 100*ABS(objval - old_objval) / objval > epsilon_LP) {
		   printf("imp rel: %.2f \n", ABS(objval - old_objval) / objval);
		   flag_solve_SP = 1;
		   old_objval = objval;
	   }
      else{*/
	   if (depth > 0 && depth % 10 == 0 && count_same_node < 2){
		   flag_solve_SP = 1;
		   old_objval = objval;

	   }
      //}
   }

   count_added = 0;
   if(flag_solve_SP == 1){

     status = CPXgetcallbacknodex (env, cbdata, wherefrom, x, 0, numcols-1);
     if ( status ) {
       fprintf(stderr, "Failed to get node solution.\n");
       goto TERMINATE;
     }

	 //if (vers != 3) { //modifiedy by ivan 5172019
		 if (vers != 3 && Flag_Integer == 1) { //modifiedy by ivan 5172019
		 Define_Core_Point();
		 Update_Core_Point(x); //modifiedy by ivan 5172019
		 //printf("would have updated\n");
	 }
	 
	//printf("Started separating\n");
     for(i=0;i<NN-1;i++){                                 //Solve (i,j) primal/dual subproblems
       for(j = i+1;j<NN;j++){
		   //printf("Separating %d and %d\n", i, j);
		   //if(MG==1) Update_CP_MW(x,i,j);
		   Check_CP_MW(x, i, j); //modifiedy by ivan 5172019
		   //sum_core=2;
		   status = NetFlow_TP(x, i, j);
	   }
	 }
	 //printf("Finished separating\n");
	 sepcut.cutnz = 0;
	 sepcut.cutind[sepcut.cutnz] = pos_eta;
     sepcut.cutval[sepcut.cutnz++] = 1;
	 lhs = x[pos_eta];
	   // printf("eta%d%d - ", i,j);
	 for (i = 0; i < NN; i++) {
		 if (fixed_zero[i] == 0){
			 for (k = 0; k < NN; k++) {
				 if (fixed_zero[k] == 0){
					 coeff_z = 0;
					 for (j = 0; j < NN; j++)
						 coeff_z += (-alpha[i][j][k] - beta[j][i][k]);
					 if (ABS(coeff_z) > 0.00001) {
						 sepcut.cutind[sepcut.cutnz] = pos_z[i][k];
						 sepcut.cutval[sepcut.cutnz++] = coeff_z;
						 lhs += coeff_z*x[pos_z[i][k]];
					 }
				 }
			 }
		 }
		 else{
			 for (k = 0; k < NN; k++) {
				 if (i != k && fixed_zero[k] == 0){
					 coeff_z = 0;
					 for (j = 0; j < NN; j++)
						 coeff_z += (-alpha[i][j][k] - beta[j][i][k]);
					 if (ABS(coeff_z) > 0.00001) {
						 sepcut.cutind[sepcut.cutnz] = pos_z[i][k];
						 sepcut.cutval[sepcut.cutnz++] = coeff_z;
						 lhs += coeff_z*x[pos_z[i][k]];
					 }
				 }
			 }

		 }
	 }
	 if(lhs < tolerance_sep){
	     count_added++;
	     status = CPXcutcallbackadd (env, cbdata, wherefrom, sepcut.cutnz, 0, 'G', sepcut.cutind, sepcut.cutval, CPX_USECUT_PURGE);
         if( status ) {
           fprintf (stderr, "Failed to add cut.\n");
           goto TERMINATE;
		 }
	   }
   }

   if(count_added > 0){
	   stop_cutgen = 0;
	   //printf("Added cuts: %d, percent decrease: %.3f \n", count_violcuts, (objval - last_objval)/objval*100);
	   //old_lower_bound = new_lower_bound;
   }
   else
	   stop_cutgen = 1;

   if(depth == 0 && count_added == 0)
      LP_lower_bound = cutinfo->nodeobjval;

   /* Tell CPLEX that cuts have been created */ 

   if(stop_cutgen == 0)
     *useraction_p = CPX_CALLBACK_SET; 

TERMINATE:

   return (status);

} /* END mycutcallback */



/* This simple routine frees up the pointer *ptr, and sets *ptr
   to NULL */

void free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */ 



int NetFlow_TP(double *sol_z, int Origin_Node, int Destin_Node)
{
   /* Declare variables and arrays for retrieving problem data and
      solution information later on. */

   int      narcs;
   int      nnodes;
   int      solstat;
   double   objval;
   double   *x     = NULL;//Arc Flows
   double   *pi    = NULL;//Shadow Prices for Flow Balance Constraints at nodes
   double   *slack = NULL;//Slack in flow balance constraints
   double   *dj    = NULL;//Reduced costs for arc flows

   CPXENVptr env = NULL;
   CPXNETptr net = NULL;
   int       status;
   int       i, j, k, m;

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);


   if ( env == NULL ) {
      char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      goto TERMINATE;
   }

   /* Create the problem. */

   net = CPXNETcreateprob (env, &status, "netex1");

if (net == NULL) {
	fprintf(stderr, "Failed to create network object.\n");
	goto TERMINATE;
}

//system("pause");
status = buildNetwork(env, net, sol_z, Origin_Node, Destin_Node);
//system("pause");
if (status) {
	fprintf(stderr, "Failed to build network problem.\n");
	goto TERMINATE;
}


/* Optimize the problem and obtain solution. */
CPXsetdblparam(env, CPX_PARAM_TILIM, 100);
//CPXsetintparam(env, CPX_PARAM_NETPPRIIND, 2);
status = CPXNETprimopt(env, net);
if (status) {
	fprintf(stderr, "Failed to optimize network.\n");
	goto TERMINATE;
}

/* get network dimensions */

narcs = CPXNETgetnumarcs(env, net);
nnodes = CPXNETgetnumnodes(env, net);

/* allocate memory for solution data */
pi = (double*)malloc(nnodes * sizeof(double));

status = CPXNETsolution(env, net, &solstat, &objval, NULL, pi, NULL, NULL);
if (solstat != CPX_STAT_OPTIMAL) {
	missed++;
	/*printf("solstat: %d \n", solstat);
	printf("Warning: sum_supply: %.4f sum_demand: %.4f \n", sum_supply_j, sum_supply_i);
	getchar();*/
	// printf("Warning did not solve status=%d for %d and %d\n", solstat, Origin_Node, Destin_Node); 
	 //buildRealNetwork(env, net, sol_z, Origin_Node, Destin_Node);
	// status = CPXNETprimopt (env, net);
	 //status = CPXNETsolution (env, net, &solstat, &objval, NULL, pi, NULL, NULL);
	 //printf("After second try, solution status is %d\n", solstat);
	 //getchar();

}
if (status) {
	fprintf(stderr, "Failed to obtain solution.\n");
	goto TERMINATE;
}

for (k = 0; k < Breakpoint_O_D; k++) {
	alpha[Origin_Node][Destin_Node][index_hub_oi[k]] = pi[k];					//giving
}
for (m = Breakpoint_O_D; m < nnodes; m++) {
	beta[Origin_Node][Destin_Node][index_hub_dj[m-Breakpoint_O_D]] = (-pi[m]); //receiving
}


TERMINATE:

/* Free memory for solution data */
free_and_null((char**)& pi);

/* Free up the problem as allocated by CPXNETcreateprob, if necessary */

if (net != NULL) {
	status = CPXNETfreeprob(env, &net);
	if (status) {
		fprintf(stderr, "CPXNETfreeprob failed, error code %d.\n", status);
	}
}

/* Free up the CPLEX environment, if necessary */

if (env != NULL) {
	status = CPXcloseCPLEX(&env);
	if (status) {
		char  errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "Could not close CPLEX environment.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
	}
}

return (status);

}  /* END main */


int buildNetwork(CPXENVptr env, CPXNETptr net, double* z_sol, int Origin_Node, int Destin_Node)
{
	int status = 0;
	int k, m;
	int NNODES_i = 0; //No. of potential 1st hub (k)
	int NNODES_j = 0; //No. of potential 2nd hub (m)
	int NNODES = 0;//NNODES_i hubs k and NNODES_j hubs m
	int NARCS;
	int count = 0;

	
	for (k = 0; k < NN; k++) {
		if (not_eligible_hub[Origin_Node][k] == 0) {						
			index_hub_oi[NNODES_i] = k;
			NNODES_i++;
			NNODES++;
		}
		if (not_eligible_hub[Destin_Node][k] == 0) {
			index_hub_dj[NNODES_j] = k; 
			NNODES_j++;
			NNODES++;
		}
	}

	Breakpoint_O_D = NNODES_i;
	NARCS = NNODES_i * NNODES_j;
   
   //int *tail = new int[NARCS], *head = new int[NARCS];
   //double *obj = new double[NARCS],*ub = new double[NARCS],*lb = new double[NARCS];
   //double *supply = new double[NNODES];

   int *tail, *head;
   double *obj, *ub, *lb;
   double *supply;

   i_vector(&tail,NARCS,"open_cplex:4");
   i_vector(&head,NARCS,"open_cplex:6");
   d_vector(&obj,NARCS,"open_cplex:2");
   d_vector(&ub,NARCS,"open_cplex:7");
   d_vector(&lb,NARCS,"open_cplex:7");
   d_vector(&supply,NNODES,"open_cplex:7");
   
   if(vers==3)sum_core=NN;

   sum_core = 1;
   ////Using Maganti-Wong Appproach to obtain Pareto-Optimal cuts
   if (MG == 1){
	   sum_supply_i = 0;
	   for (k = 0; k < NNODES_i; k++){
		   supply[k] = core[Origin_Node][index_hub_oi[k]] + sum_core*z_sol[pos_z[Origin_Node][index_hub_oi[k]]];//Supply at hub k = Zik
		   sum_supply_i += supply[k];
	   }
	   sum_supply_j = 0;
	   for (m = 0; m < NNODES_j; m++){
		   supply[NNODES_i + m] = -(core[Destin_Node][index_hub_dj[m]] + sum_core*z_sol[pos_z[Destin_Node][index_hub_dj[m]]]);//Supply at hub m = -Zjm
		   sum_supply_j += supply[NNODES_i + m];
	   }
   }
   else{
	   ////Using Papadakos Appproach to obtain Pareto-Optimal cuts
	   //sum_supply_i = 0;
	   for (k = 0; k < NNODES_i; k++){
		   supply[k] = z_sol[pos_z[Origin_Node][index_hub_oi[k]]];//Supply at hub k = Zik
		   //sum_supply_i += supply[k];
	   }
	   //sum_supply_j = 0;
	   for (m = 0; m < NNODES_j; m++){
		   supply[NNODES_i + m] = -z_sol[pos_z[Destin_Node][index_hub_dj[m]]];//Supply at hub m = -Zjm
		  //sum_supply_j += supply[NNODES_i + m];
	   }
   }
   /*if (ABS(sum_supply_j + sum_supply_i) > 0.001)
	   printf("Warning: %.2f %.2f \n", sum_supply_j, sum_supply_i);*/
   //Using Standard Appproach to obtain Optimality cuts
   //for(k=0;k<NNODES_i;k++)
	  // supply[k] = z_sol[pos_z[Origin_Node][k]];//Supply at hub k = Zik

   //for(m=0;m<NNODES_j;m++)
	  // supply[NNODES_i+m] = -z_sol[pos_z[Destin_Node][m]];//Supply at hub m = -Zjm

   for(k=0;k<NNODES_i;k++)//From all first (k) hubs to all second (m) hubs
   {
	   for(m=0;m<NNODES_j;m++)
	   {
		   tail[count]= k;//Hub k node number starts from 0
		   head[count]= NNODES_i+m;//
		   obj[count]=(W[Origin_Node][Destin_Node]*c_t[index_hub_oi[k]][index_hub_dj[m]] + W[Destin_Node][Origin_Node]*c_t[index_hub_dj[m]][index_hub_oi[k]]);
		   //obj[count]=(W[Origin_Node][Destin_Node])*(c_c[Origin_Node][k] + c_t[k][m] + c_d[m][Destin_Node]);
		   //obj[count] = (W[Origin_Node][Destin_Node])*(c_t[cand_hubs[k]][cand_hubs[m]]);
		   ub[count]=CPX_INFBOUND;
		   lb[count]=0;
		   count++;
	   }
   }

   //if ( CPXNETgetnumnodes (env, net) > 0 ) {
   //   status = CPXNETdelnodes (env, net, 0,
   //                            CPXNETgetnumnodes (env, net)-1);
   //   if ( status ) goto TERMINATE;
   //}

   /* Set optimization sense */

   status = CPXNETchgobjsen (env, net, CPX_MIN);
   if ( status ) goto TERMINATE;

   /* Add nodes to network along with their supply values,
      but without any names. */

   status = CPXNETaddnodes (env, net, NNODES, supply, NULL);
   if ( status ) goto TERMINATE;

   /* Add arcs to network along with their objective values and
      bounds, but without any names. */

   status = CPXNETaddarcs (env, net, NARCS, tail, head, lb, ub, obj, NULL);
   if ( status ) goto TERMINATE;

TERMINATE:

   free_and_null ((char **) &tail);
   free_and_null ((char **) &head);
   free_and_null ((char **) &obj);
   free_and_null ((char **) &ub);
   free_and_null ((char **) &lb);
   free_and_null ((char **) &supply);

   return (status);

}  /* END buildnetwork */

int buildRealNetwork (CPXENVptr env, CPXNETptr net, double *z_sol, int Origin_Node, int Destin_Node)
{
   int status = 0;
   int k,m;
   double sum_supply_i, sum_supply_j;
   int NNODES_i = count_cand_hubs;//No. of potential 1st hub (k)
   int NNODES_j = count_cand_hubs;//No. of potential 2nd hub (m)
   int NNODES = NNODES_i + NNODES_j;//NNODES_i hubs k and NNODES_j hubs m
   int NARCS = NNODES_i*NNODES_j;
   int count=0;
   
   //int *tail = new int[NARCS], *head = new int[NARCS];
   //double *obj = new double[NARCS],*ub = new double[NARCS],*lb = new double[NARCS];
   //double *supply = new double[NNODES];

   int *indices;
   double *obj, *ub, *lb;
   double *supply;

   d_vector(&supply,NNODES,"open_cplex:7");
   i_vector(&indices,NNODES,"open_cplex:7");
   
   ////Using the actual solution value.
   //sum_supply_i = 0;
	for (k = 0; k < NNODES_i; k++){
	 supply[k] = z_sol[pos_z[Origin_Node][cand_hubs[k]]];//Supply at hub k = Zik
	 indices[k]=k;
	  //sum_supply_i += supply[k];
	 }
	//sum_supply_j = 0;
	 for (m = 0; m < NNODES_j; m++){
	  supply[NNODES_i + m] = -(z_sol[pos_z[Destin_Node][cand_hubs[m]]]);//Supply at hub m = -Zjm
	  indices[NNODES_i + m]=NNODES_i + m;
	 //sum_supply_j += supply[NNODES_i + m];
	 }

   /* Set optimization sense */
   status = CPXNETchgsupply (env, net, NNODES, indices, supply);

  if ( status ) goto TERMINATE;

TERMINATE:

   free_and_null ((char **) &indices);
   free_and_null ((char **) &supply);

   return (status);

}  /* END buildnetwork */


//void Define_Core_Point(void)
//{
//	int i, k;
//	double precount, epsilon;
//	
//	if(hybrid==0 || hybrid==3){
//		precount=(double) (p_hubs*1.0/(count_cand_hubs));
//		epsilon = (double)(1 / (2 * count_cand_hubs));
//		//printf("Precount=%lf\n",precount);
//		for (i = 0; i < NN; i++){
//			if (fixed_zero[i] == 0){
//				for (k = 0; k < NN; k++){
//					if (i == k){
//						core[i][k] = (double)(precount - epsilon);
//						//printf("corepoint=%lf\n",core[i][k]);
//						sum_core += core[i][k];
//					}
//					else{
//						if (fixed_zero[k] == 0){
//							core[i][k] = (double) ((1 - (precount-epsilon)) / (count_cand_hubs - 1));
//							//("corepoint=%lf\n",core[i][k]);
//							sum_core += core[i][k];
//						}
//					}
//				}
//			}
//			else{
//				for (k = 0; k < NN; k++){
//					if (i != k && fixed_zero[k] == 0){
//						core[i][k] = (double) (1.0 / (count_cand_hubs));
//						//printf("corepoint=%lf\n",core[i][k]);
//						sum_core += core[i][k];
//					}
//				}
//			}
//		}
//	}
//	else{
//		sum_core = 0;
//		if(count_cand_hubs>1){
//			for (i = 0; i < NN; i++){
//				if (fixed_zero[i] == 0){
//					for (k = 0; k < NN; k++){
//						if (i == k){
//							core[i][k] = 0.5;
//							sum_core += core[i][k];
//						}
//						else{
//							if (fixed_zero[k] == 0){
//								core[i][k] = 0.5 / (count_cand_hubs - 1);
//								sum_core += core[i][k];
//							}
//						}
//					}
//				}
//				else{
//					for (k = 0; k < NN; k++){
//						if (i != k && fixed_zero[k] == 0){
//							core[i][k] = (double) (1.0 / (count_cand_hubs));
//							sum_core += core[i][k];
//						}
//					}
//				}
//			}
//		}
//		else{ //This is when there is only one candidate hub
//			for(i=0;i<NN;i++){
//				core[i][cand_hubs[0]]=1.0;
//			}
//		}
//	}
//}

void Define_Core_Point(void)
{
	int i, k;
	double precount, epsilon;
	//Initialize
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			core[i][k] = 0;
		}
	}

	if (hybrid == 0 || hybrid == 3) {
		precount = (double)(p_hubs * 1.0 / (count_cand_hubs));
		epsilon = (double)(1 / (2 * count_cand_hubs));
		//printf("Precount=%lf\n",precount);
		for (i = 0; i < NN; i++) {
			if (eli_per_com[i] > 1) { //If there is more than one candidate then you can do this
				if (fixed_zero[i] == 0) {
					for (k = 0; k < NN; k++) {
						if (i == k) {
							core[i][k] = (double)(precount - epsilon);
							//printf("corepoint=%lf\n",core[i][k]);
							sum_core += core[i][k];
						}
						else {
							if (/*fixed_zero[k] == 0 && */not_eligible_hub[i][k] == 0) {
								core[i][k] = (double)((1 - (precount - epsilon)) / (eli_per_com[i] - 1));
								//("corepoint=%lf\n",core[i][k]);
								sum_core += core[i][k];
							}
						}
					}
				}
				else {
					for (k = 0; k < NN; k++) {
						if (i != k && /*fixed_zero[k] == 0 && */not_eligible_hub[i][k] == 0) {
							core[i][k] = (double)(1.0 / (eli_per_com[i]));
							//printf("corepoint=%lf\n",core[i][k]);
							sum_core += core[i][k];
						}
					}
				}
			}
			else {  //What to do in the event that there is only one potential assignment.
				for (k = 0; k < NN; k++) {
					if (not_eligible_hub[i][k] == 0) {
						core[i][k] = 1.0;
					}
				}
			}
		}
	}
	else {
		sum_core = 0;
		if (count_cand_hubs > 1) {
			for (i = 0; i < NN; i++) {
				if (eli_per_com[i] > 1) { //If there is more than one candidate then you can do this
					if (fixed_zero[i] == 0) {
						for (k = 0; k < NN; k++) {
							if (i == k) {
								core[i][k] = 0.5;
								sum_core += core[i][k];
							}
							else {
								if (/*fixed_zero[k] == 0 && */not_eligible_hub[i][k] == 0) {
									core[i][k] = 0.5 / (eli_per_com[i] - 1);
									sum_core += core[i][k];
								}
							}
						}
					}
					else {
						for (k = 0; k < NN; k++) {
							if (i != k && /*fixed_zero[k] == 0 &&*/ not_eligible_hub[i][k] == 0) {
								core[i][k] = (double)(1.0 / (eli_per_com[i]));
								sum_core += core[i][k];
							}
						}
					}
				}
				else {  //What to do in the event that there is only one potential assignment.
					for (k = 0; k < NN; k++) {
						if (not_eligible_hub[i][k] == 0) {
							core[i][k] = 1.0;
						}
					}
				}
			}
		}
		else { //This is when there is only one candidate hub
			for (i = 0; i < NN; i++) {
				core[i][cand_hubs[0]] = 1.0;
			}
		}
	}
}


void Update_Core_Point(double *z_sol){

  int i,k,j,m;
  double sum_i, sum_j;
  double fact=0.5;		

  for (i = 0; i < NN; i++) {
	  for (j = i; j < NN; j++) {
		  sum_i = 0;
		  sum_j = 0;
		  for (k = 0; k < count_cand_hubs; k++)
			  sum_i += /*z_sol[pos_z[i][cand_hubs[k]]]*/core[i][cand_hubs[k]]/**(1-not_eligible_hub[i][cand_hubs[k]])*/;
		  for (m = 0; m < count_cand_hubs; m++)
			  sum_j += /*z_sol[pos_z[j][cand_hubs[m]]]*/ core[j][cand_hubs[m]]/** (1 - not_eligible_hub[j][cand_hubs[m]])*/;
		 /* if (ABS(sum_i - sum_j) > 0.0001)
			  printf("something wrong with core value for %d and %d: sum_i:%.2f sum_j:%.2f \n", i,j, sum_i, sum_j);*/
	  }
  }

  for (i = 0; i < NN; i++){
	  for (k = 0; k < count_cand_hubs; k++){
		//  if (z_sol[pos_z[i][cand_hubs[k]]] > 0.0001)
		//	  printf("z(%d %d): %.2f \n",i, cand_hubs[k], z_sol[pos_z[i][cand_hubs[k]]]);
		  core[i][cand_hubs[k]] = fact*core[i][cand_hubs[k]] + (1-fact)*z_sol[pos_z[i][cand_hubs[k]]];
	  }
  }

  for (i = 0; i < NN; i++) {
	  for (j = i; j < NN; j++) {
		  sum_i = 0;
		  sum_j = 0;
		  for (k = 0; k < count_cand_hubs; k++)
			  sum_i += core[i][cand_hubs[k]];
		  for (m = 0; m < count_cand_hubs; m++)
			  sum_j += core[j][cand_hubs[m]];
		  /*if (ABS(sum_i - sum_j) > 0.0001)
			  printf("something wrong with core point: sum_i:%.2f sum_j:%.2f \n", sum_i, sum_j);*/
	  }
  }


}

void Update_CP_MW(double *z_sol,int i, int j){
	int h,s;
	double acum_sup, acum_dem=0,temp=0, open=0, checkcore=0;
	sum_core=0;
	for(h=0;h<count_cand_hubs-1;h++){
		open+=z_sol[pos_z[cand_hubs[h]][cand_hubs[h]]];		
		//printf("Hub %d is open with %lf\n", cand_hubs[h],z_sol[pos_z[cand_hubs[h]][cand_hubs[h]]]);
	}	
	for(h=0;h<count_cand_hubs;h++){
		temp=z_sol[pos_z[j][cand_hubs[h]]]-z_sol[pos_z[i][cand_hubs[h]]];
		if(temp>0){
			core[j][cand_hubs[h]]=(double) (0.75/open);			
			core[i][cand_hubs[h]]=core[j][cand_hubs[h]]+temp;
		}
		else{
			if(temp<0){
				core[i][cand_hubs[h]]=(double) (0.75/open);
				/*printf("core val is %lf\n",core[i][cand_hubs[h]]);
				getchar();*/
				core[j][cand_hubs[h]]=core[i][cand_hubs[h]]-temp;
			}
			else {
				core[i][cand_hubs[h]]=(double) (0.1/count_cand_hubs);
				core[j][cand_hubs[h]]=(double) (0.1/count_cand_hubs);
			}
		}
		/****Updating the acumulations and sum_core****/
		sum_core+=core[j][cand_hubs[h]]+core[i][cand_hubs[h]];
		/*acum_dem+=core[j][cand_hubs[h]];
		acum_sup+=core[i][cand_hubs[h]];*/
	}
	/*core[i][cand_hubs[count_cand_hubs-1]]=1-acum_sup;
	core[j][cand_hubs[count_cand_hubs-1]]=1-acum_dem;*/
	/*printf("Total sum_core is %lf\n", sum_core);
	getchar();*/
	/****Printing Corepoint*******/
	/*printf("Corepooint assignments of origin with %d candidate hubs\n", count_cand_hubs);
	for(h=0;h<count_cand_hubs-1;h++){
		if(core[i][cand_hubs[h]]+z_sol[pos_z[i][cand_hubs[h]]]>.000001)printf("core[%d][%d]=%lf vs solution[%d][%d]=%lf\n",i,cand_hubs[h],core[i][cand_hubs[h]],i,cand_hubs[h],z_sol[pos_z[i][cand_hubs[h]]] );
	}

	printf("Corepooint assignments of destination\n");
	for(h=0;h<count_cand_hubs-1;h++){
		if(core[j][cand_hubs[h]]+z_sol[pos_z[j][cand_hubs[h]]]>.000001)printf("core[%d][%d]=%lf vs solution[%d][%d]=%lf\n",j,cand_hubs[h],core[j][cand_hubs[h]],j,cand_hubs[h],z_sol[pos_z[j][cand_hubs[h]]]);
	}
	getchar();*/
}

void Check_CP_MW(double *z_sol,int i, int j){
	int h,s;
	double acum,temp=0, open=0, temp2,checkcore1=0,checkcore2;
	double *values_hub;
	values_hub=create_double_vector(count_cand_hubs);
	acum=0;
	sum_core=0;
	//printf("There are %d candidate hubs\n",count_cand_hubs);
	for(h=0;h<count_cand_hubs;h++){
		/*checkcore1+=core[i][cand_hubs[h]];
		checkcore2+=core[j][cand_hubs[h]];*/
		temp=(z_sol[pos_z[j][cand_hubs[h]]]-z_sol[pos_z[i][cand_hubs[h]]]);
		if(ABS(temp)>0.00001){
			values_hub[h]=(core[i][cand_hubs[h]]-core[j][cand_hubs[h]])/temp;
			/*if(values_hub[h]<0){
				temp2=core[j][cand_hubs[h]];
				core[j][cand_hubs[h]]=core[i][cand_hubs[h]];
				core[i][cand_hubs[h]]=temp2;
				values_hub[h]=-values_hub[h];
			}*/
		}
		else values_hub[h]=0;
		acum+=values_hub[h];
		if(values_hub[h]>sum_core) sum_core=values_hub[h];
		//printf("The delta[%d] value is %lf\n", h, values_hub[h]);
	}
	//printf("Assigned total of %lf from origin %d and %lf from destination %d \n",checkcore1,i,checkcore2,j);getchar();
	
	/*core[i][cand_hubs[count_cand_hubs-1]]=1-acum_sup;
	core[j][cand_hubs[count_cand_hubs-1]]=1-acum_dem;*/
	/*printf("Total sum_core is %lf\n", sum_core);
	getchar();*/
	/****Printing Corepoint*******/
	/*if(acum>0.0001){
		printf("Corepoint assignments of origin with %d candidate hubs\n", count_cand_hubs);
		for(h=0;h<count_cand_hubs-1;h++){
			if(core[i][cand_hubs[h]]+z_sol[pos_z[i][cand_hubs[h]]]>.000001)printf("core[%d][%d]=%lf vs solution[%d][%d]=%lf\n",i,cand_hubs[h],core[i][cand_hubs[h]],i,cand_hubs[h],z_sol[pos_z[i][cand_hubs[h]]] );
		}

		printf("Corepoint assignments of destination\n");
		for(h=0;h<count_cand_hubs-1;h++){
			if(core[j][cand_hubs[h]]+z_sol[pos_z[j][cand_hubs[h]]]>.000001)printf("core[%d][%d]=%lf vs solution[%d][%d]=%lf\n",j,cand_hubs[h],core[j][cand_hubs[h]],j,cand_hubs[h],z_sol[pos_z[j][cand_hubs[h]]]);
		}
		printf("Accumulated value is %lf.\n Sum_core is %lf\n", acum, sum_core);
		getchar();
	}*/
	if(acum<0.001) sum_core=0;
	else sum_core=sum_core;
	free(values_hub);
}

int CPXPUBLIC Heur (CPXCENVptr env, void  *cbdata, int wherefrom, void *cbhandle, double  *objval_p, double *x, int *checkfeas_p,int *useraction_p)
{
   int status = 0; 
   int depth, SEQNUM;
   int       i,k,m,index, j, cols, addcuts=0,opened;
   double    etaval,currincumbent;
   double   *bestx;
   double	*temp_x;
   int flag_newinc=0;
   cols=count_cand_hubs*NN+1; //candidate hubs by the nodes plus the artificial variable
   bestx=create_double_vector(cols);
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
   *useraction_p = CPX_CALLBACK_DEFAULT;
	if(use_firstsolution==1 ){
		use_firstsolution=0; 
		etaval=0;
		for(i=0;i<cols;i++)x[i]= 0;
		for (i = 0; i < NN; i++) {
			x[pos_z[i][best_sol_assignments[i]]] = 1;
			//printf("Assigned %d to hub %d\n",i,best_assigmnent[i]);			
			for (m = 0; m < NN; m++)
				etaval += W[i][m] * (c_t[best_sol_assignments[i]][best_sol_assignments[m]]);
		}
		x[pos_eta]=etaval;
		printf("Value of etaval is %.2lf\n",etaval);	
         *objval_p =UpperBound;   //Now telling it what the new objective value is
		 Prev_incumbent=UpperBound;
		 printf("Writing the best found solution with value %lf\n",UpperBound);
		 *checkfeas_p = 1;
		  /* Tell CPLEX that a solution is being returned */
		  *useraction_p = CPX_CALLBACK_SET;
	}	
	

TERMINATE:
	free(bestx);
   return (0);
} /* END Lagrange_Heur */
