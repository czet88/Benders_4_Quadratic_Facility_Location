#include "def.h"

extern double     **c, **c_c, **c_t, **c_d, **f, **W, *O, *D, **b;
extern double     collect,transfer,distribute;
extern int        NN, Q;
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
extern double     sum_core;
extern double     old_objval;
extern int        MG;
extern double     *initial_x;
extern int        *best_assigmnent;
extern int        *open_plants;
// int Capacitated_instances;

void SSCFLP_model(void)
{
	int i,j,k,l,m, count;
    FILE     *out;
    clock_t  start,end;
	int index,index1, index11;  // indices auxiliares para rellenar las matrices
    double   cputime;
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each column ...............
	int       *matcnt; // number of non-zero element in each column... .................
	int       *matind; // associated row of each non-zelo element ......................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	//double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int *indices, *priority;
	int		nodecount;
	int cur_numcols;
	double *localx;
	int    *local_assign; //Contain the local optimal assignment

	UpperBound=MAX_DOUBLE;
	start = clock();

    objsen = 1; //min
	
	//Initialize CPLEX environment
	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"UHLPSA");
	lp = CPXcreateprob (env, &status, probname);
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}


	                                        //Define z_ik variables
    index1 = 0;  // index of columns
	numcols = NN*NN;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");

    for(i=0;i<NN;i++){
       for(k=0;k<NN;k++){
	     pos_z[i][k] = index1;
	     if(i==k){
           obj[index1] = f[i][0];
		 }
	     else{
			 obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
			 //	 obj[index1] = 0;
		 }
         ctype[index1] = 'B';
         lb[index1] = 0;
         ub[index1] = 1;
         index1++;
	   }
	}
	status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);

	
	                                             //Add assignment constraints  \sum_{k \in N} z_ik = 1
	numrows = NN;
	numnz = NN*NN;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");

    index = 0;
    index1 = 0;
    for(i=0;i<NN;i++){
      sense[index1]='E';
      rhs[index1]= 1;
	  matbeg[index1++] = index;
	  for(k=0;k<NN;k++){
		   matind[index] = pos_z[i][k];
           matval[index++] = 1;
	  }
    }
	status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if( status ) 
      fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


		                                                   //Add linking constraints  z_ik <= z_kk

	numrows = NN*(NN-1);
	numnz = 2*NN*(NN-1);
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");

    index = 0;
    index1 = 0;
	for(i=0;i<NN;i++){
      for(k=0;k<NN;k++){
		if(i!=k){
          sense[index1]='L';
          rhs[index1]= 0;
		  matbeg[index1++] = index;
		  matind[index] = pos_z[i][k];
          matval[index++] = 1;
		  matind[index] = pos_z[k][k];
          matval[index++] = -1;
		}
	  }
	}
	status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if( status ) 
      fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

								//Add exactly p_hubs
	if(hybrid==0 || hybrid==3){
		numrows = 1;
		numnz = NN;
		d_vector(&rhs,numrows,"open_cplex:2");
		c_vector(&sense,numrows,"open_cplex:3");
		i_vector(&matbeg,numrows,"open_cplex:4");
		i_vector(&matind,numnz,"open_cplex:6");
		d_vector(&matval,numnz,"open_cplex:7");
		index = 0;
		index1 = 0;
		sense[index1]='E';
		rhs[index1]= p_hubs;
		matbeg[index1++] = index;
		for(i=0;i<NN;i++){
			matind[index] = pos_z[i][i];
			matval[index++] = 1;
		}
		status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if( status ) 
		  fprintf (stderr,"CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}

	//Add capacity constraints sum(i in NN) O_i z_ik <= b*z_kk
	if (Capacitated_instances == 1){
		numrows = NN;
		numnz = NN*(NN + 1);
		d_vector(&rhs, numrows, "open_cplex:2");
		c_vector(&sense, numrows, "open_cplex:3");
		i_vector(&matbeg, numrows, "open_cplex:4");
		i_vector(&matind, numnz, "open_cplex:6");
		d_vector(&matval, numnz, "open_cplex:7");

		index = 0;
		index1 = 0;
		for (k = 0; k < NN; k++){
			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_z[k][k];
			matval[index++] = -(b[k][0] - O[k]);
			for (i = 0; i < NN; i++){
				if (i != k){
					matind[index] = pos_z[i][k];
					matval[index++] = O[i];
				}
			}
		}
		status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if (status)
			fprintf(stderr, "CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}


	//branch and bound parameters
    CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_OFF); //output display
    //CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
    CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,4); //different levels of output display
	//CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1);
    //CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
    CPXsetdblparam(env,CPX_PARAM_TILIM,7200); // time limit
    //CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
   // CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.00001); // e-optimal solution (%gap)
    //CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
    //CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
//	CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
    //CPXsetintparam(env,CPX_PARAM_HEURFREQ, 0); //heuristic frequency and intensisty 
    //CPXoptimize(env,lp);   solve a linear program
    //CPXgetbestobjval(env,lp,&value);   \\obtain the LP relaxation bound
    //printf("LP lower bound: %f   ",value);
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,238016.277+0.001); // provide an initial upper bound
	CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables

   

   /* Assure linear mappings between the presolved and original
      models */
   status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
   if ( status )  goto TERMINATE;

   /* Turn on traditional search for use with control callbacks */
   status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
   if ( status )  goto TERMINATE;

   /* Let MIP callbacks work on the original model */
   status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
     if ( status )  goto TERMINATE;

	 /* Set up to use MIP callback */

	 status = CPXsetlazyconstraintcallbackfunc (env, mycheckcallback, NULL); 
	  if ( status )  goto TERMINATE;


	CPXmipopt(env,lp);  //solve the integer program

	i=CPXgetstat(env,lp);
	if(i==101)
		printf("Optimal solution found\n");
	else if(i==102)
		printf("e-optimal solution found\n");
	else if(i==107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n",i);
  
    // out = open_file("result_CHLP_3index.txt","a+");

	
	// retrive solution values
	CPXgetmipobjval(env,lp,&UFLPrelval);
	//printf("MIP value: %.2f \n", UFLPrelval);
	

	//CPXgetbestobjval(env,lp,&value);  //best lower bound in case thew problem was not solve to optimality

    end = clock();
    cputime = (double)(end - start) / CLOCKS_PER_SEC;
   
	/*Retrieve the solutiona and calculate corresponding quadratic costs*/
	numcols=CPXgetnumcols(env,lp);
	localx=create_double_vector(numcols);
	local_assign=create_int_vector(NN);
	CPXgetmipx(env,lp,localx,0, numcols-1);  // obtain the values of the decision variables
	for(i=0;i<NN;i++){
		for(k=0;k<NN;k++){
			if (localx[pos_z[i][k]] > 0.5) {
				local_assign[i] = k;
				//printf("%d(%d) ", i+1, local_assign[i]+1);
			}
		}
	}
	//printf("\n");

	UBOptUFLP=UFLPrelval;
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			UBOptUFLP+=W[i][j]*(c_t[local_assign[i]][local_assign[j]]);
		}
	}
	UpperBound = UBOptUFLP;

	printf("Initial LB from SSCFLP: %.2f UB: %.2f open facilities: ", UFLPrelval, UpperBound);
	for (i = 0; i < NN; i++) {
		if (localx[pos_z[i][i]] > 0.5) {
			printf("%d (%.2f) ", i + 1, f[i][0]);
		}
	}
	printf("\n");
	printf("Time to solve SSCFLP: %.2f \n", cputime);

 //   index=0;
	//for (i = 0; i < NN; i++) {
	//	for (k = 0; k < NN; k++) {
	//		if (i == k) {
	//			if (initial_x[index++] > 0.5)
	//				open_plants[i] = 1;
	//			else
	//				open_plants[i] = 0;
	//		}
	//		else{
	//	 	  if (initial_x[index++] > 0.5) {
	//			best_assigmnent[i] = k;
	//			//printf("z[%d][%d] = 1 \n", i+1,k+1);
	//		  }
	//		}
	//	}
	//}

	/*for (i = 0; i < NN; i++) {
	   if (initial_x[pos_z[i][i]] > 0.5)
		  open_plants[i] = 1;
	   else
		  open_plants[i] = 0;
	}
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
		   if (initial_x[pos_z[i][k]] > 0.5){
			   best_assigmnent[i] = k;
		       break;
		   }
		}
	}*/

	//printf("Initial solution: ");
	//value = 0;

	//prevsols[countsols].preeta=0;     //Carlos's modification
	//prevsols[countsols].num_comb=0;    //Carlos's modification
	/*for (k = 0; k < NN; k++) {
		if (open_plants[k] > 0.5) {
			value += f[k][0];
			printf("%d ", k + 1);
		}
	}*/
	//	prevsols[countsols].indhub[k]=best_assigmnent[k];   //Carlos's modification
	//	prevsols[countsols].num_comb++;
	//	for (m = 0; m < NN; m++){
	//		value += W[k][m] * (c_c[k][best_assigmnent[k]] + c_t[best_assigmnent[k]][best_assigmnent[m]] + c_d[best_assigmnent[m]][m]);
	//		prevsols[countsols].preeta+=W[k][m]*c_t[best_assigmnent[k]][best_assigmnent[m]]; //calculating the preeta  //Carlos's modification
	//	}
	//}
	//prevsols[countsols].preeta-=combtol;
	//countsols++;         //Carlos's modification
	//prevsols[countsols].indhub=create_int_vector(NN);   //Carlos's modification
	//printf("\n");
	//printf("Initial Lowerbound: %lf \t Upperbound from SSCFLP: %lf\t Upperbound from optimal SSCFLP: %lf \n",UFLPrelval, UpperBound, UBOptUFLP);
	UBsearchUFLP=UpperBound;
	//getchar();

	//UpperBound = value;
	//index = 0;
	//for (i = 0; i < NN; i++) {
	//	for (k = 0; k < NN; k++) {
	//		if (best_assigmnent[i] == k)
	//			initial_x[index++] = 1;
	//		else
	//			initial_x[index++] = 0;
	//	}
	//}

     //getchar();
	 TERMINATE:

       /* Free the allocated vectors */


    if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
        fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
	}
    if ( env != NULL ) {
       status = CPXcloseCPLEX (&env);
       if ( status ) {
         char  errmsg[1024];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
       }
	}

	//free(x);

}

int CPXPUBLIC 
mycheckcallback (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               int        *useraction_p)
{
   int status = 0;



   
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
   double   epsilon_LP = 100.0;
   double *x;

   /***New modifications***/
   double locval;
   int *tempop_plants;
   int *temp_assign;
   int index;

   tempop_plants=create_int_vector(NN);
   temp_assign=create_int_vector(NN);
   x=create_double_vector(NN*NN);
   *useraction_p = CPX_CALLBACK_DEFAULT; 



   // Get current depth in BB
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
   if (status) {
	   fprintf(stdout, "Can't get depth for node.");
	   goto TERMINATE;
   }

	   
	   if (status) {
		   fprintf(stderr, "Failed to get node id.\n");
		   goto TERMINATE;
	   }

	   
		   status = CPXgetcallbacknodex (env, cbdata, wherefrom, x, 0, NN*NN-1);
		   for (i = 0; i < NN; i++) {
			   if (x[pos_z[i][i]] > 0.5)
				   tempop_plants[i] = 1;
			   else
				   tempop_plants[i] = 0;
		   }
		   for (i = 0; i < NN; i++) {
			   for (k = 0; k < NN; k++) {
				   if (x[pos_z[i][k]] > 0.5){
					   temp_assign[i] = k;
					   break;
				   }
			   }
		   }

		   status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, &objval);
		   /*printf("Initial solution: ");*/
		   locval= 0;

		   prevsols[countsols].preeta=0;     //Carlos's modification
		   prevsols[countsols].num_comb=0;    //Carlos's modification
		   for (k = 0; k < NN; k++) {
			   if (tempop_plants[k] >0.5){
				   locval += f[k][0];
				   /*printf("%d ", k + 1);*/
			   }
			   prevsols[countsols].indhub[k]=temp_assign[k];   //Carlos's modification
			   prevsols[countsols].num_comb++;
			   for (m = 0; m < NN; m++){
				   locval += W[k][m] * (c_c[k][temp_assign[k]] + c_t[temp_assign[k]][temp_assign[m]] + c_d[temp_assign[m]][m]);
				   prevsols[countsols].preeta+=W[k][m]*c_t[temp_assign[k]][temp_assign[m]]; //calculating the preeta  //Carlos's modification
			   }
		   }
		   prevsols[countsols].preeta-=combtol;
		   countsols++;         //Carlos's modification
		   prevsols[countsols].indhub=create_int_vector(NN);   //Carlos's modification
		   //printf("\n");
		   printf("MIP solution %lf---Evaluated a new upperbound from SSCFLP: %.12f \n", objval, locval);

		   if(UpperBound > locval){
			   printf("Updated new solution\n");
			   //getchar();
			   UpperBound=locval;
			   index = 0;
			   for (i = 0; i < NN; i++) {
				   best_assigmnent[i]=temp_assign[i];
				   open_plants[i]=tempop_plants[i];
				   for (k = 0; k < NN; k++) {				   
					   if (temp_assign[i] == k)
						   initial_x[index++] = 1;
					   else
						   initial_x[index++] = 0;
				   }
			   }
		   }

TERMINATE:

   return (status);

} /* END mycutcallback */



int Reassign_nodes(int *best_assigmnent1)
{
	int i, j, k, l, m, count;
	FILE     *out;
	clock_t  start, end;
	int index, index1, index11;  // indices auxiliares para rellenar las matrices
	double   cputime;
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each column ...............
	int       *matcnt; // number of non-zero element in each column... .................
	int       *matind; // associated row of each non-zelo element ......................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int *indices, *priority;
	int		nodecount, statusP;
	int cur_numcols;


	start = clock();

	objsen = 1; //min

	//Initialize CPLEX environment
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UHLPSA");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}


	//Define z_ik variables
	index1 = 0;  // index of columns
	numcols = NN*NN;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (i = 0; i<NN; i++){
		for (k = 0; k<NN; k++){
			if (i != k && open_plants[i] == 0 && open_plants[k] == 1){
				pos_z[i][k] = index1;
				obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
				ctype[index1] = 'B';
				lb[index1] = 0;
				ub[index1] = 1;
				index1++;
			}
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);


	//Add assignment constraints  \sum_{k \in N} z_ik = 1
	numrows = NN;
	numnz = NN*NN;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++){
		if (open_plants[i] == 0){
			sense[index1] = 'E';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			for (k = 0; k < NN; k++){
				if (i != k && open_plants[k] == 1){
					matind[index] = pos_z[i][k];
					matval[index++] = 1;
				}
			}
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add capacity constraints sum(i in NN) O_i z_ik <= b*z_kk

	numrows = NN;
	numnz = NN*(NN + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (k = 0; k<NN; k++){
		if (open_plants[k] == 1){
			sense[index1] = 'L';
			rhs[index1] = b[k][0] - O[k];
			matbeg[index1++] = index;
			for (i = 0; i < NN; i++){
				if (i != k && open_plants[i] == 0){
					matind[index] = pos_z[i][k];
					matval[index++] = O[i];
				}
			}
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//branch and bound parameters
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF); //output display
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4); //different levels of output display
	CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1);
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env, CPX_PARAM_TILIM, 500); // time limit
	CPXsetdblparam(env, CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//	CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, 0); //heuristic frequency and intensisty 
	//CPXoptimize(env,lp);   solve a linear program
	//CPXgetbestobjval(env,lp,&value);   \\obtain the LP relaxation bound
	//printf("LP lower bound: %f   ",value);
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,238016.277+0.001); // provide an initial upper bound
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables


	///* Turn on traditional search for use with control callbacks */
	//status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	//if ( status )  goto TERMINATE;

	CPXmipopt(env, lp);  //solve the integer program

	i = CPXgetstat(env, lp);
	/*if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);*/
	if (i == 101 || i == 102){
		// retrive solution values
		statusP = 1;
		numcols = CPXgetnumcols(env, lp);
		d_vector(&x, numcols, "open_cplex:9");
		CPXgetmipobjval(env, lp, &value);
		CPXgetmipx(env, lp, x, 0, numcols - 1);  // obtain the values of the decision variables
		for (i = 0; i < NN; i++) {
			if (open_plants[i] == 0){
				for (k = 0; k < NN; k++) {
					if (i != k && open_plants[k] == 1 && x[pos_z[i][k]] > 0.5){
						best_assigmnent1[i] = k;
						break;
					}
				}
			}
			else{
				best_assigmnent1[i] = i;
			}
		}
		free(x);
	}
	else
		statusP = 0;

TERMINATE:

	/* Free the allocated vectors */


	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	//free(x);
	return statusP;
}

