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
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600); // time limit
	CPXsetdblparam(env, CPX_PARAM_TRELIM, 14000); // B&B memory limit
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.0000001); // e-optimal solution (%gap)
											  //CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0); // add cuts or not
	//CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, 0.05);
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); //output display
	acumstatus += CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);											/* Do not use presolve */
	acumstatus += CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 					/* Turn on traditional search for use with control callbacks */
	acumstatus += CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);									 /* Let MIP callbacks work on the original model */



	return 0;
}

void Benders_BC(void)
{
	int i,j,k,l,m,r, count;
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
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int *indices, *priority;
	int		nodecount, num_z_var;
	int cur_numcols;
	double best_upper_bound;
	double best_lower_bound, coeff_z;
    CUTINFO cutinfo;
    cutinfo.x   = NULL;
    cutinfo.beg = NULL;
    cutinfo.ind = NULL; 
    cutinfo.val = NULL;
    cutinfo.rhs = NULL;

	priority = create_int_vector(NN*NN);
    indices = create_int_vector(NN*NN);
	x = create_double_vector(NN*NN+NN);
	
	start = clock();

	Define_Core_Point();
	//printf("So far added %d optimality guts\n", count_added);getchar();

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

	printf("C:");
	for (i = 0; i < NN; i++)
		printf("%d ", fixed_zero[i]);
	printf("\n O:");
	for (i = 0; i < NN; i++)
		printf("%d ", fixed_one[i]);

	printf("\n Cand_hubs:");
	for (i = 0; i < count_cand_hubs; i++)
		printf("%d ", cand_hubs[i]);


	                                        //Define z_ik variables
    index1 = 0;  // index of columns
	numcols = NN*NN;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");

	for (i = 0; i < NN; i++){
		if (fixed_zero[i] == 0){
			for (k = 0; k < NN; k++){
				if (fixed_zero[k] == 0){
					pos_z[i][k] = index1;
					if (i == k){
						obj[index1] = f[i][0];
						priority[index1] = 2;
					}
					else{
						obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
						//	 obj[index1] = 0;
						priority[index1] = 1;
					}
					ctype[index1] = 'B';
					lb[index1] = 0;
					ub[index1] = 1;
					indices[index1] = index1;
					index1++;
				}
			}
		}
		else{
			for (k = 0; k < NN; k++) {
				if (i != k && fixed_zero[k] == 0){
					pos_z[i][k] = index1;
					obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
					priority[index1] = 1;
					ctype[index1] = 'B';
					lb[index1] = 0;
					ub[index1] = 1;
					indices[index1] = index1;
					index1++;
				}
			}
		}
	}
	status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_z_var = index1;

	                                        //Define eta variable
	pos_eta = index1;
	index1 = 0;  // index of columns
	numcols = 1;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");

	//pos_eta = NN*NN + index1;
    obj[index1] = 1;
    ctype[index1] = 'C';
    lb[index1] = 0;
    ub[index1] = CPX_INFBOUND;
    index1++;

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
		  if (fixed_zero[k] == 0){
			  matind[index] = pos_z[i][k];
			  matval[index++] = 1;
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
	for (i = 0; i < NN; i++){
		for (k = 0; k < NN; k++){
			if (i != k && fixed_zero[k] == 0){
				sense[index1] = 'L';
				rhs[index1] = 0;
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

	if(hybrid==0 || hybrid==3){																//Add exactly p_hubs
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
			if(fixed_zero[i] == 0){
				matind[index] = pos_z[i][i];
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
			if (fixed_zero[k] == 0){
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
	/********************************************************** Adding initial set of optimality cuts****/
	numrows =count_added;
	numnz = count_added*(NN*NN+1);
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");
	index = 0;
    index1 = 0;
	for(r=0;r<count_added;r++){
		sense[index1]='G';
		rhs[index1]= initial_cuts[r].rhs[0];
		matbeg[index1++] = index;
		matind[index] = pos_eta;
		matval[index++] = 1;
		for (i = 0; i < NN; i++) {
			 if (fixed_zero[i] == 0){
				 for (k = 0; k < NN; k++) {
					 if (fixed_zero[k] == 0){
						 if(ABS(initial_cuts[r].origval[i*NN+k])>0.001)
						 {
							 matind[index]= pos_z[i][k];
							 matval[index++]=initial_cuts[r].origval[i*NN+k];
						 }
					 }
				 }
			 }
			 else{
				 for (k = 0; k < NN; k++) {
					 if (i != k && fixed_zero[k] == 0){
						 if (ABS(initial_cuts[r].origval[i*NN+k]) > 0.00001) {
							 matind[index]= pos_z[i][k];
							 matval[index++]=initial_cuts[r].origval[i*NN+k];
						 }
					 }
				 }

			 }
		 }

	}
	//printf("We've added %d of the old optimality cuts\n with %d constraints and %d non zeros", count_added, index1, index);getchar();
	status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if( status ) 
      fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	/***************************************************************************************************/ 
	   
											//Add fixed open hub constraints

	numrows = NN;
	numnz = NN;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (k = 0; k<NN; k++){
		if (fixed_one[k] == 1){
			sense[index1] = 'E';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			matind[index] = pos_z[k][k];
			matval[index++] = 1;
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


	//CPXwriteprob(env, lp, "BendersCHLPSA_MIP.lp", NULL);
	//branch and bound parameters
    CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
    //CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
    CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output display
	//CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1);
    //CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
    CPXsetdblparam(env,CPX_PARAM_TILIM,86400); // time limit
    CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
    CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000001); // e-optimal solution (%gap)
    //CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
    //CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
   // CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
    //CPXoptimize(env,lp);   solve a linear program
    //CPXgetbestobjval(env,lp,&value);   \\obtain the LP relaxation bound
    //printf("LP lower bound: %f   ",value);
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetintparam(env, CPX_PARAM_COVERS, 3); //strategy for adding cover cuts
	//CPXsetdblparam(env, CPX_PARAM_CUTUP, 251540.9); // provide an initial upper bound
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UpperBound + 0.001); // provide an initial upper bound
	CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	if(vers!=2) CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables

	status = CPXcopyorder(env, lp, num_z_var, indices, priority, NULL);
    if( status )
      fprintf (stderr,"CPXcopyorder failed.\n");

	 /* Set branch-and-cut parameters */

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

   cur_numcols = CPXgetnumcols (env, lp);

   cutinfo.lp = lp;
   cutinfo.numcols = cur_numcols;
   printf("Columns loaded in Cplex: %d \n", cur_numcols);
   cutinfo.x = (double *) malloc (cur_numcols * sizeof (double));
   if ( cutinfo.x == NULL ) {
      fprintf (stderr, "No memory for solution values.\n");
      goto TERMINATE;
   }

    /* Set up to use MIP callback */

     status = CPXsetusercutcallbackfunc (env, mycutcallback, &cutinfo);
	  if ( status )  goto TERMINATE;

	 status = CPXsetlazyconstraintcallbackfunc (env, mycutcallback, &cutinfo);
	  if ( status )  goto TERMINATE;

	  /* Code to use Heuristic Callback*/
	/************************************************/
	status = CPXsetheuristiccallbackfunc (env, Heur, NULL);
	/**********************************************************************************************************/

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
	CPXgetmipobjval(env,lp,&value);
	printf("Upper bound: %f   ",value);
	best_upper_bound=value;
    //fprintf(out," %.3f  ", value);
    // If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound

	CPXgetbestobjval(env,lp,&value);  //best lower bound in case thew problem was not solve to optimality
	best_lower_bound=value;
	printf("Lower bound: %f   ",value);
        
	nodecount = CPXgetnodecnt (env, lp);
	printf(" the number of BB nodes : %ld   ",nodecount);
	//add one line of code I will give you

    end = clock();
    cputime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Bender's Time: %.2f \n", cputime);


	CPXgetmipx(env,lp,x,0, cur_numcols-1);  // obtain the values of the decision variables

    index=0;

	out = open_file(output_text, "a+");
	fprintf(out, "%.2f;  %.2f; %.2f; %d; ", best_upper_bound, best_lower_bound,  cputime, nodecount);
	
	printf("Optimal set of hubs: ");
	fprintf(out, "hubs:");
	for (i = 0; i<NN; i++){
		if (x[pos_z[i][i]] > 0.5 && fixed_zero[i]==0) {
			printf("%d ", i + 1);
			fprintf(out,"%d ", i + 1);
		}
	}
	printf("\n");
	fprintf(out,";");
	fclose(out);


      
	 TERMINATE:

       /* Free the allocated vectors */

       free_and_null ((char **) &cutinfo.x);
       free_and_null ((char **) &cutinfo.beg);
       free_and_null ((char **) &cutinfo.ind);
       free_and_null ((char **) &cutinfo.val);
       free_and_null ((char **) &cutinfo.rhs);


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

	free(x);
	free(indices);
	free(priority);
	free(cand_hubs);
	free(fixed_one);
	free(fixed_zero);
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
