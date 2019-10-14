#include "def.h"
extern double** c, ** c_c, ** c_t, ** c_d, ** f, ** W, * O, * D, ** b;
extern double     collect, transfer, distribute;
extern int        NN, Q, p_hubs;
extern double     MAX_DOUBLE;
extern double     UpperBound;
extern coordinate* pts;
extern int** pos_z;
extern int        pos_eta;
extern double     old_lower_bound;
extern double*** alpha;
extern double*** beta;
extern double** core;
extern double     LP_lower_bound;
extern SCUTS      sepcut;
extern INCUTS* initial_cuts;
extern double     sum_core;
extern double     old_objval;
extern int        MG;
extern double* initial_x;
extern    int     count_same_node;
extern ZVAL* z_open;
extern ZVAL* z_closed;
extern int* cand_hubs;
extern int        count_cand_hubs;
extern int* fixed_zero;
extern int* fixed_one;
extern int        count_added;
extern COV        cover;
extern NORD* ord_nodes;
extern double     AggregatedDemand;
extern double* coeff_ES;
extern int Capacitated_instances;
//extern int        *best_assigmnent;



void Benders_root_node_heur(void)
{
	int i, j, k, l, m, count;
	FILE* out;
	clock_t  start, end, start_SP, end_SP;
	int index, index1, index11;  // indices auxiliares para rellenar las matrices
	double   cputime, cputime_SP;
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double* obj;    // objective function coefficients ..............................
	double* rhs;    // right and side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each column ...............
	int* matcnt; // number of non-zero element in each column... .................
	int* matind; // associated row of each non-zelo element ......................
	double* matval; // coefficient values fo the non-zero elements of constraints....
	double* lb;     // lower bounds of variables.....................................
	double* ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double* x;      // solution vector (double, even if the problem is integer) .....
	double* dj;
	char probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    valuef; //Value of objective function during Partial enumeration.......
	int		nodecount;
	int cur_numcols, cur_rows, count_fixed;
	double best_upper_bound;
	double best_lower_bound, coeff_z, lhs;
	double   tolerance_sep, tol_sep_cover, aa;
	int      terminate, flag_added, count_o, count_c, flag_fixed;
	double sum_O, max_O, sum_x, eval_cut;
	int is_cover, count_cover;
	int* cand_cover;
	double   EPSOBJ = 0.1;
	double   epsilon_LP;
	int      iter = 0;
	int count_cover_cuts = 0;
	int count_Fenchel_cuts = 0;
	int assign_fixed;

	/*******************************/
	cand_cover = create_int_vector(NN);
	coeff_ES = create_double_vector(NN);

	UpperBound = MAX_DOUBLE;
	MG = 1;
	count_added = 0;
	//tolerance_sep = -0.0001;
	tolerance_sep = 0.05;
	epsilon_LP = 0.005;
	tol_sep_cover = 0.001;
	cputime_SP = 0;
	initial_cuts = (INCUTS*)calloc(1000, sizeof(INCUTS));
	initial_cuts[count_added].sense = (char*)calloc(1, sizeof(char));
	initial_cuts[count_added].rhs = create_double_vector(1);
	initial_cuts[count_added].beg = create_int_vector(1);
	initial_cuts[count_added].ind = create_int_vector(NN * NN + 1);
	initial_cuts[count_added].val = create_double_vector(NN * NN + 1);
	initial_cuts[count_added].origval = create_double_vector(NN * NN);
	x = create_double_vector(NN * NN + NN);
	z_open = (ZVAL*)calloc(NN, sizeof(ZVAL));
	z_closed = (ZVAL*)calloc(NN, sizeof(ZVAL));
	cand_hubs = create_int_vector(NN);
	fixed_one = create_int_vector(NN);
	fixed_zero = create_int_vector(NN);
	count_cand_hubs = NN;
	/***Some precalculations*****/
	for (k = 0; k < NN; k++) {
		cand_hubs[k] = k;
		fixed_one[k] = 0;
		fixed_zero[k] = 0;
	}
	start = clock();
	Define_Core_Point();
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
	numcols = NN * NN;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");

	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			pos_z[i][k] = index1;
			if (i == k) {
				obj[index1] = f[i][0];
			}
			else {
				obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
			}
			lb[index1] = 0;
			ub[index1] = 1;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	//Define eta variable
	index1 = 0;  // index of columns
	numcols = 1;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");

	pos_eta = NN * NN + index1;
	obj[index1] = 1;
	lb[index1] = 0;
	ub[index1] = CPX_INFBOUND;
	index1++;

	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);

	//Add assignment constraints  \sum_{k \in N} z_ik = 1
	numrows = NN;
	numnz = NN * NN;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++) {
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (k = 0; k < NN; k++) {
			matind[index] = pos_z[i][k];
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


	//Add linking constraints  z_ik <= z_kk

	numrows = NN * (NN - 1);
	numnz = 2 * NN * (NN - 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (i != k) {
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
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add exactly p_hubs
	if (hybrid == 0 || hybrid == 3) {
		numrows = 1;
		numnz = NN;
		d_vector(&rhs, numrows, "open_cplex:2");
		c_vector(&sense, numrows, "open_cplex:3");
		i_vector(&matbeg, numrows, "open_cplex:4");
		i_vector(&matind, numnz, "open_cplex:6");
		d_vector(&matval, numnz, "open_cplex:7");
		index = 0;
		index1 = 0;
		sense[index1] = 'E';
		rhs[index1] = p_hubs;
		matbeg[index1++] = index;
		for (i = 0; i < NN; i++) {
			matind[index] = pos_z[i][i];
			matval[index++] = 1;
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


	//Add capacity constraints sum(i in NN) O_i z_ik <= b*z_kk
	if (Capacitated_instances == 1) {
		numrows = NN;
		numnz = NN * (NN + 1);
		d_vector(&rhs, numrows, "open_cplex:2");
		c_vector(&sense, numrows, "open_cplex:3");
		i_vector(&matbeg, numrows, "open_cplex:4");
		i_vector(&matind, numnz, "open_cplex:6");
		d_vector(&matval, numnz, "open_cplex:7");

		index = 0;
		index1 = 0;
		for (k = 0; k < NN; k++) {
			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_z[k][k];
			matval[index++] = -(b[k][0] - O[k]);
			for (i = 0; i < NN; i++) {
				if (i != k) {
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


	//CPXwriteprob(env, lp, "BendersCHLPSA_LP.lp", NULL);
	//CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_OFF); //output display
	//CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output display
	CPXsetdblparam(env, CPX_PARAM_TILIM, 86400); // time limit
   // CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetintparam(env,CPX_PARAM_NUMERICALEMPHASIS,1); //Numerical precision of the time
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	CPXsetintparam(env, CPX_PARAM_LPMETHOD, 4);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);

	cur_rows = CPXgetnumrows(env, lp);
	cur_numcols = CPXgetnumcols(env, lp);

	dj = create_double_vector(cur_numcols);
	d_vector(&rhs, 1, "open_cplex:2");
	c_vector(&sense, 1, "open_cplex:3");
	i_vector(&matbeg, 1, "open_cplex:4");
	i_vector(&matind, NN, "open_cplex:6");
	d_vector(&matval, NN, "open_cplex:7");
	count_fixed = 0;
	//Solve the LP relaxation for the first time
	status = CPXlpopt(env, lp);
	if (status) fprintf(stderr, "Failed to optimize LP.\n");
	CPXsolution(env, lp, &status, &value, x, NULL, NULL, dj);
	PopulateLPSupport(x);

	//Run heuristc using info from fractional solution x
	printf("Started solving heuristic\n");
	status = Construct_Feasible_Solution(x, dj);
	printf("Finished solving heuristic.\n");

	//Cutting plane algorithm for solving the root node
	do {
		/*for (i = 0; i < NN; i++) {
			for (k = 0; k < NN; k++) {
				if (x[pos_z[i][k]] > 0.001)
					printf("z(%d %d): %.2f \n", i, k, x[pos_z[i][k]]);
			}
		}*/
		//Simple variable fixing test with reduced cost coefficients
		printf("Starting elimination test\n");
		flag_fixed = 0;
		if (vers != 5) {
			assign_fixed = 0;
			for (i = 0; i < count_cand_hubs; i++) {
				if (value + dj[pos_z[cand_hubs[i]][cand_hubs[i]]] > UpperBound + 0.01) {
					//printf("Fix z[%d] = 0  obj+redcost:%.5f  UB:%.5f \n", cand_hubs[i] + 1, value + dj[pos_z[cand_hubs[i]][cand_hubs[i]]], UpperBound);
					fixed_zero[cand_hubs[i]] = 1;
					index = 0;
					index1 = 0;
					sense[index1] = 'E';
					rhs[index1] = 0;
					matbeg[index1++] = index;
					matind[index] = pos_z[cand_hubs[i]][cand_hubs[i]];
					matval[index++] = 1;
					status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
					if (status) fprintf(stderr, "CPXaddrows failed.\n");
					count_fixed++;
					flag_fixed = 1;
				}
				else {
					for (k = 0; k < NN; k++) {						
						if (value + dj[pos_z[k][cand_hubs[i]]] > UpperBound + 0.01) {
							assign_fixed++;
							index = 0;
							index1 = 0;
							sense[index1] = 'E';
							rhs[index1] = 0;
							matbeg[index1++] = index;
							matind[index] = pos_z[k][cand_hubs[i]];
							matval[index++] = 1;
							status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
						}
					}
				}
			}
			if (assign_fixed > 0) printf("Fixed %d assignment variables in this iteration\n", assign_fixed);

			if (flag_fixed == 1) { //If I fixed to 0 some hubs
				count_cand_hubs = 0;
				for (k = 0; k < NN; k++) {
					if (fixed_zero[k] == 0)
						cand_hubs[count_cand_hubs++] = k;
				}
				printf("Elimination test performed \n Fixed hubs: %d Remaining hubs: %d \n", count_fixed, count_cand_hubs);
				Define_Core_Point();
			}
		}
		printf("Finished elimination test.\n");

		count_c = 0;
		for (k = 0; k < NN; k++) {
			if (fixed_zero[k] == 0 && x[pos_z[k][k]] <= 0.2) {
				z_closed[count_c].k = k;
				z_closed[count_c++].value = f[k][0];
			}
		}
		//Partial Enumeration phase: solve LPs to try to remove potential locations

		if (iter>0 /*&& iter % 3 == 0*/ && vers != 4 && ((UpperBound - value) / UpperBound * 100< 2.0 && flag_fixed==1) ){
			printf("Entered partial enumeration phase\n");
			flag_fixed = 0;
			qsort((ZVAL*)z_closed, count_c, sizeof(z_closed[0]), Comparevalue_zc);
			cur_rows = CPXgetnumrows(env, lp);
			for (k = 0; k < NN; k++) {  			// Temporarily fix z_kk = 1 to try to permantely close it (i.e. z_k = 0 from now on)
				index = 0;
				index1 = 0;
				sense[index1] = 'E';
				rhs[index1] = 1;
				matbeg[index1++] = index;
				matind[index] = pos_z[z_closed[k].k][z_closed[k].k];
				matval[index++] = 1;
				status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
				if (status) fprintf(stderr, "CPXaddrows failed.\n");
				cur_rows++;
				//CPXwriteprob(env, lp, "BendersPE.lp", NULL);
				status = CPXlpopt(env, lp);
				if (status) fprintf(stderr, "Failed to optimize LP.\n");
				CPXgetobjval(env, lp, &valuef);
				if (status) fprintf(stderr, "Failed to getx LP.\n");
				CPXdelrows(env, lp, cur_rows - 1, cur_rows - 1);
				if (status) fprintf(stderr, "CPXdelrows failed.\n");
				cur_rows--;
				//	CPXwriteprob(env, lp, "BendersPE2.lp", NULL);
				if (valuef > UpperBound) {
					//printf("Fix z[%d] = 0 \n", z_closed[k].k + 1);
					fixed_zero[z_closed[k].k] = 1;
					index = 0;
					index1 = 0;
					sense[index1] = 'E';
					rhs[index1] = 0;
					matbeg[index1++] = index;
					matind[index] = pos_z[z_closed[k].k][z_closed[k].k];
					matval[index++] = 1;
					status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
					if (status) fprintf(stderr, "CPXaddrows failed.\n");
					cur_rows++;
					count_fixed++;
					flag_fixed = 1;
					//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
				}
			}
			printf("Finished partial enumeration phase\n");
			if (flag_fixed == 1) {
				count_cand_hubs = 0;
				for (k = 0; k < NN; k++) {
					if (fixed_zero[k] == 0)
						cand_hubs[count_cand_hubs++] = k;
				}
				printf("Partial enumeration performed \n Fixed hubs: %d Remaining hubs: %d \n", count_fixed, count_cand_hubs);
				Define_Core_Point();
				//goto EVALUATE;
			}
		}
		if (flag_fixed == 1) {
			status = CPXlpopt(env, lp);
			if (status) fprintf(stderr, "Failed to optimize LP.\n");
			CPXsolution(env, lp, &status, &value, x, NULL, NULL, dj);
		}

		//Solve the primal subproblem to find Benders cuts
		start_SP = clock();
		iter++;
		flag_added = 0;
		index = 0;
		index1 = 0;

		printf("Starting to separate Benders cuts\n");
		if (vers != 3)
			Update_Core_Point(x);
		for (i = 0; i < NN; i++) {
			for (j = i; j < NN; j++) {
				Check_CP_MW(x, i, j);
				status = NetFlow_TP(x, i, j);
			}
		}
		printf("Finished separating Benders cuts.");
		initial_cuts[count_added].sense[index1] = 'G';
		initial_cuts[count_added].rhs[index1] = 0;
		initial_cuts[count_added].beg[index1++] = index;
		initial_cuts[count_added].ind[index] = pos_eta;
		initial_cuts[count_added].val[index++] = 1;
		lhs = x[pos_eta];
		// printf("eta >= ");
		for (i = 0; i < NN; i++) {
			if (fixed_zero[i] == 0) {
				for (k = 0; k < NN; k++) {
					if (fixed_zero[k] == 0) {
						coeff_z = 0;
						for (j = 0; j < NN; j++)
							coeff_z += (-alpha[i][j][k] - beta[j][i][k]);
						if (ABS(coeff_z) > 0.00001) {
							initial_cuts[count_added].ind[index] = pos_z[i][k];
							initial_cuts[count_added].val[index++] = coeff_z;
							initial_cuts[count_added].origval[pos_z[i][k]] = coeff_z;
							lhs += coeff_z * x[pos_z[i][k]];
							//if (x[pos_z[i][k]] > 0.00001)
							 //   printf("%.2f z[%d][%d] + ", coeff_z, i + 1, k + 1);
						}
					}
				}
			}
			else {
				for (k = 0; k < NN; k++) {
					if (i != k && fixed_zero[k] == 0) {
						coeff_z = 0;
						for (j = 0; j < NN; j++)
							coeff_z += (-alpha[i][j][k] - beta[j][i][k]);
						if (ABS(coeff_z) > 0.00001) {
							initial_cuts[count_added].ind[index] = pos_z[i][k];
							initial_cuts[count_added].val[index++] = coeff_z;
							initial_cuts[count_added].origval[pos_z[i][k]] = coeff_z;
							lhs += coeff_z * x[pos_z[i][k]];
							//if (x[pos_z[i][k]] > 0.00001)
							 //   printf("%.2f z[%d][%d] + ", coeff_z, i + 1, k + 1);
						}
						else initial_cuts[count_added].origval[pos_z[i][k]] = 0;
					}
				}
			}
		}
		//  printf("\n");
		end_SP = clock();
		cputime_SP += (double)(end_SP - start_SP) / CLOCKS_PER_SEC;
		end = clock();
		//printf("Finished separating more optimality cuts took %lf secs\n",(double)(end_SP - start_SP) / CLOCKS_PER_SEC);
		cputime = (double)(end - start) / CLOCKS_PER_SEC;
		printf("iter:%d LP bound: %.2f gap: %.2f viol: %.2f viol rel %.2f time:%.2f sec SP time: %.2f (%.2f per) \n", iter, value, (UpperBound - value) / UpperBound * 100, lhs, (-lhs / value) * 100, cputime, cputime_SP, cputime_SP / cputime * 100);
		if ((-lhs / value) * 100 < tolerance_sep && ((UpperBound - value) / UpperBound) * 100 < 5.0) {
			flag_added = 0;
		}
		else flag_added = 1;

		if (lhs > -0.001)
			flag_added = 0;

		status = CPXaddrows(env, lp, 0, index1, index, initial_cuts[count_added].rhs, initial_cuts[count_added].sense, initial_cuts[count_added].beg, initial_cuts[count_added].ind, initial_cuts[count_added].val, NULL, NULL);
		if (status)  fprintf(stderr, "CPXaddrows failed.\n");
		initial_cuts[count_added].index = index;
		count_added++;
		initial_cuts[count_added].sense = (char*)calloc(1, sizeof(char));
		initial_cuts[count_added].rhs = create_double_vector(1);
		initial_cuts[count_added].beg = create_int_vector(1);
		initial_cuts[count_added].ind = create_int_vector(NN * NN + 1);
		initial_cuts[count_added].origval = create_double_vector(NN * NN);
		initial_cuts[count_added].val = create_double_vector(NN * NN + 1);

		//Solve again and if there is a different support, then find another feasible solution.
	//EVALUATE:
		status = CPXlpopt(env, lp);
		if (status) fprintf(stderr, "Failed to optimize LP.\n");
		CPXsolution(env, lp, &status, &value, x, NULL, NULL, dj);
		if (CompareLPSupport(x) >= 4) {			
			printf("Started running heuristic\n");
			status = Construct_Feasible_Solution(x, dj);
			PopulateLPSupport(x);
			printf("Finished running heuristic\n");
		}

		//Should there be a second round or not?
		if (flag_added != 0 && (UpperBound - value) / UpperBound > epsilon_LP) {
			terminate = 0;
			old_objval = value;
		}
		else {
			terminate = 1;
		}
		printf("Finished iteration %d.\n***************************************************\n", iter);

	} while (terminate != 1);

	end = clock();
	cputime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Root node LP bound: %.2f \n Root node Time: %.2f \n Num Benders cuts: %d \n  Num cover cuts: %d Num Fenchel cuts: %d \n", value, cputime, count_added, count_cover_cuts, count_Fenchel_cuts);

	count_o = 0;
	for (k = 0; k < NN; k++) {
		if (fixed_zero[k] == 0 && x[pos_z[k][k]] > 0.2) {
			z_open[count_o].k = k;
			z_open[count_o++].value = x[pos_z[k][k]];
		}
	}
	qsort((ZVAL*)z_open, count_o, sizeof(z_open[0]), Comparevalue_zo);

	count_c = 0;
	for (k = 0; k < NN; k++) {
		if (fixed_zero[k] == 0 && x[pos_z[k][k]] <= 0.2) {
			z_closed[count_c].k = k;
			z_closed[count_c++].value = f[k][0];
		}
	}
	qsort((ZVAL*)z_closed, count_c, sizeof(z_closed[0]), Comparevalue_zc);


	if(vers!=4){														//Partial Enumeration phase part I: Try to permanently open a set of hubs
		printf("Starting Partial Enumeration: \n");
		cur_rows = CPXgetnumrows(env, lp);
		for (k = 0; k < count_o; k++){  						// Temporarily fix z_kk = 0 to try to permantely open it (i.e. z_k = 1 from now on)
			index = 0;
			index1 = 0;
			sense[index1] = 'E';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_z[z_open[k].k][z_open[k].k];
			matval[index++] = 1;
			status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if (status) fprintf(stderr, "CPXaddrows failed.\n");
			cur_rows++;
			//	CPXwriteprob(env, lp, "BendersPE.lp", NULL);
			status = CPXlpopt(env, lp);
			if (status) fprintf(stderr, "Failed to optimize LP.\n");
			CPXgetobjval(env, lp, &value);
			if (status) fprintf(stderr, "Failed to getx LP.\n");
			CPXdelrows(env, lp, cur_rows - 1, cur_rows - 1);
			if (status) fprintf(stderr, "CPXdelrows failed.\n");
			cur_rows--;
			//	CPXwriteprob(env, lp, "BendersPE2.lp", NULL);
			if (value > UpperBound){
				//printf("Fix z[%d] = 1 \n", z_open[k].k+1);
				fixed_one[z_open[k].k] = 1;
				index = 0;
				index1 = 0;
				sense[index1] = 'E';
				rhs[index1] = 1;
				matbeg[index1++] = index;
				matind[index] = pos_z[z_open[k].k][z_open[k].k];
				matval[index++] = 1;
				status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
				if (status) fprintf(stderr, "CPXaddrows failed.\n");
				cur_rows++;
				count_fixed++;
				//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
			}
		}
		printf("Ended phase fix 0\n Started phase fix 1\n");
		
		//Partial Enumeration phase part II: Try to permanently close a set of hubs
		for (k = 0; k < count_c; k++){  						// Temporarily fix z_kk = 1 to try to permantely open it (i.e. z_k = 0 from now on)
			index = 0;
			index1 = 0;
			sense[index1] = 'E';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			matind[index] = pos_z[z_closed[k].k][z_closed[k].k];
			matval[index++] = 1;
			status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if (status) fprintf(stderr, "CPXaddrows failed.\n");
			cur_rows++;
			//CPXwriteprob(env, lp, "BendersPE.lp", NULL);
			status = CPXlpopt(env, lp);
			if (status) fprintf(stderr, "Failed to optimize LP.\n");
			CPXgetobjval(env, lp, &value);
			if (status) fprintf(stderr, "Failed to getx LP.\n");
			CPXdelrows(env, lp, cur_rows - 1, cur_rows - 1);
			if (status) fprintf(stderr, "CPXdelrows failed.\n");
			cur_rows--;
			//	CPXwriteprob(env, lp, "BendersPE2.lp", NULL);
			if (value > UpperBound){
				//printf("Fix z[%d] = 0 \n", z_closed[k].k+1);
				fixed_zero[z_closed[k].k] = 1;
				index = 0;
				index1 = 0;
				sense[index1] = 'E';
				rhs[index1] = 0;
				matbeg[index1++] = index;
				matind[index] = pos_z[z_closed[k].k][z_closed[k].k];
				matval[index++] = 1;
				status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
				if (status) fprintf(stderr, "CPXaddrows failed.\n");
				cur_rows++;
				count_fixed++;
				//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
			}
		}
		printf("Ended phase fix 1\n");
		
		count_cand_hubs = 0;
		for (k = 0; k < NN; k++){
			if (fixed_zero[k] == 0)
				cand_hubs[count_cand_hubs++] = k;
		}
		printf("Partial enumeration performed \n Fixed hubs: %d Remaining hubs: %d \n", count_fixed, count_cand_hubs);
	}
	status = CPXlpopt(env, lp);
	if (status) fprintf(stderr, "Failed to optimize LP.\n");
	CPXgetobjval(env, lp, &value);
	if (status) fprintf(stderr, "Failed to getx LP.\n");
	end = clock();
	cputime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Root node LP bound after PE: %.2f \n PE Time: %.2f \n Hubs fixed: %d \n", value, cputime, count_fixed);

	out = open_file(output_text, "a+");
	fprintf(out, "%.2f ; %.2f; %d; %d; %d; ", value, cputime, count_fixed,iter,count_added);
	fclose(out);

	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
      
	TERMINATE:

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
	free(dj);
	free(z_open);
	free(z_closed);
	free(cand_cover);
	free(coeff_ES);
}


void PopulateLPSupport(double* x) {
	int k;
	for (k = 0; k < NN; k++) {
		if (x[pos_z[k][k]] > 0.001) {
			FlagHubLPsupport[k] = 1;
		}
		else {
			FlagHubLPsupport[k] = 0;
		}
	}
}

int CompareLPSupport(double* x) {
	int k, retval=0;

	for (k = 0; k < NN; k++) {
		if ((x[pos_z[k][k]] < 0.001 && FlagHubLPsupport[k] == 1)) {
			retval++;
		}
	}

	return retval;
}
