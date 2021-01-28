#include "def.h"

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
	double   eta_cost = 0;
	int      iter = 0;
	int count_cover_cuts = 0;
	int count_Fenchel_cuts = 0;
	int assign_fixed, total_assign_fixed;
	int* indices, * priority;
	int *beg, *varindices, *effortlevel;
	double  *values;
	CUTINFO cutinfo;
	cutinfo.x = NULL;
	cutinfo.beg = NULL;
	cutinfo.ind = NULL;
	cutinfo.val = NULL;
	cutinfo.rhs = NULL;
	//For the bound changes
	char* btype;
	int* bind;
	double* realub;
	i_vector(&bind, 1, "open_cplex:4");
	c_vector(&btype, 1, "open_cplex:3");
	d_vector(&realub, 1, "open_cplex:3");
	realub[0] = 0;
	btype[0] = 'U';
	priority = create_int_vector(NN * NN);
	indices = create_int_vector(NN * NN);
	//Declaring the global types
	glob_numcols = NN * NN + 1;
	globvarind = create_int_vector(glob_numcols);  //Array containing all the Types for the complete model
	c_vector(&globvarctype, glob_numcols, "open_cplex:01");
	/*******************************/
	cand_cover = create_int_vector(NN);
	coeff_ES = create_double_vector(NN);
	srand(123456789);
	MG = 1;
	count_added = 0;
	//tolerance_sep = -0.0001;
	tolerance_sep = 0.50;
	epsilon_LP = 0.005;
	//epsilon_LP = 0.0000005;
	tol_sep_cover = 0.001;
	cputime_SP = 0;
	total_assign_fixed = 0;
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
				priority[index1] = 2;
			}
			else {
				obj[index1] = (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
				priority[index1] = 1;
			}
			if (initial_x[index1] == 1)
				eta_cost += obj[index1];
			lb[index1] = 0;
			ub[index1] = 1;
			indices[index1] = index1;
			globvarind[index1] = index1;
			globvarctype[index1] = 'B';
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	//The type of the continuous variable
	globvarind[index1] = index1;
	globvarctype[index1++] = 'C';
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
	printf("Initial LP bound: %.2lf\n", value);
	PopulateLPSupport(x);
	//Run heuristc using info from fractional solution x
	if (FlagHeuristic) {
		printf("Started solving heuristic\n");
		status = Construct_Feasible_Solution(x, dj);
		printf("Finished solving heuristic.\n");
	}
	//Cutting plane algorithm for solving the root node
	do {
		//Simple variable fixing test with reduced cost coefficients
		printf("Starting elimination test\n");
		flag_fixed = 0;
		assign_fixed = 0;
		if (vers != 5) {			
			for (i = 0; i < count_cand_hubs; i++) {
				if (value + dj[pos_z[cand_hubs[i]][cand_hubs[i]]] > UpperBound + 0.01) {
					//printf("Fix z[%d] = 0  obj+redcost:%.5f  UB:%.5f \n", cand_hubs[i] + 1, value + dj[pos_z[cand_hubs[i]][cand_hubs[i]]], UpperBound);
					fixed_zero[cand_hubs[i]] = 1;
					bind[0] = pos_z[cand_hubs[i]][cand_hubs[i]];
					realub[0] = 0;
					btype[0] = 'U';
					status = CPXchgbds(env, lp, 1, bind, btype, realub);
					if (status) fprintf(stderr, "CPXchange bounds failed.\n");
					count_fixed++;
					flag_fixed = 1;
					not_eligible_hub[cand_hubs[i]][cand_hubs[i]] = 1;
					eli_per_com[cand_hubs[i]]--;
					for (k = 0; k < NN; k++) {
						if (not_eligible_hub[k][cand_hubs[i]] == 0) {
							bind[0] = pos_z[k][cand_hubs[i]];
							realub[0] = 0;
							btype[0] = 'U';
							status = CPXchgbds(env, lp, 1, bind, btype, realub);
							not_eligible_hub[k][cand_hubs[i]] = 1;
							eli_per_com[k]--;
							assign_fixed++;
						}
					}
				}
				else {
					for (k = 0; k < NN; k++) {
						if (not_eligible_hub[k][cand_hubs[i]] == 0 && value + dj[pos_z[k][cand_hubs[i]]] > UpperBound + 0.01) {
							assign_fixed++;
							bind[0] = pos_z[k][cand_hubs[i]];
							realub[0] = 0;
							btype[0] = 'U';
							CPXchgbds(env, lp, 1, bind, btype, realub);
							not_eligible_hub[k][cand_hubs[i]] = 1;
							eli_per_com[k]--;
						}
					}
				}
			}
			if (assign_fixed > 0) {
				total_assign_fixed += assign_fixed;
				printf("Fixed %d assig variables in current iter, total fixed %d, (%.2f percentage) \n", assign_fixed, total_assign_fixed, (float) total_assign_fixed/(NN*NN)*100);
			}
			if (flag_fixed == 1) { //If I fixed to 0 some hubs
				count_cand_hubs = 0;
				for (k = 0; k < NN; k++) {
					if (fixed_zero[k] == 0)
						cand_hubs[count_cand_hubs++] = k;
				}
				printf("Elimination test performed \n Fixed hubs: %d Remaining hubs: %d \n", count_fixed, count_cand_hubs);
			}
		}
		printf("Finished elimination test.\n");
		//Partial Enumeration phase: solve LPs to try to remove potential locations
		if (((UpperBound - value) / UpperBound * 100 < 1.0 && flag_fixed >=0)) {
			count_c = 0;
			for (k = 0; k < count_cand_hubs; k++) {
				if (fixed_zero[cand_hubs[k]] == 0 && x[pos_z[cand_hubs[k]][cand_hubs[k]]] <= 0.2) {
					z_closed[count_c].k = cand_hubs[k];
					z_closed[count_c].value = f[cand_hubs[k]][0] * (1 - x[pos_z[cand_hubs[k]][cand_hubs[k]]]);
					for (i = 0; i < NN; i++) {
						z_closed[count_c].value += (O[i] * c_c[i][cand_hubs[k]] + D[i] * c_d[i][cand_hubs[k]]) * (1 - x[pos_z[cand_hubs[k]][cand_hubs[k]]]);
					}
					count_c++;
				}
			}
			printf("Entered partial enumeration phase to check %d of them\n", count_c);
			flag_fixed = 0;
			qsort((ZVAL*)z_closed, count_c, sizeof(z_closed[0]), Comparevalue_zc);
			cur_rows = CPXgetnumrows(env, lp);
			for (k = 0; k < count_c; k++) {  			// Temporarily fix z_kk = 1 to try to permantely close it (i.e. z_k = 0 from now on)
				printf("k:%d NLB %lf, eps: %.3f ", z_closed[k].k, z_closed[k].value, (value + dj[pos_z[cand_hubs[k]][cand_hubs[k]]])/UpperBound);
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
					printf("Fixed z[%d] = 0, NLB/UB: %.3f, (NLB - LB + RC)/UB: %.3f \n", z_closed[k].k + 1, valuef/UpperBound, (valuef - value + dj[pos_z[cand_hubs[k]][cand_hubs[k]]])/UpperBound);
					fixed_zero[z_closed[k].k] = 1;
					bind[0] = pos_z[z_closed[k].k][z_closed[k].k];
					realub[0] = 0;
					btype[0] = 'U';
					CPXchgbds(env, lp, 1, bind, btype, realub);
					count_fixed++;
					flag_fixed = 1;
					for (i = 0; i < NN; i++) {
						if (not_eligible_hub[i][z_closed[k].k] == 0) {
							not_eligible_hub[i][z_closed[k].k] = 1;
							eli_per_com[i]--;
						}
					}
					//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
				}
				else {
					printf("\n");
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
				//Define_Core_Point();
				//goto EVALUATE;
			}
		}
		if (flag_fixed == 1 || assign_fixed > 0) {
			Define_Core_Point();
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
		if (vers != 3) {
			Define_Core_Point();
			Update_Core_Point(x);
		}
		for (i = 0; i < NN-1; i++) {
			for (j = i+1; j < NN; j++) {
				//Check_CP_MW(x, i, j);				
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
		//printf("iter:%d LP bound: %.2f gap: %.2f viol: %.2f viol rel %.2f time:%.2f sec SP time: %.2f (%.2f per) \n", iter, value, (UpperBound - value) / UpperBound * 100, lhs, (-lhs / value) * 100, cputime, cputime_SP, cputime_SP / cputime * 100);
		if ((-lhs / value) * 100 < tolerance_sep && ((UpperBound - value) / UpperBound) * 100 < 10.0) {
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
		//Define_Core_Point();
		status = CPXlpopt(env, lp);
		if (status) fprintf(stderr, "Failed to optimize LP.\n");
		CPXsolution(env, lp, &status, &value, x, NULL, NULL, dj);
		printf("iter:%d LP bound: %.2f gap: %.2f viol: %.2f viol rel %.2f time:%.2f sec SP time: %.2f (%.2f per) \n", iter, value, (UpperBound - value) / UpperBound * 100, lhs, (-lhs / value) * 100, cputime, cputime_SP, cputime_SP / cputime * 100);
		if (FlagHeuristic) {
			if (CompareLPSupport(x) >= 2) {
				printf("Started running Matheuristic\n");
				status = Construct_Feasible_Solution(x, dj);
				PopulateLPSupport(NULL);
				printf("Finished running Matheuristic\n");
			}
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
	flag_fixed = 0;
	count_o = 0;
	for (k = 0; k < NN; k++) {
		if (fixed_zero[k] == 0 && x[pos_z[k][k]] > 0.8) {
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
				bind[0] = pos_z[z_open[k].k][z_open[k].k];
				btype[0] = 'B';
				realub[0] = 1;
				CPXchgbds(env, lp, 1, bind, btype, realub);				
				fixed_one[z_open[k].k] = 1;
				if (status) fprintf(stderr, "CPXaddrows failed.\n");
				count_fixed++;
				//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
			}
		}
		printf("Ended phase fix 0\n Started phase fix 1\n");
		//Partial Enumeration phase part II: Try to permanently close a set of hubs
		//printf("we will look at %d hubs\n", count_c); getchar();
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
				bind[0] = pos_z[z_closed[k].k][z_closed[k].k];
				btype[0] = 'U';
				realub[0] = 0;
				CPXchgbds(env, lp, 1, bind, btype, realub);
				flag_fixed = 1;
				count_fixed++;
				for (i = 0; i < NN; i++) {
					if (not_eligible_hub[i][z_closed[k].k] == 0) {
						not_eligible_hub[i][z_closed[k].k] = 1;
						//if(i==74)printf("Remove %d for hub %d from PE test\n", i, z_closed[k].k + 1);
						eli_per_com[i]--;
					}
					else {
						//if(i==74)printf("assigned %d for hub %d is already not eligible with value %d\n", i, z_closed[k].k+1, not_eligible_hub[i][z_closed[k].k]);
					}
				}				
				//	CPXwriteprob(env, lp, "BendersPE3.lp", NULL);
			}
		}
		printf("Ended phase fix 1\n");
		//printf("Node 74 has %d eligible hubs\n", eli_per_com[74]);
		/*for (i = 0; i < NN; i++) {
			if (not_eligible_hub[74][i] == 0) {
				printf("Hub %d is the only eligible hub\n", i);
			}
		}*/
		if (flag_fixed == 1) {
			count_cand_hubs = 0;
			for (k = 0; k < NN; k++) {
				if (fixed_zero[k] == 0)
					cand_hubs[count_cand_hubs++] = k;
			}
		}
		Define_Core_Point();
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
	fprintf(out, "%.2lf;%.2f;%.2f;%d;%d;", UpperBound, value, cputime,iter, count_fixed);
	fclose(out);
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	//Now solving the Integer Program
	 /******************************************************************************************/
	SetBranchandCutParam(env, lp);
	if (vers != 2) CPXsetintparam(env, CPX_PARAM_MIPORDIND, CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	status = CPXcopyorder(env, lp, NN*NN, indices, priority, NULL);
	cutinfo.lp = lp;
	cutinfo.numcols = cur_numcols;
	printf("Columns loaded in Cplex: %d \n", cur_numcols);
	cutinfo.x = (double*)malloc(cur_numcols * sizeof(double));
	if (cutinfo.x == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	/* Set up to use MIP callback */
	status = CPXsetusercutcallbackfunc(env, mycutcallback, &cutinfo);
	if (status)  goto TERMINATE;
	status = CPXsetlazyconstraintcallbackfunc(env, mycutcallback, &cutinfo);
	if (status)  goto TERMINATE;
	/* Code to use Heuristic Callback*/
  /************************************************/
	status = CPXsetheuristiccallbackfunc(env, Heur, NULL);
	// Add best solution from Matheuristic to MIP
	index = 0;
	eta_cost = 0;
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (best_sol_assignments[i] == k) {
				initial_x[index++] = 1;
				if (i == k)
					eta_cost += f[i][0];
				else
					eta_cost += (O[i] * c_c[i][k] + D[i] * c_d[i][k]);
			}
			else
				initial_x[index++] = 0;
		}
	}
	i_vector(&beg, 1, "open_cplex:4");
	i_vector(&effortlevel, 1, "open_cplex:4");
	i_vector(&varindices, cur_numcols, "open_cplex:6");
	d_vector(&values, cur_numcols, "open_cplex:7");
	beg[0] = 0;
	effortlevel[0] = 5;
	for (i = 0; i < cur_numcols-1; i++) {
		values[i] = (double) initial_x[i];
		//values[i] = 0;
		varindices[i] = i;
	//	if (values[i] > 0.001)
		//	printf("%.2f, %d \n", values[i], varindices[i]);
	}
	values[cur_numcols - 1] = (double) (UpperBound - eta_cost);
	//values[cur_numcols - 1] = UpperBound;
	varindices[cur_numcols - 1] = cur_numcols - 1;
	printf("glonumcols: %d cur_cols: %d (%d), %.2f, %.2f \n", glob_numcols, cur_numcols, NN*NN, UpperBound, eta_cost);
	status = CPXaddmipstarts(env, lp, 1, cur_numcols, beg, varindices, values, effortlevel, NULL);
	printf("status: %d \n", status);
	free(beg);
	free(effortlevel);
	free(varindices);
	free(values);
	CPXmipopt(env, lp);  //solve the integer program
	i = CPXgetstat(env, lp);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);
	// out = open_file("result_CHLP_3index.txt","a+");
	// retrive solution values
	CPXgetmipobjval(env, lp, &value);
	printf("Upper bound: %f   ", value);
	best_upper_bound = value;
	//fprintf(out," %.3f  ", value);
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env, lp, &value);  //best lower bound in case thew problem was not solve to optimality
	best_lower_bound = value;
	printf("Lower bound: %f   ", value);
	nodecount = CPXgetnodecnt(env, lp);
	printf(" the number of BB nodes : %ld   ", nodecount);
	//add one line of code I will give you
	end = clock();
	cputime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Bender's Time: %.2f \n", cputime);
	CPXgetmipx(env, lp, x, 0, cur_numcols - 1);  // obtain the values of the decision variables
	index = 0;
	out = open_file(output_text, "a+");
	fprintf(out, "%.2f;%.2f;%.2f;%.2lf;%d; ", best_upper_bound, best_lower_bound,  cputime, (best_upper_bound-best_lower_bound)*100/ best_upper_bound, nodecount);
	printf("Optimal set of hubs: ");
	fprintf(out, "hubs:");
	for (i = 0; i < NN; i++) {
		if (x[pos_z[i][i]] > 0.5) {
			printf("%d ", i + 1);
			fprintf(out, "%d ", i + 1);
		}
	}
	printf("eta: %.2f \n", x[cur_numcols - 1]);
	fprintf(out, ";");
	fclose(out);
TERMINATE:
	/* Free the allocated vectors */
	free_and_null((char**)& cutinfo.x);
	free_and_null((char**)& cutinfo.beg);
	free_and_null((char**)& cutinfo.ind);
	free_and_null((char**)& cutinfo.val);
	free_and_null((char**)& cutinfo.rhs);
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
	free(globvarctype);
	free(globvarind);
	free(indices);
	free(priority);
	//free(cand_hubs);
	//free(fixed_one);
	//free(fixed_zero);
}

void PopulateLPSupport(double* x) {
	int k;
	double val;
	for (k = 0; k < NN; k++) {
		if (x == NULL) val = best_sol_facilities[k];
		else val = x[pos_z[k][k]];
		if ( val > 0.001) {
			FlagHubLPsupport[k] = 1;
			//printf("Found one\n");
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

int Comparevalue_zo(const void* a, const void* b)
{
	if (((ZVAL*)a)->value < ((ZVAL*)b)->value)
		return 1;
	if (((ZVAL*)a)->value > ((ZVAL*)b)->value)
		return -1;
	return 0;
}

int Comparevalue_zc(const void* a, const void* b)
{
	if (((ZVAL*)a)->value < ((ZVAL*)b)->value)
		return 1;
	if (((ZVAL*)a)->value > ((ZVAL*)b)->value)
		return -1;
	return 0;
}

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
	CPXsetintparam(env, CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
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

int CPXPUBLIC mycutcallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	int status = 0;
	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	int      numcols = cutinfo->numcols;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	int* cutind = NULL;
	double* cutval = NULL;
	int* feas = NULL;
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
	int		Flag_Integer = 0;
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
	if (oldnodeid == cutinfo->nodeid) {
		count_same_node++;
	}
	else {
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
	if (status) {
		fprintf(stdout, "Can't get depth for node.");
		goto TERMINATE;
	}
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	// status = CPXgetcallbackgloballb (env, cbdata, wherefrom, &objval, 0, 0); 
	if (status) {
		fprintf(stdout, "Can't get objective value for node.");
		goto TERMINATE;
	}
	// We must do something in every integer solution:
	tolerance_sep = -100.0;
	//tolerance_sep = -0.01;
	flag_solve_SP = 0;
	if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS) {
		flag_solve_SP = 1;
		tolerance_sep = -0.000000000001;
		Flag_Integer = 1;
		//tolerance_sep = -0.001;
	}
	else {
		// printf("imp rel: %.6f objval: %.2f old_objval: %.2f \n", (double) ABS(objval - old_objval)/objval*100, objval, old_objval);
		 /*if (depth == 0 && 100*ABS(objval - old_objval) / objval > epsilon_LP) {
			 printf("imp rel: %.2f \n", ABS(objval - old_objval) / objval);
			 flag_solve_SP = 1;
			 old_objval = objval;
		 }
		else{*/
		if (depth > 0 && depth % 10 == 0 && count_same_node < 2) {
			flag_solve_SP = 1;
			old_objval = objval;
		}
		//}
	}
	count_added = 0;
	if (flag_solve_SP == 1) {
		status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols - 1);
		if (status) {
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
		for (i = 0; i < NN - 1; i++) {                                 //Solve (i,j) primal/dual subproblems
			for (j = i + 1; j < NN; j++) {
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
			if (fixed_zero[i] == 0) {
				for (k = 0; k < NN; k++) {
					if (fixed_zero[k] == 0) {
						coeff_z = 0;
						for (j = 0; j < NN; j++)
							coeff_z += (-alpha[i][j][k] - beta[j][i][k]);
						if (ABS(coeff_z) > 0.00001) {
							sepcut.cutind[sepcut.cutnz] = pos_z[i][k];
							sepcut.cutval[sepcut.cutnz++] = coeff_z;
							lhs += coeff_z * x[pos_z[i][k]];
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
							sepcut.cutind[sepcut.cutnz] = pos_z[i][k];
							sepcut.cutval[sepcut.cutnz++] = coeff_z;
							lhs += coeff_z * x[pos_z[i][k]];
						}
					}
				}
			}
		}
		if (lhs < tolerance_sep) {
			count_added++;
			status = CPXcutcallbackadd(env, cbdata, wherefrom, sepcut.cutnz, 0, 'G', sepcut.cutind, sepcut.cutval, CPX_USECUT_PURGE);
			if (status) {
				fprintf(stderr, "Failed to add cut.\n");
				goto TERMINATE;
			}
		}
	}
	if (count_added > 0) {
		stop_cutgen = 0;
		//printf("Added cuts: %d, percent decrease: %.3f \n", count_violcuts, (objval - last_objval)/objval*100);
		//old_lower_bound = new_lower_bound;
	}
	else
		stop_cutgen = 1;
	if (depth == 0 && count_added == 0)
		LP_lower_bound = cutinfo->nodeobjval;
	/* Tell CPLEX that cuts have been created */
	if (stop_cutgen == 0)
		*useraction_p = CPX_CALLBACK_SET;
TERMINATE:
	return (status);
} /* END mycutcallback */

void free_and_null(char** ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */

int NetFlow_TP(double* sol_z, int Origin_Node, int Destin_Node)
{
	/* Declare variables and arrays for retrieving problem data and
	   solution information later on. */
	int      narcs;
	int      nnodes;
	int      solstat;
	double   objval;
	double* x = NULL;//Arc Flows
	double* pi = NULL;//Shadow Prices for Flow Balance Constraints at nodes
	double* slack = NULL;//Slack in flow balance constraints
	double* dj = NULL;//Reduced costs for arc flows
	CPXENVptr env = NULL;
	CPXNETptr net = NULL;
	int       status;
	int       i, j, k, m;
	/* Initialize the CPLEX environment */
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[CPXMESSAGEBUFSIZE];
		fprintf(stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		goto TERMINATE;
	}
	/* Create the problem. */
	net = CPXNETcreateprob(env, &status, "netex1");
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
		beta[Origin_Node][Destin_Node][index_hub_dj[m - Breakpoint_O_D]] = (-pi[m]); //receiving
	}
TERMINATE:
	/* Free memory for solution data */
	free_and_null((char**)&pi);
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
	int* tail, * head;
	double* obj, * ub, * lb;
	double* supply;
	i_vector(&tail, NARCS, "open_cplex:4");
	i_vector(&head, NARCS, "open_cplex:6");
	d_vector(&obj, NARCS, "open_cplex:2");
	d_vector(&ub, NARCS, "open_cplex:7");
	d_vector(&lb, NARCS, "open_cplex:7");
	d_vector(&supply, NNODES, "open_cplex:7");
	if (vers == 3)sum_core = NN;
	sum_core = 1;
	////Using Maganti-Wong Appproach to obtain Pareto-Optimal cuts
	if (MG == 1) {
		sum_supply_i = 0;
		for (k = 0; k < NNODES_i; k++) {
			supply[k] = core[Origin_Node][index_hub_oi[k]] + sum_core * z_sol[pos_z[Origin_Node][index_hub_oi[k]]];//Supply at hub k = Zik
			sum_supply_i += supply[k];
		}
		sum_supply_j = 0;
		for (m = 0; m < NNODES_j; m++) {
			supply[NNODES_i + m] = -(core[Destin_Node][index_hub_dj[m]] + sum_core * z_sol[pos_z[Destin_Node][index_hub_dj[m]]]);//Supply at hub m = -Zjm
			sum_supply_j += supply[NNODES_i + m];
		}
	}
	else {
		////Using Papadakos Appproach to obtain Pareto-Optimal cuts
		//sum_supply_i = 0;
		for (k = 0; k < NNODES_i; k++) {
			supply[k] = z_sol[pos_z[Origin_Node][index_hub_oi[k]]];//Supply at hub k = Zik
			//sum_supply_i += supply[k];
		}
		//sum_supply_j = 0;
		for (m = 0; m < NNODES_j; m++) {
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
	for (k = 0; k < NNODES_i; k++)//From all first (k) hubs to all second (m) hubs
	{
		for (m = 0; m < NNODES_j; m++)
		{
			tail[count] = k;//Hub k node number starts from 0
			head[count] = NNODES_i + m;//
			obj[count] = (W[Origin_Node][Destin_Node] * c_t[index_hub_oi[k]][index_hub_dj[m]] + W[Destin_Node][Origin_Node] * c_t[index_hub_dj[m]][index_hub_oi[k]]);
			//obj[count]=(W[Origin_Node][Destin_Node])*(c_c[Origin_Node][k] + c_t[k][m] + c_d[m][Destin_Node]);
			//obj[count] = (W[Origin_Node][Destin_Node])*(c_t[cand_hubs[k]][cand_hubs[m]]);
			ub[count] = CPX_INFBOUND;
			lb[count] = 0;
			count++;
		}
	}
	//if ( CPXNETgetnumnodes (env, net) > 0 ) {
	//   status = CPXNETdelnodes (env, net, 0,
	//                            CPXNETgetnumnodes (env, net)-1);
	//   if ( status ) goto TERMINATE;
	//}
	/* Set optimization sense */
	status = CPXNETchgobjsen(env, net, CPX_MIN);
	if (status) goto TERMINATE;
	/* Add nodes to network along with their supply values,
	   but without any names. */
	status = CPXNETaddnodes(env, net, NNODES, supply, NULL);
	if (status) goto TERMINATE;
	/* Add arcs to network along with their objective values and
	   bounds, but without any names. */
	status = CPXNETaddarcs(env, net, NARCS, tail, head, lb, ub, obj, NULL);
	if (status) goto TERMINATE;
TERMINATE:
	free_and_null((char**)&tail);
	free_and_null((char**)&head);
	free_and_null((char**)&obj);
	free_and_null((char**)&ub);
	free_and_null((char**)&lb);
	free_and_null((char**)&supply);
	return (status);
}  /* END buildnetwork */

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

void Update_Core_Point(double* z_sol) {
	int i, k, j, m;
	double sum_i, sum_j;
	double fact = 0.5;
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
	for (i = 0; i < NN; i++) {
		for (k = 0; k < count_cand_hubs; k++) {
			//  if (z_sol[pos_z[i][cand_hubs[k]]] > 0.0001)
			//	  printf("z(%d %d): %.2f \n",i, cand_hubs[k], z_sol[pos_z[i][cand_hubs[k]]]);
			core[i][cand_hubs[k]] = fact * core[i][cand_hubs[k]] + (1 - fact) * z_sol[pos_z[i][cand_hubs[k]]];
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

void Update_CP_MW(double* z_sol, int i, int j) {
	int h, s;
	double acum_sup, acum_dem = 0, temp = 0, open = 0, checkcore = 0;
	sum_core = 0;
	for (h = 0; h < count_cand_hubs - 1; h++) {
		open += z_sol[pos_z[cand_hubs[h]][cand_hubs[h]]];
		//printf("Hub %d is open with %lf\n", cand_hubs[h],z_sol[pos_z[cand_hubs[h]][cand_hubs[h]]]);
	}
	for (h = 0; h < count_cand_hubs; h++) {
		temp = z_sol[pos_z[j][cand_hubs[h]]] - z_sol[pos_z[i][cand_hubs[h]]];
		if (temp > 0) {
			core[j][cand_hubs[h]] = (double)(0.75 / open);
			core[i][cand_hubs[h]] = core[j][cand_hubs[h]] + temp;
		}
		else {
			if (temp < 0) {
				core[i][cand_hubs[h]] = (double)(0.75 / open);
				/*printf("core val is %lf\n",core[i][cand_hubs[h]]);
				getchar();*/
				core[j][cand_hubs[h]] = core[i][cand_hubs[h]] - temp;
			}
			else {
				core[i][cand_hubs[h]] = (double)(0.1 / count_cand_hubs);
				core[j][cand_hubs[h]] = (double)(0.1 / count_cand_hubs);
			}
		}
		/****Updating the acumulations and sum_core****/
		sum_core += core[j][cand_hubs[h]] + core[i][cand_hubs[h]];
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

void Check_CP_MW(double* z_sol, int i, int j) {
	int h, s;
	double acum, temp = 0, open = 0, temp2, checkcore1 = 0, checkcore2;
	double* values_hub;
	values_hub = create_double_vector(count_cand_hubs);
	acum = 0;
	sum_core = 0;
	//printf("There are %d candidate hubs\n",count_cand_hubs);
	for (h = 0; h < count_cand_hubs; h++) {
		/*checkcore1+=core[i][cand_hubs[h]];
		checkcore2+=core[j][cand_hubs[h]];*/
		temp = (z_sol[pos_z[j][cand_hubs[h]]] - z_sol[pos_z[i][cand_hubs[h]]]);
		if (ABS(temp) > 0.00001) {
			values_hub[h] = (core[i][cand_hubs[h]] - core[j][cand_hubs[h]]) / temp;
			/*if(values_hub[h]<0){
				temp2=core[j][cand_hubs[h]];
				core[j][cand_hubs[h]]=core[i][cand_hubs[h]];
				core[i][cand_hubs[h]]=temp2;
				values_hub[h]=-values_hub[h];
			}*/
		}
		else values_hub[h] = 0;
		acum += values_hub[h];
		if (values_hub[h] > sum_core) sum_core = values_hub[h];
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
	if (acum < 0.001) sum_core = 0;
	else sum_core = sum_core;
	free(values_hub);
}

int CPXPUBLIC Heur(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p)
{
	int status = 0;
	int depth, SEQNUM;
	int       i, k, m, index, j, cols, addcuts = 0, opened;
	double    etaval, currincumbent;
	double* bestx;
	double* temp_x;
	int flag_newinc = 0;
	cols = count_cand_hubs * NN + 1; //candidate hubs by the nodes plus the artificial variable
	bestx = create_double_vector(cols);
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	*useraction_p = CPX_CALLBACK_DEFAULT;
	if (use_firstsolution == 1) {
		use_firstsolution = 0;
		etaval = 0;
		for (i = 0; i < cols; i++)x[i] = 0;
		for (i = 0; i < NN; i++) {
			x[pos_z[i][best_sol_assignments[i]]] = 1;
			//printf("Assigned %d to hub %d\n",i,best_assigmnent[i]);			
			for (m = 0; m < NN; m++)
				etaval += W[i][m] * (c_t[best_sol_assignments[i]][best_sol_assignments[m]]);
		}
		x[pos_eta] = etaval;
		printf("Value of etaval is %.2lf\n", etaval);
		*objval_p = UpperBound;   //Now telling it what the new objective value is
		Prev_incumbent = UpperBound;
		printf("Writing the best found solution with value %lf\n", UpperBound);
		*checkfeas_p = 1;
		/* Tell CPLEX that a solution is being returned */
		*useraction_p = CPX_CALLBACK_SET;
	}
TERMINATE:
	free(bestx);
	return (0);
} /* END Lagrange_Heur */