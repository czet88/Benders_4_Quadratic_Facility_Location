#include "def.h"

int Construct_Feasible_Solution(double *x, double *dj)
{
	int i, k, m, r, rr, count_c, status, num_modify, target_open_facilties, min_range;
	int l, p, iter, count_open, iter_max, done, target_modify, lower_limit, upper_limit;
	double rn, sum_cap, keep_perc_open, objvalue;
	ZVAL *z_cand;
	int *current_assigmnent;
	int *current_open_plants;
	double *current_capacity;
	int *initial_open_plants;
	int *which_open_plants;
	int *modified;
	current_assigmnent = create_int_vector(NN);
	current_open_plants = create_int_vector(NN);
	current_capacity = create_double_vector(NN);
	initial_open_plants = create_int_vector(NN);
	which_open_plants = create_int_vector(NN);
	modified = create_int_vector(NN);
	z_cand = (ZVAL *)calloc(NN, sizeof(ZVAL));
	iter = 0;
	iter_max = 10;
	keep_perc_open = 0.2;
	count_c = 0;
	printf("************************************\nStarting Heuristic Procedure\n");
	printf("candidate facilites in support: ");
	for (k = 0; k < NN; k++) {
		current_open_plants[k] = 0;
		if (fixed_zero[k] == 0 && x[pos_z[k][k]] > 0.001) {
			z_cand[count_c].k = k;
			z_cand[count_c++].value = x[pos_z[k][k]];
			printf("%d ", k + 1);
		}
	}
	printf("\n");
	//qsort((ZVAL *)z_cand, count_c, sizeof(z_cand[0]), Comparevalue_zc);
	// Obtain initial feasible solution by solving a MILP on the support of the current LP solution of the Benders reformulation
	status = CFLP_reduced_model(count_c, z_cand, current_assigmnent, current_open_plants);
	if (status != 0) {
		printf("something went wrong when solving reduced CFLP \n");
		return 0;
	}
	//Evaluate quadratic objective function
	objvalue = Evaluate_Quadratic_Objective(current_open_plants, current_assigmnent, current_capacity);
	if (objvalue < UpperBound) {
		UpperBound = objvalue;
		memcpy(best_sol_assignments, current_assigmnent, NN * sizeof(int));
		memcpy(best_sol_facilities, current_open_plants, NN * sizeof(int));
		printf("Initial Upperbound from Matheuristic: %.2f \n", objvalue);
	}
	if (FlagLocalSearch) {
		printf("...\nBeginning Local Search\n");
		//Disversification-Intensification strategies
		Facility_Change_Phase(current_assigmnent, current_open_plants, current_capacity, z_cand, count_c, &objvalue, dj);
		printf("UB after facility change: %.2f\n", objvalue);
		Assignment_Change_Phase(current_assigmnent, current_open_plants, current_capacity, z_cand, &objvalue);
		printf("UB after assignment change: %.2f\n", objvalue);
		printf("Open facilities: \n", objvalue);
		for (i = 0; i < NN; i++) {
			if (current_open_plants[i] == 1)
				printf("%d ", i + 1);
		}
		printf("\n");
		//Update incumbent solution
		if (objvalue < UpperBound) {
			UpperBound = objvalue;
			memcpy(best_sol_assignments, current_assigmnent, NN * sizeof(int));
			memcpy(best_sol_facilities, current_open_plants, NN * sizeof(int));
			printf("Improved Upperbound from Local Search: %.2f \n", objvalue);
		}
		memcpy(initial_open_plants, current_open_plants, NN * sizeof(int));
		printf("Finished Local Search\n...\n");
	}
	printf("Heuristic final UB: %.2f  open facilities: \n", objvalue);
	for (i = 0; i < NN; i++) {
		if (current_open_plants[i] == 1)
			printf("%d ", i + 1);
	}
	printf("\n");
	free(z_cand);
	free(current_assigmnent);
	free(current_open_plants);
	free(current_capacity);
	free(initial_open_plants);
	free(which_open_plants);
	free(modified);
	printf("Finished Heuristic Procedure\n************************************\n");
	return 1;
}

double Evaluate_Quadratic_Objective(int *current_open_plants, int *current_assigmnent, double *current_capacity)
{
	int i, k, m;
	double objvalue = 0;
	for (k = 0; k < NN; k++) {
		current_capacity[k] = b[k][0];
		if (current_open_plants[k] == 1) {
			objvalue += f[k][0];
		}
		for (m = 0; m < NN; m++)
			objvalue += W[k][m] * (c_c[k][current_assigmnent[k]] + c_t[current_assigmnent[k]][current_assigmnent[m]] + c_d[current_assigmnent[m]][m]);
	}
	for (i = 0; i < NN; i++)
		current_capacity[current_assigmnent[i]] -= O[i];
	return objvalue;
}

void Facility_Change_Phase(int *current_assigmnent, int *current_open_plants,  double *current_capacity, ZVAL *z_cand, int count_c, double *objvalue, double *dj){
	int k, flag;
	double value;
	value = *objvalue;
	for (k = 0; k < NN; k++)
		current_capacity[k] = b[k][0];
	for (k = 0; k < NN; k++)
		current_capacity[current_assigmnent[k]] -= O[k];
	flag = 1;
	while (flag) {
		flag = clients_swap_red(current_assigmnent, current_open_plants, current_capacity, &value);
		if (flag == 0) {
			flag = clients_shift_red(current_assigmnent, current_open_plants, current_capacity, &value);
			if (flag == 0 && w_p_median_constr == 0) {
				flag = open_hub_red(current_assigmnent, current_open_plants, current_capacity, &value);
				if (flag == 0 && w_p_median_constr == 0){
					flag = close_hub_red(current_assigmnent, current_open_plants, current_capacity, &value);
					if (flag == 0)
						flag = open_close_hub_red(current_assigmnent, current_open_plants, current_capacity, &value);
				}
			}
		}
	}
	*objvalue = value;
}

int Assignment_Change_Phase(int *current_assigmnent, int *current_open_plants, double *current_capacity, ZVAL *z_cand, double *objvalue){
	int i, k, m, status, flag;
	//double objvalue;
	double value;
	status = Reassign_nodes_red(current_assigmnent, current_open_plants);
	if (status == 1){
		value = 0;
		for (k = 0; k < NN; k++) {
			current_capacity[k] = b[k][0];
			if (current_open_plants[k] == 1){
				value += f[k][0];
			}
			for (m = 0; m < NN; m++)
				value += W[k][m] * (c_c[k][current_assigmnent[k]] + c_t[current_assigmnent[k]][current_assigmnent[m]] + c_d[current_assigmnent[m]][m]);
		}
		for (i = 0; i < NN; i++)
			current_capacity[current_assigmnent[i]] -= O[i];
		flag = 1;
		while (flag) {
			flag = clients_swap_red(current_assigmnent, current_open_plants, current_capacity, &value);
			if (flag == 0) {
				flag = clients_shift_red(current_assigmnent, current_open_plants, current_capacity, &value);
			}
		}
		*objvalue = value;
	}
	else
		return 0;
	return 1;
}

int clients_swap_red(int *assignment, int *current_open_plants, double *current_capacity, double *vs) {
	register i, j;
	double      MinCostChange = MAX_DOUBLE;
	double      CanMovCostChange;
	int      client1, client2;
	int      plant1, plant2;
	int      flag = 0;
	for (i = 0; i < NN; i++) {
		for (j = 0; j < NN; j++) {
			if (i != j && assignment[i] != assignment[j]) {
				if (current_capacity[assignment[i]] + O[i] - O[j] >= 0 && current_capacity[assignment[j]] + O[j] - O[i] >= 0) {
					//CanMovCostChange=suma2_F_ijkm(i, j, assignment[j], assignment[i], assignment);
					//CanMovCostChange-=suma2_F_ijkm(i, j, assignment[i], assignment[j], assignment);
					CanMovCostChange = sum2_F_ijkm(i, j, assignment[j], assignment[i], assignment);
					CanMovCostChange -= sum2_F_ijkm(i, j, assignment[i], assignment[j], assignment);
					if (CanMovCostChange < 0) {
						client1 = i;
						plant1 = assignment[i];
						client2 = j;
						plant2 = assignment[j];
						MinCostChange = CanMovCostChange;
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1)
			break;
	}
	if (MinCostChange < 0) {
		current_capacity[plant1] += O[client1] - O[client2];
		current_capacity[plant2] += O[client2] - O[client1];
		assignment[client1] = plant2;
		assignment[client2] = plant1;
		*vs += MinCostChange;
		return 1;
	}
	return 0;
}

int clients_shift_red(int *assignment, int *current_open_plants, double *current_capacity, double *vs) {
	register i, k;
	double      MinCostChange = MAX_DOUBLE;
	double      CanMovCostChange;
	int      client;
	int      plant1, plant2;
	//int      *open_plants;
	int      flag = 0;
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (current_open_plants[k] && k != assignment[i]) {
				if (current_capacity[k] - O[i] >= 0) {
					//CanMovCostChange=suma_F_ijkm(i,k, assignment)-suma_F_ijkm(i, assignment[i], assignment);
					CanMovCostChange = sum_F_ijkm(i, k, assignment) - sum_F_ijkm(i, assignment[i], assignment);
					if (b[assignment[i]][0] == current_capacity[assignment[i]] + O[i])
						CanMovCostChange -= f[assignment[i]][0];
					if (CanMovCostChange < 0) {
						client = i;
						plant1 = assignment[i];
						plant2 = k;
						MinCostChange = CanMovCostChange;
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1)
			break;
	}
	if (MinCostChange < 0) {
		current_capacity[plant1] += O[client];
		current_capacity[plant2] -= O[client];
		if (b[assignment[client]][0] == current_capacity[assignment[client]])
			current_open_plants[plant1] = 0;
		assignment[client] = plant2;
		*vs += MinCostChange;
		return 1;
	}
	return 0;
}

int CFLP_reduced_model(int count_c, ZVAL *z_cand, int *assigmnents, int *open_plants)
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
	//double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int *indices, *priority;
	int		nodecount;
	int cur_numcols;
	int    **pos_z_red;
	double *sol_x;
	int    *local_assign; //Contain the local optimal assignment
	pos_z_red = create_int_matrix(NN,count_c);
	//UpperBound = MAX_DOUBLE;
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
	numcols = count_c * NN;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (i = 0; i < NN; i++) {
		for (k = 0; k < count_c; k++) {
			pos_z_red[i][k] = index1;
			if (i == z_cand[k].k) {
				obj[index1] = f[i][0];
			}
			else {
				obj[index1] = (O[i] * c_c[i][z_cand[k].k] + D[i] * c_d[i][z_cand[k].k]);
				//	 obj[index1] = 0;
			}
			ctype[index1] = 'B';
			lb[index1] = 0;
			ub[index1] = 1;
			index1++;
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
	numnz = count_c * NN;
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
		for (k = 0; k < count_c; k++) {
			matind[index] = pos_z_red[i][k];
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
	numrows = count_c * (NN - 1);
	numnz = 2 * count_c*(NN - 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++) {
		for (k = 0; k < count_c; k++) {
			if (i != z_cand[k].k) {
				sense[index1] = 'L';
				rhs[index1] = 0;
				matbeg[index1++] = index;
				matind[index] = pos_z_red[i][k];
				matval[index++] = 1;
				matind[index] = pos_z_red[z_cand[k].k][k];
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
	if (w_p_median_constr==1) {
		numrows = 1;
		numnz = count_c;
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
		for (i = 0; i < count_c; i++) {
			matind[index] = pos_z_red[z_cand[i].k][i];
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
		numrows = count_c;
		numnz = count_c * (NN + 1);
		d_vector(&rhs, numrows, "open_cplex:2");
		c_vector(&sense, numrows, "open_cplex:3");
		i_vector(&matbeg, numrows, "open_cplex:4");
		i_vector(&matind, numnz, "open_cplex:6");
		d_vector(&matval, numnz, "open_cplex:7");
		index = 0;
		index1 = 0;
		for (k = 0; k < count_c; k++) {
			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_z_red[z_cand[k].k][k];
			matval[index++] = -(b[z_cand[k].k][0] - O[z_cand[k].k]);
			for (i = 0; i < NN; i++) {
				if (i != z_cand[k].k) {
					matind[index] = pos_z_red[i][k];
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
	//CPXwriteprob(env, lp, "reduced_model.lp", NULL);
	//branch and bound parameters
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF); //output display
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4); //different levels of output display
	//CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1);
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env, CPX_PARAM_TILIM, 1000); // time limit
	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.00001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
//	CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, 0); //heuristic frequency and intensisty 
	//CPXoptimize(env,lp);   solve a linear program
	//CPXgetbestobjval(env,lp,&value);   \\obtain the LP relaxation bound
	//printf("LP lower bound: %f   ",value);
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,238016.277+0.001); // provide an initial upper bound
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
   /* Assure linear mappings between the presolved and original
	  models */
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if (status)  goto TERMINATE;
	/* Turn on traditional search for use with control callbacks */
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if (status)  goto TERMINATE;
	/* Let MIP callbacks work on the original model */
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if (status)  goto TERMINATE;
	/* Set up to use MIP callback */
	//status = CPXsetlazyconstraintcallbackfunc(env, mycheckcallback, NULL);
	//if (status)  goto TERMINATE;
	i = CPXmipopt(env, lp);  //solve the integer program
	i = CPXgetstat(env, lp);
	if (i != 101 && i != 102 && i != 107)
		status = 0;
	else
		status = 1;
	/*if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);*/
	// out = open_file("result_CHLP_3index.txt","a+");
	// retrive solution values
	CPXgetmipobjval(env, lp, &UFLPrelval);
	//printf("MIP value: %.2f \n", UFLPrelval);
	//CPXgetbestobjval(env,lp,&value);  //best lower bound in case thew problem was not solve to optimality
	end = clock();
	cpuFacLocIni += (double)(end - start) / CLOCKS_PER_SEC;
	//printf("Time to solve SSCFLP: %.2f \n", cputime);
	/*Retrieve the solutiona and calculate corresponding quadratic costs*/
	numcols = CPXgetnumcols(env, lp);
	sol_x = create_double_vector(numcols);
	//local_assign = create_int_vector(NN);
	CPXgetmipx(env, lp, sol_x, 0, numcols - 1);  // obtain the values of the decision variables
	for (i = 0; i < NN; i++) {
		for (k = 0; k < count_c; k++) {
			if (sol_x[pos_z_red[i][k]] > 0.5) {
				assigmnents[i] = z_cand[k].k;
				//printf("%d(%d) ", i + 1, assigmnents[i] + 1);
			}
		}
	}
	UBOptUFLP = UFLPrelval;
	for (i = 0; i < NN; i++) {
		for (j = 0; j < NN; j++) {
			UBOptUFLP += W[i][j] * (c_t[assigmnents[i]][assigmnents[j]]);
		}
	}
	//printf("LB: %.2f UB: %.2f initial facilities: ", UFLPrelval, UBOptUFLP);
	for (i = 0; i < count_c; i++) {
		if (sol_x[pos_z_red[z_cand[i].k][i]] > 0.5) {
			open_plants[z_cand[i].k] = 1;
			//printf("%d (%.2f) ", z_cand[i].k + 1, f[z_cand[i].k][0]);
		}
	}
	//printf("\n");
	//printf("Initial Lowerbound: %lf \t Upperbound from SSCFLP: %lf\t Upperbound from optimal SSCFLP: %lf \n", UFLPrelval, UFLPrelval, UBOptUFLP);
	//UBsearchUFLP = UpperBound;
	free(sol_x);
	free(pos_z_red);
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
	return status;
}

int open_hub_red(int *current_assigmnent, int *current_open_plants, double *current_capacity, double *vs)
{
	int i, j, k, m;
	double   MinCostChange;
	double   valor;
	int      plant1, node;
	int      flag = 0;
	MinCostChange = *vs;
	for (i = 0; i < NN; i++) {
		if (fixed_zero[i] == 0 && current_open_plants[i] == 0 && b[i][0] - O[i] >= 0) {
			//current_open_plants[i] = 1;
			memcpy(allocation, current_assigmnent, NN * sizeof(int));
			memcpy(capacity, current_capacity, NN * sizeof(double));
			capacity[allocation[i]] += O[i];
			allocation[i] = i;
			capacity[i] -= O[i];
			for (j = 1; j < NN; j++) {
				node = costoso[i].A[j].j;
				if (capacity[i] - O[node] >= 0 && c[node][allocation[node]] > c[node][i]) {
					capacity[i] -= O[node];
					capacity[allocation[node]] += O[node];
					allocation[node] = i;
				}
				if (capacity[i] < menor_O)
					break;
			}
			valor = 0;
			for (k = 0; k < NN; k++) {
				if (capacity[k] < b[k][0])
					valor += f[k][0];
				for (m = 0; m < NN; m++) {
					//valor+=F_ijkm[k][m][allocation[k]][allocation[m]];
					valor += W[k][m] * (c_c[k][allocation[k]] + c_t[allocation[k]][allocation[m]] + c_d[allocation[m]][m]);
				}
			}
			if (valor < MinCostChange) {
				memcpy(best_allocation, allocation, NN * sizeof(int));
				memcpy(best_capacity, capacity, NN * sizeof(double));
				plant1 = i;
				MinCostChange = valor;
			}
			//current_open_plants[i] = 0;
		}
	}
	if (MinCostChange < *vs) {
		memcpy(current_assigmnent, best_allocation, NN * sizeof(int));
		memcpy(current_capacity, best_capacity, NN * sizeof(double));
		current_open_plants[plant1] = 1;
		*vs = MinCostChange;
		return 1;
	}
	return 0;
}

int close_hub_red(int *current_assigmnent, int *current_open_plants, double *current_capacity, double *vs)
{
	int i, j, k, m;
	double   MinCostChange;
	double   valor, min_val;
	int      plant1;
	int      flag;
	int      node, node2;
	MinCostChange = *vs;
	for (i = 0; i < NN; i++) {
		//if (avail_capacity[i] < b[i]) {
		if (current_open_plants[i] == 1) {
			flag = 0;
			//open_plants[i] = 0;
			//current_open_plants[i] = 0;
			memcpy(allocation, current_assigmnent, NN * sizeof(int));
			memcpy(capacity, current_capacity, NN * sizeof(double));
			for (j = 1; j < NN; j++) {
				node = costoso[i].A[j].j;
				if (capacity[node] < b[node][0] && capacity[node] - O[i] >= 0) {
					capacity[node] -= O[i];
					allocation[i] = node;
					capacity[i] = b[i][0];
					flag = 1;
					break;
				}
			}
			if (flag == 1) {
				for (j = 0; j < NN; j++) {
					node = orden_O[NN - j - 1].hub;
					if (allocation[node] == i) {
						for (k = 1; k < NN; k++) {
							node2 = costoso[node].A[k].j;
							if (capacity[node2] < b[node2][0] && capacity[node2] - O[node] >= 0 && node2 != i) {
								capacity[node2] -= O[node];
								allocation[node] = node2;
								break;
							}
						}
						if (allocation[node] == i) {
							flag = 0;
							break;
						}
					}
				}
			}
			if (flag == 1) {
				valor = 0;
				for (k = 0; k < NN; k++) {
					if (capacity[k] < b[k][0])
						valor += f[k][0];
					for (m = 0; m < NN; m++) {
						//valor+=F_ijkm[k][m][allocation[k]][allocation[m]];
						valor += W[k][m] * (c_c[k][allocation[k]] + c_t[allocation[k]][allocation[m]] + c_d[allocation[m]][m]);
					}
				}
				if (valor < MinCostChange) {
					memcpy(best_allocation, allocation, NN * sizeof(int));
					memcpy(best_capacity, capacity, NN * sizeof(double));
					plant1 = i;
					MinCostChange = valor;
				}
			}
			//current_open_plants[i] = 1;
		}
	}
	if (MinCostChange < *vs) {
		memcpy(current_assigmnent, best_allocation, NN * sizeof(int));
		memcpy(current_capacity, best_capacity, NN * sizeof(double));
		current_open_plants[plant1] = 0;
		*vs = MinCostChange;
		return 1;
	}
	return 0;
}

int open_close_hub_red(int *current_assigmnent, int *current_open_plants, double *current_capacity, double *vs)
{
	int i, j, k, m, kk;
	double   MinCostChange;
	double   valor;
	int      plant1;
	int      plant2, node1, node2, node;
	double   min_val;
	int      flag = 0;
	MinCostChange = *vs;
	for (k = 0; k < NN; k++) {
		for (m = 0; m < NN; m++) {
			if (fixed_zero[k] == 0 && current_open_plants[m] == 1 && current_open_plants[k] == 0 && b[k][0] - O[k] > 0 ) {
				flag = 0;
				current_open_plants[k] = 1;
				current_open_plants[m] = 0;
				memcpy(allocation, current_assigmnent, NN * sizeof(int));
				memcpy(capacity, current_capacity, NN * sizeof(double));
				node1 = allocation[k];
				capacity[allocation[k]] += O[k];
				allocation[k] = k;
				capacity[k] -= O[k];
				for (j = 1; j < NN; j++) {
					node = costoso[m].A[j].j;
					if (capacity[node] < b[node][0] && capacity[node] - O[m] >= 0) {
						capacity[node] -= O[m];
						allocation[m] = node;
						capacity[m] = b[m][0];
						flag = 1;
						break;
					}
				}
				if (flag == 1) {
					for (j = 0; j < NN; j++) {
						node = orden_O[NN - j - 1].hub;
						if (allocation[node] == m) {
							for (kk = 1; kk < NN; kk++) {
								node2 = costoso[node].A[kk].j;
								if (capacity[node2] < b[node2][0] && capacity[node2] - O[node] >= 0 && node2 != m) {
									capacity[node2] -= O[node];
									allocation[node] = node2;
									break;
								}
							}
							if (allocation[node] == m) {
								flag = 0;
								break;
							}
						}
					}
				}
				if (flag == 1) {
					for (j = 1; j < NN; j++) {
						node = costoso[k].A[j].j;
						if (capacity[k] - O[node] >= 0 && c[node][allocation[node]] > c[node][k]) {
							capacity[k] -= O[node];
							capacity[allocation[node]] += O[node];
							allocation[node] = k;
						}
						if (capacity[k] < menor_O)
							break;
					}
					valor = 0;
					for (i = 0; i < NN; i++) {
						if (current_open_plants[i] == 1)
							valor += f[i][0];
						for (j = 0; j < NN; j++) {
							valor += W[i][j] * (c_c[i][allocation[i]] + c_t[allocation[i]][allocation[j]] + c_d[allocation[j]][j]);
						}
					}
					if (valor < MinCostChange) {
						//memcpy(best_allocation, allocation, NN * sizeof(int));
						memcpy(current_assigmnent, allocation, NN * sizeof(int));
						memcpy(current_capacity, capacity, NN * sizeof(double));
						MinCostChange = valor;
						*vs = MinCostChange;
						return 1;
					}
				}
				current_open_plants[k] = 0;
				current_open_plants[m] = 1;
				allocation[k] = node1;
			}
		}
	}
	return 0;
}

int Reassign_nodes_red(int *best_assigmnent1, int *current_open_plants)
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
	int **pos_z_red;
	int		nodecount, statusP;
	int cur_numcols;
	pos_z_red = create_int_matrix(NN, NN);
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
	numcols = NN * NN;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (i != k && current_open_plants[i] == 0 && current_open_plants[k] == 1) {
				pos_z_red[i][k] = index1;
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
	numnz = NN * NN;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++) {
		if (current_open_plants[i] == 0) {
			sense[index1] = 'E';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			for (k = 0; k < NN; k++) {
				if (i != k && current_open_plants[k] == 1) {
					matind[index] = pos_z_red[i][k];
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
	numnz = NN * (NN + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	index = 0;
	index1 = 0;
	for (k = 0; k < NN; k++) {
		if (current_open_plants[k] == 1) {
			sense[index1] = 'L';
			rhs[index1] = b[k][0] - O[k];
			matbeg[index1++] = index;
			for (i = 0; i < NN; i++) {
				if (i != k && current_open_plants[i] == 0) {
					matind[index] = pos_z_red[i][k];
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
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); // Number of threads to use
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
	if (i == 101 || i == 102|| i==107) {
		// retrive solution values
		statusP = 1;
		numcols = CPXgetnumcols(env, lp);
		d_vector(&x, numcols, "open_cplex:9");
		CPXgetmipobjval(env, lp, &value);
		CPXgetmipx(env, lp, x, 0, numcols - 1);  // obtain the values of the decision variables
		for (i = 0; i < NN; i++) {
			if (current_open_plants[i] == 0) {
				for (k = 0; k < NN; k++) {
					if (i != k && current_open_plants[k] == 1 && x[pos_z_red[i][k]] > 0.5) {
						best_assigmnent1[i] = k;
						break;
					}
				}
			}
			else {
				best_assigmnent1[i] = i;
			}
		}
		free(x);
	}
	else 
		statusP = 0;
	end = clock();
	cpuGenAss+=(double)(end - start) / CLOCKS_PER_SEC;
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
	free(pos_z_red);
	//free(x);
	return statusP;
}

int Reassign_nodes(int* best_assigmnent1)
{
	int i, j, k, l, m, count;
	FILE* out;
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
	char probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int* indices, * priority;
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
	numcols = NN * NN;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (i != k && open_plants[i] == 0 && open_plants[k] == 1) {
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
	numnz = NN * NN;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	index = 0;
	index1 = 0;
	for (i = 0; i < NN; i++) {
		if (open_plants[i] == 0) {
			sense[index1] = 'E';
			rhs[index1] = 1;
			matbeg[index1++] = index;
			for (k = 0; k < NN; k++) {
				if (i != k && open_plants[k] == 1) {
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
	numnz = NN * (NN + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");
	index = 0;
	index1 = 0;
	for (k = 0; k < NN; k++) {
		if (open_plants[k] == 1) {
			sense[index1] = 'L';
			rhs[index1] = b[k][0] - O[k];
			matbeg[index1++] = index;
			for (i = 0; i < NN; i++) {
				if (i != k && open_plants[i] == 0) {
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
	CPXsetintparam(env, CPX_PARAM_THREADS, 1); // Number of threads to use
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
	if (i == 101 || i == 102) {
		// retrive solution values
		statusP = 1;
		numcols = CPXgetnumcols(env, lp);
		d_vector(&x, numcols, "open_cplex:9");
		CPXgetmipobjval(env, lp, &value);
		CPXgetmipx(env, lp, x, 0, numcols - 1);  // obtain the values of the decision variables
		for (i = 0; i < NN; i++) {
			if (open_plants[i] == 0) {
				for (k = 0; k < NN; k++) {
					if (i != k && open_plants[k] == 1 && x[pos_z[i][k]] > 0.5) {
						best_assigmnent1[i] = k;
						break;
					}
				}
			}
			else {
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

double sum2_F_ijkm(int i, int j, int k, int m, int* assignment) {
	int t;
	double suma = 0;
	for (t = 0; t < NN; t++) {
		if (i != t && j != t)
			suma += W[i][t] * (c_c[i][k] + c_t[k][assignment[t]] + c_d[assignment[t]][t]) + W[t][i] * (c_c[t][assignment[t]] + c_t[assignment[t]][k] + c_d[k][i]);
		else {
			if (i == t)
				suma += W[i][t] * (c_c[i][k] + c_t[k][k] + c_d[k][t]);
		}
		if (j != t && i != t)
			suma += W[j][t] * (c_c[j][m] + c_t[m][assignment[t]] + c_d[assignment[t]][t]) + W[t][j] * (c_c[t][assignment[t]] + c_t[assignment[t]][m] + c_d[m][j]);
		else {
			if (j == t)
				suma += W[j][t] * (c_c[j][m] + c_t[m][m] + c_d[m][t]);
		}
	}
	suma += W[i][j] * (c_c[i][k] + c_t[k][m] + c_d[m][j]) + W[j][i] * (c_c[j][m] + c_t[m][k] + c_d[k][i]);
	return suma;
}

double sum_F_ijkm(int i, int k, int* assignment) {
	int j;
	double suma = 0;
	for (j = 0; j < NN; j++) {
		if (i != j)
			suma += W[i][j] * (c_c[i][k] + c_t[k][assignment[j]] + c_d[assignment[j]][j]) + W[j][i] * (c_c[j][assignment[j]] + c_t[assignment[j]][k] + c_d[k][i]);
		else
			suma += W[i][j] * (c_c[i][k] + c_t[k][k] + c_d[k][j]);
	}
	return suma;
}

void PopulateLPSupport(double* x) {
	int k;
	double val;
	for (k = 0; k < NN; k++) {
		if (x == NULL) val = best_sol_facilities[k];
		else val = x[pos_z[k][k]];
		if (val > 0.001) {
			FlagHubLPsupport[k] = 1;
			//printf("Found one\n");
		}
		else {
			FlagHubLPsupport[k] = 0;
		}
	}
}

int CompareLPSupport(double* x) {
	int k, retval = 0;
	for (k = 0; k < NN; k++) {
		if ((x[pos_z[k][k]] < 0.001 && FlagHubLPsupport[k] == 1)) {
			retval++;
		}
	}
	return retval;
}
