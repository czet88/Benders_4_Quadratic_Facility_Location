#include "def.h"

int add_variables(CPXENVptr env, CPXLPptr lp) {
	int i, index1, k, numcols, status;
	double* obj;    // objective function coefficients of the design variables
	double* lb;     //lower bound of the design variables
	double* ub;     //upper bound of the design variables
	char** varname;  //char array stating the name of the variable

	//Initializing the global problem variables
	/**********************************************/
	numcols = NN * NN + 1;
	pos_z = create_int_matrix(NN, NN);
	priority = create_int_vector(NN * NN);
	indices = create_int_vector(NN * NN);
	varname = create_stringarray(numcols, 20);
	
	globvarind = create_int_vector(numcols);  //Array containing all the Types for the complete model
	c_vector(&globvarctype, numcols, "open_cplex:01");

	//Define z_ik variables
	index1 = 0;  // index of columns
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
			lb[index1] = 0;
			ub[index1] = 1;
			indices[index1] = index1;
			globvarind[index1] = index1;
			globvarctype[index1] = 'B';
			sprintf(varname[index1], "z_(%d)(%d", i + 1, k +1);
			index1++;
		}
	}
	//The type of the continuous variable eta
	globvarind[index1] = index1;
	globvarctype[index1] = 'C';
	pos_eta =index1;
	obj[index1] = 1;
	lb[index1] = 0;
	ub[index1] = CPX_INFBOUND;
	index1++;
	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	free_stringarray(&varname, numcols);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);

	return numcols;
}

int add_assignment_constr(CPXENVptr env, CPXLPptr lp) {		//Add assignment constraints  \sum_{k \in N} z_ik = 1
	int		i, k;
	int		index, index1, acumstatus = 0;  // auxiliar indices to fill in the constraint matrix
	double* rhs;    // right hand side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zero element ...................
	double* matval; // coefficient values for the non-zero elements of constraints....
	int		numrows; // number of constraints.........................................
	int		numnz;   // number of non-zero elements in the matrix ....................

	numrows = NN + 1;
	numnz = NN * (numrows - 1) + 1;
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
	if (index >= numnz || index1 >= numrows) {
		printf("Exceded estimate\n");
		//system("pause");
	}
	acumstatus = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (acumstatus)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	return acumstatus;
}

int add_linking_constr(CPXENVptr env, CPXLPptr lp) {	//Add linking constraints  z_ik <= z_kk
	int		i, k;
	int		index, index1, acumstatus = 0;  // auxiliar indices to fill in the constraint matrix
	double* rhs;    // right hand side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zero element ...................
	double* matval; // coefficient values for the non-zero elements of constraints....
	int		numrows; // number of constraints.........................................
	int		numnz;   // number of non-zero elements in the matrix ....................

	numrows = NN * (NN - 1) + 1;
	numnz = 2 * (numrows - 1) + 1;
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
	if (index >= numnz || index1 >= numrows) {
		printf("Exceded estimate\n");
		getchar();
	}
	acumstatus = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (acumstatus)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	return acumstatus;
}

int add_capacity_constr(CPXENVptr env, CPXLPptr lp) {	//Add capacity constraints sum(i in NN) O_i z_ik <= b*z_kk
	int		i, k;
	int		index, index1, acumstatus = 0;  // auxiliar indices to fill in the constraint matrix
	double* rhs;    // right hand side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zero element ...................
	double* matval; // coefficient values for the non-zero elements of constraints....
	int		numrows; // number of constraints.........................................
	int		numnz;   // number of non-zero elements in the matrix ....................

	numrows = NN + 1;
	numnz = (numrows - 1) * (NN + 1) + 1;
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
	if (index >= numnz || index1 >= numrows) {
		printf("Exceded estimate\n");
		//system("pause");
	}
	acumstatus = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (acumstatus)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	return acumstatus;
}

int add_phubs_constr(CPXENVptr env, CPXLPptr lp) {	//Add exactly p_hubs
	int		i;
	int		index, index1, acumstatus = 0;  // auxiliar indices to fill in the constraint matrix
	double* rhs;    // right hand side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zero element ...................
	double* matval; // coefficient values for the non-zero elements of constraints....
	int		numrows; // number of constraints.........................................
	int		numnz;   // number of non-zero elements in the matrix ....................

	numrows = 1 + 1;
	numnz = NN + 1;
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
	if (index >= numnz || index1 >= numrows) {
		printf("Exceded estimate\n");
		//system("pause");
	}
	acumstatus = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (acumstatus)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	return acumstatus;
}