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