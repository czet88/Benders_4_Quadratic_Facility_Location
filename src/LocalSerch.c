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
//extern int        *best_assigmnent;
extern int        *open_plants;
extern int        *allocation;
extern int        *best_allocation;
extern double     *capacity;
extern double     *avail_capacity;
extern double     *best_capacity;
extern double     menor_O;
extern PORD       *costoso;
extern SELEC      *orden_O;
extern double     AggregatedDemand;
extern int Capacitated_instances;


double Det_Iterated_Local_Search(void) {

	int i, k, m, r, rr, flag, index, num_modify;
	FILE     *out;
	clock_t  start, end;
	double objvalue, objvalue1;
	double   cputime;
	int *best_assigmnent1;
	int l, p, iter, status, count_open, iter_max, done, target_modify, lower_limit, upper_limit;
	int co, onp; // co : count open plants, onp: open new plant randomly
	double rn, sum_cap, keep_perc_open;
	int temp[10];
	int *best_open_plants;
	int *which_open_plants;
	int *modified;

	objvalue1 = MAX_DOUBLE;
	best_assigmnent1 = create_int_vector(NN);
	best_open_plants = create_int_vector(NN);
	which_open_plants = create_int_vector(NN);
	modified = create_int_vector(NN);
	memcpy(best_assigmnent1, best_assigmnent, NN* sizeof(int));
	memcpy(best_open_plants, open_plants, NN* sizeof(int));

	start = clock();
	srand(123456789);
	iter = 0;
	iter_max = 20;
	keep_perc_open = 0.2;

																//Phase I: Apply Pure local search to SSCFLP solution
		//Evaluate objective function
		objvalue = 0;
		for (k = 0; k < NN; k++) {
			avail_capacity[k] = b[k][0];
			if (open_plants[k] == 1){
				objvalue += f[k][0];
			}
			for (m = 0; m < NN; m++)
				objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
		}
		for (i = 0; i < NN; i++)
			avail_capacity[best_assigmnent1[i]] -= O[i];

		flag = 1;
		while (flag) {
			flag = clients_swap1_f(best_assigmnent1, &objvalue);
			if (flag == 0) {
				flag = clients_shift1_f(best_assigmnent1, &objvalue);
				if (flag == 0 && hybrid==1) {
					flag = open_hub(best_assigmnent1, &objvalue);
					if (flag == 0 && hybrid==1){
						flag = close_hub(best_assigmnent1, &objvalue);
						if (flag == 0)
							flag = open_close_hub(best_assigmnent1, &objvalue);
					}
				}
			}
		}
		status = Reassign_nodes(best_assigmnent1);
		objvalue = 0;
		for (k = 0; k < NN; k++) {
			avail_capacity[k] = b[k][0];
			if (open_plants[k] == 1){
				objvalue += f[k][0];
				//	printf("%d ", k + 1);
			}
			for (m = 0; m < NN; m++)
				objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
		}
		for (i = 0; i < NN; i++)
			avail_capacity[best_assigmnent1[i]] -= O[i];

		flag = 1;
		while (flag) {
			flag = clients_swap1_f(best_assigmnent1, &objvalue);
			if (flag == 0) {
				flag = clients_shift1_f(best_assigmnent1, &objvalue);
				///*if (flag == 0) {
				//	flag = open_hub(best_assigmnent1, &objvalue);
				//	if (flag == 0){
				//		flag = close_hub(best_assigmnent1, &objvalue);*/
				//		/*if (flag == 0)
				//			flag = open_close_hub(best_assigmnent1, &objvalue);*/
				//	}
				//}
			}
		}
		//Update incumbent solution
		if (objvalue < objvalue1){
			objvalue1 = objvalue;
			memcpy(best_assigmnent, best_assigmnent1, NN* sizeof(int));
			memcpy(best_open_plants, open_plants, NN* sizeof(int));
			printf("Improved Upperbound from PLS: %.2f \n", objvalue);
			count_open = 0;
			/**********************Carlos's modifications*************************/
			prevsols[countsols].preeta=0;     //Carlos's modification
			prevsols[countsols].num_comb=0;    //Carlos's modification
			for (k = 0; k < NN; k++) {
				if (best_open_plants[k] >0.5){
					count_open++;
				}
				prevsols[countsols].indhub[k]=best_assigmnent[k];   //Carlos's modification
				prevsols[countsols].num_comb++;
				for (m = 0; m < NN; m++){
					prevsols[countsols].preeta+=W[k][m]*c_t[best_assigmnent[k]][best_assigmnent[m]]; //calculating the preeta  //Carlos's modification
				}
			}
			prevsols[countsols].preeta-=combtol;
			countsols++;         //Carlos's modification
			prevsols[countsols].indhub=create_int_vector(NN);   //Carlos's modification
		}


		//Perturb current set of hubs (Variant II of the heuristic)
		// This works for problems with fixed cost for opening facilities
		if(hybrid==1){
			for (l = 0; l < NN; l++)
				open_plants[l] = best_open_plants[l];
			for (l = 0; l < NN; l++) {
				if (best_open_plants[l] == 1){
					open_plants[l] = 0;
					do{
						status = 0;
						sum_cap = 0;
						for (r = 0; r < NN; r++) {
							if (best_open_plants[r] == 0){
								rn = rand_double();
								if (rn < (double)2 / (NN - count_open) && b[r][0] - O[r] >= 0){
									open_plants[r] = 1;
									sum_cap += b[r][0];
								}
								else
									open_plants[r] = 0;
							}
							else{
								if (l != r)
									open_plants[r] = 1;
							}
						}
						if (AggregatedDemand <= sum_cap)
							status = Reassign_nodes(best_assigmnent1);
					} while (status != 1);

					//Evaluate objective function
					objvalue = 0;
					for (k = 0; k < NN; k++) {
						avail_capacity[k] = b[k][0];
						if (open_plants[k] == 1){
							objvalue += f[k][0];
						}
						for (m = 0; m < NN; m++)
							objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
					}
					for (i = 0; i < NN; i++)
						avail_capacity[best_assigmnent1[i]] -= O[i];

					flag = 1;
					while (flag) {
						flag = clients_swap1_f(best_assigmnent1, &objvalue);
						if (flag == 0) {
							flag = clients_shift1_f(best_assigmnent1, &objvalue);
							if (flag == 0 && hybrid==1) {
								flag = open_hub(best_assigmnent1, &objvalue);
								if (flag == 0 && hybrid==1){
									flag = close_hub(best_assigmnent1, &objvalue);
									if (flag == 0)
										flag = open_close_hub(best_assigmnent1, &objvalue);
								}
							}
						}
					}
					status = Reassign_nodes(best_assigmnent1);
					objvalue = 0;
					for (k = 0; k < NN; k++) {
						avail_capacity[k] = b[k][0];
						if (open_plants[k] == 1){
							objvalue += f[k][0];
							//	printf("%d ", k + 1);
						}
						for (m = 0; m < NN; m++)
							objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
					}
					for (i = 0; i < NN; i++)
						avail_capacity[best_assigmnent1[i]] -= O[i];

					flag = 1;
					while (flag) {
						flag = clients_swap1_f(best_assigmnent1, &objvalue);
						if (flag == 0) {
							flag = clients_shift1_f(best_assigmnent1, &objvalue);
							//if (flag == 0) {
							//	flag = open_hub(best_assigmnent1, &objvalue);
							//	if (flag == 0){
							//		flag = close_hub(best_assigmnent1, &objvalue);
							//		/*if (flag == 0)
							//			flag = open_close_hub(best_assigmnent1, &objvalue);*/
							//	}
							//}
						}
					}
					//Update incumbent solution
					if (objvalue < objvalue1){
						objvalue1 = objvalue;
						memcpy(best_assigmnent, best_assigmnent1, NN* sizeof(int));
						memcpy(best_open_plants, open_plants, NN* sizeof(int));
						printf("Improved Upperbound from ILS: %.2f \n", objvalue);
						count_open = 0;
						/**********************Carlos's modifications*************************/
						prevsols[countsols].preeta=0;     //Carlos's modification
						prevsols[countsols].num_comb=0;    //Carlos's modification
						for (k = 0; k < NN; k++) {
							if (best_open_plants[k] >0.5){
								count_open++;
							}
							prevsols[countsols].indhub[k]=best_assigmnent[k];   //Carlos's modification
							prevsols[countsols].num_comb++;
							for (m = 0; m < NN; m++){
								prevsols[countsols].preeta+=W[k][m]*c_t[best_assigmnent[k]][best_assigmnent[m]]; //calculating the preeta  //Carlos's modification
							}
						}
						prevsols[countsols].preeta-=combtol;
						countsols++;         //Carlos's modification
						prevsols[countsols].indhub=create_int_vector(NN);   //Carlos's modification
					}

				}
			}
		}
		else {
		
			for (iter = 0; iter < iter_max; iter++) {
				for (l = 0; l < NN; l++)
					open_plants[l] = best_open_plants[l];
				// This works for p-median variants
				count_open = 0;
				sum_cap = 0;
				for (l = 0; l < NN; l++) {
					if (open_plants[l] == 1) {
						which_open_plants[count_open] = l;
						modified[count_open++] = 0;
						sum_cap += b[l][0];
					}
				}
				target_modify = (int)ceil(keep_perc_open*count_open);
				if (target_modify <= 1){
					lower_limit = 1;
					upper_limit = 2;
				}
				else{
					lower_limit = target_modify - 1;
					upper_limit = target_modify + 1;
				}
				num_modify = getrandom(lower_limit, upper_limit);
				//num_modify = (int)ceil(keep_perc_open*count_open);
				printf("open: %d p: %d nummodify: %d \n",count_open,p_hubs, num_modify);
				for (l = 0; l < num_modify; l++) {
					done = 0;
					do {
						r = getrandom(0, count_open-1);
						if (modified[r] == 0) {
							done = 1;
							modified[r] = 1;
							sum_cap -= b[which_open_plants[r]][0];
							status = 0;
							open_plants[which_open_plants[r]] = 0;
							//printf("r:%d which: %d \n", r, which_open_plants[r]);
							do {
								rr = getrandom(0, NN-1);
								if (open_plants[rr] == 0 && which_open_plants[r] !=rr) {
									open_plants[rr] = 1;
									sum_cap += b[rr][0];
									status = 1;
								//	printf("r:%d which: %d rr:%d \n", r, which_open_plants[r], rr);
								}
							} while (status != 1);
						}
					} while (done != 1);
				}
				//printf("flag1\n");
				//feasible = 0;
				if (AggregatedDemand <= sum_cap)
					status = Reassign_nodes(best_assigmnent1);
				if (status == 1) {
					//Evaluate objective function
					objvalue = 0;
					for (k = 0; k < NN; k++) {
						avail_capacity[k] = b[k][0];
						if (open_plants[k] == 1) {
							objvalue += f[k][0];
						}
						for (m = 0; m < NN; m++)
							objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
					}
					for (i = 0; i < NN; i++)
						avail_capacity[best_assigmnent1[i]] -= O[i];

					flag = 1;
					while (flag) {
						flag = clients_swap1_f(best_assigmnent1, &objvalue);
						if (flag == 0) {
							flag = clients_shift1_f(best_assigmnent1, &objvalue);
							//if (flag == 0 && hybrid == 1) {
						//		flag = open_hub(best_assigmnent1, &objvalue);
							//	if (flag == 0 && hybrid == 1) {
						//			flag = close_hub(best_assigmnent1, &objvalue);
									if (flag == 0)
										flag = open_close_hub(best_assigmnent1, &objvalue);
								}
						//	}
						//}
					}
				//	printf("flag2\n");
					status = Reassign_nodes(best_assigmnent1);
					objvalue = 0;
					for (k = 0; k < NN; k++) {
						avail_capacity[k] = b[k][0];
						if (open_plants[k] == 1) {
							objvalue += f[k][0];
							//	printf("%d ", k + 1);
						}
						for (m = 0; m < NN; m++)
							objvalue += W[k][m] * (c_c[k][best_assigmnent1[k]] + c_t[best_assigmnent1[k]][best_assigmnent1[m]] + c_d[best_assigmnent1[m]][m]);
					}
					for (i = 0; i < NN; i++)
						avail_capacity[best_assigmnent1[i]] -= O[i];
					//printf("flag3\n");
					flag = 1;
					while (flag) {
						flag = clients_swap1_f(best_assigmnent1, &objvalue);
						if (flag == 0) {
							flag = clients_shift1_f(best_assigmnent1, &objvalue);
							//if (flag == 0) {
							//	flag = open_hub(best_assigmnent1, &objvalue);
							//	if (flag == 0){
							//		flag = close_hub(best_assigmnent1, &objvalue);
							//		/*if (flag == 0)
							//			flag = open_close_hub(best_assigmnent1, &objvalue);*/
							//	}
							//}
						}
					}
				//	printf("flag4\n");
					//Update incumbent solution
					if (objvalue < objvalue1) {
						objvalue1 = objvalue;
						memcpy(best_assigmnent, best_assigmnent1, NN * sizeof(int));
						memcpy(best_open_plants, open_plants, NN * sizeof(int));
						printf("Improved Upperbound from ILS: %.2f \n", objvalue);
						count_open = 0;
						/**********************Carlos's modifications*************************/
						prevsols[countsols].preeta = 0;     //Carlos's modification
						prevsols[countsols].num_comb = 0;    //Carlos's modification
						for (k = 0; k < NN; k++) {
							if (best_open_plants[k] > 0.5) {
								count_open++;
							}
							prevsols[countsols].indhub[k] = best_assigmnent[k];   //Carlos's modification
							prevsols[countsols].num_comb++;
							for (m = 0; m < NN; m++) {
								prevsols[countsols].preeta += W[k][m] * c_t[best_assigmnent[k]][best_assigmnent[m]]; //calculating the preeta  //Carlos's modification
							}
						}
						prevsols[countsols].preeta -= combtol;
						countsols++;         //Carlos's modification
						prevsols[countsols].indhub = create_int_vector(NN);   //Carlos's modification
						
					}

				}
			}
		}


	UpperBound = objvalue1;
	index = 0;
	for (i = 0; i < NN; i++) {
		for (k = 0; k < NN; k++) {
			if (best_assigmnent[i] == k)
				initial_x[index++] = 1;
			else
				initial_x[index++] = 0;
		}
	}




	//printf("Upperbound from Local Search: %.12f \n", objvalue);
	//
	//printf("Improved solution: ");
	//objvalue = 0;

	//for (k = 0; k < NN; k++) {
	//	if (open_plants[k] == 1){
	//		//co++;	// added
	//		//temp[co] = k;	// new added
	//		objvalue += f[k];
	//		printf("%d ",k+1);
	//	}
	//	for (m = 0; m < NN; m++)
	//		objvalue += W[k][m] * (c_c[k][best_assigmnent[k]] + c_t[best_assigmnent[k]][best_assigmnent[m]] + c_d[best_assigmnent[m]][m]);
	//}
	//printf("\n");

	//for (k = 0; k < NN; k++)
	//	printf("z[%d][%d]=1 \n", k + 1, best_assigmnent[k]+1);
	//
	//printf(" objvalue: %.12f %.2f \n", UpperBound, objvalue);

	end = clock();
	cputime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Time to solve LS: %.2f \n", cputime);

	/*out = open_file("results_ILS.txt", "a+");
	fprintf(out, " %.3f\n", UpperBound);
	fclose(out);*/

	free(best_assigmnent1);
	free(best_open_plants);
	free(which_open_plants);
	free(modified);

	return UpperBound;

} // end of Local_Search() function

int clients_swap1_f(int *assignment, double *vs) {

	register i, j;
	double      MinCostChange = MAX_DOUBLE;
	double      CanMovCostChange;
	int      client1, client2;
	int      plant1, plant2;
	int      flag = 0;

	for (i = 0; i<NN; i++) {
		for (j = 0; j<NN; j++) {
			if (i != j&&assignment[i] != assignment[j]) {
				if (avail_capacity[assignment[i]] + O[i] - O[j] >= 0 && avail_capacity[assignment[j]] + O[j] - O[i] >= 0) {
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

	if (MinCostChange<0) {
		avail_capacity[plant1] += O[client1] - O[client2];
		avail_capacity[plant2] += O[client2] - O[client1];
		assignment[client1] = plant2;
		assignment[client2] = plant1;
		*vs += MinCostChange;
		return 1;
	}
	return 0;
}


int clients_shift1_f(int *assignment, double *vs) {

	register i, k;
	double      MinCostChange = MAX_DOUBLE;
	double      CanMovCostChange;
	int      client;
	int      plant1, plant2;
	//int      *open_plants;
	int      flag = 0;

//	open_plants = row_create(NN);
	/*for (k = 0; k<NN; k++) {
		if (avail_capacity[k]<b[k])
			open_plants[k] = 1;
		else
			open_plants[k] = 0;
	}*/

	for (i = 0; i<NN; i++) {
		for (k = 0; k<NN; k++) {
			if (open_plants[k] && k != assignment[i]) {
				if (avail_capacity[k] - O[i] >= 0) {
					//CanMovCostChange=suma_F_ijkm(i,k, assignment)-suma_F_ijkm(i, assignment[i], assignment);
					CanMovCostChange = sum_F_ijkm(i, k, assignment) - sum_F_ijkm(i, assignment[i], assignment);
					if (b[assignment[i]][0] == avail_capacity[assignment[i]] + O[i] )
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

	if (MinCostChange<0) {
		avail_capacity[plant1] += O[client];
		avail_capacity[plant2] -= O[client];
		if (b[assignment[client]][0] == avail_capacity[assignment[client]])
			open_plants[plant1] = 0;
		assignment[client] = plant2;
		*vs += MinCostChange;
		return 1;
	}
	return 0;
}


int open_hub(int *assignment, double *vs)
{
	int i, j, k, m;
	double   MinCostChange;
	double   valor;
	int      plant1, node;
	int      flag = 0;

	MinCostChange = *vs;

	for (i = 0; i<NN; i++) {
		if (open_plants[i] == 0 && b[i][0] - O[i] >= 0) {
			//open_plants[i] = 1;
			memcpy(allocation, assignment, NN* sizeof(int));
			memcpy(capacity, avail_capacity, NN* sizeof (double));
			capacity[allocation[i]] += O[i];
			allocation[i] = i;
			capacity[i] -= O[i];
			for (j = 1; j<NN; j++) {
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
			for (k = 0; k<NN; k++) {
				if (capacity[k] < b[k][0])
					valor += f[k][0];
				for (m = 0; m<NN; m++){
					//valor+=F_ijkm[k][m][allocation[k]][allocation[m]];
					valor += W[k][m] * (c_c[k][allocation[k]] + c_t[allocation[k]][allocation[m]] + c_d[allocation[m]][m]);
				}
			}

			if (valor < MinCostChange) {
				memcpy(best_allocation, allocation, NN* sizeof(int));
				memcpy(best_capacity, capacity, NN* sizeof (double));
				plant1 = i;
				MinCostChange = valor;
			}
			//open_plants[i] = 0;
		}
	}

	/*for (k = 0; k < NN; k++) {
		if (open_plants[k] == 1){
			printf("%d ", k + 1);
		}
	}
	printf("\n");*/

	if (MinCostChange < *vs) {
		memcpy(assignment, best_allocation, NN* sizeof(int));
		memcpy(avail_capacity, best_capacity, NN* sizeof (double));
		open_plants[plant1] = 1;
		*vs = MinCostChange;
		/*printf("set of hubs: ");
		for (i = 0; i < NN; i++) {
			if (open_plants[i] == 1)
				printf("%d ", i + 1);
		}
		printf("\n");*/
		return 1;
	}

	return 0;

}



int open_close_hub(int *assignment, double *vs)
{
	int i, j, k, m, kk;
	double   MinCostChange;
	double   valor;
	int      plant1;
	int      plant2, node1, node2, node;
	double   min_val;
	int      flag = 0;

	MinCostChange = *vs;

	for (k = 0; k<NN; k++) {
	  for (m = 0; m<NN; m++) {
		  if (open_plants[k] == 0 && b[k][0] - O[k] > 0 && open_plants[m] == 1) {
			flag = 0;
			open_plants[k] = 1;
			open_plants[m] = 0;
			memcpy(allocation, assignment, NN* sizeof(int));
			memcpy(capacity, avail_capacity, NN* sizeof (double));
			node1 = allocation[k];
			capacity[allocation[k]] += O[k];
			allocation[k] = k;
			capacity[k] -= O[k];
			for (j = 1; j<NN; j++) {
				node = costoso[m].A[j].j;
				if (capacity[node] < b[node][0] && capacity[node] - O[m] >= 0){
					capacity[node] -= O[m];
					allocation[m] = node;
					capacity[m] = b[m][0];
					flag = 1;
					break;
				}
			}
			if (flag == 1){
				for (j = 0; j<NN; j++) {
					node = orden_O[NN - j - 1].hub;
					if (allocation[node] == m) {
						for (kk = 1; kk<NN; kk++){
							node2 = costoso[node].A[kk].j;
							if (capacity[node2] < b[node2][0] && capacity[node2] - O[node] >= 0 && node2 != m){
								capacity[node2] -= O[node];
								allocation[node] = node2;
								break;
							}
						}
						if (allocation[node] == m){
							flag = 0;
							break;
						}
					}
				}
			}
			if (flag == 1){
				for (j = 1; j<NN; j++) {
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
					if (open_plants[i] == 1)
						valor += f[i][0];
					for (j = 0; j < NN; j++) {
						valor += W[i][j] * (c_c[i][allocation[i]] + c_t[allocation[i]][allocation[j]] + c_d[allocation[j]][j]);
					}
				}
				if (valor < MinCostChange) {
					memcpy(best_allocation, allocation, NN* sizeof(int));
					memcpy(assignment, best_allocation, NN* sizeof(int));
					memcpy(avail_capacity, capacity, NN* sizeof (double));
					MinCostChange = valor;
					*vs = MinCostChange;
					/*printf("set of hubs: ");
					for (i = 0; i < NN; i++) {
						if (open_plants[i] == 1)
							printf("%d ", i + 1);
					}
					printf("\n");
					for (k = 0; k < NN; k++)
						avail_capacity[k] = b[k];
					for (k = 0; k < NN; k++){
						avail_capacity[assignment[k]] -= O[k];
						printf("z[%d][%d]=1 \n", k + 1, assignment[k] + 1);
					}
					for (k = 0; k < NN; k++)
						printf("cap %d: (%.2f, %.2f) \n", k + 1, avail_capacity[k], b[k]);*/

					return 1;
				}
			}
			open_plants[k] = 0;
			open_plants[m] = 1;
			allocation[k] = node1;
		}
	  }
	}

	/*if (MinCostChange < *vs) {
		memcpy(assignment, best_allocation, NN* sizeof(int));
		open_plants[plant1] = 1;
		*vs = MinCostChange;
		return 1;
	}*/

	return 0;

}




int close_hub(int *assignment, double *vs)
{
	int i, j, k, m;
	double   MinCostChange;
	double   valor, min_val;
	int      plant1;
	int      flag;
	int      node, node2;


	MinCostChange = *vs;

	/*for (k = 0; k < NN; k++) {
		if (open_plants[k] == 1){
			printf("%d ", k + 1);
		}
	}
	printf("\n");*/

	for (i = 0; i < NN; i++) {
		//if (avail_capacity[i] < b[i]) {
		if (open_plants[i] == 1) {
			flag = 0;
			//open_plants[i] = 0;
			memcpy(allocation, assignment, NN* sizeof(int));
			memcpy(capacity, avail_capacity, NN* sizeof (double));
			for (j = 1; j<NN; j++) {
				node = costoso[i].A[j].j;
				if (capacity[node] < b[node][0] && capacity[node] - O[i] >= 0){
					capacity[node] -= O[i];
					allocation[i] = node;
					capacity[i] = b[i][0];
					flag = 1;
					break;
				}
			}
			if (flag == 1){
				for (j = 0; j<NN; j++) {
					node = orden_O[NN - j - 1].hub;
					if (allocation[node] == i) {
						for (k = 1; k<NN; k++){
							node2 = costoso[node].A[k].j;
							if (capacity[node2] < b[node2][0] && capacity[node2] - O[node] >= 0 && node2 != i){
								capacity[node2] -= O[node];
								allocation[node] = node2;
								break;
							}
						}
						if (allocation[node] == i){
							flag = 0;
							break;
						}
					}
				}
			}
			if (flag == 1){
				valor = 0;
				for (k = 0; k<NN; k++) {
					if (capacity[k] < b[k][0])
						valor += f[k][0];
					for (m = 0; m<NN; m++){
						//valor+=F_ijkm[k][m][allocation[k]][allocation[m]];
						valor += W[k][m] * (c_c[k][allocation[k]] + c_t[allocation[k]][allocation[m]] + c_d[allocation[m]][m]);
					}
				}
				if (valor < MinCostChange) {
					memcpy(best_allocation, allocation, NN* sizeof (int));
					memcpy(best_capacity, capacity, NN* sizeof (double));
					plant1 = i;
					MinCostChange = valor;
				}
			}
		   // open_plants[i] = 1;
		}
	}

	/*for (k = 0; k < NN; k++) {
		if (open_plants[k] == 1){
			printf("%d ", k + 1);
		}
	}
	printf("\n");*/

	if (MinCostChange < *vs) {
		memcpy(assignment, best_allocation, NN* sizeof(int));
		memcpy(avail_capacity, best_capacity, NN* sizeof (double));
		open_plants[plant1] = 0;
		*vs = MinCostChange;
		/*printf("set of hubs: ");
		for (i = 0; i < NN; i++) {
			if (open_plants[i] == 1)
				printf("%d ", i + 1);
		}
		printf("\n");
		for (k = 0; k < NN; k++)
			avail_capacity[k] = b[k];
		for (k = 0; k < NN; k++){
			avail_capacity[assignment[k]] -= O[k];
			printf("z[%d][%d]=1 \n", k + 1, assignment[k] + 1);
		}
		for (k = 0; k < NN; k++)
			printf("cap %d: (%.2f, %.2f) \n", k + 1, avail_capacity[k], b[k]);*/
		return 1;
	}

	return 0;

}

double sum2_F_ijkm(int i, int j, int k, int m, int *assignment) {

	int t;
	double suma = 0;

	for (t = 0; t<NN; t++) {
		if (i != t&&j != t)
			suma += W[i][t] * (c_c[i][k] + c_t[k][assignment[t]] + c_d[assignment[t]][t]) + W[t][i] * (c_c[t][assignment[t]] + c_t[assignment[t]][k] + c_d[k][i]);
		else {
			if (i == t)
				suma += W[i][t] * (c_c[i][k] + c_t[k][k] + c_d[k][t]);
		}
		if (j != t&&i != t)
			suma += W[j][t] * (c_c[j][m] + c_t[m][assignment[t]] + c_d[assignment[t]][t]) + W[t][j] * (c_c[t][assignment[t]] + c_t[assignment[t]][m] + c_d[m][j]);
		else {
			if (j == t)
				suma += W[j][t] * (c_c[j][m] + c_t[m][m] + c_d[m][t]);
		}
	}
	suma += W[i][j] * (c_c[i][k] + c_t[k][m] + c_d[m][j]) + W[j][i] * (c_c[j][m] + c_t[m][k] + c_d[k][i]);

	return suma;
}


double sum_F_ijkm(int i, int k, int *assignment) {

	int j;
	double suma = 0;

	for (j = 0; j<NN; j++) {
		if (i != j)
			suma += W[i][j] * (c_c[i][k] + c_t[k][assignment[j]] + c_d[assignment[j]][j]) + W[j][i] * (c_c[j][assignment[j]] + c_t[assignment[j]][k] + c_d[k][i]);
		else
			suma += W[i][j] * (c_c[i][k] + c_t[k][k] + c_d[k][j]);
	}

	return suma;
}
