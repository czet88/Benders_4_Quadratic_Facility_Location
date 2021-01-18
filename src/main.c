#include "def.h"

 double     **c, **c_c, **c_t, **c_d, **f, **W, *O, *D, **b;
 double     collect,transfer,distribute, AD;
 int        NN, Q, p_hubs;
 double     MAX_DOUBLE = 10000000000;
 double     UpperBound;
 coordinate *pts;
 int        **pos_z;
 int        pos_eta;
 double     old_lower_bound;
 double     ***alpha;
 double     ***beta;
 double     **core;
 double     LP_lower_bound;
 SCUTS      sepcut;
 INCUTS     *initial_cuts;
 double     sum_core;
 int        MG;
 double     *initial_x;
 int        *allocation;
 int        *best_allocation;
 int        *open_plants;
 double     *capacity;
 double     *best_capacity;
 double     *avail_capacity;
 double     menor_O;
 PORD       *costoso;
 SELEC      *orden_O;
 double     AggregatedDemand;
 ZVAL       *z_open;
 ZVAL       *z_closed;
 int        *cand_hubs;
 int        count_cand_hubs;
 int        *fixed_zero;
 int        *fixed_one;
 int        count_added;
 COV        cover;
 NORD       *ord_nodes;
 /*item       *mochila;
 stype      mochilaC;
 stype      mochilaZ;*/
 double     *multiple;
 double        *coeff_ES;
 ORD        *ord_O;
 ORD        *ord_D;
 double     sum_supply_i, sum_supply_j;
 
 

 double     old_objval = 0;
 int        count_same_node = 0;

 int main(int argc, char* argv[])
 {
	 int i, num_inst, count_opt, pp;
	 FILE* ini;
	 FILE* out;

	 clock_t  start, end;
	 char instance[10];
	 double act_gap, cputime_heur, cputime;
	 int trans_fact;
	 double coll, trans, distr;
	 double UB_heur;
	 char fake_text[20];

	 coll = 0;
	 trans = 0;
	 distr = 0;
	 combtol = 0;
	 count_opt = 0;
	 old_objval = 0;

	 //Default heuristic solution parameters
	 FlagHeuristic = 1;
	 FlagLocalSearch = 1;
	 FlagIteratedLocalSearch = 1;

	 //Read command line arguments
	 if (argc < 3) {
		 printf("Error: Input file and output file not specified \n");
		 exit(8);
	 }
	 //read output file name
	 sprintf(output_text, argv[2]);

	 // If we're given a heuristic parameter file, then read it
	 if (argc >= 4) read_heur_cl_param(argv[3]);
	 printf("FlagHeur %d,\nFlagLS %d,\nFlagILS %d\n", FlagHeuristic, FlagLocalSearch, FlagIteratedLocalSearch);

	 ini = open_file(argv[1], "r");
	 //ini=open_file ("25ll.txt","r");
	 fscanf(ini, "%d", &num_inst);
	 fscanf(ini, "%s", &fake_text);

	 /*Printing the time and date these are being executed*/
	 /************************************************/
	 t = time(NULL);
	 tm = localtime(&t);
	 out = open_file(output_text, "a+");
	 fprintf(out, "\n %s\nVers;instance;APset;Cap;p;hybrid;heurParam;UBPre;LBPre;CPU_Pre;Num_Iter;Num_fixed;UB;LB;time_BC;GAP;BBnodes;Hubs;CPU_all;Missed; CPUFLP;CPUGAss\n", asctime(tm));
	 fclose(out);
	 for (i = 0; i < num_inst; i++) {
		 use_firstsolution = 1;
		 missed = 0;
		 fscanf(ini, "%d", &vers);
		 fscanf(ini, "%d", &APset);
		 fscanf(ini, "%d", &Capacitated_instances);
		 fscanf(ini, "%d", &p_hubs);
		 fscanf(ini, "%d", &hybrid);
		 fscanf(ini, "%s", &instance);
		 //fscanf(ini,"%d %lf %lf %lf",&pp,&coll,&trans,&distr);
		 read_instance(instance, 5 - 4 * APset, coll, trans, distr, APset);
		 printf(" %s  alpha= %0.2f \n", instance, trans);
		 //Obtain initial solution from a SSCFLP
		 start = clock();
		 UpperBound = MAX_DOUBLE;
		 // if(vers!=0) SSCFLP_model();
		 // if(vers!=1) UB_heur = Det_Iterated_Local_Search(); //Improve solution with a local search
		 // end = clock();
		 // cputime_heur = (double)(end - start) / CLOCKS_PER_SEC;
		 int totalFlag = FlagHeuristic + FlagLocalSearch + FlagIteratedLocalSearch;
		 out = open_file(output_text, "a+");
		 fprintf(out, "%d;%s;%d;%d;%d;%d;%d;", vers, instance, APset, Capacitated_instances, p_hubs, hybrid, totalFlag);
		 fclose(out);
		 //Solve root node
		 Benders_root_node_heur();
		 //Benders_root_node();
		 //Benders_BC();							//Solve Benders reformulation with branch and cut algorithm
		 end = clock();
		 cputime = (double)(end - start) / CLOCKS_PER_SEC;
		 printf("\nFInal CPU time: %.2lf\n", cputime);
		 out = open_file(output_text, "a+");
		 fprintf(out, "%.2f; %d; %.2lf; %.2lf\n", cputime, missed, cpuFacLocIni, cpuGenAss);
		 fclose(out);
		 //CHLPSA_model();

		 //Original_Model();  //Solve four-index formulation directly with CPLEX
		 free_memory();
		 printf("Finished solving\n");
		 //getchar();
	 }
	 fclose(ini);

	 return 0;
 }

void read_heur_cl_param(const char* name) {
	/* Modify global parameters */
	if (strcmp(name, "NoHeur") == 0) {
		FlagHeuristic = 0;
		FlagLocalSearch = 0;
		FlagIteratedLocalSearch = 0;
	}

	if (strcmp(name, "Heur") == 0) {
		FlagLocalSearch = 0;
		FlagIteratedLocalSearch = 0;
	}

	if (strcmp(name, "LS") == 0) {
		FlagIteratedLocalSearch = 0;
	}

	if (strcmp(name, "ILS") == 0) {
		FlagHeuristic = 1;
		FlagLocalSearch = 1;
		FlagIteratedLocalSearch = 1;
	}
}


