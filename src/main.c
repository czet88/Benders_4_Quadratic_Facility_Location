#include "def.h"

/* Print a usage message to stderr and abort. */
static void
usage(const char* progname)
{
	fprintf(stderr, "Usage: %s input_file_name output_file_name [Heuristic_Level] \n", progname);
	fprintf(stderr,
		" By default, the code will use the maximum heuristic level (2)\n"
		" to find a feasible solution at the beginning of the solving \n"
		" process to be able to perform variable elimination"
		" Supported options are:\n"
		"  0       Do not use any heuristic method\n"
		"  1       Use the matheuristic based on \n"
		"          ignoring the quadratic costs\n"
		"  2       Use the matheuristic along with\n"
		"          a local search algorithm\n");
	exit(2);
}

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
	 char input_text[50];
	 coll = 0;
	 trans = 0;
	 distr = 0;
	 combtol = 0;
	 count_opt = 0;
	 old_objval = 0;
	 //Default heuristic solution parameters
	 FlagHeuristic = 1;
	 FlagLocalSearch = 1;
	 vers = -1; //Version -1 contains all the bells and whistles
	 use_firstsolution = 1; //Use the first solution found as warm start for the MIP
	 missed = 0;

	 //Read command line arguments
	 if (argc < 3) {
		 usage(argv[0]);
	 }
	 //read output file name
	 sprintf(output_text, "Results/Raw_Results/");
	 strcat(output_text, argv[2]);
	 // If we're given a heuristic parameter file, then read it
	 if (argc >= 4 && !read_heur_cl_param(argv[3])) {
		 usage(argv[0]);
	 }
	 int totalFlag = FlagHeuristic + FlagLocalSearch;

	 //read input file name
	 sprintf(input_text, "Inputs/");
	 strcat(input_text, argv[1]);
	 printf("Input: %s\nOutput: %s\n", input_text, output_text);
	 //Write the header of the output file
	 out = open_file(output_text, "a+");
	 fprintf(out, "instance;APset;Cap;p;p_median_constr;fixed_costs;heurParam;UBPre;LBPre;CPU_Pre;Num_Iter;Num_fixed;UB;LB;time_BC;GAP;BBnodes;Hubs;CPU_all;Missed; CPUFLP;CPUGAss\n");
	 fclose(out);
	 ini = open_file(input_text, "r");
	 fscanf(ini, "%d", &num_inst);
	 for (i = 0; i < num_inst; i++) {
		 fscanf(ini, "%d", &APset); //{0,1} flag to denote f reading from AP data set
		 fscanf(ini, "%d", &Capacitated_instances); //{0,1} flag to denote if we ignore or not the capacity constraints
		 fscanf(ini, "%d", &p_hubs);//Maximum number of hubs to be opened
		 fscanf(ini, "%d", &w_p_median_constr);// Consider a p_median constraint
		 fscanf(ini, "%d", &w_fixed_costs);// Consider fixed costs
		 fscanf(ini, "%s", &instance); //Name of the instance
		 printf("++++++++++++++++++++++++++++++++++\nSolving %s\n", instance);
		 //Read the data from the instance file
		 read_instance(instance, 5 - 4 * APset, coll, trans, distr, APset);
		 UpperBound = MAX_DOUBLE;
		 //Write the header names into the output file
		 out = open_file(output_text, "a+");
		 fprintf(out, "%s;%d;%d;%d;%d;%d;%d;", instance, APset, Capacitated_instances, p_hubs, w_p_median_constr, w_fixed_costs, totalFlag);
		 fclose(out);
		 start = clock();
		 Benders_framework();
		 end = clock();
		 cputime = (double)(end - start) / CLOCKS_PER_SEC;
		 printf("Final CPU time: %.2lf\n", cputime);
		 out = open_file(output_text, "a+");
		 fprintf(out, "%.2f; %d; %.2lf; %.2lf\n", cputime, missed, cpuFacLocIni, cpuGenAss);
		 fclose(out);
		 free_memory();
		 printf("++++++++++++++++++++++++++++++++++\nFinished solving %s\n", instance);
	 }
	 fclose(ini);
	 return 0;
 }

bool read_heur_cl_param(const char* name) {	
	/* Modify global parameters */
	if (strcmp(name, "0") == 0) {
		FlagHeuristic = 0;
		FlagLocalSearch = 0;
		return true;
	}
	if (strcmp(name, "1") == 0) {
		FlagLocalSearch = 0;
		return true;
	}
	if (strcmp(name, "2") == 0) {
		FlagHeuristic = 1;
		FlagLocalSearch = 1;
		return true;
	}
	return false;
}


