LIB = /encs/pkg/cplex-12.9.0/root/cplex
CXX=gcc -O2
CPLXST=-DIL_STD  -L$(LIB)/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -ldl
LIB2=-I$(LIB)/include/ilcplex
OBJ= Matheuristic.o main.o ReadData.o solve_LP_Benders_Heur.o
BCHLPSA : $(OBJ)
	$(CXX) -I$(LIB)/include/ilcplex $(OBJ) -o BCHLPSA $(CPLXST)

Matheuristic.o : Matheuristic.c def.h 
	$(CXX) -c Matheuristic.c $< $(LIB2)
main.o : main.c def.h 
	$(CXX) -c main.c $< $(LIB2)
ReadData.o : ReadData.c def.h 
	$(CXX) -c ReadData.c $< $(LIB2)
solve_LP_Benders_Heur.o : solve_LP_Benders_Heur.c def.h
	$(CXX) -c solve_LP_Benders_Heur.c $< $(LIB2)

clean :
	rm BCHLPSA $(OBJ)

