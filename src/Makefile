LIB = /encs/pkg/cplex-12.9.0/root/cplex
CXX=gcc -O2
CPLXST=-DIL_STD  -L$(LIB)/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -ldl
LIB2=-I$(LIB)/include/ilcplex
OBJ= combo.o LocalSerch.o Matheuristic.o main.o ReadData.o solve_UFLP.o solve_Benders_model.o solve_LP_Benders.o solve_LP_Benders_Heur.o
BCHLPSA : $(OBJ)
	$(CXX) -I$(LIB)/include/ilcplex $(OBJ) -o BCHLPSA $(CPLXST)

combo.o: combo.c def.h combo.h
	$(CXX) -c combo.c $< $(LIB2)
LocalSerch.o : LocalSerch.c def.h combo.h
	$(CXX) -c LocalSerch.c $< $(LIB2)
Matheuristic.o : Matheuristic.c def.h combo.h
	$(CXX) -c Matheuristic.c $< $(LIB2)
main.o : main.c def.h combo.h
	$(CXX) -c main.c $< $(LIB2)
ReadData.o : ReadData.c def.h combo.h
	$(CXX) -c ReadData.c $< $(LIB2)
solve_UFLP.o : solve_UFLP.c def.h combo.h
	$(CXX) -c solve_UFLP.c $< $(LIB2)
solve_Benders_model.o : solve_Benders_model.c def.h combo.h
	$(CXX) -c solve_Benders_model.c $< $(LIB2)
solve_LP_Benders.o : solve_LP_Benders.c def.h
	$(CXX) -c solve_LP_Benders.c $< $(LIB2)
solve_LP_Benders_Heur.o : solve_LP_Benders_Heur.c def.h
	$(CXX) -c solve_LP_Benders_Heur.c $< $(LIB2)

clean :
	rm BCHLPSA $(OBJ)

