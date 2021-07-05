LIB = /encs/pkg/cplex-12.9.0/root/cplex
CXX=gcc -O2
CPLXST=-DIL_STD  -L$(LIB)/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -ldl
LIB2=-I$(LIB)/include/ilcplex
OBJ= Benders_framework.o main.o Matheuristic.o Model_info.o ReadData.o
B4QFL : $(OBJ)
	$(CXX) -I$(LIB)/include/ilcplex $(OBJ) -o B4QFL $(CPLXST)

Benders_framework.o : Benders_framework.c def.h
	$(CXX) -c Benders_framework.c $< $(LIB2)
main.o : main.c def.h 
	$(CXX) -c main.c $< $(LIB2)
Matheuristic.o : Matheuristic.c def.h
	$(CXX) -c Matheuristic.c $< $(LIB2)
Model_info.o : Model_info.c def.h
	$(CXX) -c Model_info.c $< $(LIB2)
ReadData.o : ReadData.c def.h 
	$(CXX) -c ReadData.c $< $(LIB2)

clean :
	rm B4QFL $(OBJ)

