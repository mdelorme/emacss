CPP=g++
CFLAGS=-g -O3 -pedantic -Wall #-fastmath
LM=-lm #-std=c99


CLUSTER=input.cpp info.cpp params.cpp cluster.cpp output.cpp
DYN=dynamics.cpp
SE=stellar_evolution.cpp 
DYN_OBJ=$(DYN:.cpp=.o)
CLUSTER_OBJ=$(CLUSTER:.cpp=.o)
SE_OBJ=$(SE:.cpp=.o)

all: emacss tester

emacss: emacss.cpp $(DYN_OBJ) $(CLUSTER_OBJ) $(SE_OBJ)
	$(CPP) $(CFLAGS) emacss.cpp $(DYN_OBJ) $(CLUSTER_OBJ) $(LM) $(SE_OBJ) -o $@

tester: tester.cpp $(DYN_OBJ) $(CLUSTER_OBJ) $(SE_OBJ) 
	$(CPP) $(CFLAGS) tester.cpp $(CLUSTER_OBJ) $(DYN_OBJ) $(LM) $(SE_OBJ) -o $@

clean:
	rm *.o

emacss.o: emacss.cpp
	$(CPP) -g -c $^ -o $@

input.o: emacss_dir/input.cpp
	$(CPP) -g -c $^ -o $@

info.o: emacss_dir/info.cpp
	$(CPP) -g -c $^ -o $@

params.o: emacss_dir/params.cpp
	$(CPP) -g -c $^ -o $@

cluster.o: emacss_dir/cluster.cpp
	$(CPP) -g -c $^ -o $@

output.o: emacss_dir/output.cpp
	$(CPP) -g -c $^ -o $@

dynamics.o: dynamics/dynamics.cpp
	$(CPP) -g -c $^ -o $@

stellar_evolution.o: stevo/stellar_evolution.cpp
	$(CPP) -g -c $^ -o $@
