# Compile type
DEBUGFLAG = -g3
RELEASEFLAG = -O3 -march=native

# Compiler handle
CXXFLAGS=  -std=c++0x -I./includes -D__LINUX 
# CXXFLAGS += $(DEBUGFLAG)
CXXFLAGS += $(RELEASEFLAG)

LDFLAG= -L./includes -pthread -Wl,--no-as-needed

# Compiled files
OBJ= ./src/BrownianDynamicEngine.o ./src/BrownianDynamics3D.o ./src/BrownianDynamicsWithEP.o ./src/Main.o 

# Compiler commands
test:$(OBJ)
	g++ -fopenmp -o $@ $^ $(LDFLAG)
test_static:$(OBJ)
	g++ -static -fopenmp -o $@ $^ $(LDFLAG) -lquadmath -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
	mv ./test_static ./Executive\ Package/
# here the linking flag is from https://gcc.gnu.org/ml/gcc-help/2010-05/msg00029.html
%.o : Scripts/%.cpp
	g++ -c $(CXXFLAGS) $^
clean:   
	rm -r -f core* ./Executive\ Package/test_static test */*.o *.il *~ */*~
clear:
	rm -rf ./traj/*.dat ./*.dat
