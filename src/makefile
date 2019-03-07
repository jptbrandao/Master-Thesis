# Make file test


# Variables to control makefile

CXX = g++
CXXFlags = -Wall -g -std=c++11
#MosekHEADER = ~/projects/github/mestrado/API/mosek/8/tools/platform/osx64x86/h/mosek.h
#MosekSHARED = ~/projects/github/mestrado/API/mosek/8/tools/platform/osx64x86/bin -Wl,-headerpad,128
#LIBS = -lmosek64

# Targets need to bring executables up to date

main: main.o weaklyPoly.o Data.o Solution.o commonFuncs.o experiments.o newqs.o searchSpace.o
	$(CXX) $(CXXFlags) -o main main.o weaklyPoly.o Data.o Solution.o commonFuncs.o experiments.o newqs.o searchSpace.o
	@- make clean
	@- ./main

main.o: main.cpp 
	$(CXX) $(CXXFlags) -c main.cpp

weaklyPoly.o: weaklypoly.cpp weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp
	$(CXX) $(CXXFlags) -c weaklypoly.cpp weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp

Data.o: DataStructures/Data.hpp DataStructures/Data.cpp DataStructures/Solution.hpp
	$(CXX) $(CXXFlags) -c DataStructures/Data.cpp DataStructures/Data.hpp DataStructures/Solution.hpp

Solution.o: DataStructures/Solution.hpp DataStructures/Solution.cpp
	$(CXX) $(CXXFlags) -c DataStructures/Solution.cpp DataStructures/Solution.hpp

commonFuncs.o: DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.cpp commonFuncs.hpp
	$(CXX) $(CXXFlags) -c commonFuncs.cpp DataStructures/Solution.hpp DataStructures/Data.hpp 

experiments.o: weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp newQS.hpp
	$(CXX) $(CXXFlags) -c experiments.cpp weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp newQS.hpp

newqs.o: newQS.hpp newQS.cpp DataStructures/Solution.hpp DataStructures/Data.hpp SearchSpace.hpp
	$(CXX) $(CXXFlags) -c newQS.cpp newQS.hpp DataStructures/Solution.hpp DataStructures/Data.hpp SearchSpace.hpp

searchSpace.o: SearchSpace.cpp SearchSpace.hpp
	$(CXX) $(CXXFlags) -c SearchSpace.cpp SearchSpace.hpp

clean:
	@- rm -f *.o
	@- rm -f *.gch
	@- rm -f *.gch*
	@- rm -f subprocedureTests/*.o
	@- rm -f subprocedureTests/*.gch
	@- rm -f subprocedureTests/*.gch*
	@- rm -f DataStructures/*.o
	@- rm -f DataStructures/*.gch
	@- rm -f DataStructures/*.gch*
	@- rm -f testFiles/*.o
	@- rm -f testFiles/*.gch
