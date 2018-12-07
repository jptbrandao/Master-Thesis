# Make file test


# Variables to control makefile

CXX = g++
CXXFlags = -Wall -g -std=c++11
MosekHEADER = ~/projects/github/mestrado/API/mosek/8/tools/platform/osx64x86/h/mosek.h
MosekSHARED = ~/projects/github/mestrado/API/mosek/8/tools/platform/osx64x86/bin -Wl,-headerpad,128
LIBS = -lmosek64

# Targets need to bring executables up to date

main: main.o quicksort.o weaklyPoly.o testQuickSort.o Data.o Solution.o testSolver.o scopt-ext.o commonFuncs.o testweakly.o experiments.o
	$(CXX) $(CXXFlags) -L../../../APis/mosek/bin -Wl,-headerpad,128 -o main main.o quicksort.o weaklyPoly.o testQuickSort.o Data.o Solution.o testSolver.o scopt-ext.o commonFuncs.o testweakly.o experiments.o -lmosek64
	install_name_tool -change libmosek64.8.1.dylib ~/projects/APIs/mosek/bin/libmosek64.8.1.dylib $@
	@- make clean
	@- ./main

main.o: main.cpp quicksort.hpp testFiles/TestQuickSort.hpp subprocedureTests/greedy2Test.hpp testFiles/testWeakly.hpp
	$(CXX) $(CXXFlags) -c main.cpp

quicksort.o: quicksort.cpp quicksort.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp
	$(CXX) $(CXXFlags) -c quicksort.cpp quicksort.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp

weaklyPoly.o: weaklypoly.cpp weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp
	$(CXX) $(CXXFlags) -c weaklypoly.cpp weaklypoly.hpp DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.hpp

testQuickSort.o: testFiles/TestQuickSort.cpp testFiles/TestQuickSort.hpp DataStructures/Data.hpp
	$(CXX) $(CXXFlags) -c testFiles/TestQuickSort.cpp testFiles/TestQuickSort.hpp

testweakly.o: testFiles/testWeakly.cpp testFiles/testWeakly.hpp DataStructures/Data.hpp 
	$(CXX) $(CXXFlags) -c testFiles/testWeakly.cpp testFiles/testWeakly.hpp

Data.o: DataStructures/Data.hpp DataStructures/Data.cpp DataStructures/Solution.hpp
	$(CXX) $(CXXFlags) -c DataStructures/Data.cpp DataStructures/Data.hpp DataStructures/Solution.hpp

Solution.o: DataStructures/Solution.hpp DataStructures/Solution.cpp
	$(CXX) $(CXXFlags) -c DataStructures/Solution.cpp DataStructures/Solution.hpp

testSolver.o: testSolver.hpp testSolver.cpp scopt-ext.h
	$(CXX) $(CXXFlags) -c testSolver.cpp testSolver.hpp scopt-ext.h 

scopt-ext.o: scopt-ext.h scopt-ext.cpp
	$(CXX) $(CXXFlags) -c scopt-ext.cpp scopt-ext.h mosek.h 

#greedy2Test.o: subprocedureTests/greedy2Test.hpp DataStructures/Data.hpp DataStructures/Solution.hpp quicksort.hpp testSolver.hpp
#	$(CXX) $(CXXFlags) -c subprocedureTests/greedy2Test.cpp subprocedureTests/greedy2Test.hpp DataStructures/Data.hpp DataStructures/Solution.hpp quicksort.hpp testSolver.hpp

commonFuncs.o: DataStructures/Solution.hpp DataStructures/Data.hpp commonFuncs.cpp commonFuncs.hpp
	$(CXX) $(CXXFlags) -c commonFuncs.cpp DataStructures/Solution.hpp DataStructures/Data.hpp 

experiments.o: weaklypoly.hpp quicksort.hpp testSolver.hpp DataStructures/Solution.hpp DataStructures/Data.hpp
	$(CXX) $(CXXFlags) -c experiments.cpp weaklypoly.hpp quicksort.hpp testSolver.hpp DataStructures/Solution.hpp DataStructures/Data.hpp

clean:
	@- rm -f *.o
	@- rm -f *.gch
	@- rm -f subprocedureTests/*.o
	@- rm -f subprocedureTests/*.gch
	@- rm -f DataStructures/*.o
	@- rm -f DataStructures/*.gch
	@- rm -f testFiles/*.o
	@- rm -f testFiles/*.gch
