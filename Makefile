CFLAGS = -Wall -std=c++17
INC      = /Library/gurobi811/mac64/include/
CPPLIB   = -L/Library/gurobi811/mac64/lib -lgurobi_c++ -lgurobi81

all: build

build: NFoldSolve

clean:
	rm -f main.o nfold.o objVals.txt times.txt

NFoldSolve: main.o nfold.o
	g++ $(CFLAGS) -o NFoldSolve main.o nfold.o -I$(INC) $(CPPLIB)

main.o: main.cc
	g++ $(CFLAGS) -c main.cc -I$(INC)

nfold.o: nfold.cc nfold.h
	g++ $(CFLAGS) -c nfold.cc -I$(INC)
