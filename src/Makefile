CXXFLAGS=-Wall -pedantic -std=c++11 -DNDEBUG -O3
LDFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 -ljellyfish-2.0 -lpthread -lz

all: clean bt

bt: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o bt
# DO NOT DELETE
