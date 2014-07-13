CXXFLAGS=-Wall -pedantic -std=c++11 -DNDEBUG -O3
LDFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 -ljellyfish-2.0

bt: main.o Query.o Kmers.o BloomTree.o BF.o util.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o bt
# DO NOT DELETE
