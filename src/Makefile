#CXXFLAGS= -std=c++11 -O3 $(shell pkg-config --cflags jellyfish-2.0) -M -MF $@.dep
#-I/opt/local/stow/sdsl-lite/include
CXXFLAGS= -std=c++11 -O3 $(shell pkg-config --cflags jellyfish-2.0 sdsl-lite)
LDFLAGS=$(shell pkg-config --libs sdsl-lite) $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz $(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L/-Wl,-rpath,/')
#LDFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 $(shell pkg-config --libs jellyfish-2.0 sdsl-lite) -lpthread -lz $(shell pkg-config --libs-only-L jellyfish-2.0 sdsl-lite | sed -e 's/-L/-Wl,-rpath,/g') -fsanitize=address
#LDFLAGS=-lsdsl -ldivsufsort -ldivsufsort64 $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz
#CXXFLAGS=-Wall -v -pedantic -std=c++11 -DNDEBUG -O3
#LDFLAGS=-static -lsdsl -ldivsufsort -ldivsufsort64 -ljellyfish-2.0 -lpthread -lz

#all: clean bt

bt: main.o Build.o Query.o Kmers.o BloomTree.o BF.o util.o Count.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o bt
# DO NOT DELETE
