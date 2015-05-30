## Requirements

You will need:

* g++
* CMake
* zlib
* the Succinct Data Structure Library, https://github.com/simongog/sdsl-lite
* Jellyfish 2.2.0, https://github.com/gmarcais/Jellyfish

## Build instructions

The command

    make

will build the binary `bt` in the `src/` directory.  You can specify header
paths in INCLUDES and library paths in LIBS, e.g. `make INCLUDES=... LIBS=...`.

## Example build for Linux

These instructions work on a blank Ubuntu 14.04 machine.

You will need the following Ubuntu packages:

    g++ git zlib1g-dev cmake

Then, install SDSL:

    cd ~/
    git clone https://github.com/simongog/sdsl-lite.git
    cd sdsl-lite && ./install.sh
    cd

Install Jellyfish:

    curl -L -O https://github.com/gmarcais/Jellyfish/releases/download/v2.2.0/jellyfish-2.2.0.tar.gz
    tar xzf jellyfish-2.2.0.tar.gz
    cd jellyfish-2.2.0
    ./configure && make

And build bloomtree:

    cd
    git clone https://github.com/ctb/bloomtree.git -b aws-build
    cd bloomtree/src
    make INCLUDES="-I ~/include/ -I ~/jellyfish-2.2.0/include" \
         LIBS="-L ~/lib -L ~/jellyfish-2.2.0/.libs"

