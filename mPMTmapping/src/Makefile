CXXFLAGS=`root-config --cflags --libs`

all: make_mPMTmap

make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
