CXXFLAGS=`root-config --cflags --libs`

all: make_mPMTmap make_onePositionMap Fitter

make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
Fitter: Fitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/Fitter Fitter.cc
