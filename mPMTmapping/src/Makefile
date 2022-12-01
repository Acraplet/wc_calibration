CXXFLAGS=`root-config --cflags --libs`
LIBS := -L ../lib -lChisq
INCLUDES := ../chisq

all: make_mPMTmap make_onePositionMap AttenuationFitter Fitter make_referenceMap

make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
make_referenceMap: make_referenceMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_referenceMap make_referenceMap.cc
Fitter: Fitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
AttenuationFitter: AttenuationFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I ../chisq
