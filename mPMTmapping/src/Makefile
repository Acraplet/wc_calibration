CXXFLAGS=`root-config --cflags --libs`
LIBS := -L ../lib -lChisq
INCLUDES := ../chisq -I ../include

all: make_mPMTmap AttenuationFitter Fitter makeRef_maxNbHitAtPos ScatteringFitter make_onePositionMap

make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
makeRef_maxNbHitAtPos: makeRef_maxNbHitAtPos.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/makeRef_maxNbHitAtPos makeRef_maxNbHitAtPos.cc
Fitter: Fitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
AttenuationFitter: AttenuationFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
ScatteringFitter: ScatteringFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
