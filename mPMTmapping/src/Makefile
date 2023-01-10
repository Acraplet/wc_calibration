CXXFLAGS=`root-config --cflags --libs`
LIBS := -L ../lib -lChisq
INCLUDES := ../chisq -I ../include

all: BinnedFitter Fitter
       #	make_mPMTmap AttenuationFitter make_oneBin_reference Fitter makeRef_maxNbHitAtPos ScatteringFitter make_onePositionMap

make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
make_oneBin_reference: make_oneBin_reference.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_reference make_oneBin_reference.cc -I$(INCLUDES) $(LIBS)
makeRef_maxNbHitAtPos: makeRef_maxNbHitAtPos.cc
	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/makeRef_maxNbHitAtPos makeRef_maxNbHitAtPos.cc
Fitter: Fitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
BinnedFitter: BinnedFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
AttenuationFitter: AttenuationFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
ScatteringFitter: ScatteringFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
