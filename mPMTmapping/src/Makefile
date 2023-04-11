CXXFLAGS=$(shell $(ROOTSYS)/bin/root-config --cflags --libs)
ifneq (,$(findstring -std=c,$(CXXFLAGS)))
$(info Using ROOT C++ standard)
$(info $(CXXFLAGS))
else
$(info ROOT C++ standard not found. Default to c++11)
CXXFLAGS := $(CXXFLAGS) -std=c++11
$(info $(CXXFLAGS))
endif
LIBS := -L ../lib -lChisq
INCLUDES := ../chisq -I ../include

ifdef WCSIMDIR
all: ScatteringFitter2DBinned BinnedFitter Fitter make_mPMTmap AttenuationFitter make_oneBin_reference make_oneBin_test Fitter makeRef_maxNbHitAtPos ScatteringFitter make_onePositionMap ScatteringFitterBinned TestFitter WCSIM_TreeConvert
else
all: ScatteringFitter2DBinned BinnedFitter Fitter make_mPMTmap AttenuationFitter make_oneBin_reference make_oneBin_test Fitter makeRef_maxNbHitAtPos ScatteringFitter make_onePositionMap ScatteringFitterBinned TestFitter
endif

plot_Scattering_Absorption_ExtrapolateR: plot_Scattering_Absorption_ExtrapolateR.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
make_mPMTmap: make_mPMTmap.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
make_onePositionMap: make_onePositionMap.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
make_oneBin_reference: make_oneBin_reference.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
make_oneBin_test: make_oneBin_test.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
makeRef_maxNbHitAtPos: makeRef_maxNbHitAtPos.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
Fitter: Fitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
BinnedFitter: BinnedFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ -ggdb3 $(LIBS) -I$(INCLUDES)
AttenuationFitter: AttenuationFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
ScatteringFitter: ScatteringFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
ScatteringFitterBinned: ScatteringFitterBinned.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
ScatteringFitter2DBinned: ScatteringFitter2DBinned.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
TestFitter: TestFitter.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
WCSIM_TreeConvert: WCSIM_TreeConvert.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES) -L$(WCSIMDIR) -lWCSimRoot -I$(WCSIMDIR)/include



#CXXFLAGS=`root-config --cflags --libs` -std=c++17
#LIBS := -L ../lib -lChisq
#INCLUDES := ../chisq -I ../include
#
#all: AttenuationFitterBinned AttenuationFitter make_reference_abs
#	#ScatteringFitter2DBinned BinnedFitter Fitter make_mPMTmap AttenuationFitter make_oneBin_reference make_oneBin_test Fitter makeRef_maxNbHitAtPos ScatteringFitter make_onePositionMap
##ScatteringFitter2DBinned
#
#make_reference_abs: make_reference_abs.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_reference_abs make_reference_abs.cc -I$(INCLUDES) $(LIBS)
#plot_Scattering_Absorption_ExtrapolateR: plot_Scattering_Absorption_ExtrapolateR.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#make_mPMTmap: make_mPMTmap.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_mPMTmap make_mPMTmap.cc
#make_onePositionMap: make_onePositionMap.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
#make_oneBin_reference: make_oneBin_reference.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_reference make_oneBin_reference.cc -I$(INCLUDES) $(LIBS)
#make_oneBin_noAttenuation_reference: make_oneBin_noAttenuation_reference.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_noAttenuation_reference make_oneBin_noAttenuation_reference.cc -I$(INCLUDES) $(LIBS)
#make_oneBin_test: make_oneBin_test.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_test make_oneBin_test.cc -I$(INCLUDES) $(LIBS)
#makeRef_maxNbHitAtPos: makeRef_maxNbHitAtPos.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/makeRef_maxNbHitAtPos makeRef_maxNbHitAtPos.cc
#Fitter: Fitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
#BinnedFitter: BinnedFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ -ggdb3 $(LIBS) -I$(INCLUDES)
#AttenuationFitter: AttenuationFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#AttenuationFitterBinned: AttenuationFitterBinned.cc
#	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitterBinned.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#ScatteringFitter: ScatteringFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
##ScatteringFitterBinned: ScatteringFitterBinned.cc
##	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#ScatteringFitter2DBinned: ScatteringFitter2DBinned.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
