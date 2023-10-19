CXXFLAGS=$(shell $(ROOTSYS)/bin/root-config --cflags)
ifneq (,$(findstring -std=c,$(CXXFLAGS)))
$(info Using ROOT C++ standard)
$(info $(CXXFLAGS))
else
$(info ROOT C++ standard not found. Default to c++11)
CXXFLAGS := $(CXXFLAGS) -std=c++11
$(info $(CXXFLAGS))
endif
LIBS := $(shell $(ROOTSYS)/bin/root-config --libs) -L ../lib -lChisq
INCLUDES := -I ../chisq -I ../include


all: make_mPMTmap_binned splitReference_intoMaps AttenuationFitterPMT-basedBinned 


make_mPMTmap_binned: make_mPMTmap_binned.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -g -o ../bin/$@ $(LIBS) $(INCLUDES)
splitReference_intoMaps: splitReference_intoMaps.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) $(INCLUDES)
AttenuationFitterPMT-basedBinned: AttenuationFitterPMT-basedBinned.cc
	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) $(INCLUDES)


#ifdef WCSIMDIR
#INCLUDES := $(INCLUDES) -I $(WCSIMDIR)/include
#all: ScatteringFitter2DBinned AttenuationFitterBinned make_mPMTmap makeRef_maxNbHitAtPos WCSIM_TreeConvert TestFitter 
#else
#TestFitter ScatteringFitter2DBinned AttenuationFitterBinned make_mPMTmap makeRef_maxNbHitAtPos TestFitter
#endif

#these are the old codes that aren't used in the binned approach, they are now in the mPMTmapping/Archive/OldCodes folder if one is interested
#BinnedFitter Fitter 
#AttenuationFitter make_oneBin_reference make_oneBin_test Fitter 
#ScatteringFitter make_onePositionMap ScatteringFitterBinned

#plot_Scattering_Absorption_ExtrapolateR: plot_Scattering_Absorption_ExtrapolateR.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#make_mPMTmap: make_mPMTmap.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
#makeRef_maxNbHitAtPos: makeRef_maxNbHitAtPos.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
#AttenuationFitterBinned: AttenuationFitterBinned.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#ScatteringFitter2DBinned: ScatteringFitter2DBinned.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#TestFitter: TestFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#WCSIM_TreeConvert: WCSIM_TreeConvert.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES) -L$(WCSIMDIR) -lWCSimRoot -I$(WCSIMDIR)/include
#


#And here the older codes for 1D and simplistic spline analaysis 
#make_onePositionMap: make_onePositionMap.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_onePositionMap make_onePositionMap.cc
#make_oneBin_reference: make_oneBin_reference.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_reference make_oneBin_reference.cc -I$(INCLUDES) $(LIBS)
#make_oneBin_test: make_oneBin_test.cc
#	g++ -Wall -Wextra $(CXXFLAGS) -o ../bin/make_oneBin_test make_oneBin_test.cc -I$(INCLUDES) $(LIBS)
#
#Fitter: Fitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ $(LIBS) -I$(INCLUDES)
#BinnedFitter: BinnedFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o ../bin/$@ -ggdb3 $(LIBS) -I$(INCLUDES)
#AttenuationFitter: AttenuationFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) AttenuationFitter.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#ScatteringFitter: ScatteringFitter.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
#ScatteringFitterBinned: ScatteringFitterBinned.cc
#	g++ -Wall -Wextra $(CXXFLAGS) $@.cc -o../bin/$@ $(LIBS) -I$(INCLUDES)
