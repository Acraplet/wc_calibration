CXXFLAGS=$(shell $(ROOTSYS)/bin/root-config --cflags --glibs) -fPIC -I../include
ifneq (,$(findstring -std=c,$(CXXFLAGS)))
$(info Using ROOT C++ standard)
$(info $(CXXFLAGS))
else
$(info ROOT C++ standard not found. Default to c++11)
CXXFLAGS := $(CXXFLAGS) -std=c++11
$(info $(CXXFLAGS))
endif
LIBDIR=../lib
OBJECTS := chisq.o
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS :=  $(shell root-config --glibs)


%.o:%.cc
	g++ -c $(CXXFLAGS) $<

all: libChisq

libChisq: libChisq.a $(LIBDIR)/libChisq.so

libChisq.a: $(OBJECTS)
	ar rcsv $@ $^

$(LIBDIR)/libChisq.so: libChisq.a $(OBJECTS)
	g++ -shared -Wl,-export-dynamic $(ROOTCFLAGS) $^ -o $@ $(ROOTLIBS) 


clean:
	rm -f ../lib/chisq.so
	rm -f *.o lib*
