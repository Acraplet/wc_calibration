CXXFLAGS=`root-config --cflags --glibs` -fPIC -I../include
LIBDIR=../lib
OBJECTS := chisq.o
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS :=  $(shell root-config --cflags)


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
