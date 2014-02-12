ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib -lEG
ROOTGLIBS     = $(shell root-config --glibs) -Lrootlib -lEG

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -g -O2

ARCH         := $(shell root-config --arch)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  -j 4
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = suppress
endif

CXXFLAGS      += $(ROOTCFLAGS) -I$(ROOFITSYS)/include/ -I$(CMSSW_BASE)/src

NGLIBS         = $(ROOTGLIBS) 
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit -lboost_regex
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./

OUTLIB           = ./lib/

$(OUTLIB)HZZ4lBase.o: $(INCLUDEDIR)/HZZ4lBase.C
	@mkdir -p lib/
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)HZZ4lBase.o $<

$(OUTLIB)HZZ4l.o: $(INCLUDEDIR)/HZZ4l.C $(OUTLIB)HZZ4lBase.o $(OUTLIB)XSecReader.o $(OUTLIB)PUReweight.o
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)HZZ4l.o $<

run_HZZ4l: $(INCLUDEDIR)run_HZZ4l.cpp $(OUTLIB)HZZ4l.o
	$(CXX) $(CXXFLAGS) -o run_HZZ4l $(OUTLIB)/*.o $(GLIBS) $<

run_HZZ4l_CR: $(INCLUDEDIR)run_HZZ4l_CR.cpp $(OUTLIB)HZZ4l.o
	$(CXX) $(CXXFLAGS) -o run_HZZ4l_CR $(OUTLIB)/*.o $(GLIBS) $<

$(OUTLIB)XSecReader.o: $(INCLUDEDIR)/../Plotter/root_lib/XSecReader.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)XSecReader.o $<

$(OUTLIB)PUReweight.o: $(INCLUDEDIR)/../../src/PUReweight.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)PUReweight.o $<

saveN:	$(INCLUDEDIR)saveNormValues.cpp $(OUTLIB)XSecReader.o
	$(CXX) $(CXXFLAGS) -o saveNormValues $(OUTLIB)/*.o $(GLIBS) $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f $(OUTLIB)*.so
	rm -f run_HZZ4l

all: run_HZZ4l run_HZZ4l_CR
