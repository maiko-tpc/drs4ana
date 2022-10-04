GCC = g++ -std=c++0x
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
DEBUG = -Wall
all: anadrs.cpp
	$(GCC) anadrs.cpp $(ROOTFLAGS) $(DEBUG) $(ROOTLIBS) -o anadrs
clean:                                                                          
	rm -f anadrs
