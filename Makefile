ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --libs`
DEBUG = -Wall

all: anadrs.cpp
	c++ $(ROOTLIBS) $(ROOTFLAGS) $(DEBUG) -o anadrs anadrs.cpp

clean:
	rm -f anadrs *.o
