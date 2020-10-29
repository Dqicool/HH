PWD = $(shell pwd)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = -L/usr/lib/root -lGenVector -lEve -lEG -lGeom -lGed -lRGL -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic
ROOTINC   = -I/mnt/NVME/HH/src

CXX = g++

convert: src/convert.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/convert.cpp -o build/convert $(ROOTLIBS) $(ROOTINC)

presel: src/presel.cpp src/presel.h src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/presel.cpp -o build/presel $(ROOTLIBS) $(ROOTINC)

selection: src/selection.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/selection.cpp -o build/selection $(ROOTLIBS) $(ROOTINC)

draw: src/draw.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/draw.cpp -o build/draw $(ROOTLIBS) $(ROOTINC)

stack: src/stack.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/stack.cpp -o build/stack $(ROOTLIBS) $(ROOTINC)

test: src/test.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/test.cpp -o build/test $(ROOTLIBS) $(ROOTINC)

test2: src/test2.cpp src/genAna.h src/presel.h
	$(CXX) $(ROOTFLAGS) -g src/test2.cpp -o build/test2 $(ROOTLIBS) $(ROOTINC)