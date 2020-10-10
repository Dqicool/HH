PWD = $(shell pwd)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = -L/usr/lib/root -lGenVector -lEve -lEG -lGeom -lGed -lRGL -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic
ROOTINC   = -I/mnt/SSD/VBS-PHE/src

CXX = g++

convert: src/convert.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/convert.cpp -o build/convert $(ROOTLIBS) $(ROOTINC)

selection: src/selection.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/selection.cpp -o build/selection $(ROOTLIBS) $(ROOTINC)

draw: src/draw.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/draw.cpp -o build/draw $(ROOTLIBS) $(ROOTINC)