PWD = $(shell pwd)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = -L/usr/lib/root -lGenVector -lEve -lEG -lGeom -lGed -lRGL -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic
PWDINC   = -I/mnt/NVME/HH/src

CUDA_CC = /opt/cuda/bin/nvcc
CUDA_FLAG = -arch compute_75
CUDA_LIB = -L/opt/cuda/lib64 -lcudart
CUDA_INC = -I/opt/cuda/include

CXX = c++

convert: src/convert.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/convert.cpp -o build/convert $(ROOTLIBS) $(PWDINC)

presel: src/presel.cpp src/presel.h src/genAna.h
	$(CXX) $(ROOTFLAGS) -g  src/presel.cpp -o build/presel $(ROOTLIBS) $(PWDINC)

selection: build/selection.o build/GPUScaleB.o  
	$(CXX) $(ROOTFLAGS) -g build/selection.o build/GPUScaleB.o -o build/selection $(ROOTLIBS) $(CUDA_LIB) 
	rm -rf build/selection.o build/GPUScaleB.o

draw: src/draw.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/draw.cpp -o build/draw $(ROOTLIBS) $(PWDINC)

stack: src/stack.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/stack.cpp -o build/stack $(ROOTLIBS) $(PWDINC)

test: src/test.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/test.cpp -o build/test $(ROOTLIBS) $(PWDINC)

test2: src/test2.cpp src/genAna.h src/presel.h
	$(CXX) $(ROOTFLAGS) -g src/test2.cpp -o build/test2 $(ROOTLIBS) $(PWDINC)

build/selection.o: src/selection.h src/genAna.h
	$(CXX) $(ROOTFLAGS) -c -g src/selection.cpp -o build/selection.o $(PWDINC)

build/GPUScaleB.o: src/GPUScaleB.cu src/GPUScaleB.cuh
	$(CUDA_CC) -O3 -c src/GPUScaleB.cu -o build/GPUScaleB.o $(CUDA_FLAG) $(CUDA_INC) $(PWDINC)

fit: src/fit.cpp src/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/fit.cpp -o build/fit $(ROOTLIBS) $(PWDINC)

clean:
	rm -rf build/*