CFLAGS=-O3 -std=c++0x -larpack 

all: scimap 

scimap: main.cpp
	g++ -o main main.cpp Configuration.cpp PDBLoader.cpp LRMSD.cpp DPESScIMAP.cpp -I../../ann_1.1.2/include -L../../ann_1.1.2/lib -lANN $(CFLAGS) 

clean:
	rm *.o
