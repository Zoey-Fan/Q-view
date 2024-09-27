HEADERS =-I./ -I/usr/local/include/ -I./glm -I./tinyCG
LIBS =-L./ -L/usr/local/lib -L./glm -L./tinyCG

Mergeindex:  mergeindex.cpp
	g++ $(HEADERS) $(LIBS) -std=c++17 -fopenmp -Wall -O3 -o Mergeindex mergeindex.cpp  ./tinyCG/threeCGAL.cpp
clean:
	rm  Mergeindex
