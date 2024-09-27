# HEADERS =-I./ -I/usr/local/include/ -I./glm -I./tinyCG 
# LIBS =-L./ -L/usr/local/lib -L/opt/GDAL/lib -L./glm -L./tinyCG 
HEADERS =-I./ -I/usr/local/include/ -I./glm -I./tinyCG 
LIBS =-L./ -L/usr/local/lib -L./glm -L./tinyCG 

# Build3dindex:  build3dindex.cpp
# 	g++ $(HEADERS) $(LIBS) -std=c++17 -fopenmp -o Build3dindex build3dindex.cpp  ./tinyCG/threeCGAL.cpp -lgdal -lgmp -Wunused-result
# clean:
# 	rm  Build3dindex
Build3dindex:  build3dindex.cpp
	g++ $(HEADERS) $(LIBS) -std=c++17 -fopenmp -o Build3dindex build3dindex.cpp  ./tinyCG/threeCGAL.cpp -lgmp -Wunused-result
clean:
	rm  Build3dindex
