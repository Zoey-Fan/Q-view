
HEADERS =-I./ -I/usr/local/include -I/usr/include
LIBS =-L./ -L/usr/local/lib


# Viewshed3d:  Q-curve.cpp
# 	mpicxx $(HEADERS) $(LIBS) -std=c++17 -fopenmp -Wall -O3 -g -o Q-CURVE Q-curve.cpp -lproj -lpng 
# clean:
# 	rm  Q-CURVE


Viewshed3d:  Q-view.cpp
	mpicxx $(HEADERS) $(LIBS) -std=c++17 -fopenmp -Wall -O3 -g -o Q-VIEW Q-view.cpp -lproj -lpng 
clean:
	rm  Q-VIEW
