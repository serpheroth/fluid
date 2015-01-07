all: clean water

water: water.o grid.o util.o sph.o
	g++ water.o grid.o util.o sph.o -framework OpenGL -framework Cocoa -framework IOkit -framework CoreVideo -lglfw3 -lGLEW -std=c++11 -O2 -o water

water.o: grid.h water.h sph.h water.cpp
	g++ water.cpp -Wall -std=c++11 -O2 -c 

grid.o: grid.h grid.cpp
	g++ grid.cpp -std=c++11 -O2 -c

util.o: util.h util.cpp
	g++ util.cpp -std=c++11 -O2 -c

sph.o: sph.h sph.cpp
	g++ sph.cpp -std=c++11 -O2 -c

clean:
	rm -rf water water.o  grid.o util.o  sph.o
