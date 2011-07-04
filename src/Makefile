OPT= -L/usr/local/lib -I/usr/local/include -fopenmp -Wall -g 
OPT += -D_APPLE 
LIB =  -framework OPENGL -framework GLUT -lpng


all:
	gcc $(OPT) -std=c99 main.c -c
	gcc $(LIB) -std=c99 *.o -o nbody

clean:
	rm -vf *.o
	rm -vf nbody
