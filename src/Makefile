UNAME := $(shell uname)
CC=gcc
OPENGL=0
ifeq ($(UNAME), Linux)
	OPT = -L/usr/lib -I/usr/include -Wall  -g
	LIB =  -lm
ifeq ($(OPENGL), 1)
	OPT+= -DOPENGL 
	LIB+= -lglut 
endif
endif
ifeq ($(UNAME), Darwin)
	OPT = -L/usr/local/lib -I/usr/local/include  -Wall -g 
	OPT+= -D_APPLE  
	LIB =  
ifeq ($(OPENGL), 1)
	OPT+= -DOPENGL 
	LIB+= -framework OPENGL -framework GLUT
endif
endif

all:
	$(CC) $(OPT) -std=c99 main.c -c
	$(CC) $(OPT) -std=c99 problem.c -c
	$(CC) $(OPT) -std=c99 tree.c -c
	$(CC) $(OPT) -std=c99 particle.c -c
	$(CC) $(OPT) -std=c99 gravity.c -c
	$(CC) $(OPT) -std=c99 integrator.c -c
	$(CC) $(OPT) -std=c99 zpr.c -c
	$(CC) $(OPT) -std=c99 boundaries.c -c
	$(CC) $(OPT) -std=c99 opengl.c -c
	$(CC) $(LIB) -std=c99 *.o -o nbody

clean:
	rm -vf *.o
	rm -vf nbody
