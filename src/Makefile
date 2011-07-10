UNAME := $(shell uname)
CC=gcc
ifeq ($(UNAME), Linux)
	OPT = -L/usr/lib -I/usr/include -Wall  -g
	LIB =  -lm
endif
ifeq ($(UNAME), Darwin)
	OPT = -L/usr/local/lib -I/usr/local/include  -Wall -g 
	OPT+= -D_APPLE  
	LIB =  
endif

# Nothing needs to be changed below

ifeq ($(OPENGL), 1)
ifeq ($(UNAME), Linux)
	OPT+= -DOPENGL 
	LIB+= -lglut 
endif
ifeq ($(UNAME), Darwin)
	OPT+= -DOPENGL 
	LIB+= -framework OPENGL -framework GLUT
endif
endif

# Create precompiler definitions from links
OPT+= -D$(shell basename `readlink gravity.c` '.c' | tr '[a-z]' '[A-Z]')
OPT+= -D$(shell basename `readlink integrator.c` '.c' | tr '[a-z]' '[A-Z]')
OPT+= -D$(shell basename `readlink boundaries.c` '.c' | tr '[a-z]' '[A-Z]')
OPT+= -D$(shell basename `readlink collisions.c` '.c' | tr '[a-z]' '[A-Z]')

all:
	$(CC) $(OPT) -std=c99 main.c -c
	$(CC) $(OPT) -std=c99 problem.c -c
	$(CC) $(OPT) -std=c99 tree.c -c
	$(CC) $(OPT) -std=c99 particle.c -c
	$(CC) $(OPT) -std=c99 gravity.c -c
	$(CC) $(OPT) -std=c99 integrator.c -c
	$(CC) $(OPT) -std=c99 boundaries.c -c
	$(CC) $(OPT) -std=c99 output.c -c
	$(CC) $(OPT) -std=c99 collisions.c -c
	$(CC) $(OPT) -std=c99 zpr.c -c
	$(CC) $(OPT) -std=c99 opengl.c -c
	$(CC) $(LIB) -std=c99 *.o -o nbody

clean:
	rm -vf *.o
	rm -vf nbody
