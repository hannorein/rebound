UNAME=Darwin
ifeq ($(UNAME), Linux)
OPT= -L/home/rein/localinstalls_linux/lib -I/home/rein/localinstalls_linux/include  -O3 -Wall  -g 
LIB=  -lm -lglut  -lconfig
endif
ifeq ($(UNAME), Darwin)
OPT= -L/usr/local/lib -I/usr/local/include  -Wall -g 
OPT += -D_APPLE  
OPT += -DOPENGL 
LIB =  -framework OPENGL -framework GLUT
endif

all:
	gcc $(OPT) -std=c99 main.c -c
	gcc $(OPT) -std=c99 particle.c -c
	gcc $(OPT) -std=c99 gravity.c -c
	gcc $(OPT) -std=c99 integrator.c -c
	gcc $(OPT) -std=c99 zpr.c -c
	gcc $(OPT) -std=c99 opengl.c -c
	gcc $(LIB) -std=c99 *.o -o nbody

clean:
	rm -vf *.o
	rm -vf nbody
