
CC= gcc 
OPT= -g -DPYTHON -DCOLLISIONS_NONE

all: 
	 $(CC) $(OPT) -c -fPIC -Wall -Werror ../src/integrator_ias15.c
	 $(CC) $(OPT) -c -fPIC -Wall -Werror ../src/gravity_direct.c
	 $(CC) $(OPT) -c -fPIC -I../src/ -Wall -Werror integrator_lib.c
	 $(CC) $(OPT) -shared -o libias15.so integrator_ias15.o gravity_direct.o integrator_lib.o


