
CC= gcc 
OPT= -O3 -DCOLLISIONS_NONE

all: 
	$(CC) $(OPT) -c -fPIC -Wall -Werror ../src/integrator_ias15.c
	$(CC) $(OPT) -c -fPIC -Wall -Werror ../src/gravity_direct.c
	$(CC) $(OPT) -c -fPIC -Wall -Werror ../src/tools.c
	$(CC) $(OPT) -c -fPIC -Wall -Werror ../src/particle.c
	$(CC) $(OPT) -c -fPIC -I../src/ -Wall -Werror libias15.c
	$(CC) $(OPT) -shared -o libias15.so particle.o tools.o integrator_ias15.o gravity_direct.o libias15.o

	@echo ""	
	@echo "The shared library libias15.so has been created."

