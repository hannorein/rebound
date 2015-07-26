export OPENGL=1

all: librebound
	@echo ""
	@echo "Compiling problem file..."
	@echo ""
	cc -I../../src/ -L. problem.c -lrebound -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."
	@echo ""

librebound:
	$(MAKE) -C ../../src/ librebound
	@cp ../../src/librebound.so .

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf rebound
