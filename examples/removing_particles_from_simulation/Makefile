export OPENGL=0
export SERVER=0
include ../../src/Makefile.defs

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
all: problem.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""

clean:
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_console.html -o rebound_console.html
