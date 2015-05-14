
mercury: mercury6.tar compilemercury

mercury6.tar:
	wget http://star.arm.ac.uk/~jec/mercury/mercury6.tar

compilemercury:
	rm -rf mercury6
	mkdir mercury6
	tar xvf mercury6.tar -C mercury6
	cd mercury6 && gfortran -o mercury mercury6_2.for

swifter: swifter.tar.gz compileswifter

swifter.tar.gz:
	wget http://www.boulder.swri.edu/swifter/swifter.tar.gz

compileswifter:
	rm -rf swifter
	tar xvfz swifter.tar.gz
	patch swifter/tool/tool_follow.f90 swifter.patch
	cp swifter.Makefile.Defines swifter/Makefile.Defines
	cd swifter && make clean && make
