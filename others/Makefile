
all: mercury6.tar compile

mercury6.tar:
	wget http://star.arm.ac.uk/~jec/mercury/mercury6.tar

compile:
	rm -rf mercury6
	mkdir mercury6
	tar xvf mercury6.tar -C mercury6
	cd mercury6 && gfortran -o mercury mercury6_2.for

