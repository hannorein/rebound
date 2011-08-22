#!/bin/bash
export OS=`uname`
for example in `ls -d */`
	do echo $example
	cd $example
	make && ./nbody
	cd ..
done
