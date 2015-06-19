#!/bin/bash
set -o pipefail
echo "Running tests...\n"
for dir in ./*/
do
	echo ${dir}
	pushd ${dir} &> /dev/null
	make > /dev/null
	if [ $? -eq 1 ]; then
		echo "Error. Did not compile."
	fi
	popd  &> /dev/null
done
