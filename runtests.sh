#!/bin/bash
function cleanup {
	cd ..
  rm -r build_test
}

mkdir build_test
cd build_test
cmake ../
make
if [ -f "./lqcdng" ]
then
	echo "Compile successful !"
else
	echo "Compile NOT successful !"
	exit
fi
echo "Performing test run..."
echo "Calling './lqcdng -s 6 -t 4 -b 5.00 -n 10 -k 1 -e 10' ..."
./lqcdng -s 6 -t 4 -b 5.00 -n 10 -k 1 -e 10
echo "Cleaning up..."
cleanup
cd ..
