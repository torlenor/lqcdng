#!/bin/bash
function cleanup {
	cd ..
  rm -r build_test
}

mkdir build_test
cd build_test
cmake ../
make
if [ -f "./lqcdng_puregauge" ]
then
	echo "Compile successful !"
else
	echo "Compile NOT successful !"
  cleanup
	exit
fi
echo "Performing test run..."
echo "Calling './lqcdng_puregauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 10 -c -m' ..."
./lqcdng_puregauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 10 -c -m
echo "Cleaning up..."
cleanup
cd ..
