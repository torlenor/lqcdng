#!/bin/bash
function cleanup {
	cd ..
  rm -r build_test
}

mkdir build_test
cd build_test
cmake ../
make
if [ -f "./puresu2gauge" -a -f "./puresu3gauge" ]
then
	echo "Compile successful !"
else
	echo "Compile NOT successful !"
  cleanup
	exit
fi
echo "Performing test run..."
echo "Calling './puresu2gauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 20 -c -m' ..."
./puresu3gauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 40 -c -m
echo "Calling './puresu3gauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 20 -c -m' ..."
./puresu2gauge -s 8 -t 4 -b 5.00 -n 10 -k 5 -e 40 -c -m
echo "Cleaning up..."
cleanup
cd ..
