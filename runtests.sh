#!/bin/bash

su2params="-s 8 -t 4 -b 5.00 -n 20 -k 5 -e 40 -c -m"
su3params="-s 8 -t 4 -b 6.00 -n 20 -k 5 -e 40 -c -m"

function cleanup {
	cd ..
  rm -r build_test
}

mkdir build_test
cd build_test
cmake ../
make -j 2
if [ -f "./puresu2gauge" -a -f "./puresu3gauge" ]
then
	echo "Compile successful !"
else
	echo "Compile NOT successful !"
  cleanup
	exit
fi
echo "Performing test run..."
echo "Calling './puresu2gauge ${su2params}' ..."
./puresu2gauge $su3params
echo "Calling './puresu3gauge ${su3params}' ..."
./puresu3gauge $su2params 
echo "Cleaning up..."
cleanup
cd ..
