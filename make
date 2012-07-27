#!/bin/bash
optimize=0
debug=0
while (( "$#" ))
do
if [ "$1" == "debug" ]
then
debug=1
fi
if [ "$1" == "optimize" ]
then
    optimize=1
fi
shift
done
g++ -c lodepng.cpp
ar rvs lodepng.a lodepng.o
if [ $debug -eq 1 ]
then
g++ -Wall Nuclei21.cpp -g -o Nuclei21 lodepng.a
else
    if [ $optimize -eq 1 ]
    then
g++ -Wall Nuclei21.cpp -O2 -o Nuclei21 lodepng.a
else
g++ -Wall Nuclei21.cpp -o Nuclei21 lodepng.a
fi
fi
g++ pos_lapl_init.cpp -o pos_lapl_init
