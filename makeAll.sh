#!/bin/bash


declare -a array=(rates generateLines cellfinder mklos mkspec anal mkALLabs)

base="./funcs/"
for i in ${array[@]}; do
    loc=$base$i
    echo 
    echo $i
    cd $loc 
    if ls *.f 1> /dev/null 2>&1; then
        touch *.f
    else
        touch *.c
        make clean 
    fi
    make
    cd ../..
done




