#!/bin/bash

clear

rm -f ../database/*.dat ../plots/*.png

echo -e "\ncompiling... \n"

make

echo -e "running... "

time ./execut

out_scalar=$(sed -n '/#define out_scalar / p' params.h | cut -f3 -d' ')
sed -i '/out_scalar = / s/.*/out_scalar = '$out_scalar'/' scalar.plt

out_vector=$(sed -n '/#define out_vector / p' params.h | cut -f3 -d' ')
sed -i '/out_vector = / s/.*/out_vector = '$out_vector'/' vector.plt

out_tensor=$(sed -n '/#define out_tensor / p' params.h | cut -f3 -d' ')
sed -i '/out_tensor = / s/.*/out_tensor = '$out_tensor'/' tensor.plt


key=$(sed -n '/#define key_dynamics / p' params.h | cut -f3 -d' ')
sed -i '/key = / s/.*/key = '$key'/' scalar.plt
sed -i '/key = / s/.*/key = '$key'/' vector.plt
sed -i '/key = / s/.*/key = '$key'/' tensor.plt

Nf=$(sed -n '/#define Nf / p' params.h | cut -f3 -d' ')
sed -i '/Nf = / s/.*/Nf = '$Nf'/' scalar.plt
sed -i '/Nf = / s/.*/Nf = '$Nf'/' vector.plt
sed -i '/Nf = / s/.*/Nf = '$Nf'/' tensor.plt

gnuplot scalar.plt
gnuplot vector.plt
gnuplot tensor.plt

key_eog=$(sed -n '/#define key_eog / p' params.h | cut -f3 -d' ')

rm -f execut

echo -e "\nfinished !\n"

if (( $(echo "$key_eog == 1" | bc) )); then

	eog ../plots/*
fi

