#!/bin/bash
xccc -Wall matrix.c -o matrix;
for i in 120 480
do
	xcbatch 1 matrix $i 100;
done


for j in 1920 4096 6000
do
	xcbatch 1 matrix $j 10;
done

exit 0;
