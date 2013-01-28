#!/bin/bash
xccc -Wall matrix.c -o matrix;


for i in 120 480 1920 4096 6000 9000 12000 15000
do
	xcbatch 36 matrix $i 10;
done
