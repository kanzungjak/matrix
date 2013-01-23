#!/bin/bash
# Zweimal laufen lassen!
xccc -Wall ../matrix.c -o matrix;

p=4 # 1,4,9,16,25
n=120 # 120,480,1920
# Wiederholungsmessungen z.B. 10, 20, 40, 80, 100 ->
# Bei uns besser 8, 16, 40, 80, 96
count=0

while [ $count -le 8 ]
	do
		count=$[$count+1]
		xcbatch $p matrix $n;
	done

exit 0
