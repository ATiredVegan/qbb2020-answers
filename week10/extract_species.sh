#!/bin/bash
touch bin_identities.bar
for FILENAME in bins/*
do 
	line=$(head -n 1 $FILENAME);
	line="${line:1}"
	echo "$line"
	grep "$line" week13_data/KRAKEN/assembly.kraken >>bin_identities.txt
done 