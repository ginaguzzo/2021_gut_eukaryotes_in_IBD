#!/bin/bash

# Find average read length of each fastq file and save to text file:
for file in $(ls *gz) 
do 
    gzip -dc "$file" | \
    head -n 10000 | \
    awk '{ if (NR%4==2){count++; bases += length}} END{printf "%3.0f\n", bases/count}' \
    >> average_read_lengths.txt
done
