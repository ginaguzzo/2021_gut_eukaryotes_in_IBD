#!/bin/bash
# Calculate average in text file with script:

count=0;
total=0; 

for i in $( awk '{ print $1; }' average_read_lengths.txt )
   do 
     total=$(echo $total+$i | bc )
     ((count++))
   done
echo "scale=2; $total / $count" | bc
