#!/bin/bash
#combine all maf files into one
#$1 - dir to maf files
#$2 - the output single maf file
for i in $1/*/*.gz; do
  # cp $i $2/$(basename $i);
  gzcat i >> $2
  echo $(basename $i);
done