#!/bin/bash

for i in $(ls *.pdb)
do
  name=$(echo $i | awk -F. '{print $1}')
  natoms=$(grep -c "AT.*M" $i)
  echo $natoms > "$name".xyz
  echo "converted from $i" >> "$name".xyz
  grep "AT.*M" $i | awk '{print $11 "\t" $6 "\t" $7 "\t" $8}' >> "$name".xyz
done
