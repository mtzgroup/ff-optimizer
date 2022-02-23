#!/bin/bash

for i in $(seq 1 30)
do
  rm -rf 1_opt/opt_"$i"*
  rm -rf 1_opt/valid_"$i"*
  rm -rf 2_sampling/"$i"_*
done
rm ObjectiveFunction.png
rm ParameterChange.png
rm -rf 1_opt/forcefield
rm -rf 1_opt/targets
rm -rf 1_opt/result
mv 1_opt/opt_0.in temp
rm -rf 1_opt/opt_0*
mv temp 1_opt/opt_0.in
