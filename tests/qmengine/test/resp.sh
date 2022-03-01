#!/bin/bash

ml TeraChem/2021.11-intel-2017.8.262-CUDA-11.4.1
for i in $(ls *pdb)
do
  name=$(echo $i | awk -F. '{print $1}')
  sed -i "s/.*coordinates.*/coordinates $i/g" tc_resp.in
  terachem tc_resp.in > tc_resp_"$name".out
  mv scr."$name"/esp.xyz esp_"$name".xyz
done
