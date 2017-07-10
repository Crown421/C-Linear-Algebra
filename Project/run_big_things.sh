#!/bin/bash
#shopt -s nullglob

for file in ../Data/randomOrthClust*.txt;
do
  echo $(basename "$file" .txt)
  # ./gmres1 $(basename "$file" .txt)
  # ./gmres1SAFE $(basename "$file" .txt) SAFE
  # ./gmres1O3 $(basename "$file" .txt) O3
  #
   ./gmres2 $(basename "$file" .txt)
  #./gmres2SAFE $(basename "$file" .txt) SAFE
  # ./gmres2O3 $(basename "$file" .txt) O3

  # ./gaussj $(basename "$file" .txt)
  # ./gaussjSAFE $(basename "$file" .txt) SAFE
  # ./gaussjO3 $(basename "$file" .txt) O3
  #
  # ./LUsolve $(basename "$file" .txt)
  # ./LUsolveSAFE $(basename "$file" .txt) SAFE
  # ./LUsolveO3 $(basename "$file" .txt) O3

done
