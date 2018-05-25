#!/bin/bash
rm result.txt

intra=0

dossier=data/Litt_Real/

T=48
n=60

for sym in 2 3 4 ; do
  for id in {1..20}; do
    for met in -1 -4 -7 -2 -3 ; do	
      ./mf $met $dossier $n $T 1 3 $sym 0 $intra $id
    done
    printf "\n" >> result.txt
  done
  printf "\n" >> result.txt    
done
done

printf "\n" >> result.txt
printf "\n" >> result.txt
