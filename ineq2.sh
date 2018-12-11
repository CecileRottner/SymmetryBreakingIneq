#!/bin/bash
rm result.txt

intra=0


dossier=data/smaller/

for y in 10 ; do





dossier=data/Litt_Real/



T=48
n=60

for T in 48 96 ; do
  for sym in 2 3 4 ; do
    for id in {1..20}; do
      for met in $(( -1 * $y )) $(( -3 * $y )) ; do	
        ./mf $met $dossier $n $T 1 3 $sym 0 $intra $id
      done
      printf "\n" >> result.txt
    done
    printf "\n" >> result.txt    
  done
done


printf "\n" >> result.txt
printf "\n" >> result.txt

done

