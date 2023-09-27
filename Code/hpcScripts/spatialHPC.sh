#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=1gb:ompthreads=8
#PBS -J 1-100

## So I need 
seeds=()
for s in $(seq 100 100 10000)
do
    seeds+=("$s")
done

# echo "${args[@]}"

r0=(100)
K0=(1000)
I0=(0.5)

args=()
s=0

for r in "${r0[@]}"
do
    for k in "${K0[@]}"
    do
        for int in "${I0[@]}"
        do
            for s in "${seeds[@]}"
            do
                args+=("$s" "$r" "$k" "$int")
            done
        done
    done
done

cp $HOME/MTaNa/Code/spatialMTaNaHPC.cpp $TMPDIR
g++ spatialMTaNaHPC.cpp -fopenmp

./a.out ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 4)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 3)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 2)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 1)]}
