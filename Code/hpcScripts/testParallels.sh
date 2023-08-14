#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=1gb
#PBS -J 1-24

## Just a test to see if the parallel code works on HPC
## randomly chose to do 24 runs

## So I need 
seeds=()
for s in $(seq 100 100 2400)
do
    seeds+=("$s")
done

# echo "${args[@]}"

r0=(10)
K0=(10)
I0=(0.1)

args=()
s=0

for r in "${r0[@]}"
do
    for k in "${K0[@]}"
    do
        for int in "${I0[@]}"
        do
            for i in $(seq 1 24)
            do
                args+=("${seeds[s]}" "$r" "$k" "$int")
                s=$((s+1))
            done
        done
    done
done

cp $HOME/MTaNa/Code/MTaNaHPC.cpp $TMPDIR
g++ MTaNaHPC.cpp -fopenmp

./a.out ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 4)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 3)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 2)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 4 - 1)]}
