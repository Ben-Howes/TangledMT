#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-405

## How many different runs will I need to do?
## 3 different r0 (intrinsic growth) values 0.1 1 10
## 3 different K0 (varrying capacity) values 1 10 100
## 3 different I0 (intraspecific competition) values 0.1 1 10
## 3 P0 (pOff slope) values 0.1 1 10

## So that is 3*3*3*3 = 81
## But I want to run each of these for 5 communities
## So 81*5 = 405

## So I need 405
seeds=()
for s in $(seq 100 100 40500)
do
    seeds+=("$s")
done

# echo "${args[@]}"

r0=(0.1 1 10)
K0=(1 10 100)
I0=(0.1 1 10)
P0=(0.1 1 10)

args=()
s=0

for r in "${r0[@]}"
do
    for k in "${K0[@]}"
    do
        for int in "${I0[@]}"
        do
            for p in "${P0[@]}"
            do
                for i in $(seq 1 5)
                do
                    args+=("${seeds[s]}" "$r" "$k" "$int" "$p")
                    s=$((s+1))
                done
            done
        done
    done
done

cp $HOME/MTaNa/Code/MTaNaHPC.cpp $TMPDIR
g++ MTaNaHPC.cpp

./a.out ${args[$( expr "$PBS_ARRAY_INDEX" '*' 5 - 5)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 5 - 4)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 5 - 3)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 5 - 2)]} ${args[$( expr "$PBS_ARRAY_INDEX" '*' 5 - 1)]}
