#!/bin/bash
cd sample-input/benchmarks/c5g7
for n in {1..4}
do
    for i in {1..20}
    do
        echo $i
        python segmentation_c5g7.py $n
    done
done
