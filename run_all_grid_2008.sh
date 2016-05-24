#!/bin/sh

need sgegrid

NUM_RUNS=50

for i in {1..8}; do
  qsub -t 1-$NUM_RUNS:1 variable_length_ga.sh ~/workspace/wsc2008/Set0${i}MetaData 2008-variable-length-ga${i} variable-length-ga.params;
done
