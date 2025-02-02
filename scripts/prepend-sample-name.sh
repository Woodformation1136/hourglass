#!/bin/bash

PROJECTIONS=$1
SAMPLE_NAME=$2
OUTPUT=$3

# write head line
head -n 1 $PROJECTIONS > $OUTPUT

# append the rest
awk -v sample=$SAMPLE_NAME 'IF NR > 1 {print sample":"$0}' $PROJECTIONS >> $OUTPUT
