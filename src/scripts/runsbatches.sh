#!/bin/bash

# runs all batch files in the current directory, use with caution

for i in *batch*; do sbatch $i; done
