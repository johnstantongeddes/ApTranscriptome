#!/bin/bash

export PATH=/home/ubuntu/software/Sailfish-0.6.3-Linux_x86-64/bin:$PATH
export LD_LIBRARY_PATH=/home/ubuntu/software/Sailfish-0.6.3-Linux_x86-64/lib:$LD_LIBRARY_PATH

# script to bootstrap gene expression quantification using sailfish
parallel Rscript sailfish_quant.R {1} ::: {1..20}


