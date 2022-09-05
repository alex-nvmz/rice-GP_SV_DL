#!/bin/bash

# Go to folder where bash script is. Coherence with relative paths
cd $(dirname $0)

# 1. Parameters ---------------------------------------------------------------------------------

# Define Parameters file and Script
IN_PAR="./utility/pars_GWAS.csv"
IN_SCRIPT="./09b_GWAS.R"
# Suffix
NOHUP_SUF="GWAS"

# Depends on the machine
N_CORES=4

# 2. Script -------------------------------------------------------------------------------------

# Get number of parameter combinations
N_PAR=$(awk 'NR > 1' $IN_PAR | wc -l)

# For every param combination
for PAR in $(seq $N_PAR)
do
    nohup Rscript $IN_SCRIPT $PAR > "./nohup_logs/nohup_${NOHUP_SUF}_${PAR}.out" &

    # Execute up to N_CORES jobs in parallel
    # While number of jobs running (r) in the background (p, just the ID) is >= than N_CORES; sleep
    while [[ $(jobs -pr | wc -l) -ge $N_CORES ]]
    do
        # Sleep 1 second
        sleep 1
    done
done
