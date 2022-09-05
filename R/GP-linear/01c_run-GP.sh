#!/bin/bash

# Go to folder where bash script is. Coherence with relative paths
cd $(dirname $0)

# Help function
help()
{
    echo "Usage: $(basename $0) [ -p | --par ]
              [ -s | --script ]
              [ -f | --suf  ]
              [ -c | --cores ]"
    exit 2
}

# 1. Parse arguments ---------------------------------------------------------------------

# IN_PAR (-p | --par): path to parameters file
# IN_SCRIPT (-s | --script): path to script to run
# NOHUP_SUF (-f | --suf): suffix for the nohup logs
# N_CORES (-c | --cores): number of cores to use

SHORT=p:,s:,f:,c:,h
LONG=par:,script:,suf:,cores:,help
OPTS=$(getopt --alternative --name run.sh --options $SHORT --longoptions $LONG -- "$@")

# If no arguments, help
VALID_ARGUMENTS=$#
if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  help
fi

eval set -- "$OPTS"

while [ : ] ; do
    case "$1" in
        -p | --par)
            IN_PAR="$2"
            shift 2
            ;;
        -s | --script)
            IN_SCRIPT="$2"
            shift 2
            ;;
        -f | --suf)
            NOHUP_SUF="$2"
            shift 2
            ;;
        -c | --cores)
            N_CORES="$2"
            shift 2
            ;;
        -h | --help)
            help
            ;;
        --)
            shift;
            break
            ;;
        *)
            echo "Unexpected option $1"
            ;;
    esac
done

# 2. Run parallel jobs -------------------------------------------------------------------

# Get number of parameter combinations
N_PAR=$(awk 'NR > 1' $IN_PAR | wc -l)

# For every param combination
for PAR_ROW in $(seq $N_PAR)
do
    # Rscript is given Par_file and Par_file_row
    nohup Rscript $IN_SCRIPT $(basename $IN_PAR) $PAR_ROW > "./nohup_logs/nohup_${NOHUP_SUF}_${PAR_ROW}.out" &

    # Execute up to N_CORES jobs in parallel
    # While number of jobs running (r) in the background (p, just the ID) is >= than N_CORES; sleep
    while [[ $(jobs -pr | wc -l) -ge $N_CORES ]]
    do
        # Sleep 1 second
        sleep 1
    done
done
