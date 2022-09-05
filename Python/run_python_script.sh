#!/bin/bash

# Go to folder where bash script is. Coherence with relative paths
cd $(dirname $0)

# Help function
help()
{
    echo "Usage: $(basename $0) [ -p | --par ]
              [ -s | --script ]
              [ -f | --suf  ]
              [ -i | --init ]"
    exit 2
}

# 1. Parse arguments ---------------------------------------------------------------------

# IN_PAR (-p | --par): path to parameters file
# IN_SCRIPT (-s | --script): path to script to run
# NOHUP_SUF (-f | --suf): suffix for the nohup logs
# INITIAL_ROW (-i | --init): integer indicating first row of the parameters file

SHORT=p:,s:,f:,i:,h
LONG=par:,script:,suf:,init:,help
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
        -i | --init)
            INITIAL_ROW="$2"
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

# 2. Run job -------------------------------------------------------------------

# Get number of parameter combinations
N_PAR=$(awk 'NR > 1' $IN_PAR | wc -l)
# Python indexing starts at 0

# For every param combination (row numbers from initial row to number of par combs)
for PAR_ROW in $(seq $INITIAL_ROW $(($N_PAR - 1)))
do
    # python is given Par_file and Par_file_row
    nohup python $IN_SCRIPT $(basename $IN_PAR) $PAR_ROW > "./nohup_logs/nohup_${NOHUP_SUF}_${PAR_ROW}.out" &

    # Execute jobs 1 by 1
    # While number of jobs running (r) in the background (p, just the ID) is >= than 1; sleep
    while [[ $(jobs -pr | wc -l) -ge 1 ]]
    do
        # Sleep 1 second
        sleep 1
    done
done
