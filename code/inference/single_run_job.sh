SIM_NAME=$1
SIM_ID=$2
USE_FEATURES=$3
USE_RJ=$4
MAX_SUBRANGE_SPLIT_SIZE=$5
N_PROC=$6

echo "input_prefix=\"${SIM_NAME}\"; sim_id=${SIM_ID}; use_features=${USE_DIST}; use_add=false; use_rj=${USE_RJ}; max_subrange_split_size=${MAX_SUBRANGE_SPLIT_SIZE}; n_proc=${N_PROC}; source(\"run_FIG.Rev\")" | rb-mrm
