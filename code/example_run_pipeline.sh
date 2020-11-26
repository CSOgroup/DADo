#!/usr/bin/bash


#  ./example_run_pipeline.sh example1000_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAluad_norm_luad.R 1 2 3 4 5fc 5corr 6 7 8 9 10 11 12 13 # tab 10

script_name="$0"

start_time=$(date -R)    

scriptFolder="."  # where the script files are located

set -e

args=( "$@" )
args_len=${#args[@]}

settingF=${args[0]}

i=1
while [[ $i -lt args_len ]]; do
    scriptF="$(ls $scriptFolder | grep -P ^${args[$i]}_.+R$)"
    cmd="Rscript $scriptFolder/$scriptF $settingF"
    echo "> $cmd"
    $cmd
    ((i++))
done

echo "*** DONE - $script_name"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0



