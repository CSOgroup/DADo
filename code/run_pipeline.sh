#!/usr/bin/bash

# running like:
#  ./run_pipeline.sh <run_file> <step#>

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

