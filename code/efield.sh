#!/bin/bash

start_time=$(date +%s)

# Run the efield.py script with parameters
python3 ./code/efield.py $1 $2
$efield_temp = $1

# Get the size of the temporary output file
file_size=$(stat -c%s ".results/$efield_temp.csv")

# Maximum file size in bytes (100MB)
max_size=$((50 * 1024 * 1024))
if [ $file_size -le $max_size ]; then
    mv .results/$efield_temp.csv .results/efield.csv
else
    split -b $max_size .results/$efield_temp.csv .results/efield_p_
    for f in .results/efield_p_*; do
        mv -- "$f" "$f.csv"
    done
    
    rm .results/efield_temp.csv
fi


end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Execution time: $execution_time seconds"