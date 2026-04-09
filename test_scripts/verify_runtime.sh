#!/bin/bash
# expect the first argument to be the column name : contain what platform(OS) is running the script
expected_column="$1"

WD="test_scripts/test_data"
input_file="${WD}/expected_runtime.tsv"
reported_file="time_log.tsv"
selected_columns_file="${WD}/selected_columns.tsv"

# Get the column index (1-based) of the expected column name
col_index=$(head -1 "$input_file" | tr '\t' '\n' | awk -v col="$expected_column" '{if ($0 == col) print NR}')

# Check if column was found
if [[ -z "$col_index" ]]; then
  echo "Column '$expected_column' not found in $input_file"
  exit 1
fi

cut -f1,2,"$col_index" "$input_file" > "${selected_columns_file}"
final_file="${WD}/combined_with_reported.tsv"

# assuming the reported file and expected file have the same order of commands
cut -f2 "$reported_file" > /tmp/reported_column.tsv
paste "$selected_columns_file" /tmp/reported_column.tsv > "$final_file"


fail_count=0

# Skip header
while IFS=$'\t' read -r command threshold expected reported; do
    allowed=$(echo "$expected + $threshold" | bc -l)
    is_exceed=$(echo "$reported > $allowed" | bc -l)

    if [ "$is_exceed" = "1" ]; then
        diff=$(echo "$reported - $expected" | bc -l)
        echo "❌ $command exceeded the allowed runtime usage."
        echo "Expected: $expected S, Threshold: $threshold S, Reported: $reported S, Difference: $diff S"
        ((fail_count++))
    else
        echo "✅ $command passed the runtime check."
        diff=$(echo "$reported - $expected" | bc -l)
        echo "Expected: $expected S, Threshold: $threshold S, Reported: $reported S, Difference: $diff S"
    fi
done < <(tail -n +2 "$final_file")

if [ "$fail_count" -eq 0 ]; then
    echo "✅ All runtime checks passed."
else
    echo "❌ $fail_count checks failed."
    exit 1
fi
