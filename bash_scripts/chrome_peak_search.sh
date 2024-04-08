#!/bin/bash
# This is the function for manually check the existence the file name of some Chromatic peaks XCMS cannot fetch

# The arguements are /the/path/to/specific/mzXML/data/ basePeakMz retentionTime 
dir="$1"
basePeakMz="$2"
retentionTime="$3"

# Create the output .txt file
output_file="matches.txt"

# Check if directory is right
if [[ -z "$dir" ]]
then
  echo "Directory argument is missing. Usage: /the/path/to/specific/mzXML/data/"
  exit 1
fi

# Check if the M/Z and RT is given
if [[ -z "$basePeakMz" ]] || [[ -z "$retentionTime" ]]
then
  echo "One or more arguments are missing. Usage: basePeakMz retentionTime after your directory"
  exit 1
fi

# Loop through all .mzXML files in the given directory
for file in "$dir"/*.mzXML
do
  echo "Processing $file"

  # grep the lines that contain the specified 'basePeakMz' and 'retentionTime'
  matches=$(grep -P "basePeakMz=\"$basePeakMz\d*\".*retentionTime=\"$retentionTime\d*S\"" "$file")

  if [[ -z "$matches" ]]
  then
    echo "$file: no match" >> "$output_file"
  else
    echo "$file: $matches" >> "$output_file"
  fi
done

echo "Processing complete. Results are in $output_file."
