#!/bin/bash

# This function is for automatically transvert the inventory file name to ftp server path of IBDMDB HMP2 files
# use this after you cut the inventory.csv file use IBDMDB_csv_cut.sh
echo "Usage: ./this_bash.sh <input_file> <output_file>"
sleep 1

# Set here the input and ootput file path
input_file=$1
output_file=$2

# The Inputs are changing every different bath, should check or set as aurgument
# However, I am lazy
# Sh*t
url_prefix="https://ibdmdb.org/tunnel/static/HMP2/Metabolites/1723/HILIC-pos/"
url_suffix=".gz"

# Create the output file
> "$output_file"

# Read each line from the input file
while IFS= read -r line
do
  # Print the full URL and append it to the output file
  echo "${url_prefix}${line}${url_suffix}" >> "$output_file"
done < "$input_file" 
