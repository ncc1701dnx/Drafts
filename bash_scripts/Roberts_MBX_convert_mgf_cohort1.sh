#!bin/bash

# Define the input directory
INPUT_DIR="/nfs/data/metabolomics_chenrui/cohort1/massive.ucsd.edu/MSV000082094/updates/2019-01-10_rhmills_712042dc/raw/"

# Define the output directory
OUTPUT_DIR="/nfs/data/Roberts_MBX_results/MGFs_cohort1"

# Create the output directory if it doesn't exist
sudo mkdir -p $OUTPUT_DIR

# Loop through each .mzXML file in the input directory and convert it
for file in $INPUT_DIR*.mzXML; do
    # Extract the filename without the path and extension
    filename=$(basename -- "$file")
    base="${filename%.*}"

    # Convert the .mzXML file to .mgf format using msconvert (via Docker)
    sudo docker run -it --rm -e WINEDEBUG=-all -v $INPUT_DIR:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/$filename --mgf --filter "msLevel 2" -o /data/$OUTPUT_DIR/
done
