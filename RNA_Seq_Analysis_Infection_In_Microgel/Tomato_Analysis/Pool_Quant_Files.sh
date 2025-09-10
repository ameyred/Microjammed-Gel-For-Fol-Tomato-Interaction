#!/bin/bash

# Source directory where sample folders are
source_dir="/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/All_Quant_Files"

# Destination directory to collect renamed .sf files
dest=/Users/vidhasrivastava/Amey_Lab/Microgel_Project/Big_Experiment/Tomato_Analysis/All_Quant_SF_Files
mkdir -p "$dest"

# Loop through each subdirectory
for dir in "$source_dir"/*/; do
  if [[ -f "${dir}quant.sf" ]]; then
    folder_name=$(basename "$dir")
    sample_name="$folder_name"  # No renaming needed

    echo "Copying ${folder_name}/quant.sf as ${sample_name}.sf"
    cp "${dir}quant.sf" "${dest}/${sample_name}.sf"
  else
    echo "No quant.sf in $dir"
  fi
done
