#!/bin/bash

#First read in the config file and store the values in a dictionary

config=$1

declare -A config_vals

while IFS=$'\t' read -r -a myArray
do
	config_vals["${myArray[0]}"]="${myArray[1]}"
done < $config

#run mosaicism script

./mosaicism $config

# run python script to generate images, arguments are input file and output_directory

python mosaicism_image_generator_script.py -i "${config_vals[OUTPUT_FILE_PATH]}" -o "${config_vals[IMAGE_DIRECTORY]}"
