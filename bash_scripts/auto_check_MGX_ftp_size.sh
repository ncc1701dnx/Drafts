#!/bin/bash
input_file=$1
# This is used for check the file existence and size of each downloadable file in MGX_ftp.txt url file
# Make sure output.txt is empty before we start
if [ -f "ftp_size.txt" ]; then
	echo "You already have MGX_ftp_size.txt file in this folder, please delete!"
	sleep 2
else
	if [ -f "$input_file" ]; then
		> ftp_size.txt
		while read -r url; do # Read the ftp url file line by line
			curl_info=$(curl --head "$url")
			size=$(echo "$curl_info" | grep 'Content-Length:' | awk '{print $2}') # Extract the size from curl --head information

			if [ ! -z "$size" ]; then # write the size line after the url line
				echo "$url - $size" >> ftp_size.txt
			else
				echo "$url - size not found" >> ftp_size.txt
			fi
		done < "$input_file" # Input site of the input ftp url file (MGX_ftp.txt)
	else
		echo "No list file for wget found, please make sure you put it into the same folder with this .sh file"
		sleep 2
	fi
fi
