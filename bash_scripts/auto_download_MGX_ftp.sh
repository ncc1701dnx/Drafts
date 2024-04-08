#!/bin/bash
#make sure the MGX_ftp.txt file is in the same directory i.e. same path with this ssh file
echo "Make sure you have MGX_ftp.txt in the same folder"

sleep 2

if [ -f "MGX_ftp.txt" ]; then
    mkdir MGX_files
    wget -i MGX_ftp.txt -P MGX_files/
else
    echo "No list file for wget found, please make sure you put it into the same folder with this .sh file"
fi
