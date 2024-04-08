#bin!bash!
# What you need to do with this bash is to copy the following codes and paste when you are in the folder of .gz files and run

mkdir rawfiles
for file in *.gz ; do gunzip -c "$file" > rawfiles/"${file%.*}" ; done
