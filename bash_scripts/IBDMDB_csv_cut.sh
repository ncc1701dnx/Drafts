#bin!bash!

# This is used for cut inventory.csv file of HMP2 on IBDMDB to only the filename that we need

inventory=$1
filenames=$2
cut -d ',' -f 3 "$inventory" | tail -n +2 >"$filenames"

