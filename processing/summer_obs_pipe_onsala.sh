#!/bin/bash

# $1 - The directory of the data
# $2 - The name of the target directory
# $3 - The date of the observaton in YYYY-MM-DD
# $4 - The name of the target
usage="$(basename "$0") [-h] [-d n] [-f n] [-t n] [-n n] -- process the data on the Swedish LOFAR station into filterbank files

where:
    -h  show this help text
    -d  the directory of the data sess_sid[d]
    -f  the name of the target directory scan_[f]
    -t  the date of the observation
    -n  the name of the target"

while getopts ':hd:f:t:n:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    d) directory=$OPTARG ;;
    f) targetdir=$OPTARG ;;
    t) obsdate=$OPTARG   ;;
    n) targetname=$OPTARG ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

start=`date`

echo "Starting at" $start

# Make sure necessary directories exist
mkdir -p /datax2/brendan/$obsdate
mkdir -p /datax2/brendan/$obsdate/$targetname

# Move to data directory
cd /datax/Projects/proj21/sess_sid$directory/scan_$targetdir/SE607*/

# Copy the data that isn't .zst
echo "Copying"
cp  ./*.log /datax2/brendan/$obsdate/$targetname
cp  ./*.yml /datax2/brendan/$obsdate/$targetname
cp  ./*.h /datax2/brendan/$obsdate/$targetname

# Extract the .zst data to the output folder
echo "Extracting"

# Assign the filenames
file0=$(ls *.zst | grep "udp_SE607_16070")
file0=${file0::-4}
file1=$(ls *.zst | grep "udp_SE607_16071")
file1=${file1::-4}
file2=$(ls *.zst | grep "udp_SE607_16072")
file2=${file2::-4}
file3=$(ls *.zst | grep "udp_SE607_16073")
file3=${file3::-4}

# Extract them to the output
/home/brendan/zstd/build/cmake/programs/zstd -d ./$file0.zst -o /datax2/brendan/$obsdate/$targetname/$file0
/home/brendan/zstd/build/cmake/programs/zstd -d ./$file1.zst -o /datax2/brendan/$obsdate/$targetname/$file1
/home/brendan/zstd/build/cmake/programs/zstd -d ./$file2.zst -o /datax2/brendan/$obsdate/$targetname/$file2
/home/brendan/zstd/build/cmake/programs/zstd -d ./$file3.zst -o /datax2/brendan/$obsdate/$targetname/$file3

# Move to that directory - everything else will be done here so it's just easier
cd /datax2/brendan/$obsdate/$targetname

# Convert the data into GUPPI raw files
echo "Converting to GUPPI raw files"

# Find the filename
file0=${file0:16}

singularity exec /datax2/obs/singularity/lofar-upm_latest.simg bash -c "lofar_udp_extractor -p 30 -M GUPPI -I ./*.h -S 1 -b 0,411 -i ./udp_SE607_1607[[port]].$file0 -o ./$targetname.[[iter]].raw -m 4096"

# Remove the extracted files
echo "Deleting extracted files"
rm -r udp_SE607*

# Channelise the .raw files
echo "Running rawspec"

rawspec -f 65536,8 -t 2,16 $targetname
echo "$targetname.rawspec.0000.fil is the high spectral resolution product"
echo "$targetname.rawspec.0001.fil is the high time resolution product"

## Not happening until rawspec gets back working
# Remove .raw files to save space
echo "Deleting .raw files"
rm ./*.raw

echo " "
echo "Started at " $start
echo "Ended at " `date`
