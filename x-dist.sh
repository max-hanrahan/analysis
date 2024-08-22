#!/bin/bash

FILENAME=$1
nbin=100

xlo=$(awk 'NR==6 {print $1}' "$FILENAME")
xhi=$(awk 'NR==6 {print $2}' "$FILENAME")
ylo=$(awk 'NR==7 {print $1}' "$FILENAME")
yhi=$(awk 'NR==7 {print $2}' "$FILENAME")
zlo=$(awk 'NR==8 {print $1}' "$FILENAME")
zhi=$(awk 'NR==8 {print $2}' "$FILENAME")

dx=$(awk "BEGIN {print ($xhi - $xlo) / $nbin}")

totconf=$(($(grep "TIMESTEP" "$FILENAME" | wc -l) / 2))

# Calculate the normalization factor
normalization_factor=$(awk "BEGIN {print $totconf * ($yhi - $ylo)*($zhi - $zlo) * $dx}")

# Determine the total number of data lines
total_data_lines=$(wc -l < "$FILENAME")

# Calculate the number of lines to skip the first half
lines_to_skip=$((total_data_lines / 2 + total_data_lines % 2))

# Calculate the number of data lines to process
lines_to_process=$((total_data_lines - lines_to_skip))

# Create the output file and add the header
echo -e "#\tbin\tall\tt1\tt2" > $FILENAME.dat

# Tail the file to get the last half of the data lines and append the output to x-dist.dat
tail -n "$lines_to_process" "$FILENAME" | awk -v xlo="$xlo" -v dx="$dx" -v nbin="$nbin" -v normalization_factor="$normalization_factor" '
    NR>9 {val=$4; type=$3; bin=int((val - xlo) / dx);
    if (bin >= 0 && bin < nbin) {n[bin][type]++}}
    END {for (i=0; i<nbin; i++) print i, (n[i][1]+0+n[i][2]+0) / normalization_factor, (n[i][1]+0) / normalization_factor, (n[i][2]+0) / normalization_factor}' >> $FILENAME.dat

