#!/bin/bash

FILENAME="../equil_traj.dump"
nbin=20

ylo=$(awk 'NR==6 {print $1}' "$FILENAME")
yhi=$(awk 'NR==6 {print $2}' "$FILENAME")
dy=$(awk "BEGIN {print ($yhi - $ylo) / $nbin}")

totconf=$(($(grep "TIMESTEP" "$FILENAME" | wc -l) / 2))

# Calculate the normalization factor
normalization_factor=$(awk "BEGIN {print $totconf * ($yhi - $ylo)**2 * $dy}")

# Determine the total number of data lines
total_data_lines=$(wc -l < "$FILENAME")

# Calculate the number of lines to skip the first half
lines_to_skip=$((total_data_lines / 2 + total_data_lines % 2))

# Calculate the number of data lines to process
lines_to_process=$((total_data_lines - lines_to_skip))

# Create the output file and add the header
echo -e "#\tbin\tall\tt1\tt2" > y-dist.dat

# Tail the file to get the last half of the data lines and append the output to y-dist.dat
tail -n "$lines_to_process" "$FILENAME" | awk -v ylo="$ylo" -v dy="$dy" -v nbin="$nbin" -v normalization_factor="$normalization_factor" '
    NR>9 {val=$5; type=$3; bin=int((val - ylo) / dy);
    if (bin >= 0 && bin < nbin) {n[bin][type]++}}
    END {for (i=0; i<nbin; i++) print i, (n[i][1]+0+n[i][2]+0) / normalization_factor, (n[i][1]+0) / normalization_factor, (n[i][2]+0) / normalization_factor}' >> y-dist.dat

