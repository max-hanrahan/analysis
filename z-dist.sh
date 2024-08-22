#!/bin/bash

FILENAME="../equil_traj.dump"
nbin=20

zlo=$(awk 'NR==6 {print $1}' "$FILENAME")
zhi=$(awk 'NR==6 {print $2}' "$FILENAME")
dz=$(awk "BEGIN {print ($zhi - $zlo) / $nbin}")

totconf=$(($(grep "TIMESTEP" "$FILENAME" | wc -l) / 2))

# Calculate the normalization factor
normalization_factor=$(awk "BEGIN {print $totconf * ($zhi - $zlo)**2 * $dz}")

# Determine the total number of data lines
total_data_lines=$(wc -l < "$FILENAME")

# Calculate the number of lines to skip the first half
lines_to_skip=$((total_data_lines / 2 + total_data_lines % 2))

# Calculate the number of data lines to process
lines_to_process=$((total_data_lines - lines_to_skip))

# Create the output file and add the header
echo -e "#\tbin\tall\tt1\tt2" > z-dist.dat

# Tail the file to get the last half of the data lines and append the output to z-dist.dat
tail -n "$lines_to_process" "$FILENAME" | awk -v zlo="$zlo" -v dz="$dz" -v nbin="$nbin" -v normalization_factor="$normalization_factor" '
    NR>9 {val=$6; type=$3; bin=int((val - zlo) / dz);
    if (bin >= 0 && bin < nbin) {n[bin][type]++}}
    END {for (i=0; i<nbin; i++) print i, (n[i][1]+0+n[i][2]+0) / normalization_factor, (n[i][1]+0) / normalization_factor, (n[i][2]+0) / normalization_factor}' >> z-dist.dat

