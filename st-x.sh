#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <dump-file> <nconf>"
    exit 1
fi

# Define the filename
FILENAME="$1"
totconf="$2"

nbin=20
nline=7609

# Extract xlo and xhi from the file
xlo=$(awk 'NR==6 {print $1}' "$FILENAME")
xhi=$(awk 'NR==6 {print $2}' "$FILENAME")

# Calculate dx
dx=$(awk -v xlo="$xlo" -v xhi="$xhi" -v nbin="$nbin" 'BEGIN {print (xhi - xlo) / nbin}')

# Output the calculated values
echo $xlo $xhi $dx

# Process each surface-tension file


for x in $(seq 1 $totconf); do
  k=$((nline * x))
  head -n $k "$FILENAME" | tail -n $nline | awk -v xlo="$xlo" -v dx="$dx" -v nbin="$nbin" 'NR>9 {
    for (i=0; i<nbin; i++) {
      if ($1 > xlo && $1 < (xlo + (dx * (i + 1)))) {
        tension[i]+= (0.5*($5+$6) - $4)
        #tension[i] += $2
        count[i]++
      }
    }
  }
  END {
    for (i=0; i<nbin; i++) {
      print (xlo + (dx * i)), tension[i]
    }
  }' > surface-tension-bin-$x.dat
  echo $x
done

rm -f *avg*

# Run the averaging script
sh ~mhassan01/analysis/shell_script/average.sh surface-tension-bin-*.dat

# Move the output to the final file
mv avg.dat surface-tension-bin-avg.dat


