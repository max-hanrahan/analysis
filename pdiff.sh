#!/bin/sh 
awk 'BEGIN {maxbin=-1; dx=0.1} NR>8004500 && NR%16009>9 {delta=-$4+0.5*($5+$6); bin=int($1/dx); if (bin>maxbin) maxbin=bin; deltap[bin]+=delta; count[bin]++} END {for (n=0; n<=maxbin; n++) if (count[n]>0) print dx*n, deltap[n]/count[n]}' stress.dump > pdiff.dat
