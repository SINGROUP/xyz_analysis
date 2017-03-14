#!/bin/bash
#
# Script for multiplying an xyz structure over its periodic
# boundaries.
#
# Eero Holmstrom, 2013

# usage
if [[ -z $7 ]];
then
    echo "Usage: $(basename $0) [nx1 nx2 ny1 ny2 nz1 nz2] [input xyz file]";
    exit 1;
fi

# assign input parameters
nx1=$1
nx2=$2
ny1=$3
ny2=$4
nz1=$5
nz2=$6
infile=$7

# get box size
xsize=$(awk '{n++}(n==2){print $(NF-2); exit}' $infile)
ysize=$(awk '{n++}(n==2){print $(NF-1); exit}' $infile)
zsize=$(awk '{n++}(n==2){print $NF; exit}' $infile)

# strip the file of its header
rm -f $infile.coords
tail -$(($(wc $infile | awk '{print $1}')-2)) $infile > $infile.coords
rm -f final.coords

# loop over all dimensions, multiplying the cell as you go
for (( ix=$nx1; ix<=$nx2; ix++ ))
do
    for (( iy=$ny1; iy<=$ny2; iy++ ))
    do
	for (( iz=$nz1; iz<=$nz2; iz++ ))
	do

	    # do the multiplication for the block (ix, iy, iz)
	    awk -v xsize=$xsize -v ysize=$ysize -v zsize=$zsize -v ix=$ix -v iy=$iy -v iz=$iz '{print $1, $2+ix*xsize, $3+iy*ysize, $4+iz*zsize, $5;}' $infile.coords >> final.coords
	    
	done
    done
done

#
# turn final.coords into an xyz file
#

# create header first
totalnatoms=$(wc final.coords | awk '{print $1}')
totalxsize=$(awk -v xsize=$xsize -v nx1=$nx1 -v nx2=$nx2 'BEGIN{print (nx2-nx1+1)*xsize; exit}')
totalysize=$(awk -v ysize=$ysize -v ny1=$ny1 -v ny2=$ny2 'BEGIN{print (ny2-ny1+1)*ysize; exit}')
totalzsize=$(awk -v zsize=$zsize -v nz1=$nz1 -v nz2=$nz2 'BEGIN{print (nz2-nz1+1)*zsize; exit}')

echo $totalnatoms
echo "Frame number 1 1 fs boxsize $totalxsize $totalysize $totalzsize"

# then cat the position data
cat final.coords

# clean up
rm -f final.coords $infile.coords

exit 0;

