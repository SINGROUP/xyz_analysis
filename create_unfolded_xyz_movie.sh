#!/bin/bash
#
# Take a periodic xyz movie file and turn it into an unfolded xyz
# movie file.
#
# Eero Holmstrom, 2014
#

# usage
if [[ -z "$1" ]];
then
    echo "Usage: $(basename $0) [movie.xyz]"
    exit 1;
fi

# assign input parameters
infile=$1

# get number of ions
nions=$(awk '{print $1; exit}' $infile)

echo ""
echo "Found $nions ions."
echo "Starting loop over movie file."
echo ""

# clean up previous runs
rm -f "unfolded_"$infile

# loop over frames, outputting the unfolded coordinates of each ion at each frame
awk -v nions=$nions 'BEGIN{ionline=nions; cframe=0;};

{ionline++};

# print out the number of atoms at header of each frame
(ionline==nions+1){print $1};

/Frame/{print $0; ionline=0; cframe++; header=$0; xsize=$(NF-2); ysize=$(NF-1); zsize=$NF};

(ionline >= 1 && ionline <= nions){

thisx=$2; thisy=$3; thisz=$4;

if(cframe>1) {
# now, if the displacement from the previous frame was more than half
# of the simulation cell size, you are sure to have witnessed the
# crossing of a periodic boundary. take this into account by keeping
# track of how many times each ion has crossed a boundary and in which
# direction.
dispprevx=thisx-xprev[ionline];
dispprevy=thisy-yprev[ionline];
dispprevz=thisz-zprev[ionline];

if(dispprevx > xsize/2.0){ xboundarycrossed[ionline]--;};
if(dispprevx < -1.0*xsize/2.0){ xboundarycrossed[ionline]++;};
if(dispprevy > ysize/2.0){ yboundarycrossed[ionline]--;};
if(dispprevy < -1.0*ysize/2.0){ yboundarycrossed[ionline]++;};
if(dispprevz > zsize/2.0){ zboundarycrossed[ionline]--;};
if(dispprevz < -1.0*zsize/2.0){ zboundarycrossed[ionline]++;};

# find unfolded position of this atom
thisx=thisx + xsize*xboundarycrossed[ionline];
thisy=thisy + ysize*yboundarycrossed[ionline];
thisz=thisz + zsize*zboundarycrossed[ionline];

};

# output the unfolded coordinates of this ion at this frame
print $1, thisx, thisy, thisz, $5;

# save current position of ion for comparison with next frame to come
xprev[ionline]=$2;
yprev[ionline]=$3;
zprev[ionline]=$4;

}' $infile > "unfolded_"$infile

echo "Done! Exiting.";
echo ""

exit 0;
