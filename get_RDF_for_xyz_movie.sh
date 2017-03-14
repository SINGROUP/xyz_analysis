#!/bin/bash
#
# Compute g(r) from a given xyz movie file for a given element
# pair. If no element pair is given, calculate total g(r) over all
# atoms.
#
# Eero Holmstrom, 2015
#

# usage
if [[ -z "$5" ]]
then
    echo "Usage: $(basename $0) [xyz file] [dr] [rmin] [rmax] [skip this many steps from beginning] [e1] [e2]"
    exit 1;
fi

# assign input arguments
infile=$1
dr=$2
rmin=$3
rmax=$4
nskip=$5

# get number of atoms
natoms=$(awk '{print $1; exit}' $infile)

#
# if e1 and e2 are given by the user, then compute the g(r) for e1-e2
#
if [[ -n "$7" ]];
then
 
e1=$6
e2=$7

echo ""
echo "Computing g(r) for $e1-$e2..."

#
# calculate partial e1-e2 g(r)
#
awk -v nskip=$nskip -v e1=$e1 -v e2=$e2 -v natoms=$natoms -v dr=$dr -v rmin=$rmin -v rmax=$rmax 'BEGIN{n=0; nframesprocessed=0; nframestotal=0; pi=3.14159265359;};

# keep track of line number within frame
{n++;};

# get size of cell for later calculating total density
(n==2){xsize=$(NF-2); ysize=$(NF-1); zsize=$NF;};

# read in atomic coordinates for e1 and e2 for this frame
(n>2){

if ($1==e1) {i++; xe1[i]=$2; ye1[i]=$3; ze1[i]=$4; };
if ($1==e2) {j++; xe2[j]=$2; ye2[j]=$3; ze2[j]=$4; };

};

# compute the RDF once you reach the end of this frame
(n==natoms+2){

nframestotal++;

if(nframestotal > nskip){

ne1=i;
ne2=j;

# calculate relevant total density of particles
totalrho=ne2/(xsize*ysize*zsize);

#
# calculate g = g(r)
#

rcounter=0;

# loop over r
for(r=rmin+dr; r<=rmax; r=r+dr){

rcounter++;
rarray[rcounter] = r;
garray[rcounter] = 0;

# loop over all atom pairs
for (i=1; i<=ne1; i++){
for (j=1; j<=ne2; j++){

# calculate distance between atoms i and j, taking periodic boundary
# conditions into account
distx=xe1[i]-xe2[j];
disty=ye1[i]-ye2[j];
distz=ze1[i]-ze2[j];

if(distx < 0.0){distx=-1*distx;};
if(disty < 0.0){disty=-1*disty;};
if(distz < 0.0){distz=-1*distz;};

if(distx > xsize/2){distx=xsize-distx;}
if(disty > ysize/2){disty=ysize-disty;}
if(distz > zsize/2){distz=zsize-distz;}

thisr = sqrt( distx^2 + disty^2 + distz^2 );

if(thisr >= r-dr/2 && thisr < r+dr/2){garray[rcounter]++;}

};
};

# calculate density in this radial shell
Vr=4.0/3.0*pi*((r+dr/2)^3-(r-dr/2)^3);

# normalize density in this shell by total density and number of atoms used
garray[rcounter] = garray[rcounter]/Vr;
garray[rcounter] = garray[rcounter]/totalrho/ne1;

};

# update cumulative mean of RDF
for(i=1; i<=rcounter; i++){

totalrarray[i] = rarray[i];
totalgarray[i] = totalgarray[i] + garray[i];

}

# update relevant counters
nframesprocessed++;

};

# set relevant counters to zero for starting next frame
n=0; i=0; j=0;

};

# print out total cumulative mean of RDF
END{

for(i=1; i<=rcounter; i++){print totalrarray[i], totalgarray[i]/nframesprocessed};

}' $infile > this.rdf

mv this.rdf "RDF."$e1"_"$e2".dat"
echo "Done. Results output to" "RDF."$e1"_"$e2".dat"
echo ""

else

echo ""
echo "Computing total g(r) over all atoms..."

#
# calculate total g(r) over all atoms
#
awk -v nskip=$nskip -v e1=$e1 -v e2=$e2 -v natoms=$natoms -v dr=$dr -v rmin=$rmin -v rmax=$rmax 'BEGIN{n=0; i=0; nframesprocessed=0; nframestotal=0; pi=3.14159265359;};

# keep track of line number within frame
{n++;};

# get size of cell for later calculating total density
(n==2){xsize=$(NF-2); ysize=$(NF-1); zsize=$NF; totalrho=natoms/(xsize*ysize*zsize);};

# read in atomic coordinates for this frame
(n>2){ {i++; x[i]=$2; y[i]=$3; z[i]=$4;} };

# compute the RDF once you reach the end of this frame
(n==natoms+2){

nframestotal++;

if(nframestotal > nskip){

#
# calculate g = g(r)
#

rcounter=0;

# loop over r
for(r=rmin+dr; r<=rmax; r=r+dr){

rcounter++;
rarray[rcounter] = r;
garray[rcounter] = 0;

# loop over all atom pairs
for (i=1; i<=natoms; i++){
for (j=1; j<=natoms; j++){

if(i!=j){

# calculate distance between atoms i and j, taking periodic boundary
# conditions into account
distx=x[i]-x[j];
disty=y[i]-y[j];
distz=z[i]-z[j];

if(distx < 0.0){distx=-1*distx;};
if(disty < 0.0){disty=-1*disty;};
if(distz < 0.0){distz=-1*distz;};

if(distx > xsize/2){distx=xsize-distx;}
if(disty > ysize/2){disty=ysize-disty;}
if(distz > zsize/2){distz=zsize-distz;}

thisr = sqrt( distx^2 + disty^2 + distz^2 );

if(thisr >= r-dr/2 && thisr < r+dr/2){garray[rcounter]++;}

};

};
};

# calculate density in this radial shell
Vr=4.0/3.0*pi*((r+dr/2)^3-(r-dr/2)^3);

# normalize density in this shell by total density and number of atoms used
garray[rcounter] = garray[rcounter]/Vr;
garray[rcounter] = garray[rcounter]/totalrho/natoms;

};

# update cumulative mean of RDF
for(i=1; i<=rcounter; i++){

totalrarray[i] = rarray[i];
totalgarray[i] = totalgarray[i] + garray[i];

}

# update relevant counters
nframesprocessed++;

};

# set relevant counters to zero for starting next frame
n=0; i=0;

};

# print out total cumulative mean of RDF
END{

for(i=1; i<=rcounter; i++){print totalrarray[i], totalgarray[i]/nframesprocessed};

}' $infile > this.rdf

mv this.rdf "RDF.total.dat"
echo "Done. Results output to" "RDF.total.dat"
echo ""

fi

exit 0;

