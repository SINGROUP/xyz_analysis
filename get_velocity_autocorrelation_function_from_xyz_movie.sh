#!/bin/bash
#
# get_velocity_autocorrelation_function_from_xyz_movie.sh
#
# Using finite displacements, compute the velocity autocorrelation
# function (VAF) of a given xyz trajectory file. NB! Assumes unfolded
# coordinates for the atoms. Adopting convention of Dickey (1969) for
# the VAF.
#
# Eero Holmstrom, 2015
#

# usage
if [[ -z "$6" ]]
then
    echo "Usage: $(basename $0) [trajectory.xyz] [time step in simulation in fs] [first time origin in time steps] [step between origins] [last time origin in time steps] [maximum time steps per origin]"
    exit 1;
fi

# assign input variables
moviefile=$1
timestep=$2
firstt0=$3
dt0=$4
lastt0=$5
tmax=$6

# clean up previous runs
rm -f VAF.* total.VAF.dat

# get number of ions
nions=$(awk '{print $1; exit}' $moviefile)

#
# loop over time origins, producing a separate file with mean(sum v(t)*v(0)) / mean(sum v(t)^2) vs. time for each origin, where the mean is taken over atoms.
#
i=0;
for ((t0=firstt0; t0<=lastt0; t0=t0+dt0))
do

# compute VAF for this time origin
i=$(($i+1));
awk -v nions=$nions -v timestep=$timestep -v t0=$t0 -v tmax=$tmax 'BEGIN{frame=0; globalline=0; localline=-999; dotprodcumlower=0.0; dotprodcumupper=0.0;}

# keep track of line number
{globalline++;};

# when we quit the header, we are in the first frame
(globalline==2){frame=1; localline=0}

# whenever we are in a frame, keep track of the line number within the frame
(globalline>=3){localline++;}

# for the first frame, just save the particle positions
(frame==t0 && localline >= 1 && localline<=nions){x0[localline]=$2; y0[localline]=$3; z0[localline]=$4; xtprev[localline]=$2; ytprev[localline]=$3; ztprev[localline]=$4;};

# for the second frame, calculate v(0) for each particle and save the result.
(frame==t0+1 && localline >= 1 && localline<=nions){

x1[localline]=$2; y1[localline]=$3; z1[localline]=$4;

# calculate v0
vtx0[localline]=(x1[localline]-x0[localline]) / timestep;
vty0[localline]=(y1[localline]-y0[localline]) / timestep;
vtz0[localline]=(z1[localline]-z0[localline]) / timestep;

# calculate dot product of v0 and v0 for each ion
vt0vt0[localline] = vtx0[localline]^2 + vty0[localline]^2 + vtz0[localline]^2;

# compute the denominator of Dickey (1969), eq. 13
dotprodcumlower = dotprodcumlower + vt0vt0[localline];

};

# from the second frame onwards, calculate Cvv(t) = sum_1_to_nions [ dotprod(v_i(t), v_i(0)) ] / sum_1_to_nions [ dotprod(v_i(0), v_i(0)) ] for each step
(frame>=t0+1 && localline >= 1 && localline<=nions){

# get positions for this frame
xt[localline]=$2; yt[localline]=$3; zt[localline]=$4;

# calculate velocities for this frame from finite displacements
vtx[localline]=(xt[localline]-xtprev[localline]) / timestep;
vty[localline]=(yt[localline]-ytprev[localline]) / timestep;
vtz[localline]=(zt[localline]-ztprev[localline]) / timestep;

# for each ion, calculate dot product of current velocity vector with initial velocity vector. update cumulative dot product ratio.
dotprodcumupper = dotprodcumupper + (vtx[localline]*vtx0[localline] + vty[localline]*vty0[localline] + vtz[localline]*vtz0[localline]);

# save positions for calculating velocities for next frame
xtprev[localline]=$2; ytprev[localline]=$3; ztprev[localline]=$4;

};

# update counters, print out results.
(localline==nions+1 && frame >= t0+1){print (frame-(t0+1))*timestep, dotprodcumupper / dotprodcumlower; dotprodcumupper=0.0;}
(localline==nions+1){frame++; localline=-1}

# exit as designated
(frame >= t0+tmax){exit;}' $moviefile > VAF.$i

done

#
# loop over VAFs, creating an ensemble-averaged total VAF
#

# clean up from previous runs
rm -f total.VAF.dat

# get number of lines in an individual RDF file
for f in VAF.*
do
    nt=$(wc $f | awk '{print $1; exit}')
    break
done

tcounter=1
while [ $tcounter -le $nt ]
do

#
# print this value of t to the new row to be created
#
thist=$(awk -v tcounter=$tcounter 'BEGIN{n=0}{n++}(n==tcounter){print $1; exit}' $f)
echo $thist >> total.VAF.dat

# loop over VAF.* files
for g in VAF.*
do

# get the VAF(t) of this data set at this t
thisVAF=$(awk -v tcounter=$tcounter 'BEGIN{n=0}{n++}(n==tcounter){print $NF}' $g)

# add it to the current row in total.VAF.dat
awk -v tcounter=$tcounter -v thisVAF=$thisVAF 'BEGIN{n=0}{n++; if (n==tcounter) {print $0, thisVAF;} else {print $0}}' total.VAF.dat > vaf.dat
mv vaf.dat total.VAF.dat

done

tcounter=$(($tcounter+1))
done

#
# finally, calculate the average of the fields 2...NF for each row in total.VAF.dat
#
awk '{r=$1; meanvaf=0; for (i=2; i<=NF; i++){meanvaf=meanvaf+$i}; meanvaf=meanvaf / (NF-1); print r, meanvaf }' total.VAF.dat > mean.dat
mv mean.dat total.VAF.dat

exit 0;
