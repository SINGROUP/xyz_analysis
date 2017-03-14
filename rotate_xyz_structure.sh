#!/bin/bash
#
# Rotate a given xyz structure by the user-given Euler angles. If less
# than three angles are given, do a completely random rotation.
#
# Eero Holmstrom, 2013
#

# usage
if [[ -z "$1" ]];
then
    echo "Usage: $(basename $0) [input xyz file] <phi> <v> <psi>";
    exit 1;
fi

# assign input xyz file path
infile=$1

# define pi
pi=3.14159265

# create random Euler angles if necessary
if [[ -z "$2" || -z "$3" || -z "$4" ]];
then
    echo "NB: Three angles not input, performing a completely random rotation!"
    # Create random Euler angles. Formulas are obtained from analytically inverting Eqs 13.38 from Brannon (2002). 
    # There is no arccos() in awk, so we need to use arctan() instead.
    f=$(awk -v pi=$pi -v r=$RANDOM 'BEGIN{r=r/32767.0; print (2.0*r-1.0)*pi}')
    v=$(awk -v pi=$pi -v r=$RANDOM 'BEGIN{r=r/32767.0; x=1.0-2.0*r; print 2.0*atan2( sqrt(1.0-x*x), 1.0+x)}')
    p=$(awk -v pi=$pi -v r=$RANDOM 'BEGIN{r=r/32767.0; print (2.0*r-1.0)*pi}')
else
    f=$2
    v=$3
    p=$4
    echo "Doing a rotation for ($f, $v, $p) radians."
fi

# go through the xyz structure, rotating each point separately
awk -v f=$f -v v=$v -v p=$p 'BEGIN{n=0;}{n++;}(n<=2){print $0}(n>2){x0=$2; y0=$3; z0=$4; x1=(cos(f)*cos(p)-sin(f)*cos(v)*sin(p))*x0 - (cos(f)*sin(p)+sin(f)*cos(v)*cos(p))*y0 + sin(f)*sin(v)*z0; y1=(sin(f)*cos(p)+cos(f)*cos(v)*sin(p))*x0 + (cos(f)*cos(v)*cos(p)-sin(f)*sin(p))*y0 - sin(v)*cos(f)*z0; z1=sin(v)*sin(p)*x0 + sin(v)*cos(p)*y0 + cos(v)*z0; print $1, x1, y1, z1, $5}' $infile > $infile.euler.rotated

exit 0;
