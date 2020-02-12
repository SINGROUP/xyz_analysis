#!/usr/bin/perl

# reads xyz file (incl. MD trajectory) and wraps all atoms outside the box back
# output: xyz2.xyz, 2nd-generation xyz format
# Change of basis, 2x2 matrix, Ref: xyz2vasp.pl by Stephan Blankenburg blankenburg@phys.upb.de
# yang@mpie.de, 18.05.2016
# lei.yang@aalto.fi, 12.02.2020

print "put in \$alat first \n";
@alat[0]=[9.6688995361,         0.0000000000,         0.0000000000];
@alat[1]=[3.2229630756,         9.1159271198,         0.0000000000];
@alat[2]=[0.0000000000,         0.0000000000,        24.8682003021];
#@alat[0]=[4, 0, 0]; 
#@alat[1]=[2, 4, 0]; 
#@alat[2]=[0, 0, 1];
#for ($k=0;$k<3;$k++){
#$l[$k]=sqrt($alat[$k][0]**2+$alat[$k][1]**2+$alat[$k][2]**2);
#}
#print "$l[1]\n";

open(IN,$ARGV[0]);
@in=<IN>;
$Nlines=$.;
close(IN);
#print "$Nlines\n";

($Nmol) = split(" ",$in[0]); # number of atoms
$Nrun=$Nlines/($Nmol+2);   # number of stemps
print "$Nmol, $Nrun\n";

# read data
for ($step=0;$step<$Nrun;$step++) {
  for ($mol=0;$mol<$Nmol;$mol++) {
    $line = ($Nmol+2)*$step+$mol+2;
    ($new[$step][$mol][3],$new[$step][$mol][0],$new[$step][$mol][1],$new[$step][$mol][2]) = split(" ",$in[$line]);
		# converting
    $detM=$alat[0][0]*$alat[1][1]*$alat[2][2]+$alat[0][1]*$alat[1][2]*$alat[2][0]+$alat[0][2]*$alat[1][0]*$alat[2][1]-$alat[2][0]*$alat[1][1]*$alat[0][2]-$alat[2][1]*$alat[1][2]*$alat[0][0]-$alat[2][2]*$alat[1][0]*$alat[0][1];
    $xi[0]=($alat[1][1]*$alat[2][2]-$alat[2][1]*$alat[1][2])/$detM;
    $yi[0]=($alat[2][1]*$alat[0][2]-$alat[0][1]*$alat[2][2])/$detM;
    $zi[0]=($alat[0][1]*$alat[1][2]-$alat[1][1]*$alat[0][2])/$detM;
    $xi[1]=($alat[2][0]*$alat[1][2]-$alat[1][0]*$alat[2][2])/$detM;
    $yi[1]=($alat[0][0]*$alat[2][2]-$alat[2][0]*$alat[0][2])/$detM;
    $zi[1]=($alat[1][0]*$alat[0][2]-$alat[0][0]*$alat[1][2])/$detM;
    $xi[2]=($alat[1][0]*$alat[2][1]-$alat[2][0]*$alat[1][1])/$detM;
    $yi[2]=($alat[2][0]*$alat[0][1]-$alat[0][0]*$alat[2][1])/$detM;
    $zi[2]=($alat[0][0]*$alat[1][1]-$alat[1][0]*$alat[0][1])/$detM;

		$a[0]=$new[$step][$mol][0]*$xi[0]+$new[$step][$mol][1]*$xi[1]+$new[$step][$mol][2]*$xi[2];
		$a[1]=$new[$step][$mol][0]*$yi[0]+$new[$step][$mol][1]*$yi[1]+$new[$step][$mol][2]*$yi[2];
		$a[2]=$new[$step][$mol][0]*$zi[0]+$new[$step][$mol][1]*$zi[1]+$new[$step][$mol][2]*$zi[2];
#		print "$a[0] $a[1] $a[2]\n";
    # processing
    for ($i=0;$i<10;$i++) { 
       for ($k=0;$k<3;$k++) {
#					 $proj[$k]=($new[$step][$mol][0]*$alat[$k][0]+$new[$step][$mol][1]*$alat[$k][1]+$new[$step][$mol][2]*$alat[$k][2]) / $l[$k];
           if ($a[$k] < 0 ){
#						  print "$step $mol $k $a[$k] \n";
#							$exess[$k]=sqrt($new[$step][$mol][0]**2+$new[$step][$mol][1]**2+$new[$step][$mol][2]**2)/$proj[$k];
              for ($j=0;$j<3;$j++) {
								$new[$step][$mol][$j]=$new[$step][$mol][$j]+$alat[$k][$j]; 
#								              print "$new[$step][$mol][$j] \n";
              }
              $a[$k]=$a[$k]+1;
           }
					 elsif ($a[$k] > 1) {
              for ($j=0;$j<3;$j++) {
								$new[$step][$mol][$j]=$new[$step][$mol][$j]-$alat[$k][$j];
              }
              $a[$k]=$a[$k]-1;
					 }
       }
    }
  }
}



# write output data in new xyz
open(OUT,">$ARGV[0]2.xyz");
for ($step=0;$step<$Nrun;$step++) {
  print OUT $in[0];
  print OUT "Lattice=\"$alat[0][0] $alat[0][1] $alat[0][2] $alat[1][0] $alat[1][1] $alat[1][2] $alat[2][0] $alat[2][1] $alat[2][2]\" \n";
  for ($mol=0;$mol<$Nmol;$mol++) {
    print OUT "$new[$step][$mol][3] $new[$step][$mol][0] $new[$step][$mol][1] $new[$step][$mol][2] \n";
  }
#  print OUT "\n";
}
close(OUT);
