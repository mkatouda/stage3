#! /bin/sh
ID='$Id: gmstoresp.sh 6001 2005-06-23 15:51:30Z ckarney $'

usage="
$ID
(c) 2004, Sarnoff Corporation, Princeton, NJ, USA

This shell script allows Gamess to be used instead of Gaussian for
assigning RESP charges using Antechamber.  The input is a Gamess output
file and this is converted to an antechamber file (with atom types, atom
positions, bond types and partial charges).  Partial charges are
calculated using the RESP model and Gamess must be run in such a way as
to compute the electrostatic potential at the points required by RESP.
Typical input to Gamess should therefore include

 \$ELPOT IEPOT=1 WHERE=PDC OUTPUT=BOTH \$END
 \$PDC PTSEL=CONNOLLY CONSTR=NONE PTDENS=1.68017 \$END

1.68017 = 6*bohr^2 gives a density of fitting points (6 per A^2) which
matches the antechamber's recommended input file for Gaussian.  A
suitable input file can be created with gjftogms -r.

Run as a filter, thus

    gms methanol 2>&1 | tee methanol.log | $0 [-d] [-h] > methanol.ac

Optional argument -d retains the temporary directory (under ${TMPDIR:-/tmp})
for debugging.

Optional argument -h prints this message.

gamesstoac depends upon OpenBabel, Antechamber and RESP.  These must be
installed and be in your PATH.

For more info see:
   http://amber.scripps.edu
   http://amber.scripps.edu/antechamber/antechamber.html
   http://openbabel.sourceforge.net
   http://www.msg.ameslab.gov/GAMESS/GAMESS.html
"
DEBUG=
while getopts dh c; do
    case $c in
	d ) DEBUG=y;;
	h ) echo "usage: $usage"; exit;;
	* ) echo "usage: $usage" 1>&2; exit 1;;
    esac
done
shift `expr $OPTIND - 1`

# If more than zero arguments passed
if [ $# -ne 0 ]; then
   echo "usage: $usage" 1>&2
   exit 1
fi

if [ -z "$DEBUG" ]; then
   # If not debugging, set trap to cleanup temp directory
   TEMP=
   trap 'trap "" 0; test "$TEMP" && rm -rf "$TEMP"; exit 1' 1 2 3 9 15
   trap		   'test "$TEMP" && rm -rf "$TEMP"'	       0
fi

# Create temporary directory under system tmp directory
TEMP=`mktemp -d ${TMPDIR:-/tmp}/gtoaXXXXXXXX`

# If problems...
if [ $? -ne 0 ]; then
   echo "$0: Can't create temp directory, exiting..." 1>&2
   exit 1
fi

[ "$DEBUG" ] && echo "Saving intermediate files in $TEMP" 1>&2

# Need to do this because antechamber and resp write files into the
# current directory and we don't want these to interfer with other
# invocations of these programs.
cd $TEMP

#
# Create an electrostatic potential (esp) file from a HyperChem HIN file
# and from a gamess log file.
#
MakeESPFile() {

   # Find and output number of atoms and esp fitting points.
   cat $1 | awk '{
      if (/TOTAL NUMBER OF ATOMS/)
	 natoms = $6;
      if (/NUMBER OF POINTS SELECTED FOR FITTING/)
	 nesp = $8;
   }
   END {
      printf "%5d%5d%5d\n", natoms, nesp, 0;
   }'

   # Find and output atom positions
   cat $2 | awk '
      BEGIN {
	 # Bohr radius in angstroms, 2002 CODATA value
	 bohr = 0.5291772108;
      }
      {
	 # Look for atom keyword and store x, y, z position
	 if ($1 == "atom")
	    # See resp.f for antechamber file formatting
	    printf "%33.7E%16.7E%16.7E\n",
	       $8 / bohr, $9 / bohr, $10 / bohr;
      }'


   # Find and output electrostatic potential information
   cat $1 | awk '
      BEGIN {
	 TRUE = 1;
	 FALSE = 0;
	 nesp = 0;
	 read_esp = FALSE;
      }
      {
	 if (/NUMBER OF POINTS SELECTED FOR FITTING/) {

	    # If the line marks the beginning of the electrostatic potential
	    # section, initialize the esp count and
	    # start reading the esp positions and values.
	    # By initializing the count each time we hit the beginning of
	    # and esp table, we are assured of reading the last esp
	    # table in the file.

	    nesp = 0;
	    read_esp = TRUE;

	 }
	 else if (/END OF PROPERTY EVALUATION/) {

	    # If we are past the electrostatic potiential section,
	    # stop reading the esp positions and values.

	    read_esp = FALSE;
	 }
	 else if (read_esp) {

	    # If we are supposed to be reading esp positions and values,
	    # remember positions and values for later.  Pick out
	    # position [$2,$3,$4] and potention $7.  However $1 and $2
	    # sometimes run into each other.  So use substr to peal of
	    # the index.

	    $1 = substr($1, 6);

	    esp[nesp]  = $7;
	    esp_x[nesp] = $2;
	    esp_y[nesp] = $3;
	    esp_z[nesp] = $4;
	    nesp++;
	 }
      }
      END {

	 # Write potential positions and values

	 for (i = 0; i < nesp; i++) {
	    # See resp.f for antechamber file formatting
	    printf "%17.7E%16.7E%16.7E%16.7E\n", \
		      esp[i], esp_x[i], esp_y[i], esp_z[i];
	 }
      }'
}

#
# Beginning of main shell script
#

GAMESSLOG="tmp.log"
AC="tmp.ac"
ACOUT="tmp.acout"
HIN="tmp.hin"

ESP="tmp.esp"
RESPIN1="tmp1.respin"
RESPIN2="tmp2.respin"
RESPOUT1="tmp1.respout"
RESPOUT2="tmp2.respout"
QOUT1="tmp1.qout"
QOUT2="tmp2.qout"

# Copy stdin to a file
cat > $GAMESSLOG

# Do a sanity check on the input file.  Does it contain the data we need
# to calculate resp charges?

if grep 'ELECTROSTATIC POTENTIAL' $GAMESSLOG > /dev/null; then :; else
   echo Error: Input file cannot be converted to Antechamber format. 1>&2
   echo Input file must contain electrostatic potential data 1>&2
   echo for RESP charge calculations. 1>&2
   exit 1
fi

# Run babel to get a HyperChem file from a gamess log file
babel -igam $GAMESSLOG -ohin $HIN | grep -v "^$" 1>&2

# Get total charge from gamess file
charge=`awk '{if (/CHARGE OF MOLECULE/) print $5}' $GAMESSLOG | tail -1`

# Run antechamber on HyperChem file to get an antechamber ac file
antechamber -i $HIN -fi hin -o $AC -fo ac -nc $charge | grep -v "^$" 1>&2

# Get electrostatic potentials and atom positions from gamess file
MakeESPFile $GAMESSLOG $HIN > $ESP

# Sanity check on esp file.  Need to extract integers mimicing Fortran's 2i5.
if [ `cat $ESP | wc -l` -ne `head -1 $ESP |
     awk '{print substr($0,1,5) + substr($0,6,5) + 1}'` ]; then
   echo "Miscount on number of atoms or fitting points" 1>&2
   exit 1
fi

# Create files that resp can read.  Run it twice because
# calculating RESP charges is a two stage process.
# The second call to respgen will generate table of symmetries.
respgen -i $AC -o $RESPIN1 -f resp1 | grep -v "^$" 1>&2
respgen -i $AC -o $RESPIN2 -f resp2 | grep -v "^$" 1>&2

# Calculate RESP charges
resp -O -i $RESPIN1 -o $RESPOUT1 -e $ESP -t $QOUT1	     | grep -v "^$" 1>&2
resp -O -i $RESPIN2 -o $RESPOUT2 -e $ESP -q $QOUT1 -t $QOUT2 | grep -v "^$" 1>&2

# Run antechamber again to copy the RESP charges back into the AC file
# Use -j 0 to turn off bond typing (this was already done).
antechamber -i $AC -fi ac -o $ACOUT -fo ac -c rc -cf $QOUT2 -j 0 |
   grep -v "^$" 1>&2

# Send to stdout
cat $ACOUT
