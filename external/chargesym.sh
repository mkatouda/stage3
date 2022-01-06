#! /bin/sh
ID='$Id: chargesym.sh 6001 2005-06-23 15:51:30Z ckarney $'

usage="
$ID
(c) 2004, Sarnoff Corporation, Princeton, NJ, USA

This shell script symmetrizes the partial charges found in an input
Antechamber file.  New Partial charges are calculated by: (1) calling
respgen on an input Antechamber file to identify charge groups,
(2) summing original partial charges across all group members (atoms),
(3) assigning the mean partial charge of the group to each group member
(atom).

Run as a filter, thus

    gms methanol 2>&1 | tee methanol.log | gmstobcc | $0 > methanol.ac

Optional argument -d retains the temporary directory (under ${TMPDIR:-/tmp})
for debugging.

Optional argument -h prints this message.

symmetrize depends upon respgen and Antechamber.  These must be
installed and be in your PATH.  Furthermore environment variable
AMBERHOME must be defined for Antechamber.

For more info see:
   http://amber.scripps.edu
   http://amber.scripps.edu/antechamber/antechamber.html
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
   trap            'test "$TEMP" && rm -rf "$TEMP"'            0
fi

# Create temporary directory under system tmp directory
TEMP=`mktemp -d ${TMPDIR:-/tmp}/csymXXXXXXXX`

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
# This function takes the an antechamber charge file as $1 and
# a respgen files as $2.  Sums the charges in a charge
# group, calculates the average and distributes the result across the
# group.
#
SymmetrizeCharges() {

   awk -v charges="$1" '
    BEGIN {

    # Initialize some variables
	nStartLine = 0;
	i = 1;

    # Split the charges into an array
	nChargeArray = split(charges, ChargeArray);
    }
    {
    # Marker preceding the table
	if (/&end/) {

    # Start reading in symmetries four lines after the &end marker
	    nStartLine = NR + 4;
	}
	else if (nStartLine > 0 && NR >= nStartLine && NF > 0) {

    # Store symmetry info into array.  Add a reference for self.
	    SymmetryArray[i] = $2;
	    RefsArray[i] = 1;

    # If this atom a secondary atom in a charge group.
	    if (SymmetryArray[i] != 0) {

    # Add the charge of this atom charge to the charge of the primary atom
		ChargeArray[SymmetryArray[i]] += ChargeArray[i];

    # Add a reference to the primary atom
		RefsArray[SymmetryArray[i]] += 1;

    # Mark this atom charge as zero for now.
		ChargeArray[i] = 0;
	    }
	    i++;
	}
    }
    END {

    # Loop for each atom.
	for (i = 1; i <= nChargeArray ; i++) {

    # If this atom is the primary atom in charge group.
	    if (RefsArray[i] > 1) {

    # Calculate the group average charge.
		ChargeArray[i] = ChargeArray[i] / RefsArray[i];
	    }
	    else if (SymmetryArray[i] != 0) {

    # If this atom is a secondary atom in charge group,
    # distribute the group average to this atom.
		ChargeArray[i] = ChargeArray[SymmetryArray[i]];
	    }

    # Print the charge for this atom
	    printf "%10.6f", ChargeArray[i];
	    if (i % 8 == 0)
		printf "\n";
	}
    }'
}

#
# Beginning of main shell script
#

AC="tmp.ac"
ACOUT="tmp.acout"
RESPIN="tmp.respin"
QOUT1="tmp.qout1"
QOUT2="tmp.qout2"

# Copy stdin to a file
cat > $AC

# Run antechamber to get the unsymmetrized charges
antechamber -i $AC -fi ac -o $ACOUT -fo ac -c wc -cf $QOUT1 -j 0 |
   grep -v "^$" 1>&2

# Run respgen to get symmetries
respgen -i $AC -o $RESPIN -f resp | grep -v "^$" 1>&2

# Symmetrize the charges
SymmetrizeCharges "`cat $QOUT1`" < $RESPIN > $QOUT2

# Run antechamber to insert the symmetrized charges
antechamber -i $AC -fi ac -o $ACOUT -fo ac -c rc -cf $QOUT2 -j 0 |
   grep -v "^$" 1>&2

# Send to stdout
cat $ACOUT
