=====
Usage
=====

STaGE is used from the command line, by executing::

    stage.py [arguments]

One of two possible input argument is required, namely ``-s/--smiles [SMILES]`` or
``-i/--inputfile [INPUTFILE]``. A base name for the output is also required. This is
specified by the ``-o/--outputfile [OUTPUTFILE]`` option.

The possible arguments are:

-h, --help                                  Show a help message and exit.
-i, --inputfile INPUTFILE                   A molecule file that can be read by
    Open Babel. This cannot be combined with the ``-s/--smiles`` option.
-s, --smiles SMILES                         A SMILES string describing the input molecule.
    This cannot be combined with the ``-i/--inputfile`` option.
-o, --outputfile OUTPUTFILE                 The base name of the output. File extensions and force
    field extensions are added to the OUTPUTFILE.
-k, --keep_ligand_name                      Do not set the name of the ligand to "LIG", but keep the
    name in the input file.
-b, --box_buffer BOX_BUFFER                 Generate a dodecahedron shaped periodic box around the
    molecule with BOX_BUFFER as minimum distance (in nm) between the molecule and the edge
    of the box. Set to 0 to disable solvation. Default: 1.0
-w, --water WATER                            The solvent model to use. This solvent model will be specified in
    the topology files and the periodic box around the molecule will be filled with water
    molecules.
-p, --ph PH                                 Protonate the molecule according to this pH using Open Babel. This does
    not always give correct results. It is safer to provide correctly protonated input
    files.
-r, --retain_charges                        Do not generate force field specific charges. This requires an
    input file in mol2 format with charges specified. The charges in the input file will be
    used for all force fields.
-q, --charge_method CHARGE_METHOD           Use the specified charge method for all force
    fields. Some CHARGE_METHODs are provided by basic STaGE functionality:

      * ``am1bcc`` - AM1-BCC is the default charge model when using GAFF, and is based on AM1
        charges with an applied bond charge correction. The charges are assigned using
        ANTECHAMBER.
      * ``am1bcc-pol`` - These are AM1-BCC charges with extra polarization applied, using the
        modifications in the ``ffOptimization/bbc-pol-parm.dat`` file. The charges are assigned
        using ANTECHAMBER.
      * ``mmff94`` - Charges used in the MMFF94 force field. The charges are assigned using
        Open Babel.
      * ``eem`` - The Electronegativity Equalization Method are charges similar to B3LYP/6-31G*.
        The charges are assigned using Open Babel.

    There are also plugins providing the following charge models:

      * ``cm1a`` - A class IV charge model based on AM1 wave functions.
        The charges are assigned using AMSOL.
      * ``cm3a`` - CM3A charges are similar to CM1A, but developed using a larger training set.
        The charges are assigned using AMSOL.
      * ``sm54a`` - Solvation model SM5 using class IV charges from AM1.
        The charges are assigned using AMSOL.
      * ``b3lyp/pcm`` - These charges are calculated using the B3LYP functional method using a
        polarizable water model (c-PCM) and a cc-pV(T+d)Z basis set. The QM calculations are
        performed using GAMESS/US. RESP charges are applied using gmstoresp.sh.

-f, --charge_multiplier CHARGE_MULTIPLIER   Multiply partial charges with this factor.
    Can only be used in combination with ``-q/--charge_method``.
-x, --calibration CALIBRATION               Modify van der Waals parameters according to the
    provided calibration file.
-c, --mergecoordinates MERGECOORDINATES     Combine the molecule coordinate files with
    an already existing coordinate file (.pdb or .gro). This can be useful for combining
    a ligand with a protein (if the coordinates match), e.g., for setting up input files for
    binding free energy calculations. If a .gro file of a protein is provided and there exists
    a corresponding (same base name) .top file that topology file will be used for the protein,
    otherwise a new topology file is generated.
--forcefields FORCEFIELDS                   A comma-separated string listing the force fields to
    generate parameters for. The available force fields are determined by the force field plugins.
    By default parameters for all forcefields are generated; gaff,cgenff,opls.
-v, --verbose                               Give verbose output of the process. This provides output from most of the
    programs that are executed and also more output from STaGE itself. This is recommended if
    the output is unexpected or if STaGE crashes.

------------
Output files
------------

A number of files are generated upon successful execution, e.g.

OUTPUTFILE.gro
    A GROMACS coordinate file of the molecule

OUTPUTFILE.mol2
    If there was not already a .mol2 for the molecule it will be generated.

posre_OUTPUTFILE.itp
    A GROMACS position restraints file that can be used by all force fields.

There are also files generated for each force field:

OUTPUTFILE_FORCEFIELD/OUTPUTFILE.top
    A GROMACS topology file.

OUTPUTFILE_FORCEFIELD/OUTPUTFILE.itp
    This file contains the force field parameters required for the molecule.

OUTPUTFILE_FORCEFIELD/index.ndx
    A GROMACS index file for selecting atoms and groups of atoms.

OUTPUTFILE_FORCEFIELD/box.gro
    If the molecule is prepared for being solvated this .gro file contains the size
    of the solvent box.

OUTPUTFILE_FORCEFIELD/solvated.gro
    If the molecule is solvated this is the file containing the molecule and the
    solvent.

OUTPUTFILE_FORCEFIELD/solvated_ionised.gro
    If the molecule is solvated and counter ions are required to neutralise the system
    this is the file containing them solvated system with counter ions.

----------------------
Force field parameters
----------------------

STaGE comes with a few additions to the standard GROMACS force fields, namely the TIP3P-M25,
TIP3P-MOD and OPC water models. In order to use these water models they must be found by GROMACS.
This can be done by either copying the contents of the forcefields/ directory to the location of
your GROMACS force fields. An alternative is to set the GMXLIB environment variable to
point at the forcefields/ directory, e.g.

``export GMXLIB=[location of your STaGE installation]/forcefields``

--------
Examples
--------

``stage.py -s 'CCO' -o ethanol``

Generates parameters for ethanol (based on SMILES input) for all available force fields using the
default charge methods for the force fields.

``stage.py -s 'CCO' -o ethanol -q b3lyp/pcm --forcefields gaff``

Generates parameters for ethanol (based on SMILES input) for GAFF using the B3LYP/PCM charge method, which
runs GAMESS/US to calculate the charges and can take a long time.

``stage.py -s 'CCO' -o ethanol -q CM1A -f 1.14 --forcefields opls``

Generates parameters for ethanol (based on SMILES input) for OPLS-AA using the CM1A charge method, which
runs AMSOL to calculate the charges, which are then multiplied by 1.14.

``stage.py -i ethanol.mol2 -o ethanol -r``

Generates parameters for ethanol using a mol2 file as input and retaining all the partial charges that
were specified in the mol2 file.

``stage.py -s 'CC(=O)O' -o acetic_acid``

Generates parameters for acetic acid (fully protonated by default) for all available force fields.

``stage.py -s 'CC(=O)O' -o acetic_acid -p 7.4``

Generates parameters for acetic acid at pH 7.4, i.e. the deprotonated form, for all available force fields.

``stage.py -s 'CC(=O)O' -o acetic_acid -p 7.4 -w opc -d 1.2``

Generates parameters for acetic acid at pH 7.4 and solvates the system in OPC water with a minimum distance
of 1.2 nm from the solute to the edge of the periodic box. Since the system is charged counter ions (1 Na+)
will also be added. Parameters will be generated for all force fields.