# STaGE3 (Small molecule Topology GEnerator for python 3) 

STaGE3 is the automatic GROMACS Topology Generation tool of organic molecules using the GAFF, OPLS-AA, and CGenFF force fields.  
STaGE3 is the python 3 port of STaGE (https://gitlab.com/gromacs/stage).  
If you use STaGE3, please cite original paper: Lundborg M., Lindahl E. Automatic GROMACS TopologyGeneration and Comparisons of 
Force Fields for Solvation Free Energy Calculations. J. Phys. Chem. B. 2014, DOI: 10.1021/jp505332p

## Licence
GNU GPLv3

## Required software

STaGE3 also depends on a number of other softwares, some of which can should be installed seperately.  
They are:  
- Required  
1. Antechamber (part of Ambertools, available after filling in form at http://ambermd.org/AmberTools-get.html )  
2. acpype (https://alanwilter.github.io/acpype/)
3. MATCH (available after filling in form at http://brooks.chem.lsa.umich.edu/index.php?page=registerSoftware&subdir=articles/resources/software&link=%3Ca%20href=%22downloads/MATCH_RELEASE.tar.gz%22%3EVersion%201.000%3C/a%3E&name=MATCH ) 
4. openbabel (http://openbabel.org/wiki/Main_Page)  
5. gromacs (https://www.gromacs.org/) 
- Optional  
6. Amsol (available after filling in form at http://t1.chem.umn.edu/license/form-user.html )  
7. GAMESS/US (http://www.msg.ameslab.gov/gamess/License_Agreement.html)  
8. Gaussian16 (https://gaussian.com/gaussian16/)  

## Installation (Required)

- Create conda virtual environment  
<pre>
conda create -n py38-stage python=3.8  
conda activate py38-stage  
</pre>

- Run the following command to install required conda packages  
<pre>
conda install -c conda-forge numpy ambertools acpype openbabel  
</pre>

- Resiter and download MATCH program from offical web cite and extract archive file
<pre>
tar xzvf MATCH_RELEASE.tar.gz
</pre>

- Set up enviroment variables for MATCH
<pre>
export PerlChemistry=/path/to/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/path/to/MATCH_RELEASE/MATCH
export PATH=${PATH}:${MATCH}/scripts
</pre>

## Command usage

<pre>
stage.py [-h] [-i INPUTFILE] [-s SMILES] -o OUTPUTFILE [-k] [-b BOX_TYPE] [-d BOX_BUFFER] [-w WATER] [-p PH] [-r] [-q {am1bcc,am1bcc-pol,mmff94,eem,qeq,qtpie}] [-f CHARGE_MULTIPLIER] [-x CALIBRATION] [-c MERGECOORDINATES] [--forcefields FORCEFIELDS] [--ffprotein FFPROTEIN] [--pname PNAME] [--nname NNAME] [-v]
</pre>

<pre>
optional arguments:  
  -h, --help            show this help message and exit  
  -i INPUTFILE, --inputfile INPUTFILE  
                        Input file name.  
  -s SMILES, --smiles SMILES  
                        Use the specified smiles string as input instead of an input file (must be inside quotes).  
  -o OUTPUTFILE, --outputfile OUTPUTFILE  
                        Name of the output files (file extensions will be appended).  
  -k, --keep_ligand_name  
                        Do not rename the ligand in the output files. When doing e.g. solvation or binding free energy it is convenient to always call the ligand the same thing - in this case "LIG". If this option is set the ligand name will not be changed to "LIG". If you need to assign parameters to e.g. co-factors it is good to keep their names to tell
                        them apart from ligands.  
  -b BOX_TYPE, --box_type BOX_TYPE  
                        Buffer from the solute to the edge of the dodecahedron shaped solvent box. Set to 0 to disable solvation (and ionisation). Default: dodecahedron  
  -d BOX_BUFFER, --box_buffer BOX_BUFFER  
                        Buffer from the solute to the edge of the solvent box. Set to 0 to disable solvation (and ionisation). Default: 1.0  
  -w WATER, --water WATER  
                        Solvent model to use in topology files. If not specified the solvent will not be specified in the topology. Suggested water models are: "opc", "spce", "tip4pew", "spc" or "tip3p".  
  -p PH, --ph PH        Protonate the molecule according to this pH (float). This does not always give correct results. It is safer to provide correctly protonated input files.  
  -r, --retain_charges  Keep the mol2 charges.  
  -q {am1bcc,am1bcc-pol,mmff94,eem,qeq,qtpie}, --charge_method {am1bcc,am1bcc-pol,mmff94,eem,qeq,qtpie}  
                        Use the specified charge method for all force fields. am1bcc: AM1 with bond charge correction (antechamber), am1bcc-pol: STaGE's own more polarized bond charge correction (antechamber), mmff94: MMFF94 (Open Babel), eem: electronegativity equalization method (Open Babel), qeq: Assign QEq (charge equilibration) partial charges (Rappe
                        and Goddard, 1991) (Open Babel), qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)  
  -f CHARGE_MULTIPLIER, --charge_multiplier CHARGE_MULTIPLIER  
                        Multiply partial charges with this factor. Can only be used in combination with --charge_method.  
  -x CALIBRATION, --calibration CALIBRATION  
                        Modify van der Waals parameters according to specified calibration file.  
  -c MERGECOORDINATES, --mergecoordinates MERGECOORDINATES  
                        Merge the created coordinates file (.gro) with an already existing coordinate file (.pdb or .gro), e.g. for combining ligand coordinates with protein coordinates. The generated topology will contain both the ligand and the protein. If a .gro file of the protein is provided and there exists a corresponding .top file that toplogy file
                        will be used for the protein, otherwise a new topology file is generated.  
  --forcefields FORCEFIELDS  
                        Force fields to generate parameters for, specified as a comma-separated string without spaces. Default: gaff,cgenff  
  --ffprotein FFPROTEIN  
                        Force field of protein.  
  --pname PNAME         Name of the positive counter ion in Solvent. Default: NA  
  --nname NNAME         Name of the negative counter ion in Solvent. Default: CL  
  -v, --verbose         Verbose output.  
</pre>

## Exmaples of usage

Generates parameters for ethanol using mol file as input for all available force fields using the default charge methods.  

<pre>
/path/to/stage.py -i ethanol.mol -o ethanol
</pre>

Generates parameters for ethanol using a mol2 file as input and retaining all the partial charges that were specified in the mol2 file.  

<pre>
/path/to/stage.py -i ethanol.mol2 -o ethanol -r
</pre>

Generates parameters for water solvated ethanol using mol file as input.
GAFF with B3LYP/PCM charge method, which runs GAMESS/US to calculate the charges and can take a long time, is used.
Mininum distance of 1.2 nm from the solute to the edge of the dodecahedron periodic box. TIP3P water model is used.  

<pre>
/path/to/stage.py -i ethanol.mol -o ethanol_solvated --forcefields gaff -q b3lyp/pcm -w tip3p -b dodecahedron  -d 1.2
</pre>

Generates parameters for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.
CGenFF with AM1-BCC charge method and Charmm27 force fields are used for ligand and protein, respectively.
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.
TIP3P water model is used and the system is neutralized adding charged counter ions (K+) or (Cl-).  

<pre>
/path/to/stage.py -i ligand.mol -c protein.pdb -o protein_ligand_solvated --forcefields cgenff -q am1bcc --ffprotein charmm27 -w tip3p -b cubic -d 1.0 --pname K --nname CL
</pre>

## Author
Michio Katouda (katouda@rist.or.jp)  
