# STaGE3 (Small molecule Topology GEnerator for python 3) 

STaGE3 is the automatic GROMACS Topology Generation tool of small organic molecules using the GAFF and CGenFF force fields with receptor system such as proteins or macro molecules.  
STaGE3 is the python 3 fork of STaGE (https://gitlab.com/gromacs/stage) with some update improving usability.    
If you use STaGE3, please cite original paper: Lundborg M., Lindahl E. Automatic GROMACS TopologyGeneration and Comparisons of Force Fields for Solvation Free Energy Calculations. J. Phys. Chem. B. 2014, DOI: 10.1021/jp505332p

## Licence
GNU GPLv3

## Required software

STaGE3 also depends on a number of other softwares, some of which can should be installed seperately. They are:  

- Required  
1. Antechamber (part of Ambertools, available after filling in form at http://ambermd.org/AmberTools-get.html )  
2. acpype (https://alanwilter.github.io/acpype/)
3. openbabel (http://openbabel.org/wiki/Main_Page)  
4. gromacs (https://www.gromacs.org/) 

- Optional  
5. MATCH (available after filling in form at http://brooks.chem.lsa.umich.edu/index.php?page=registerSoftware&subdir=articles/resources/software&link=%3Ca%20href=%22downloads/MATCH_RELEASE.tar.gz%22%3EVersion%201.000%3C/a%3E&name=MATCH ) 
6. Amsol (available after filling in form at http://t1.chem.umn.edu/license/form-user.html )  
7. Gaussian 16 or Gaussian 09 (https://gaussian.com/gaussian16/)  

## Installation (Required)

- Create conda virtual environment  
<pre>
conda create -n py38-stage3 python=3.8  
conda activate py38-stage3  
</pre>

- Run the following command to install required conda packages  
<pre>
conda install -c conda-forge numpy ambertools acpype openbabel gromacs  
</pre>

- Install STaGE3 from github  
<pre>
pip install git+https://github.com/mkatouda/stage3.git
</pre>

- Install STaGE3 from local repository  
<pre>
git clone https://github.com/mkatouda/stage3.git
cd stage3
pip install .
</pre>

## Installation (Optional)

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
usage: stage3 [-h] [-i INP] [-l LIGAND] [-s SMILES] [-o OUTPUT] [--ffligand FFLIGAND]  
              [--ffprotein FFPROTEIN] [-x CALIBRATION] [-k] [-p PH] [-r]  
              [-q CHARGE_METHOD] [-f CHARGE_MULTIPLIER] [-c MERGECOORDINATES]  
              [-t MERGETOPOLOGY] [-b BOX_TYPE] [-d BOX_BUFFER] [-w WATER]  
              [--conc CONC] [--pname PNAME] [--nname NNAME] [-v]
</pre>

<pre>
optional arguments:  
  -h, --help            show this help message and exit  
  -i INP, --inp INP     yaml style input file, overwriting argument values (default: None)  
  -l LIGAND, --ligand LIGAND  
                        ligand (PDBQT, MOL, SDF, MOL2, PDB) (default: None)  
  -s SMILES, --smiles SMILES  
                        Use the specified smiles string as input instead of an input file  
                        (must be inside quotes). (default: None)  
  -o OUTPUT, --output OUTPUT  
                        Name of the output files (file extensions will be appended).   
                        (default: None)  
  --ffligand FFLIGAND   Force fields to generate parameters for, specified as  
                        a comma-separated string without spaces. (default: gaff)  
  --ffprotein FFPROTEIN  
                        Force field of protein. (default: None)  
  -x CALIBRATION, --calibration CALIBRATION  
                        Modify van der Waals parameters according to specified calibration  
                        file. (default: None)
  -k, --keep_ligand_name  
                        Do not rename the ligand in the output files. When doing  
                        e.g. solvation or binding free energy it is convenient  
                        to always call the ligand the same thing - in this case "LIG".    
                        If this option is set the ligand name will not be changed to "LIG".  
                        If you need to assign parameters to e.g. co-factors it is good  
                        to keep their names to tell them apart from ligands. (default: False)  
  -p PH, --ph PH        Protonate the molecule according to this pH (float). 
                        This does not always give correct results.  
                        It is safer to provide correctly protonated input files. (default: None)  
  -r, --retain_charges  Keep the mol2 charges. (default: False)  
  -q CHARGE_METHOD, --charge_method CHARGE_METHOD  
                        Use the specified charge method for all force fields.am1bcc: AM1 with bond charge correction (antechamber)  
                        am1bcc-pol: STaGE's own more polarized bond charge correction (antechamber)  
                        mmff94: MMFF94 (Open Babel)  
                        eem: electronegativity equalization method (Open Babel)  
                        qeq: Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991) (Open Babel)  
                        qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)  
                        gaussian/hf: Hatree-Fock/6-31G(d) basis set followed by RESP (Gaussian)  
                         (default: am1bcc)  
  -q CHARGE_METHOD, --charge_method CHARGE_METHOD  
                        Use the specified charge method for all force fields. (default: am1bcc)  
  -f CHARGE_MULTIPLIER, --charge_multiplier CHARGE_MULTIPLIER  
                        Multiply partial charges with this factor. Can only be used  
                        in combination with --charge_method. (default: 1.0)  
  -c MERGECOORDINATES, --mergecoordinates MERGECOORDINATES  
                        Merge the created coordinates file (.gro) with an already  
                        existing coordinate file (.pdb or .gro), e.g. for combining ligand  
                        coordinates with protein coordinates. The generated  
                        topology will contain both the ligand and the protein.  
                        If a .gro file of the protein is provided and there exists  
                        a corresponding .top file that toplogy file will be used for  
                        the protein, otherwise a new topology file is generated. (default: None)  
  -t MERGETOPOLOGY, --mergetopology MERGETOPOLOGY  
                        Merge the created topology file (.top) with an already existing  
                        topology file. Must be used in combination with --mergecoordinates  
                        with a .gro file of the protein. (default: None)  
  -b BOX_TYPE, --box_type BOX_TYPE  
                        Buffer from the solute to the edge of the dodecahedron shaped solvent  
                        box. Set to 0 to disable solvation (and ionisation).  
                        Default: dodecahedron (default: dodecahedron)  
  -d BOX_BUFFER, --box_buffer BOX_BUFFER  
                        Buffer from the solute to the edge of the solvent box.  
                        Set to 0 to disable solvation (and ionisation). (default: 1.0)
  -w WATER, --water WATER  
                        Solvent model to use in topology files. If not specified  
                        the solvent will not be specified in the topology.  
                        Suggested water models are: "opc", "spce", "tip4pew", "spc" or  
                        "tip3p". (default: None)  
  --conc CONC           Specify salt concentration (mol/liter). (default: 0.0)  
  --pname PNAME         Name of the positive counter ion in Solvent. (default: NA)  
  --nname NNAME         Name of the negative counter ion in Solvent. (default: CL)  
  -v, --verbose         Verbose output. (default: False)  

</pre>

## Exmaples of command line usage

### Simple small molecule

Generates MD input files for ethanol using mol file as input for all available force fields using the default charge methods.  

<pre>
stage3 -l ethanol.mol -o ethanol
</pre>

### Simple small molecule given user determied partial charges

Generates MD input files for ethanol using a mol2 file as input and retaining all the partial charges that were specified in the mol2 file.  

<pre>
stage3 -i ethanol.mol2 -o ethanol -r
</pre>

### Small molecule with B3LYP/PCM charge in water

Generates MD input files for water solvated ethanol using mol file as input.  
GAFF with Hatree-Fock 6-31G* charge method, which runs Gaussian to calculate the charges and can take a long time, is used.  
Mininum distance of 1.2 nm from the solute to the edge of the dodecahedron periodic box. TIP3P water model is used.  

<pre>
stage3 -l ethanol.mol -o ethanol_solvated --ffligand gaff -q gaussian/hf -w tip3p -b dodecahedron -d 1.2
</pre>

### Protein-ligand binding system

Generates MD input files for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.  
GAFF2 with AM1-BCC charge method and AMBER99SB-ILDN force fields are used for ligand and protein, respectively.  
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.  
TIP3P water model is used and the system is neutralized adding charged counter ions (K+) or (Cl-).  

<pre>
stage3 -l ligand.mol -c protein.pdb -o protein_ligand_solvated --ffligand gaff2 -q am1bcc --ffprotein charmm27 \
       -w tip3p -b cubic -d 1.0 --pname K --nname CL
</pre>

## Exmaples of yaml input usage

### Protein-ligand binding system

Generates MD input files for water solvated protein-ligand complex using protein pdb file and ligand mol file as input.  
GAFF2 with AM1-BCC charge method and AMBER99SB-ILDN force fields are used for ligand and protein, respectively.  
Mininum distance of 1.0 nm from the solute to the edge of the cubic periodic box.  

Prepare input yaml file input.yml:

<pre>
ligand: './jz4.mol'
mergecoordinates: './3HTB_protein.pdb'
output: '3HTB-jz4-wat'
ffligand: 'gaff'
ffpritein: 'amber99sb-ildn'
charge_method: 'am1bcc'
water: 'tip3p'
box_type: 'dodecahedron'
box_buff: 1.0
conc: 0.1
pname: 'NA'
nname: 'CL'
verbose: True
</pre>

Then, run stage3 in command line:

<pre>
stage3 -i input.yml
</pre>

Keywards of yaml file are the same in the name of command line options.  
See above explation of command line options.  

## Author
Michio Katouda (katouda@rist.or.jp)  
