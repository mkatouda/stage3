#!/usr/bin/env python

"""
    
    STaGE3 VERSION 0.1.0

    STaGE3 is the automatic GROMACS Topology Generation tool of organic 
    molecules using the GAFF, OPLS-AA, and CGenFF force fields.
    STaGE3 is the python 3 fork of STaGE (https://gitlab.com/gromacs/stage).
    If you use STaGE3, please cite original paper: Lundborg M., Lindahl E.
    Automatic GROMACS TopologyGeneration and Comparisons of Force Fields
    for Solvation Free Energy Calculations. J. Phys. Chem. B. 2014,
    DOI: 10.1021/jp505332p

    STaGE3 can use ACPYPE (calling AnteChamber) and MATCH to generate
    topologies for GAFF, OPLS and CGENFF and compatible charges.
    OPLS atom types and charges, from ACPYPE, are improved.
    The range of parameter assignment tools might be extended in the future.

    STaGE3 was forked by Michio Katouda
    Copyright (c) 2021-2022, 
    Check out https://github.com/mkatouda/stage3 for more information.

    Originl STaGE was written by Magnus Lundborg
    Copyright (c) 2013-2015, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.

    The following applications must be available in the list of directories of
    executable programs (the PATH environment variable):

    acpype
    antechamber
    amsol7.1
    tleap
    babel
    MATCH.pl
    gmx pdb2gmx
    gmx editconf
    gmx solvate (or genbox)
    gmx grompp
    gmx genion
    gmx make_ndx
    gmx genrestr
    rungms
    gmstoresp.sh
    g16

    The following environment variables must be set

    AMBERHOME

    The following environment variables may be set when using optional programs

    MATCH
    PerlChemistry

    It might be good to set the GMXLIB environment variable to point to the ../top/ directory of
    the GROMACS installation.


    Requirements: Python version >= 3.0
                  Open Babel version >= 2.3
                  GROMACS (pdb2gmx, make_ndx, genrestr, editconf, solvate (genbox), grompp and genion)
                  ACPYPE (for GAFF and OPLS)
                  AnteChamber (for GAFF and OPLS)

    Optional:     MATCH (for CGENFF)
                  Gaussian (for HF/6-31G* RESP partial charges)
                  GAMESS (for B3LYP/PCM partial charges: NYI)
                  gmstoresp.sh


    ACPYPE
    Sousa da Silva, A. W., Vranken, W. ACPYPE - AnteChamber PYthon Parser
    interfacE. BMC Res. Notes 2012, 5: 367.

    AnteChamber
    1.  Wang, J., Wang, W., Kollman P. A., Case, D. A. Automatic atom type and
        bond type perception in molecular mechanical calculations.
        J. Mol. Graph. and Model. 2006, 25: 247-260.
    2.  Wang, J., Wolf, R. M., Caldwell, J. W., Kollman, P. A., Case, D. A.
        Development and testing of a general AMBER force field.
        J. Comput. Chem. 2004, 25: 1157-1174.

    GAMESS/US
    Schmidt, M. W.; Baldridge, K. K.; Boatz, J. A.; Elbert, S. T.; Gordon,
    M. S.; Jensen, J. H.; Koseki, S.; Matsunaga, N.; Nguyen, K. A.; Su, S.;
    Windus, T. L.; Dupuis, M.; Montgomery, J. A. General atomic and
    molecular electronic structure system. J. Comput. Chem. 1993,
    14: 1347-1363.

    gmstoresp.sh
    (c) 2004, Sarnoff Corporation, Princeton, NJ, USA

    GROMACS
    1.  Pronk, S., Pall, S., Schulz, R., Larsson, P., Bjelkmar, P., Apostolov, R.,
        Shirts, M. R., Smith, J. C., Kasson, P. M., van der Spoel, D., Hess, B., Lindahl, E.
        GROMACS 4.5: a high-throughput and highly parallel open source molecular
        simulation toolkit. Bioinformatics 2013, 29: 845-854.
    2.  Hess, B., Kutzner, C., van der Spoel, D., Lindahl, E. GROMACS 4: Algorithms
        for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation.
        J. Chem. Theory Comput. 2008, 4: 435-447.
    http://www.gromacs.org

    MATCH
    Yesselman, J. D., Price, D. J., Knight, J. L. Brooks, C. L. 3rd. MATCH: an
    atom-typing toolset for molecular mechanics force fields.
    J. Comput. Chem. 2012, 33: 189-202.

    Open Babel
    O' Boyle, M. O., Banck, M., James, C. A., Morley, C., Vandermeersch, T.,
    Hutchison, G. R. Open Babel: An open checmical toolbox.
    J. Cheminform. 2011, 3: 33.

"""

import sys
import os
import argparse
import shutil
import traceback
import inspect
from glob import glob

import yaml

from .config import ffligandsString, standardChargeMethods, chargeMethodsHelp, waterModels, boxTypes
from .ForceFieldPlugin import ForceFieldPlugin
from .GaffForceFieldPlugin import GaffForceFieldPlugin
from .CgenffForceFieldPlugin import CgenffForceFieldPlugin
#from .OplsForceFieldPlugin import OplsForceFieldPlugin
from .ChargePlugin import ChargePlugin
#from .GamessChargePlugin import GamessChargePlugin
from .GaussianChargePlugin import GaussianChargePlugin
from .util import babelConvert, renameAtoms, mol2RenameToLig, getNetChargeOfMol2, makeRestraintsRun, getChargeOfTopology
from .util import generateCharges, calibrateVdW, mergeCoordinateFiles, copyItp, modproteinItp, splitTopologyToItp, mergeTopologyFiles
from .util import hydrogens2VirtualSites, generateLinearVirtualSites, solvateSystem, neutraliseSystem, makeIndexRun, convertGmx2Amb


def _loadPlugins(path, name, baseclass):

    loadedPlugins = []

    plugins = glob(os.path.join(path, name))
    print(os.getcwd(), 'plugins: ', plugins)
    for i in range(len(plugins)):
        #plugins[i] = '.plugins.' + os.path.basename(plugins[i])[:-3]
        plugins[i] = 'from .plugins import ' + os.path.basename(plugins[i])[:-3]
        print(plugins[i])

    for plugin in plugins:
        print(plugin)

        try:
            p = __import__(plugin, fromlist = ['nonsense']) # fromlist must be a non-empty list
        except (ImportError,NotImplementedError):
            print('ImportError')
            continue

        for cls in dir(p):                         # Loop over all objects in the module's namespace.
            try:
                cls=getattr(p,cls)
                if (inspect.isclass(cls)              # It should be a class.
                    and inspect.getmodule(cls) == p   # Make sure it was defined in module, not just imported.
                    and issubclass(cls, baseclass)    # It should be a subclass of Converter.
                    and cls.__name__ != baseclass.__name__): # It should not be the Converter baseclass itself.
                        instance = cls()              # Create an object of the plugin class.
                        if not instance in loadedPlugins:
                            loadedPlugins.append(instance)
            except Exception:
                print('Error instantiating plugin:', plugin)
                traceback.print_exc()

    print('loadedPlugins:', loadedPlugins)

    print('baseclass:', baseclass)
    if baseclass == 'ForceFieldPlugin':
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order, c.forceFieldName))
    elif baseclass == 'ChargePlugin':
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order, c.programName))
    else:
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order))

    return loadedPlugins

#def getChargeMethod():
#    progDir = os.path.dirname(__file__)
#    standardChargeMethods = ['am1bcc', 'am1bcc-pol', 'mmff94', 'eem', 'qeq', 'qtpie']
#    chargeMethodsHelp = ['am1bcc: AM1 with bond charge correction (antechamber)',
#                         'am1bcc-pol: STaGE\'s own more polarized bond charge correction (antechamber)',
#                         'mmff94: MMFF94 (Open Babel)',
#                         'eem: electronegativity equalization method (Open Babel)',
#                         'qeq: Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991) (Open Babel)',
#                         'qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)']
#
#    chargePlugins = _loadPlugins(os.path.join(progDir, 'plugins'), '*ChargePlugin.py', ChargePlugin)
#    chargePluginNameList = []
#    chargePluginsHelp = []
#    for pl in chargePlugins:
#        for method, description in pl.alternativeChargeMethods.items():
#            chargePluginNameList.append(method)
#            chargeMethodsHelp.append('%s: %s (%s)' % (method, description, pl.programName))
#
#    chargeMethodsList = standardChargeMethods + chargePluginNameList
#    print('chargeMethodsList: ', chargeMethodsList)
#
#    return chargeMethodsList

def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        #formatter_class=argparse.RawTextHelpFormatter,
        formatter_class=customHelpFormatter,
        description='STaGE is a tool for generating GROMACS topologies\n'
        'of small molecules for different force fields (currently GAFF, GAFF2,\n'
        'and CGenFF). It uses many external programs to perform its tasks,\n'
        'so make sure you read the documentation to properly cite the programs if\n'
        'you use STaGE in a publication and also cite:\n'
        'Lundborg M., Lindahl E. Automatic GROMACS Topology Generation and Comparisons\n'
        'of Force Fields for Solvation Free Energy Calculations. J. Phys. Chem. B. 2014,\n'
        'DOI: 10.1021/jp505332p'
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values',
    )
    parser.add_argument(
        '-l', '--ligand', type=str,
        help="ligand (PDBQT, MOL, SDF, MOL2, PDB)"
    )
    parser.add_argument(
        '-s', '--smiles', type=str,
        help = 'Use the specified smiles string as input\n'
        'instead of an input file (must be inside quotes).'
    )
    parser.add_argument(
        '-o', '--output', type=str,
        help = 'Name of the output files (file extensions will be appended).'
    )
    parser.add_argument(
        '--ffligand', type=str, default='gaff',
        help = 'Force fields to generate parameters for, specified\n'
        'as a comma-separated string without spaces:\n'
        + ', '.join(ffligandsString)
    )
    parser.add_argument(
        '--ffprotein', type=str,
        help = 'Force field of protein.'
    )
    parser.add_argument(
        '-x', '--calibration', type=str,
        help = 'Modify van der Waals parameters according to specified\n'
        'calibration file.'
    )
    parser.add_argument(
        '-k', '--keep_ligand_name', action='store_true',
        help = 'Do not rename the ligand in the output files.\n'
        'When doing e.g. solvation or binding free energy\n'
        'it is convenient to always call the ligand the\n'
        'same thing - in this case "LIG". If this option\n'
        'is set the ligand name will not be changed to "LIG".\n'
        'If you need to assign parameters to e.g. co-factors\n'
        'it is good to keep their names to tell them apart\n'
        'from ligands.'
    )
    parser.add_argument(
        '-p', '--ph', type=float,
        help = 'Protonate the molecule according to this pH (float).\n'
        'This does not always give correct results. It is safer\n'
        'to provide correctly protonated input files.'
    )
    #parser.add_argument(
    #    '-y', '--virtualhydrogens', action='store_true',
    #    help = 'Turn hydrogens into virtual interaction sites to allow longer '
    #    'timesteps (experimental).'
    #)
    parser.add_argument(
        '-r', '--retain_charges', action='store_true',
        help = 'Keep the mol2 charges.'
    )
    parser.add_argument(
        '-q', '--charge_method', type=str, default='am1bcc',
        help = 'Use the specified charge method for all force fields:\n'
        + '\n'.join(chargeMethodsHelp) + '\n'
    )
    parser.add_argument(
        '-f', '--charge_multiplier', type=float, default=1.0,
        help = 'Multiply partial charges with this factor. Can only be used\n'
        'in combination with --charge_method.'
    )
    parser.add_argument(
        '-c', '--mergecoordinates', type=str,
        help = 'Merge the created coordinates file (.gro) with an\n'
        'already existing coordinate file (.pdb or .gro), '
        'e.g. for combining \n'
        'ligand coordinates with protein coordinates. The generated topology\n'
        'will contain both the ligand and the protein. If a .gro file of the\n'
        'protein is provided and there exists a corresponding .top file that\n'
        'toplogy file will be used for the protein, otherwise a new topology\n'
        'file is generated.'
    )
    parser.add_argument(
        '-t', '--mergetopology',
        help = 'Merge the created topology file (.top) with an\n'
        'already existing topology file.\n'
        'Must be used in combination with --mergecoordinates\n'
        'with a .gro file of the protein.'
    )
    parser.add_argument(
        '-b', '--box_type', type=str, default='dodecahedron',
        help = 'Type of simulation box: ' 
        + ', '.join(boxTypes)
    )
    parser.add_argument(
        '-d', '--box_buffer', type=float, default=1.0,
        help = 'Buffer from the solute to the edge of the\n'
        'solvent box. Set to 0 to disable solvation (and ionisation).'
    )
    parser.add_argument(
        '-w', '--water', type=str,
        help = 'Solvent model to use in topology files. If not \n'
        'specified the solvent will not be specified in \n'
        'the topology. Suggested water models are: \n'
        + ', '.join(waterModels)
    )
    parser.add_argument(
        '--conc', type=float, default=0.0,
        help = 'Specify salt concentration (mol/liter).'
    )
    parser.add_argument(
        '--pname', type=str, default='NA',
        help = 'Name of the positive counter ion in Solvent.'
    )
    parser.add_argument(
        '--nname', type=str, default='CL',
        help = 'Name of the negative counter ion in Solvent.'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help = 'Verbose output.'
    )
    args = parser.parse_args()

    #if not args.inp and not args.smiles:
    #    parser.error('No molecule input specified. Either -i or -s '
    #    'must be provided.')

    #if not args.output:
    #    parser.error('No output name specified. -o must be specified. '
    #    'File extensions will be appended to the name')

    #if args.inp and args.smiles:
    #    parser.error('--input and --smiles parameters cannot be used together.')

    #if args.retain_charges and args.charge_method:
    #    parser.error('--retain_charges and --charge_method parameters cannot be '
    #                 'used together.')

    #if args.output == args.inp:
    #    parser.error('Output name and input file may not be identical.')

    #if args.inp and not os.path.exists(args.inp):
    #    parser.error('The input file does not exists.')

    return args 

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    del conf_def['inp']
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def stage3_run(ligand, smiles, output, ffligand, ffprotein, calibration, 
               keep_ligand_name, ph, retain_charges, charge_method, charge_multiplier,
               mergecoordinates, mergetopology, box_type, box_buffer,
               water, conc, pname, nname, verbose=False):

    if mergecoordinates:
        proteinCoords = os.path.abspath(mergecoordinates)
        print('proteinCoords: ', proteinCoords)
        if os.path.splitext(proteinCoords)[1].lower() not in ('.gro', '.pdb'):
            print('The coordinates file to merge with must have a '
                  '.gro or .pdb extension')
            sys.exit(1)
    else:
        proteinCoords = None

    if mergetopology:
        if not proteinCoords or os.path.splitext(proteinCoords)[1].lower() != '.gro':
            print('A .gro coordinate file must be specified with '
                  '--mergecoordinates if --mergetopology is used.')
            sys.exit(1)
        else:
            proteinTopology = os.path.abspath(mergetopology)
    else:
        proteinTopology = None

    # chosenForcefields contains the actual force fields that were chosen by the user.
    # However, if generating OPLS parameters GAFF parameters must first be generated,
    # but these can safely be removed later on.
    chosenForcefields = ffligand.lower().split(',')
    forcefields = list(chosenForcefields)

    converters = []
    for forcefield in forcefields:
        if forcefield not in ffligandsString:
            print('Force field %s is not available.' % forcefield)
            sys.exit(1)
        elif 'gaff' in forcefield:
            converter = GaffForceFieldPlugin(forcefield)
#        elif forcefield == 'opls': # Warning: currently not available!
#            converter = OplsForceFieldPlugin()
        elif forcefield == 'cgenff':
            converter = CgenForceFieldPlugin()
        converters.append(converter)
    print('converters:', converters)

    if 'opls' in forcefields and not 'gaff' in forcefields:
        forcefields.append('gaff')

    if smiles:
        iBase, iFormat = None, 'smiles'
        outputFile = os.path.abspath(output)
    else:
        iBase, iFormat = os.path.splitext(ligand)
        iFormat = iFormat.lower()
        outputFile = os.path.abspath(output)
        #outputFile = iBase

    if ligand:
        inputFile = os.path.abspath(ligand)
    else:
        inputFile = None

    if verbose:
        print('iBase: ', iBase, 'iFormat:', iFormat)
        print('inputFile:', inputFile)
        print('outputFile:', outputFile)

    if ',' in outputFile or (inputFile and ',' in inputFile):
        print('The output and/or input files must not contain a comma (",")')
        sys.exit(1)

    if iFormat != '.mol2' or ph:
        if smiles:
            inMolecules = babelConvert(outputFile = outputFile,
                                       smiles = smiles,
                                       pH = ph,
                                       verbose = verbose)
        else:
            inMolecules = babelConvert(inputFile = inputFile, outputFile =
                                       outputFile, pH = ph,
                                       verbose = verbose)
    else:
        if inputFile != outputFile+'.mol2':
            shutil.copy(inputFile, outputFile+'.mol2')
        inMolecules = [outputFile+'.mol2']

    if verbose:
        print('inMolecules:', inMolecules)

    if not inMolecules:
        print('Could not convert molecule to mol2')
        sys.exit(1)

    if calibration and not os.path.exists(calibration):
        print('Calibration file not found. Continuing without calibrating vdW parameters.')
        calibration = None

    count = 0

    for inMolecule in inMolecules:

        if verbose:
            print('inMolecule: ', inMolecule)

        # Warning! This routine has some bugs for making wrong input mol2 file for acpype run.
        if 'cgenff' in forcefields:
            renameAtoms(inMolecule)

        if not keep_ligand_name:
            mol2RenameToLig(inMolecule)

        ligandCoords = os.path.splitext(inMolecule)[0] + '.gro'
        if verbose:
            print('ligandCoords: ', ligandCoords)

        if not os.path.exists(ligandCoords):
            try:
                babelConvert(inMolecule, ligandCoords, verbose = verbose)
            except Exception:
                print('Cannot create .gro file:')
                traceback.print_exc()

        count += 1

        netCharge = getNetChargeOfMol2(inMolecule)

        #MK
        if verbose:
            print('netCharge: ', netCharge)

        try:
            makeRestraintsRun(ligandCoords, verbose = verbose)
        except Exception:
            print('Cannot make restraints:')
            traceback.print_exc()

        totCharge = None
        ffSpecificCoords = None

        if abs(netCharge) > 0.01:
            print('Warning: The molecule has a net charge (%s). If decoupled for free energy calculations '
                  'corrections must be applied. See e.g.\n'
                  'Kastenholz and Hunenberger. J. Chem. Phys. 2006, 124, 124106-27' % netCharge)

        if charge_method:
            charge_method = charge_method.lower()
            print('charge_method: ', charge_method)
            try:
                if charge_method in standardChargeMethods:
                    generateCharges(inMolecule, charge_method, netCharge,
                                    multiplier = charge_multiplier, verbose = verbose)
                elif 'gaussian' in charge_method:
                    GaussianChargePlugin().generateCharges(inMolecule, charge_method, netCharge,
                                                         multiplier = charge_multiplier,
                                                         verbose = verbose)
                #elif 'gamess' in charge_method:
                #    GamessChargePlugin().generateCharges(inMolecule, charge_method, netCharge,
                #                                         multiplier = charge_multiplier,
                #                                         verbose = verbose)
                #else:
                #    for pl in chargePlugins:
                #        if charge_method in pl.alternativeChargeMethods:
                #            pl.generateCharges(inMolecule, charge_method, netCharge,
                #                               multiplier = charge_multiplier,
                #                               verbose = verbose)
                #            break
            except Exception:
                print('Cannot generate charges:')
                traceback.print_exc()

        outputFile = os.path.splitext(inMolecule)[0]
        if verbose:
            print('outputFile: ', outputFile)

        # The opls directory is automotically created when generating a GAFF topology.
        # If there was no opls directory before and no OPLS topology is generated make sure
        # that this directory is removed at a later stage.
        if os.path.exists(outputFile + '_opls'):
            oplsDidExist = True
        else:
            oplsDidExist = False

        # The gaff directory is automotically created when generating a OPLS topology.
        # If there was no gaff directory before and no GAFF topology is generated make sure
        # that this directory is removed at a later stage.
        if os.path.exists(outputFile + '_gaff') or os.path.exists(outputFile + '_gaff2'):
            gaffDidExist = True
        else:
            gaffDidExist = False

        for converter in converters:
            if not converter.forceFieldName in forcefields:
                continue

            if ffprotein is None:
                if 'gaff' in converter.forceFieldName:
                    ffprotein = 'amber99sb-ildn'
                elif converter.forceFieldName == 'cgenff':
                    ffprotein = 'charmm27'
            else:
                ffprotein = ffprotein

            if verbose:
                print(ffprotein, 'ffprotein: ', ffprotein)

            if verbose:
                print('Generating %s parameters' % converter.forceFieldName)
            try:
                converter.generate(inputFile = inMolecule, output = outputFile,
                                   keepMol2Charges = retain_charges or charge_method,
                                   netCharge = netCharge,
                                   verbose = verbose)
            except Exception:
                print('Error running generator for ' + converter.forceFieldName)
                traceback.print_exc()

            try:
                converter.convert2Gromacs(outputFile, verbose)
            except Exception:
                print('Error converting %s topology to GROMACS format.' % converter.forceFieldName)
                traceback.print_exc()
                continue

            if calibration:
                try:
                    calibrateVdW(outputFile, calibration, converter.forceFieldName,
                                 water, charge_method, verbose = verbose)
                except Exception:
                    print('Error calibrating van der Waals parameters.')
                    traceback.print_exc()

            ffDir = outputFile + '_%s' % converter.forceFieldName
            if verbose:
                print('ffDir: ', ffDir)

            if os.path.exists(ffDir):
                outputFileBaseName = os.path.basename(outputFile)
                if verbose:
                    print('Generating %s topology' % converter.forceFieldName)

                ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')
                #if os.path.exists(ffSpecificCoords):
                #    os.remove(ffSpecificCoords)
                #ffSpecificCoords = None

                if proteinCoords:
                    print('proteinCoords: ', proteinCoords)
                    if os.path.splitext(proteinCoords)[1].lower() != '.gro':
                        try:
                            ffProteinTopology, ffProteinCoords, ffProteinPosre = converter.coordsToTopology(outputFile, proteinCoords, ffprotein, verbose)
                        except Exception:
                            print('Error generating protein topology.')
                            traceback.print_exc()

                    else:
                        ffProteinCoords = proteinCoords
                        ffProteinTopology = proteinTopology
                        ffProteinPosre = None
                        copyItp(ffProteinTopology, ffDir, verbose = verbose)
                        modproteinItp(ffProteinTopology, ffDir, verbose = verbose)

                    #ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')

                    if verbose:
                        print('Merging coordinates:', ligandCoords, 'with', ffProteinCoords, 'to', ffSpecificCoords)
                    try:
                        mergeCoordinateFiles(ligandCoords, ffProteinCoords, ffSpecificCoords)
                    except Exception:
                        print('Cannot merge coordinate files:')
                        traceback.print_exc()
                        continue

                try:
                    ligandTop = converter.genTop(outputFile, ffprotein, water,
                                                 verbose = verbose)
                except Exception:
                    print('Cannot generate topology for %s:' % converter.forceFieldName)
                    traceback.print_exc()
                    continue

                try:
                    converter.fixAssignment(outputFile, retain_charges or charge_method,
                                            netCharge, spreadCharges = False, verbose = verbose)
                except Exception:
                    print('Cannot fix Assignments:')
                    traceback.print_exc()

                #if conf.virtualhydrogens:
                    #hydrogens2VirtualSites(ffDir, outputFileBaseName, verbose = verbose)

                #generateLinearVirtualSites(ffDir, outputFileBaseName, verbose = verbose)

                #ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')
                if not os.path.exists(ffSpecificCoords):
                    ffSpecificCoords = None

                if proteinCoords:
                    outTop = os.path.join(ffDir, outputFileBaseName + '.top')
                    if verbose:
                        print('OutTop:', outTop)

                    try:
                        if not ffProteinTopology or not os.path.exists(ffProteinTopology):
                            ffProteinTopology = os.path.splitext(ffProteinCoords)[0] + '.top'
                            if not os.path.exists(ffProteinTopology):
                                ffProteinTopology = None

                        if ffProteinTopology:
                            splitTopologyToItp(ffProteinTopology, verbose = verbose)
                            mergeTopologyFiles(ligandTop, ffProteinTopology, outTop, 
                                               verbose = verbose)

                            # Only calculate the total charge of the system once.
                            if totCharge == None:
                                totCharge = netCharge + getChargeOfTopology(ffProteinTopology)

                    except Exception:
                        print('Cannot merge files:')
                        traceback.print_exc()
                        continue

                if water:
                    if box_buffer > 0:
                        solvateSystem(ffSpecificCoords, ffDir, outputFileBaseName, water,
                        box_type, box_buffer, verbose = verbose)

                    if totCharge == None:
                        totCharge = netCharge

                    if totCharge or conc > 0.0:
                        neutraliseSystem(ffDir, outputFileBaseName, conc, 
                                         pname, nname, verbose = verbose)

                try:
                    makeIndexRun(ffDir, outputFileBaseName, verbose = verbose)
                except Exception as e:
                    print('Cannot make index file:')
                    traceback.print_exc()

                try:
                    convertGmx2Amb(ffDir, outputFileBaseName, verbose = verbose)
                except Exception as e:
                    print('Cannot make Amber topology file:')
                    traceback.print_exc()

            # If we are not generating an opls topology remove the automatically generated
            # opls directory.
            if 'gaff' in converter.forceFieldName and 'opls' not in chosenForcefields and not oplsDidExist:
                try:
                    shutil.rmtree(outputFile + '_opls')
                except Exception:
                    print('Error removing OPLS directory')
                    traceback.print_exc()

            # If we are not generating a gaff topology remove the automatically generated
            # gaff directory.
            elif converter.forceFieldName == 'opls' and 'gaff' not in chosenForcefields and not gaffDidExist:
                try:
                    shutil.rmtree(outputFile + '_gaff')
                except Exception:
                    print('Error removing GAFF directory')
                    traceback.print_exc()
            elif converter.forceFieldName == 'opls' and 'gaff2' not in chosenForcefields and not gaffDidExist:
                try:
                    shutil.rmtree(outputFile + '_gaff2')
                except Exception:
                    print('Error removing GAFF2 directory')
                    traceback.print_exc()
           
            try:
                converter.finalClean(outputFile)
            except Exception:
                print('Cannot clean up for %s' % converter.forceFieldName)
                traceback.print_exc()

    print('Finished processing molecule')

def stage3_main(conf):
    ligand = conf['ligand']
    smiles = conf['smiles']
    output = conf['output']
    ffligand = conf['ffligand']
    ffprotein = conf['ffprotein']
    calibration = conf['calibration']
    keep_ligand_name = conf['keep_ligand_name']
    ph = conf['ph']
    retain_charges = conf['retain_charges']
    charge_method = conf['charge_method']
    charge_method = conf['charge_method']
    charge_multiplier = conf['charge_multiplier']
    mergecoordinates = conf['mergecoordinates']
    mergetopology = conf['mergetopology']
    box_type = conf['box_type']
    box_buffer = conf['box_buffer']
    water = conf['water']
    conc = conf['conc']
    pname = conf['pname']
    nname = conf['nname']
    verbose = conf['verbose']

    stage3_run(ligand, smiles, output, ffligand, ffprotein, calibration, 
               keep_ligand_name, ph, retain_charges, charge_method, charge_multiplier,
               mergecoordinates, mergetopology, box_type, box_buffer,
               water, conc, pname, nname, verbose)

def main():
    args = get_parser()
    if args.verbose: print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    stage3_main(conf)

if __name__ == "__main__":
    main()
