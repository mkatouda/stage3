#!/usr/bin/env python

"""
    A wrapper for generating GROMACS topology files.
    Can use ACPYPE (calling AnteChamber) and MATCH to generate
    topologies for GAFF, OPLS and CGENFF and compatible charges.
    OPLS atom types and charges, from ACPYPE, are improved.
    The range of parameter assignment tools might be extended in the future.

                      VERSION 0.9.9

    Written by Magnus Lundborg
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

    The following environment variables must be set

    AMBERHOME
    MATCH
    PerlChemistry

    It might be good to set the GMXLIB environment variable to point to the ../top/ directory of
    the GROMACS installation.


    Requirements: Python version >= 3.0
                  Open Babel version >= 2.3
                  GROMACS (pdb2gmx, make_ndx, genrestr, editconf, solvate (genbox), grompp and genion)
                  ACPYPE (for GAFF and OPLS)
                  AnteChamber (for GAFF and OPLS)
                  MATCH (for CGENFF)
                  GAMESS (for B3LYP/PCM partial charges)
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
import time
import traceback
import inspect

from glob import glob
from itertools import product, chain
from operator import itemgetter
from ForceFieldPlugin import ForceFieldPlugin
from ChargePlugin import ChargePlugin

from util import babelConvert, renameAtoms, mol2RenameToLig, getNetChargeOfMol2, makeRestraintsRun, getChargeOfTopology
from util import generateCharges, calibrateVdW, mergeCoordinateFiles, copyItp, modproteinItp, splitTopologyToItp, mergeTopologyFiles
from util import hydrogens2VirtualSites, generateLinearVirtualSites, solvateSystem, neutraliseSystem, makeIndexRun

def _loadPlugins(path, name, baseclass):

    loadedPlugins = []

    plugins = glob(os.path.join(path, name))
    print('plugins: ', plugins)
    for i in range(len(plugins)):
        plugins[i] = 'plugins.' + os.path.basename(plugins[i])[:-3]

    for plugin in plugins:

        try:
            p = __import__(plugin, fromlist = ['nonsense']) # fromlist must be a non-empty list
        except (ImportError,NotImplementedError):
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

    if baseclass == 'ForceFieldPlugin':
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order, c.forceFieldName))
    elif baseclass == 'ChargePlugin':
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order, c.programName))
    else:
        loadedPlugins = sorted(loadedPlugins, key=lambda c: (c.order))

    return loadedPlugins

if __name__ == '__main__':

    progDir = os.path.dirname(__file__)

    converters = _loadPlugins(os.path.join(progDir, 'plugins'), '*ForceFieldPlugin.py', ForceFieldPlugin)
    #print('converters: ', converters)

    forcefieldsString = ','.join(str(conv.forceFieldName) for conv in converters)
    print('forcefieldsString: ', forcefieldsString)


    standardChargeMethods = ['am1bcc', 'am1bcc-pol', 'mmff94', 'eem', 'qeq', 'qtpie']
    chargeMethodsHelp = ['am1bcc: AM1 with bond charge correction (antechamber)',
                         'am1bcc-pol: STaGE\'s own more polarized bond charge correction (antechamber)',
                         'mmff94: MMFF94 (Open Babel)',
                         'eem: electronegativity equalization method (Open Babel)',
                         'qeq: Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991) (Open Babel)',
                         'qtpie: Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007) (Open Babel)']

    chargePlugins = _loadPlugins(os.path.join(progDir, 'plugins'), '*ChargePlugin.py', ChargePlugin)
    chargePluginNameList = []
    chargePluginsHelp = []
    for pl in chargePlugins:
        for method, description in pl.alternativeChargeMethods.items():
            chargePluginNameList.append(method)
            chargeMethodsHelp.append('%s: %s (%s)' % (method, description, pl.programName))

    chargeMethodsList = standardChargeMethods + chargePluginNameList
    print('chargeMethodsList: ', chargeMethodsList)

    parser = argparse.ArgumentParser(description='STaGE is a tool for generating GROMACS topologies '
                                     'of small molecules for different force fields (currently GAFF, OPLS-AA '
                                     'and CGenFF). It uses many external programs to perform its tasks, '
                                     'so make sure you read the documentation to properly cite the programs if '
                                     'you use STaGE in a publication and also cite: '
                                     'Lundborg M., Lindahl E. Automatic GROMACS Topology Generation and Comparisons '
                                     'of Force Fields for Solvation Free Energy Calculations. J. Phys. Chem. B. 2014, '
                                     'DOI: 10.1021/jp505332p')
    parser.add_argument('-i', '--inputfile',
                        help = 'Input file name.')
    parser.add_argument('-s', '--smiles',
                        help = 'Use the specified smiles string as input '
                        'instead of an input file (must be inside quotes).')
    parser.add_argument('-o', '--outputfile',
                        help = 'Name of the output files (file extensions '
                        'will be appended).', required=True)
    parser.add_argument('-k', '--keep_ligand_name',
                        action = 'store_true',
                        help = 'Do not rename the ligand in the output files. '
                               'When doing e.g. solvation or binding free energy '
                               'it is convenient to always call the ligand the '
                               'same thing - in this case "LIG". If this option '
                               'is set the ligand name will not be changed to "LIG". '
                               'If you need to assign parameters to e.g. co-factors '
                               'it is good to keep their names to tell them apart '
                               'from ligands.')
    parser.add_argument('-b', '--box_type',
                        default = 'dodecahedron',
                        help = 'Buffer from the solute to the edge of the '
                               'dodecahedron shaped solvent box. Set to 0 '
                               'to disable solvation (and ionisation).\n'
                               'Default: dodecahedron')
    parser.add_argument('-d', '--box_buffer',
                        type = float,
                        default = 1.0,
                        help = 'Buffer from the solute to the edge of the '
                               'solvent box. Set to 0 '
                               'to disable solvation (and ionisation).\n'
                               'Default: 1.0')
    parser.add_argument('-w', '--water',
                        help = 'Solvent model to use in topology files. If not '
                               'specified the solvent will not be specified in '
                               'the topology. Suggested water models are: '
                               '"opc", "spce", "tip4pew", "spc" or "tip3p".')
    parser.add_argument('-p', '--ph',
                        type = float,
                        help = 'Protonate the molecule according to this pH (float). '
                               'This does not always give correct results. It is safer '
                               'to provide correctly protonated input files.')
    parser.add_argument('-r', '--retain_charges',
                        action = 'store_true',
                        help = 'Keep the mol2 charges.')
    parser.add_argument('-q', '--charge_method',
                        choices = standardChargeMethods + chargePluginNameList,
                        help = 'Use the specified charge method for all force fields. ' + ', '.join(chargeMethodsHelp))
    parser.add_argument('-f', '--charge_multiplier',
                        type = float,
                        default = 1.0,
                        help = 'Multiply partial charges with this factor. Can only be used '
                               'in combination with --charge_method.')
    parser.add_argument('-x', '--calibration',
                        help = 'Modify van der Waals parameters according to specified '
                               'calibration file.')
    parser.add_argument('-c', '--mergecoordinates',
                        help = 'Merge the created coordinates file (.gro) with an '
                               'already existing coordinate file (.pdb or .gro), '
                               'e.g. for combining '
                               'ligand coordinates with protein coordinates. The generated topology '
                               'will contain both the ligand and the protein. If a .gro file of the '
                               'protein is provided and there exists a corresponding .top file that '
                               'toplogy file will be used for the protein, otherwise a new topology '
                               'file is generated.')
    parser.add_argument('-t', '--mergetopology',
                        help = 'Merge the created topology file (.top) with an '
                               'already existing topology file. '
                               'Must be used in combination with --mergecoordinates '
                               'with a .gro file of the protein.')
    #parser.add_argument('-y', '--virtualhydrogens',
                        #action = 'store_true',
                        #help = 'Turn hydrogens into virtual interaction sites to allow longer '
                               #'timesteps (experimental).')
    parser.add_argument('--forcefields',
                        default = forcefieldsString,
                        help = 'Force fields to generate parameters for, specified as a '
                               'comma-separated string without spaces. Default: %s' % forcefieldsString)
    parser.add_argument('--ffprotein',
                        default = None,
                        help = 'Force field of protein.')
    parser.add_argument('--pname',
                        default = 'NA',
                        help = 'Name of the positive counter ion in Solvent. Default: NA')
    parser.add_argument('--nname',
                        default = 'CL',
                        help = 'Name of the negative counter ion in Solvent. Default: CL')
    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        help = 'Verbose output.')

    args = parser.parse_args()

    count = 0

    if not args.inputfile and not args.smiles:
        parser.error('No molecule input specified. Either -i or -s '
        'must be provided.')

    if not args.outputfile:
        parser.error('No output name specified. -o must be specified. '
        'File extensions will be appended to the name')

    if args.inputfile and args.smiles:
        parser.error('--input and --smiles parameters cannot be used together.')

    if args.retain_charges and args.charge_method:
        parser.error('--retain_charges and --charge_method parameters cannot be '
                     'used together.')

    if args.outputfile == args.inputfile:
        parser.error('Output name and input file may not be identical.')

    if args.inputfile and not os.path.exists(args.inputfile):
        parser.error('The input file does not exists.')

    if args.mergecoordinates:
        proteinCoords = os.path.abspath(args.mergecoordinates)
        print('proteinCoords: ', proteinCoords)
        if os.path.splitext(proteinCoords)[1].lower() not in ('.gro', '.pdb'):
            parser.error('The coordinates file to merge with must have a '
            '.gro or .pdb extension')
    else:
        proteinCoords = None

    # chosenForcefields contains the actual force fields that were chosen by the user.
    # However, if generating OPLS parameters GAFF parameters must first be generated,
    # but these can safely be removed later on.
    chosenForcefields = args.forcefields.split(',')
    forcefields = list(chosenForcefields)

    for forcefield in forcefields:
        if forcefield not in forcefieldsString:
            parser.error('Force field %s is not available.' % forcefield)

    if 'opls' in forcefields and not 'gaff' in forcefields:
        forcefields.append('gaff')

    if args.mergetopology:
        if not proteinCoords or os.path.splitext(proteinCoords)[1].lower() != '.gro':
            parser.error('A .gro coordinate file must be specified with '
                         '--mergecoordinates if --mergetopology is used.')
        else:
            proteinTopology = os.path.abspath(args.mergetopology)
    else:
        proteinTopology = None

    if args.smiles:
        iBase, iFormat = None, 'smiles'
        outputFile = os.path.abspath(args.outputfile)
    else:
        iBase, iFormat = os.path.splitext(args.inputfile)
        iFormat = iFormat.lower()
        outputFile = os.path.abspath(args.outputfile)
        #outputFile = iBase

    if args.inputfile:
        inputFile = os.path.abspath(args.inputfile)
    else:
        inputFile = None

    if args.verbose:
        print('iBase: ', iBase, 'iFormat:', iFormat)
        print('inputFile:', inputFile)
        print('outputFile:', outputFile)

    if ',' in outputFile or (inputFile and ',' in inputFile):
        parser.error('The output and/or input files must not contain a comma (",")')

    if iFormat != '.mol2' or args.ph:
        if args.smiles:
            inMolecules = babelConvert(outputFile = outputFile,
                                       smiles = args.smiles,
                                       pH = args.ph,
                                       verbose = args.verbose)
        else:
            inMolecules = babelConvert(inputFile = inputFile, outputFile =
                                       outputFile, pH = args.ph,
                                       verbose = args.verbose)
    else:
        if inputFile != outputFile+'.mol2':
            shutil.copy(inputFile, outputFile+'.mol2')
        inMolecules = [outputFile+'.mol2']

    if args.verbose:
        print('inMolecules:', inMolecules)

    if not inMolecules:
        print('Could not convert molecule to mol2')
        sys.exit(1)

    if args.calibration and not os.path.exists(args.calibration):
        print('Calibration file not found. Continuing without calibrating vdW parameters.')
        args.calibration = None

    for inMolecule in inMolecules:

        if args.verbose:
            print('inMolecule: ', inMolecule)

        # Warning! This routine has some bugs for making wrong input mol2 file for acpype run.
        if 'cgenff' in forcefields:
            renameAtoms(inMolecule)

        if not args.keep_ligand_name:
            mol2RenameToLig(inMolecule)

        ligandCoords = os.path.splitext(inMolecule)[0] + '.gro'
        if args.verbose:
            print('ligandCoords: ', ligandCoords)

        if not os.path.exists(ligandCoords):
            try:
                babelConvert(inMolecule, ligandCoords, verbose = args.verbose)
            except Exception:
                print('Cannot create .gro file:')
                traceback.print_exc()

        count += 1

        netCharge = getNetChargeOfMol2(inMolecule)

        #MK
        if args.verbose:
            print('netCharge: ', netCharge)

        try:
            makeRestraintsRun(ligandCoords, verbose = args.verbose)
        except Exception:
            print('Cannot make restraints:')
            traceback.print_exc()

        totCharge = None
        ffSpecificCoords = None

        if abs(netCharge) > 0.01:
            print('Warning: The molecule has a net charge (%s). If decoupled for free energy calculations '
                  'corrections must be applied. See e.g.\n'
                  'Kastenholz and Hunenberger. J. Chem. Phys. 2006, 124, 124106-27' % netCharge)

        if args.charge_method:
            print('charge_method: ', args.charge_method)
            """
            try:
                if args.charge_method in standardChargeMethods:
                    generateCharges(inMolecule, args.charge_method, netCharge,
                                    multiplier = args.charge_multiplier, verbose = args.verbose)
                else:
                    for pl in chargePlugins:
                        if args.charge_method in pl.alternativeChargeMethods:
                            pl.generateCharges(inMolecule, args.charge_method, netCharge,
                                               multiplier = args.charge_multiplier,
                                               verbose = args.verbose)
                            break
            except Exception:
                print('Cannot generate charges:')
                traceback.print_exc()
            """

        outputFile = os.path.splitext(inMolecule)[0]
        if args.verbose:
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
        if os.path.exists(outputFile + '_gaff'):
            gaffDidExist = True
        else:
            gaffDidExist = False
        for converter in converters:
            if not converter.forceFieldName in forcefields:
                continue

            if args.ffprotein is None:
                if converter.forceFieldName == 'gaff':
                    ffprotein = 'amber99sb-ildn'
                    #ffprotein = 'amber14sb'
                elif converter.forceFieldName == 'cgenff':
                    ffprotein = 'charmm27'
                    #ffprotein = 'charmm36-jul2021'
            else:
                ffprotein = args.ffprotein

            if args.verbose:
                print(args.ffprotein, 'ffprotein: ', ffprotein)

            if args.verbose:
                print('Generating %s parameters' % converter.forceFieldName)
            try:
                converter.generate(inputFile = inMolecule, output = outputFile,
                                   keepMol2Charges = args.retain_charges or args.charge_method,
                                   netCharge = netCharge,
                                   verbose = args.verbose)
            except Exception:
                print('Error running generator for ' + converter.forceFieldName)
                traceback.print_exc()

            try:
                converter.convert2Gromacs(outputFile, args.verbose)
            except Exception:
                print('Error converting %s topology to GROMACS format.' % converter.forceFieldName)
                traceback.print_exc()
                continue

            if args.calibration:
                try:
                    calibrateVdW(outputFile, args.calibration, converter.forceFieldName,
                                 args.water, args.charge_method, verbose = args.verbose)
                except Exception:
                    print('Error calibrating van der Waals parameters.')
                    traceback.print_exc()

            ffDir = outputFile + '_%s' % converter.forceFieldName
            if args.verbose:
                print('ffDir: ', ffDir)

            if os.path.exists(ffDir):
                outputFileBaseName = os.path.basename(outputFile)
                if args.verbose:
                    print('Generating %s topology' % converter.forceFieldName)

                ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')
                #if os.path.exists(ffSpecificCoords):
                #    os.remove(ffSpecificCoords)
                #ffSpecificCoords = None

                if proteinCoords:
                    print('proteinCoords: ', proteinCoords)
                    if os.path.splitext(proteinCoords)[1].lower() != '.gro':
                        try:
                            ffProteinTopology, ffProteinCoords, ffProteinPosre = converter.coordsToTopology(outputFile, proteinCoords, ffprotein, args.verbose)
                        except Exception:
                            print('Error generating protein topology.')
                            traceback.print_exc()

                    else:
                        ffProteinCoords = proteinCoords
                        ffProteinTopology = proteinTopology
                        ffProteinPosre = None
                        copyItp(ffProteinTopology, ffDir, verbose = args.verbose)
                        modproteinItp(ffProteinTopology, ffDir, verbose = args.verbose)

                    #ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')

                    if args.verbose:
                        print('Merging coordinates: ', ligandCoords, 'with', ffProteinCoords, 'to', ffSpecificCoords)
                    try:
                        mergeCoordinateFiles(ligandCoords, ffProteinCoords, ffSpecificCoords)
                    except Exception:
                        print('Cannot merge coordinate files:')
                        traceback.print_exc()
                        continue

                try:
                    ligandTop = converter.genTop(outputFile, ffprotein, args.water,
                                                 verbose = args.verbose)
                except Exception:
                    print('Cannot generate topology for %s:' % converter.forceFieldName)
                    traceback.print_exc()
                    continue

                try:
                    converter.fixAssignment(outputFile, args.retain_charges or args.charge_method,
                                            netCharge, spreadCharges = False, verbose = args.verbose)
                except Exception:
                    print('Cannot fix Assignments:')
                    traceback.print_exc()

                #if args.virtualhydrogens:
                    #hydrogens2VirtualSites(ffDir, outputFileBaseName, verbose = args.verbose)

                #generateLinearVirtualSites(ffDir, outputFileBaseName, verbose = args.verbose)

                #ffSpecificCoords = os.path.join(ffDir, outputFileBaseName + '.gro')
                if not os.path.exists(ffSpecificCoords):
                    ffSpecificCoords = None

                if proteinCoords:
                    outTop = os.path.join(ffDir, outputFileBaseName + '.top')
                    if args.verbose:
                        print('OutTop:', outTop)

                    try:
                        if not ffProteinTopology or not os.path.exists(ffProteinTopology):
                            ffProteinTopology = os.path.splitext(ffProteinCoords)[0] + '.top'
                            if not os.path.exists(ffProteinTopology):
                                ffProteinTopology = None

                        if ffProteinTopology:
                            splitTopologyToItp(ffProteinTopology, verbose = args.verbose)
                            mergeTopologyFiles(ligandTop, ffProteinTopology, outTop, 
                                               verbose = args.verbose)

                            # Only calculate the total charge of the system once.
                            if totCharge == None:
                                totCharge = netCharge + getChargeOfTopology(ffProteinTopology)

                    except Exception:
                        print('Cannot merge files:')
                        traceback.print_exc()
                        continue

                if args.water:
                    if args.box_buffer > 0:
                        solvateSystem(ffSpecificCoords, ffDir, outputFileBaseName, args.water,
                        args.box_type, args.box_buffer, verbose = args.verbose)

                    if totCharge == None:
                        totCharge = netCharge

                    if totCharge:
                        neutraliseSystem(ffDir, outputFileBaseName, totCharge, 
                                         args.pname, args.nname, verbose = args.verbose)

                try:
                    makeIndexRun(ffDir, outputFileBaseName, verbose = args.verbose)
                except Exception as e:
                    print('Cannot make index file:')
                    traceback.print_exc()

            # If we are not generating an opls topology remove the automatically generated
            # opls directory.
            if converter.forceFieldName == 'gaff' and 'opls' not in chosenForcefields and not oplsDidExist:
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

            try:
                converter.finalClean(outputFile)
            except Exception:
                print('Cannot clean up for %s' % converter.forceFieldName)
                traceback.print_exc()

    print('Finished processing molecule')

    sys.exit(0)
