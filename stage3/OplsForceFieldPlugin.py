"""
    This file is part of STaGE - a wrapper for generating
    GROMACS topology files.

    Written by Magnus Lundborg
    Copyright (c) 2013-2015, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

import copy
import os
import subprocess
import numpy
from datetime import datetime
from itertools import product
from glob import glob

from .ForceFieldPlugin import ForceFieldPlugin
from .molScripts.molecule import MolSystem
from .molScripts.forcefield import Forcefield
from .util import getCommandOutput, isInPath, findForceFieldsDir, is_exe

progDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
obabel_atom_typer_loc = os.path.join(progDir, 'external', 'obabel_atom_typer', 'bin', 'obabel_atom_typer')

if not is_exe(obabel_atom_typer_loc):
    raise ImportError('obabel_atom_typer is not available. Cannot set OPLS atom types')

if not isInPath('gmx'):
    if isInPath('gmx_seq'):
        gmxSuffix = '_seq'
    elif isInPath('gmx_mpi'):
        gmxSuffix = '_mpi'
else:
    gmxSuffix = ''

def _pruneHydrogenAlts(heavyAtomAlts, hydrogenAlts, maxHydrogenDiff = 9):

    for i in reversed(range(len(hydrogenAlts))):
        alt = hydrogenAlts[i]
        hydrogenOplsNr = int(alt.name.split('_')[-1])
        for heavyAtomAlt in heavyAtomAlts:
            heavyAtomOplsNr = int(heavyAtomAlt.name.split('_')[-1])
            if abs(hydrogenOplsNr - heavyAtomOplsNr) <= maxHydrogenDiff:
                break
        else:
            del(hydrogenAlts[i])

    return hydrogenAlts

def _balanceCharges(origAtoms, alts, netCharge, equivalentAtomsCnt):
    """ This function tests all combinations of alternative atoms to find the
    most reasonable partial charges. It can be slow for very big molecules. """

    origNetCharge = sum([a.charge for a in origAtoms])

    minDiff = abs(origNetCharge - netCharge)
    bestAtoms = list(origAtoms)

    #print 'Original charge: %s, origDiff: %s' % (origNetCharge, minDiff)

    # Go through the list and find the atom combination with net charge closest
    # to the requested net charge.
    for alt in product(*alts):

        charges = numpy.array([a.charge if a != () else 0 for a in alt])

        charge = numpy.sum(charges * equivalentAtomsCnt)

        diff = abs(charge - netCharge)

        #print charge, diff

        if diff < minDiff:
            minDiff = diff
            bestAtoms = list(alt)
            #print 'Found new best:', [a.name if a != () else "" for a in bestAtoms], diff
            if diff <= 0.00001:
                break

    return (bestAtoms, minDiff)

class OplsForceFieldPlugin(ForceFieldPlugin):

    def __init__(self):

        self.forceFieldName = "opls"
        self.order = 9

    def _fixCoordsFileResidues(self, coordsFile, verbose = False):
        """ Replace residue and atom names to match the force field nomenclature
        Returns the name of the new coordinate file """

        replace = {'CAL': 'CA', 'CLA': 'CL', 'SOD': 'NA', 'POT': 'K', 'CES': 'CS', 'HIE': 'HISE', 'HID': 'HISD'}

        with open(coordsFile) as f:
            lines = f.readlines()

        modified = False
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname = line[17:21].strip()
                new = replace.get(resname)
                if new:
                    if verbose:
                        print('Changing')
                        print(line)
                    atomname = line[12:16].strip()
                    modified = True
                    line = '%s%-4s%s' % (line[:17], new, line[21:])
                    new = replace.get(atomname)
                    if new:
                        line = '%s%-4s%s' % (line[:12], new, line[16:])
                    lines[i] = line

                    if verbose:
                        print('to')
                        print(line)

        if modified:
            newCoordsFile = os.path.basename(coordsFile)
            with open(newCoordsFile, 'w') as f:
                for line in lines:
                    f.write(line)
            newCoordsFile = os.path.abspath(newCoordsFile)

            return newCoordsFile

        return None

    def fixAssignment(self, output, retainCharges = False, netCharge = 0, spreadCharges = False, verbose = False):
        """ Set OPLS atom types based on ACPYPE OPLS atom typing and the specified alternative types.
        If retainCharges == False the charges will be set according to the atom type.
        If spreadCharges == True non-integer charges will be spread over all atoms. N.B. non-integer charges
        indicate that the atom types are not correct. """

        outputFileBaseName = os.path.basename(output)

        oplsFileName = os.path.join(output + '_opls', outputFileBaseName + '.itp')
        gaffFileName = os.path.join(output + '_gaff', outputFileBaseName + '.itp')

        if not os.path.exists(gaffFileName):
            print('A GAFF topology is required to make an OPLS topology.')
            return

        forcefieldsDir = findForceFieldsDir('oplsaa')

        forcefieldFileName = os.path.join(forcefieldsDir, 'oplsaa.ff', 'ffnonbonded.itp')

        molFile = output + '.mol2'

        if not os.path.exists(forcefieldFileName):
            print('Cannot find forcefield %s' % forcefieldFileName)
            return

        if verbose:
            print('Setting OPLS atom types.')

        ff = Forcefield()

        ff.readGromacsItpFile(forcefieldFileName)
        forcefieldFileName = os.path.join(forcefieldsDir, 'oplsaa.ff', 'ffbonded.itp')
        ff.readGromacsItpFile(forcefieldFileName)

        with open(gaffFileName) as oplsFile:
            lines = oplsFile.readlines()

        mode = None

        nAtoms = 0

        bonds = []
        equivalentAtoms = []
        equivalentAtomsCnt = []

        for line in lines:
            if not mode:
                if line.find('[ atoms ]') != -1 or line.find('[atoms]') != -1:
                    mode = 'atoms'
                elif line.find('[ bonds ]') != -1 or line.find('[bonds]') != -1:
                    mode = 'bonds'

            # Read the atoms of the molecule to find out which atom types are used.
            elif mode == 'atoms':
                if line[0] == '[':
                    if line.find('[ bonds ]') != -1 or line.find('[bonds]') != -1:
                        mode = 'bonds'
                    else:
                        mode = None
                    continue
                elif line[0] == ';':
                    continue

                parts = line.split()
                if len(parts) > 7:
                    nAtoms += 1
                    bonds.append([])
                    equivalentAtoms.append(None)
                    equivalentAtomsCnt.append(1)

            elif mode == 'bonds':
                if line[0] == '[':
                    break
                elif line[0] == ';':
                    continue

                parts = line.split()
                if len(parts) > 5:
                    a1 = int(parts[0]) - 1
                    a2 = int(parts[1]) - 1
                    bonds[a1].append(a2)
                    bonds[a2].append(a1)

        alts = [[] for i in range(nAtoms)]

        with open(os.path.join(os.path.dirname(__file__), '..', '..', 'smarts_atom_types/opls.par')) as smartsFile:
            smartsExprs = smartsFile.readlines()
            for line in smartsExprs:
                if line[0] == '*':
                    continue

                parts=line.split('|')
                if len(parts) < 4:
                    continue

                for i, part in enumerate(parts):
                    parts[i] = part.strip()
                code = parts[2]
                smarts = parts[3]

                obabel_atom_typer_cmd = [obabel_atom_typer_loc, '%s' % smarts, molFile]
                #print smarts
                try:
                    result = getCommandOutput(obabel_atom_typer_cmd)
                except subprocess.CalledProcessError as e:
                    print('Failed running', ' '.join(obabel_atom_typer_cmd))
                    try:
                        print(result, '\n')
                    except Exception:
                        pass
                #print result

                matches = result.split()
                if matches and '====' in matches[0]:
                    #print result
                    continue

                for atomIndex in matches:
                    try:
                        atomIndex = int(atomIndex) - 1
                    except ValueError:
                        continue
                    atom = ff.findAtom(code)
                    #print atomIndex, atom.name
                    if atom and atom not in alts[atomIndex]:
                        alts[atomIndex].append(atom)

        atoms = []

        for i, alt in enumerate(alts):
            if not alt:
                print('Cannot set all OPLS atom types.')
                for i, val in enumerate(alts):
                    print(atoms[i].name, ': ', [a.name for a in val])
                raise Exception('Cancelling setting OPLS atom types.')

            atom = alt[0]
            atoms.append(copy.copy(atom))


        if retainCharges:
            chargesMolFile = output + '_bcc_gaff.mol2'
            if not os.path.exists(chargesMolFile):
                chargesMolFile = output + '.mol2'
            with open(chargesMolFile) as f:
                contents = f.readlines()
            molsys = MolSystem()
            molsys.readMol2(contents)

            # TODO: Check that we can rely on using the first molecule (probably does not work with multiple
            # molecules per file) and that the atom indexing is correct.
            if molsys.molecules:
                mol = molsys.molecules[0]
                for i, bccAtom in enumerate(mol.atoms):
                    atoms[i].charge = bccAtom.charge
        else:
            for i, atom in enumerate(atoms):
                if atom.elementNr != None and atom.elementNr != 1:
                    hydrogens = [a for a in bonds[i] if atoms[a].elementNr and atoms[a].elementNr == 1]
                    if hydrogens and equivalentAtomsCnt[hydrogens[0]] == 1:
                        for h in hydrogens[1:]:
                            equivalentAtoms[h] = hydrogens[0]
                            equivalentAtomsCnt[h] = 0
                            equivalentAtomsCnt[hydrogens[0]] += 1
                            alts[h] = [()]

        currNetCharge = sum([a.charge for a in atoms])

        if verbose:
            print('The charge of the molecule using the assigned OPLS atoms is: '
                  '%.3f. The net charge of molecule should be %d.' % (currNetCharge, netCharge))

        # If the current charge of the molecule does not agree with the requested
        # net charge try to fix the charges by changing the atom types.
        if abs(currNetCharge - netCharge) > 0.01:
            if not retainCharges:
                if verbose:
                    print('Trying to reduce the charge difference by changing atom types.')
                    print('Possible atom alternatives:')
                    for i, val in enumerate(alts):
                        if val != [()]:
                            print(atoms[i].name, ': ', [a.name for a in val])

                (newAtoms, minDiff) = _balanceCharges(atoms, alts, netCharge, equivalentAtomsCnt)

                if newAtoms != atoms:
                    atoms = [None] * len(atoms)
                    for i in range(len(atoms)):
                        if equivalentAtoms[i] == None:
                            atoms[i] = copy.copy(newAtoms[i])
                        else:
                            atoms[i] = copy.copy(newAtoms[equivalentAtoms[i]])
                    currNetCharge = sum([a.charge for a in atoms])

                if verbose:
                    print('Best alternative:')
                    for atom in atoms:
                        print(atom.name)

            if abs(currNetCharge - netCharge) > 0.01:
                if spreadCharges:
                    print('WARNING: Net OPLS charge is not an integer! '
                        'Current charge: %.3f (expected %d)' % (currNetCharge, netCharge))
                    print('The charge difference will be spread over all atoms.')

                    diff = netCharge - currNetCharge
                    diffPerAtom = diff / len(atoms)
                    print('Charge modifiction per atom: %f' % diffPerAtom)
                    for a in atoms:
                        a.charge += diffPerAtom

                else:
                    print('WARNING: The OPLS charges are not corrected. Charge is: %.3f (expected %d).' % \
                    (currNetCharge, netCharge))

        # Write the new data to the file.
        mode = None
        with open(oplsFileName, 'w') as oplsFile:
            prog = os.path.basename(__file__)
            time = datetime.today()

            oplsFile.write('; %s created by %s on %s (bonded data partially from AcPype)\n' % (os.path.basename(oplsFileName), prog, time))
            started = False
            qTot = 0
            for line in lines:
                if line[0] == '[':
                    started = True
                    mode = None
                elif line[0] == ';' and mode != 'atomtypes':
                    if started:
                        oplsFile.write(line)
                    continue

                if not mode:
                    if line.find('[ atomtypes ]') != -1 or line.find('[atomtypes]') != -1:
                        mode = 'atomtypes'
                        continue

                    elif line.find('[ atoms ]') != -1 or line.find('[atoms]') != -1:
                        mode = 'atoms'
                        cnt = 0

                    elif line.find('[ bonds ]') != -1 or line.find('[bonds]') != -1:
                        mode = 'bonds'

                    elif line.find('[ angles ]') != -1 or line.find('[angles]') != -1:
                        mode = 'angles'

                    elif line.find('[ dihedrals ]') != -1 or line.find('[dihedrals]') != -1:
                        mode = 'dihedrals'

                    if started:
                        oplsFile.write(line)

                    continue

                elif mode == 'atomtypes':
                    continue

                elif mode == 'atoms':
                    parts = line.split()
                    if len(parts) < 8:
                        oplsFile.write(line)
                        continue

                    parts[1] = atoms[cnt].name
                    parts[6] = '%.6f' % atoms[cnt].charge
                    parts[7] = '%.5f' % atoms[cnt].mass
                    qTot += atoms[cnt].charge
                    atomType = atoms[cnt].kind

                    line = '  %-4s %-8s %-3s %-5s %-5s %-4s %10s %10s     ; %f   %s' % (parts[0],
                        parts[1], parts[2], parts[3], parts[4], parts[5], parts[6],
                        parts[7], qTot, atomType)

                    line += '\n'

                    cnt += 1

                elif mode == 'bonds':
                    parts = line.split()
                    if len(parts) < 5:
                        oplsFile.write(line)
                        continue

                    atom1 = atoms[int(parts[0]) - 1]
                    atom1Type = atom1.kind
                    atom2 = atoms[int(parts[1]) - 1]
                    atom2Type = atom2.kind

                    bond = ff.findBond(atom1, atom2)
                    if not bond:
                        line = line[:-1]
                        line += '  !Warning using GAFF parameters!\n'

                        print('Warning using GAFF parameters for OPLS bond %s - %s!' % (atom1Type, atom2Type))
                        if line.count(';') > 1:
                            line = line.replace(';', '', 1)
                    else:
                        line = bond.getGromacsOutput(includeComment = False)
                        newParts = line.split()
                        newParts[0] = parts[0]
                        newParts[1] = parts[1]
                        line = '    ' + '   '.join(newParts) + ' ; %s - %s\n' % (atom1Type, atom2Type)

                elif mode == 'angles':
                    parts = line.split()
                    if len(parts) < 6:
                        oplsFile.write(line)
                        continue

                    atom1 = atoms[int(parts[0]) - 1]
                    atom1Type = atom1.kind
                    atom2 = atoms[int(parts[1]) - 1]
                    atom2Type = atom2.kind
                    atom3 = atoms[int(parts[2]) - 1]
                    atom3Type = atom3.kind

                    angle = ff.findAngle(atom1, atom2, atom3)
                    if not angle:
                        line = line[:-1]
                        line += '  !Warning using GAFF parameters!\n'

                        print('Warning using GAFF parameters for OPLS angle %s - %s - %s!' % (atom1Type, atom2Type, atom3Type))
                        if line.count(';') > 1:
                            line = line.replace(';', '', 1)
                    else:
                        ## 180 degree angles have too high force constants
                        #if abs(angle.angle - 180.00) < 0.00001 and angle.force > 350:
                            #angle.angle = 179.85
                            #angle.force = 350
                        line = angle.getGromacsOutput(includeComment = False)
                        newParts = line.split()
                        newParts[0] = parts[0]
                        newParts[1] = parts[1]
                        newParts[2] = parts[2]
                        line = '    ' + '   '.join(newParts) + ' ; %s - %s - %s\n' % (atom1Type, atom2Type, atom3Type)

                elif mode == 'dihedrals':
                    parts = line.split()
                    if len(parts) < 6:
                        oplsFile.write(line)
                        continue

                    atom1 = atoms[int(parts[0]) - 1]
                    atom1Type = atom1.kind
                    atom2 = atoms[int(parts[1]) - 1]
                    atom2Type = atom2.kind
                    atom3 = atoms[int(parts[2]) - 1]
                    atom3Type = atom3.kind
                    atom4 = atoms[int(parts[3]) - 1]
                    atom4Type = atom4.kind
                    func = int(parts[4])

                    if func == 1:
                        line = '; ' + line
                        oplsFile.write(line)
                        continue

                    dih = ff.findDihedral(atom1, atom2, atom3, atom4, func)
                    if not dih:
                        line = line[:-1]
                        line += '  !Warning using GAFF parameters!\n'

                        #if func == 1:
                            #print('Warning using GAFF parameters for OPLS improper %s - %s - %s - %s!' % (atom1Type, atom2Type, atom3Type, atom4Type))
                        #else:
                        print('Warning using GAFF parameters for OPLS dihedral %s - %s - %s - %s!' % (atom1Type, atom2Type, atom3Type, atom4Type))
                        if line.count(';') > 1:
                            line = line.replace(';', '', 1)
                    else:
                        line = dih.getGromacsOutput(includeComment = False)
                        newParts = line.split()
                        newParts[0] = parts[0]
                        newParts[1] = parts[1]
                        newParts[2] = parts[2]
                        newParts[3] = parts[3]
                        line = '    ' + '   '.join(newParts) + ' ; %s - %s - %s - %s\n' % (atom1Type, atom2Type, atom3Type, atom4Type)


                oplsFile.write(line)

    def genTop(self, output, forcefield = 'oplsaa', solvent = None, verbose = False):
        """ Generate a GROMACS topology file. If solvent is specified a line
        is added to include the parameters for that solvent. """

        outputFileBaseName = os.path.basename(output)
        outputDir = os.path.dirname(output)
        topologyDir = output + '_%s' % self.forceFieldName
        topologyFileName = os.path.join(topologyDir, outputFileBaseName + '.top')
        #ligRestrFile = os.path.join('..', 'posre_' + os.path.splitext(outputFileBaseName)[0] + '.itp')
        ligRestrFile = 'posre_' + os.path.splitext(outputFileBaseName)[0] + '.itp'

        shutil.copy(os.path.join(outputDir, ligRestrFile), os.path.join(topologyDir, ligRestrFile))

        forcefieldsDir = findForceFieldsDir(forcefield)
        forcefieldDir = forcefield + '.ff'

        if solvent and not os.path.exists(os.path.join(forcefieldsDir, forcefieldDir, solvent + '.itp')):
            if forcefieldsDir == os.path.join(progDir, 'forcefields'):
                print('Cannot find solvent (%s.itp) in %s' % (solvent, os.path.join(forcefieldsDir, forcefieldDir)))
                return

            forcefieldsDir = os.path.join(progDir, 'forcefields')

            if not os.path.exists(os.path.join(forcefieldsDir, forcefieldDir, solvent + '.itp')):
                print('Cannot find solvent (%s.itp) in %s' % (solvent, os.path.join(forcefieldsDir, forcefieldDir)))
                return

        time = datetime.today()

        name = outputFileBaseName

        itpFileName = os.path.join(topologyDir, outputFileBaseName + '.itp')
        with open(itpFileName) as f:
            lines = f.readlines()
        with open(itpFileName, 'w') as f:

            inMoleculeSection = False

            for line in lines:
                if inMoleculeSection:
                    if line.startswith('#'):
                        f.write(line)
                        continue
                    if inMoleculeSection:
                        if line.startswith('['):
                            inMoleculeSection = False
                            f.write(line)
                        continue

                if line.find('[ moleculetype ]') != -1 or line.find('[moleculetype]') != -1:
                    f.write(line)
                    f.write('; Name            nrexcl\n')
                    f.write(' %s 3\n\n' % name)

                    inMoleculeSection = True
                    continue

                f.write(line)
            #f.write('\n\n#ifdef POSRES || POSRES_LIG\n')
            f.write('\n\n#ifdef POSRES_LIG\n')
            f.write('#include "%s"\n' % ligRestrFile)
            f.write('#endif\n')


        forcefieldFile = os.path.join(forcefieldDir, 'forcefield.itp')
        if not os.path.isfile(os.path.join(forcefieldsDir, forcefieldFile)):
            print('Cannot find forcefield %s' % os.path.join(forcefieldsDir, forcefieldFile))
            return

        itpFiles = glob(os.path.join(topologyDir, outputFileBaseName + '*.itp'))
        sortedItpFiles = []
        for f in itpFiles:
            if f.find('bonded.itp') != -1:
                sortedItpFiles.insert(0, f)
            else:
                sortedItpFiles.append(f)

        with open(topologyFileName, 'w') as f:

            f.write('; Topology file generated by %s %s\n' % (os.path.join(os.path.basename(__file__), '..'), time))
            f.write('\n')
            f.write('#include "%s"\n' % forcefieldFile)

            for itpFile in sortedItpFiles:
                #itpDir = os.path.dirname(itpFile)
                #if itpDir == topologyDir:
                    #f.write('#include "%s"\n' % os.path.basename(itpFile))
                #elif itpDir == outputDir:
                    #f.write('#include "%s"\n' % os.path.join('..', os.path.basename(itpFile)))
                #else:
                f.write('#include "%s"\n' % os.path.basename(itpFile))

            if solvent:
                solventFile = os.path.join(forcefieldDir, solvent + '.itp')
                f.write('#include "%s"\n' % solventFile)
                f.write('#include "%s"\n' % os.path.join(forcefieldDir, 'ions.itp'))

            f.write('\n')
            f.write('[ system ]\n')
            f.write('%s\n' % name)
            f.write('\n')
            f.write('[ molecules ]\n')
            l = max(len('Compound') + 3, len(name) + 1)
            compStr = 'Compound'.ljust(l)
            f.write('; %s nmols\n' % compStr)
            f.write('%s 1\n' % name.ljust(l))

        return topologyFileName

    def coordsToTopology(self, output, coordsFile,  forcefield = 'oplsaa', verbose = False):
        """ Generate topology and .gro coordinates from a pdb file.
        The topology is generated using pdb2gmx.
        Return .top file and .gro file """

        outputFileBaseName = os.path.basename(output)
        outputDir = output + '_%s' % self.forceFieldName
        currDir = os.getcwd()

        coordsFilePath = os.path.dirname(coordsFile)
        coordsBaseName = os.path.basename(coordsFile)
        coordsBaseName_wo_ext = os.path.splitext(coordsBaseName)[0]

        topolFile = coordsBaseName_wo_ext + '.top'
        groFile = coordsBaseName_wo_ext + '.gro'
        restraintsFile = 'posre_' + coordsBaseName_wo_ext + '.itp'

        #if os.path.exists(topolFile) and os.path.exists(groFile) and os.path.exits(restraintsFile):
        #    return topolFile, groFile, restraintsFile

        os.chdir(outputDir)       

        pdb2gmxCommand = ['gmx'+gmxSuffix, 'pdb2gmx', '-f', coordsFile,
                        '-o', groFile,
                        '-p', topolFile,
                        '-i', restraintsFile,
                        '-ff', forcefield,
                        '-water', 'none',
                        '-ignh']

        if verbose:
            print(' '.join(pdb2gmxCommand))

        try:
            result = getCommandOutput(pdb2gmxCommand)
            if verbose:
                print(result)
        except subprocess.CalledProcessError as e:
            print('Failed running %s' % ' '.join(pdb2gmxCommand))
            try:
                print(result)
                print('\n')
            except Exception:
                pass
            print(e)
            print('Fixing residue names of %s for %s' % (coordsFile, self.forceFieldName))

            # Try renaming some residues and atoms to match the force field - see if it helps.
            newCoordsFile = self._fixCoordsFileResidues(coordsFile, verbose)
            if newCoordsFile:
                print('Trying again')
                pdb2gmxCommand = ['gmx'+gmxSuffix, 'pdb2gmx', '-f', newCoordsFile,
                                '-o', groFile,
                                '-p', topolFile,
                                '-i', restraintsFile,
                                '-ff', forcefield,
                                '-water', 'none',
                                '-ignh']

                if verbose:
                    print(' '.join(pdb2gmxCommand))

                try:
                    result = getCommandOutput(pdb2gmxCommand)
                    if verbose:
                        print(result)
                    print('pdb2gmx command successful after fixing residue names. Continuing as normal')
                except subprocess.CalledProcessError as e:
                    print('Failed running %s' % ' '.join(pdb2gmxCommand))
                    try:
                        print(result)
                        print('\n')
                    except Exception:
                        pass
                    print(e)

        topolFile = os.path.abspath(topolFile)
        groFile = os.path.abspath(groFile)
        restraintsFile = os.path.abspath(restraintsFile)
        itpFiles = glob('*.itp')

        if len(itpFiles) > 0:
            for itpFile in itpFiles:
                if not 'posre_' in itpFile:
                    with open(itpFile) as f:
                        lines = f.readlines()

                    with open(itpFile, 'w') as f:
                        for line in lines:
                            stripped = line.strip()
                            if stripped == '#ifdef POSRES':
                                #f.write('#ifdef POSRES || POSRES_PROT\n')
                                f.write('#ifdef POSRES_PROT\n')
                                continue
                            f.write(line)

        if os.path.exists(topolFile):
            with open(topolFile) as f:
                lines = f.readlines()

            with open(topolFile, 'w') as f:
                for line in lines:
                    stripped = line.strip()
                    if stripped == '#ifdef POSRES':
                        #f.write('#ifdef POSRES || POSRES_PROT\n')
                        f.write('#ifdef POSRES_PROT\n')
                        continue
                    f.write(line)
        else:
            topolFile = None
        if not os.path.exists(groFile):
            groFile = None
        if not os.path.exists(restraintsFile):
            restraintsFile = None

        return topolFile, groFile, restraintsFile

    def finalClean(self, output):

        gmxbk_files = glob('./#*', recursive=True)

        for gmxbk_file in gmxbk_files:
            if os.path.isfile(gmxbk_file):
                os.remove(gmxbk_file)
