"""
    This file is part of STaGE - a wrapper for generating
    GROMACS topology files.

    Written by Magnus Lundborg
    Copyright (c) 2013-2015, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

import os
import shutil
import subprocess
from datetime import datetime
from glob import glob

from ForceFieldPlugin import ForceFieldPlugin
from util import getCommandOutput, isInPath, findPath, findForceFieldsDir

if not isInPath('gmx'):
    if isInPath('gmx_seq'):
        gmxSuffix = '_seq'
    elif isInPath('gmx_mpi'):
        gmxSuffix = '_mpi'
else:
    gmxSuffix = ''

class GaffForceFieldPlugin(ForceFieldPlugin):

    def __init__(self):

        self.forceFieldName = "gaff"
        self.order = 1

    def _postGenerateClean(self, output, removeFiles = True):
        """ Clean up after an AcPype run. Files are copied to forcefield
        subdirectories. If removeFiles is False the original files will
        be retained - otherwise they are deleted. """

        outputFileBaseName = os.path.basename(output)

        gaffDir = output + '_%s' % self.forceFieldName
        oplsDir = output + '_opls'

        if removeFiles:
            command = shutil.move
        else:
            command = shutil.copy

        try:
            shutil.rmtree(output + '.acpype')
            if removeFiles:
                for ext in ['_AC.inpcrd', '_AC.lib', '_AC.prmtop', '_user_gaff.mol2', '_AC.frcmod']:
                    f = output + ext
                    if os.path.isfile(f):
                        os.remove(f)

            if not os.path.exists(gaffDir):
                os.mkdir(gaffDir)

            for ext in ['.itp', '.top']:
                f = output + '_GMX' + ext
                if os.path.isfile(f):
                    command(f, os.path.join(gaffDir, outputFileBaseName + ext))

            f = output + '_GMX' + '.gro'
            if os.path.isfile(f):
                command(f, output + '.gro')

            if not os.path.exists(oplsDir):
                os.mkdir(oplsDir)

            for ext in ['.itp', '.top']:
                f = output + '_GMX_OPLS' + ext

                if os.path.isfile(f):
                    command(f, os.path.join(oplsDir, outputFileBaseName + ext))

            # Fix the name in the .itp files - do not include the whole path of the
            # file.
            for fname in [os.path.join(gaffDir, outputFileBaseName + '.itp'),
                        os.path.join(oplsDir, outputFileBaseName + '.itp')]:
                with open(fname) as f:
                    lines = f.readlines()

                with open(fname, 'w') as f:
                    for line in lines:
                        l = line.replace(output + '_GMX', output)
                        l = l.replace(output, outputFileBaseName)
                        f.write(l)

        except shutil.Error:
            print('Cannot clean up')

        except os.EX_NOPERM:
            print('Insufficient permissions for cleaning up')

    def _fixCoordsFileResidues(self, outputDir, outputFileBaseName, coordsFile, verbose = False):
        """ Replace residue and atom names to match the force field nomenclature
        Returns the name of the new coordinate file """

        replace = {'CAL': 'CA', 'CLA': 'CL', 'SOD': 'NA', 'POT': 'K', 'CES': 'CS'}

        with open(coordsFile) as f:
            lines = f.readlines()

        modified = False
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname = line[17:20].strip()
                new = replace.get(resname)
                if new:
                    if verbose:
                        print('Changing')
                        print(line)
                    atomname = line[12:16].strip()
                    modified = True
                    line = '%s%3s%s' % (line[:17], new, line[20:])
                    new = replace.get(atomname)
                    if new:
                        line = '%s%-4s%s' % (line[:12], new, line[16:])
                        lines[i] = line

                    if verbose:
                        print('to')
                        print(line)

        if modified:
            coordsBaseName = os.path.basename(coordsFile)
            newCoordsFile = os.path.join(outputDir, outputFileBaseName + '_' + self.forceFieldName + '_' + coordsBaseName)
            with open(newCoordsFile, 'w') as f:
                for line in lines:
                    f.write(line)

            return newCoordsFile

        return None

    def generate(self, inputFile, output, keepMol2Charges = False, netCharge = None,
                 verbose = False):
        """ Run AcPype on inputFile. If keepMol2Charges is False partial charges
        will be assigned according to the forcefield. For OPLS the charges are
        first set to the same as for GAFF - they can be changed later on. The total charge
        of the molecule will be set to netCharge (if it is specified) or 0. If
        verbose is True more detailed output will be given. """

        oldFiles = glob(output + '_AC.*')
        if os.path.isfile(output + '_bcc_gaff.mol2'):
            oldFiles.append(output + '_bcc_gaff.mol2')
        for f in oldFiles:
            os.remove(f)

        acpypeCmd = ['acpype', '-i', inputFile, '-b', output, '-o', 'gmx']
        if keepMol2Charges:
            acpypeCmd += ['-c', 'user']
            netCharge = None

        if netCharge != None:
            acpypeCmd += ['-n', '%d' % netCharge]

        if verbose:
            print(' '.join(acpypeCmd))

        result = getCommandOutput(acpypeCmd)
        if verbose:
            print(result)

        self._postGenerateClean(output)


    def genTop(self, output, forcefield = 'amber99sb-ildn', solvent = None, verbose = False):
        """ Generate a GROMACS topology file. If solvent is specified a line
        is added to include the parameters for that solvent. """

        outputFileBaseName = os.path.basename(output)
        outputDir = os.path.dirname(output)
        topologyDir = output + '_%s' % self.forceFieldName
        topologyFileName = os.path.join(topologyDir, outputFileBaseName + '.top')
        progDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
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

    def coordsToTopology(self, output, coordsFile, forcefield = 'amber99sb-ildn',
                         verbose = False):
        """ Generate topology and .gro coordinates from a pdb file.
        The topology is generated using pdb2gmx.
        Return .top file and .gro file """

        outputFileBaseName = os.path.basename(output)
        outputDir = output + '_%s' % self.forceFieldName

        coordsFilePath = os.path.dirname(coordsFile)
        coordsBaseName = os.path.basename(coordsFile)
        coordsBaseName_wo_ext = os.path.splitext(coordsBaseName)[0]

        topolFile = coordsBaseName_wo_ext + '.top'
        groFile = coordsBaseName_wo_ext + '.gro'
        restraintsFile = 'posre_' + coordsBaseName_wo_ext + '.itp'

        #if os.path.exists(topolFile) and os.path.exists(groFile) and os.path.exits(restraintsFile):
        #    return topolFile, groFile, restraintsFile

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
            newCoordsFile = self._fixCoordsFileResidues(outputDir, outputFileBaseName, coordsFile, verbose)
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

        topolFile = os.path.join(outputDir, coordsBaseName_wo_ext + '.top')
        groFile = os.path.join(outputDir, coordsBaseName_wo_ext + '.gro')
        restraintsFile = os.path.join(outputDir, 'posre_' + coordsBaseName_wo_ext + '.itp')
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

        outputFileBaseName = os.path.basename(output)

        pkl_file = os.path.abspath(outputFileBaseName + '.pkl')
        if os.path.isfile(pkl_file):
            os.remove(pkl_file)

        gmxbk_files = glob('./#*', recursive=True)

        for gmxbk_file in gmxbk_files:
            if os.path.isfile(gmxbk_file):
                os.remove(gmxbk_file)
