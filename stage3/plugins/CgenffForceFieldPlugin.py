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
from util import getCommandOutput, getGroLigandName, setRtpLigandName, replaceRtfCharges, isInPath, findForceFieldsDir

if not isInPath('gmx'):
    if isInPath('gmx_seq'):
        gmxSuffix = '_seq'
    elif isInPath('gmx_mpi'):
        gmxSuffix = '_mpi'
else:
    gmxSuffix = ''

class CgenffForceFieldPlugin(ForceFieldPlugin):

    def __init__(self):

        self.forceFieldName = "cgenff"
        self.order = 9

    def _postGenerateClean(self, output, removeFiles = True):
        """ Clean up after running MATCH. If removeFiles == False the files will not
        be removed and instead of moving files they will be copied. """

        outputDir = os.path.abspath(os.path.dirname(output))
        outputFileBaseName = os.path.basename(output)
        currDir = os.getcwd()

        if removeFiles:
            command = shutil.move
        else:
            command = shutil.copy

        try:
            for ext in ['.rtf', '.prm']:
                f = os.path.abspath(outputFileBaseName + ext)
                if f != output + ext and os.path.isfile(f):
                    if os.path.isfile(output + ext):
                        os.remove(output + ext)
                    command(f, outputDir)

            f = os.path.abspath('top_' + outputFileBaseName + '.rtf')
            if os.path.isfile(f):
                if f != os.path.join(outputDir, f):
                    if os.path.isfile(os.path.join(outputDir, f)):
                        os.remove(os.path.join(outputDir, f))
                    command(f, outputDir)

        except shutil.Error:
            print('Cannot clean up')

        except os.EX_NOPERM:
            print('Insufficient permissions for cleaning up')

    def _postConvert2GromacsClean(self, output, removeFiles = True):
        """ Clean up after converting from Charmm to GROMACS """

        outputDir = os.path.abspath(os.path.dirname(output))
        outputFileBaseName = os.path.basename(output)
        forcefieldDir = os.path.join(output + '_%s' % self.forceFieldName, 'forcefield.ff')

        f = output + '.gro'
        if os.path.isfile(f):
            shutil.copy(f, os.path.join(output + '_%s' % self.forceFieldName, 
                                        outputFileBaseName + '.gro'))

        if removeFiles:
            command = shutil.move
        else:
            command = shutil.copy

        ligName = getGroLigandName(output + '.gro')
        if not os.path.exists(forcefieldDir):
            os.mkdir(forcefieldDir)
        for f in ['ffbonded.itp', 'ffnonbonded.itp', 'forcefield.doc',
                'forcefield.itp', 'aminoacids.rtp', 'atomtypes.atp']:
            fFull = os.path.join(output + '_%s' % self.forceFieldName, f)
            if os.path.isfile(fFull):
                if f == 'aminoacids.rtp':
                    setRtpLigandName(fFull, ligName)

                outputFile = os.path.join(output + '_%s' % self.forceFieldName, forcefieldDir, f)
                if os.path.isfile(outputFile):
                    os.remove(outputFile)
                command(fFull, outputFile)

        if removeFiles:
            f = os.path.join(outputDir, 'top_' + outputFileBaseName + '.rtf')
            if os.path.isfile(f):
                os.remove(f)
            f = os.path.join(outputDir, outputFileBaseName + '.pkl')
            if os.path.isfile(f):
                os.remove(f)
            f = os.path.join(outputDir, outputFileBaseName + '.prm')
            if os.path.isfile(f):
                os.remove(f)
            f = os.path.join(outputDir, outputFileBaseName + '.rtf')
            if os.path.isfile(f):
                os.remove(f)

    def _fixCoordsFileResidues(self, coordsFile, verbose = False):
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
            newCoordsFile = os.path.basename(coordsFile)
            with open(newCoordsFile, 'w') as f:
                for line in lines:
                    f.write(line)
            newCoordsFile = os.path.abspath(newCoordsFile)

            return newCoordsFile

        return None


    def generate(self, inputFile, output, keepMol2Charges = False, netCharge = None,
                 verbose = False):
        """ Run MATCH to generate cgenff parameters. """

        forcefield = 'top_all36_cgenff_new'
        #forcefield = 'top_all36_cgenff'
        matchCmd = ['MATCH.pl', '-Forcefield', forcefield, inputFile]

        if verbose:
            print(' '.join(matchCmd))

        #try:
        result = getCommandOutput(matchCmd)
        if verbose:
            print(result)

        if keepMol2Charges:
            rtfFile = os.path.splitext(inputFile)[0] + '.rtf'
            replaceRtfCharges(rtfFile, inputFile)

        self._postGenerateClean(output)


    def convert2Gromacs(self, output, verbose = False):
        """ Run a script (by David van der Spoel) to convert a charmm topology to
        GROMACS format """

        convCmd = [os.path.join(os.path.dirname(__file__), '..', 'charmm2gromacs-pvm.py')]
        convCmd += [output + '.rtf', output + '.prm', output + '_%s' % self.forceFieldName]

        if verbose:
            print(' '.join(convCmd))

        try:
            result = getCommandOutput(convCmd)
            if verbose:
                print(result)
        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(convCmd))
            try:
                print(result)
                print('\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')


        except shutil.Error as e:
            print('Cannot clean up:')
            print(e)

        except os.EX_NOPERM:
            print('Insufficient permissions for cleaning up')

        self._postConvert2GromacsClean(output)


    def genTop(self, output, forcefield = 'charmm27', solvent = None, verbose = False):
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

        currDir = os.getcwd()
        os.chdir(topologyDir)
        
        forcefieldArg = os.path.join('forcefield')
        pdb2gmxCommand = ['gmx'+gmxSuffix, 'pdb2gmx', '-f', output + '.gro', '-o', 'temp.gro',
                          '-p', topologyFileName,
                          #'-i', os.path.join(topologyDir, outputFileBaseName + '_restr.itp')
                          '-ff', forcefieldArg]

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
            print('Trying to continue anyhow.')

        os.chdir(currDir)
        try:
            os.remove(os.path.join(topologyDir, 'temp.gro'))
        except Exception as e:
            print(e)

        ffnonbondedFileName = os.path.join(topologyDir, 'forcefield.ff', 'ffnonbonded.itp')
        ffbondedFileName = os.path.join(topologyDir, 'forcefield.ff', 'ffbonded.itp')

        with open(ffnonbondedFileName) as ffnonbondedFile:
            ffnonbondedContents = ffnonbondedFile.readlines()
        with open(ffbondedFileName) as ffbondedFile:
            ffbondedContents = ffbondedFile.readlines()

        with open(topologyFileName) as topFile:
            lines = topFile.readlines()
        with open(topologyFileName, 'w') as topFile:
            headerWritten = False
            itpFileIncluded = False

            itpFileName = os.path.join(topologyDir, outputFileBaseName + '.itp')
            with open(itpFileName, 'w') as itpFile:

                forcefieldFile = os.path.join(forcefieldDir, 'forcefield.itp')
                if not os.path.isfile(os.path.join(forcefieldsDir, forcefieldFile)):
                    print('Cannot find forcefield %s.' % forcefieldFile)
                    return

                inMoleculeSection = 0
                for line in lines:
                    if inMoleculeSection > 1:
                        inMoleculeSection -= 1
                        continue
                    elif inMoleculeSection:
                        if line.startswith('#'):
                            #topFile.write(line)
                            continue
                        elif line.find('[ system ]') != -1:
                            inMoleculeSection = 0
                            topFile.write('\n')
                        else:
                            itpFile.write(line)
                            continue

                    # This is the end of the file
                    if line.find('[ system ]') != -1 or line.find('[system]') != -1:
                        topFile.write('[ system ]\n')
                        topFile.write('%s\n' % name)
                        topFile.write('\n')
                        topFile.write('[ molecules ]\n')
                        l = max(len('Compound') + 3, len(name) + 1)
                        compStr = 'Compound'.ljust(l)
                        topFile.write('; %s nmols\n' % compStr)
                        topFile.write('%s 1\n' % name.ljust(l))
                        break

                    elif line.find('[ moleculetype ]') != -1 or line.find('[moleculetype]') != -1:
                        if not itpFileIncluded:
                            topFile.write('#include "%s"\n\n' % os.path.basename(itpFileName))
                            itpFileIncluded = True

                        #topFile.write('#include "%s"\n\n' % outputFileBaseName)
                        itpFile.write(line)
                        itpFile.write('; Name            nrexcl\n')
                        itpFile.write('%s\t3\n' % name)
                        inMoleculeSection = 3


                    elif not headerWritten and line.find('#include') != -1:
                        topFile.write('#include "%s"\n' % forcefieldFile)
                        if not itpFileIncluded:
                            topFile.write('#include "%s"\n\n' % os.path.basename(itpFileName))
                            itpFileIncluded = True

                        # Add the contents of the nonbonded and bonded files directly in the itp file.
                        # This is not exactly elegant, but simplifies file management for e.g. running
                        # Copernicus.
                        for ffline in ffnonbondedContents:
                            itpFile.write(ffline)
                        itpFile.write('\n')

                        inDihedrals = False
                        for ffline in ffbondedContents:
                            # Dihedral parameters for triple bonds must be fixed. Currently only fixing C1T1.
                            if inDihedrals:
                                if '[' in ffline:
                                    inDihedrals = False

                                elif 'C1T1 C1T1' in ' '.join(ffline.split()):
                                    parts = ffline.split()
                                    if len(parts) >= 8:
                                        ffline = '%s\t%s\t%s\t%s\t%s\t%s\t0.00\t%s\n' % (parts[0],
                                        parts[1], parts[2], parts[3], parts[4], parts[5], parts[7])
                            elif '[ dihedraltypes ]' in ffline or '[dihedraltypes]' in ffline:
                                inDihedrals = True

                            itpFile.write(ffline)

                        headerWritten = True
                        if solvent:
                            solventFile = os.path.join(forcefieldDir, solvent + '.itp')
                            topFile.write('#include "%s"\n' % solventFile)
                            topFile.write('#include "%s"\n' % os.path.join(forcefieldDir, 'ions.itp'))

                    else:
                        topFile.write(line)

                #itpFile.write('\n\n#ifdef POSRES || POSRES_LIG\n')
                itpFile.write('\n\n#ifdef POSRES_LIG\n')
                itpFile.write('#include "%s"\n' % ligRestrFile)
                itpFile.write('#endif\n')

        return topologyFileName

    def coordsToTopology(self, output, coordsFile, forcefield = 'charmm27', verbose = False):
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

        os.chdir(currDir)

        return topolFile, groFile, restraintsFile

    def finalClean(self, output):

        forcefieldDir = os.path.join(output + '_%s' % self.forceFieldName, 'forcefield.ff')

        shutil.rmtree(forcefieldDir)

        posre_file = os.path.join(output + '_%s' % self.forceFieldName, 'posre.itp')
        if os.path.isfile(posre_file):
            os.remove(posre_file)

        gmxbk_files = glob('./#*', recursive=True)

        for gmxbk_file in gmxbk_files:
            if os.path.isfile(gmxbk_file):
                os.remove(gmxbk_file)
