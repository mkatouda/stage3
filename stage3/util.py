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
import subprocess
import shutil
import time
from glob import glob
from numpy import array, empty, average, absolute, power
from math import sqrt

from .molScripts.molecule import MolSystem, Atom

# This is for Python < 2.7
try:
    subprocess.check_output
    check_output = subprocess.check_output   
except:
    def check_output(*popenargs, **kwargs):
        """Run command with arguments and return its output as a byte string.
        Backported from Python 2.7 as it's implemented as pure python on stdlib.
        >>> check_output(['/usr/bin/python', '--version'])
        Python 2.6.2
        """
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            error = subprocess.CalledProcessError(retcode, cmd)
            error.output = output
            raise error
        return output

def getCommandOutput(cmd, inp = None):
    """ Executes a system command and returns the output. """

    # FIXME: Should check that the stack trace to limit access.
    if not inp:
        output = ''.join(check_output(cmd, stderr = subprocess.STDOUT,
                                      text=True))
    else:
        p = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.STDOUT, text=True)

        output = ''.join(p.communicate(input = inp)[0])

    return output

def is_exe(fpath):
    if not fpath:
        return False
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def isInPath(program):
    """ Determines if an executable is available by its full name or in the
    system path for executables (PATH). """

    fpath = os.path.split(program)[0]
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return False

if not isInPath('gmx'):
    if isInPath('gmx_seq'):
        gmxSuffix = '_seq'
    elif isInPath('gmx_mpi'):
        gmxSuffix = '_mpi'
else:
    gmxSuffix = ''

def findPath(program):
    """ Finds the full path to an executable by searching through the system path for
    executables (PATH). """

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def babelConvert(inputFile = None, outputFile = None, smiles = None, pH = None,
                 doMulti = False, verbose = False):
    """ Convert the input file (or a smiles molecule) to a mol2 (or another format based
    on the extension of outputFile) file. If pH is given the molecule will be protonated
    according to that pH (not fool proof).
    If doMulti is True multiple molecules in the file will be processed and
    split into separate files. If verbose is True more detailed output will be
    given. """

    if not outputFile:
        outputFile = 'temp'

    outputFileBaseName = os.path.basename(outputFile)
    parts = outputFileBaseName.split('.')
    if len(parts) == 1:
        molFile = outputFile + '.mol2'
    else:
        molFile = outputFile

    if isInPath('obabel'):
        binaryPath = 'obabel'
    else:
        progDir = os.path.dirname(__file__)
        #binaryPath = os.path.join(progDir, '..', 'external', 'openbabel', 'build', 'bin', 'obabel')
        binaryPath = os.path.join(progDir, 'external', 'openbabel', 'build', 'bin', 'obabel')
        if not isInPath(binaryPath):
            print('Cannot find openbabel binary')
            return []

    if not smiles:
        iFormat = os.path.splitext(inputFile)[1]
        if iFormat == '.smi':
            babelCmd = [binaryPath, '-ismi', inputFile, '--gen3D']
        else:
            babelCmd = [binaryPath, '-i' + iFormat[1:], inputFile]
    else:
        print(smiles)
        babelCmd = [binaryPath, "-:" + smiles, '--gen3D']
    if doMulti:
        babelCmd += ['--separate', '-m', '-v[OH2]']

    # Use mmff94 charges to get a good net charge but sometimes fails bug, thus change to gasteiger.
    #babelCmd += ['--partialcharge', 'mmff94']
    babelCmd += ['--partialcharge', 'gasteiger']
    oFormat = os.path.splitext(molFile)[1]
    babelCmd += ['-o' + oFormat[1:], '-O', molFile]

    if verbose:
        print(' '.join(babelCmd))

    try:
        #result = getCommandOutput(babelCmd, smiles)
        result = getCommandOutput(babelCmd)
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(babelCmd))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    if doMulti:
        molFiles = glob('%s[0-9].mol2' % outputFile)

        try:
            with open(list(molFiles)[0]):
                pass
        except IOError as e:
            return []
    else:
        molFiles = [molFile]

        try:
            with open(molFile):
                pass

        except IOError as e:
            return []

    # For some reason protonation according to pH cannot be done at the same
    # time as conversion from SMILES. So, for consistency it is always done
    # in a separate step
    if pH != None:

        # The default phmodel.txt is not completely correct. To use a modified version
        # change the directory to the one of this script.
        currDir = os.getcwd()
        progDir = os.path.dirname(__file__)
        os.chdir(progDir)

        try:
            for molFile in molFiles:
                babelCmd =  [binaryPath, molFile, '-p %.2f' % pH, molFile]
                # Use mmff94 charges to get a good net charge but sometimes fails bug, thus change to gasteiger.
                #babelCmd += ['--partialcharge', 'mmff94']
                babelCmd += ['--partialcharge', 'gasteiger']
                if verbose:
                    print(' '.join(babelCmd))

                result = getCommandOutput(babelCmd)
                if verbose:
                    print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(babelCmd))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        finally:
            os.chdir(currDir)

    # Convert to pdb (in case it is useful).
    # The mol2 charges are not always set correctly if changing
    # protonation states so also make a new mol2 using the mmff94
    # charge method.
    """
    try:
        for molFile in molFiles:
            pdbFile = os.path.splitext(molFile)[0] + '.pdb'
            babelCmd =  [binaryPath, '-imol', molFile, '-omol', '-O', pdbFile]
            if verbose:
                print(' '.join(babelCmd))

    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(babelCmd))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')
    """

    return molFiles

def renameAtoms(mol2File):

    with open(mol2File) as f:
        contents = f.readlines()
    molsys = MolSystem()
    molsys.readMol2(contents)

    molsys.renameAtoms()

    with open(mol2File, 'w') as f:
        molsys.writeMol2(f)

def getNetChargeOfMol2(fileName):
    """ Returns the net charge of a mol2 file (only one molecule per file) """

    totCharge = 0
    with open(fileName) as f:
        inAtom = False
        for line in f.readlines():
            line = line.strip()
            if inAtom:
                parts = line.split()
                if len(parts) < 9:
                    inAtom = False
                    break
                lineCharge = float(parts[8])
                if lineCharge:
                    totCharge += lineCharge
            else:
                if len(line) > 4 and line[0] == '@' and line[-4:] == 'ATOM':
                    inAtom = True

    return round(totCharge)

def getChargeOfTopology(topology):
    """ Reads a topology file and returns the net charge. Included itp files are
    read as well.
    FIXME: The net charge is only read from qtot comment of the last atom in each file. """

    netCharge = 0

    topologyDir = os.path.dirname(topology)

    with open(topology) as f:
        topLines = f.readlines()

    itpFiles = [topology]

    for line in topLines:
        if '#include' in line and '.itp' in line and not '.ff/' in line:
            itpFile = line.split()[-1]
            itpFile = itpFile.strip('"')
            itpFiles.append(itpFile)

    for itpFile in itpFiles:
        if not os.path.exists(itpFile):
            itpFile = os.path.join(topologyDir, itpFile)
        with open(itpFile) as f:
            itpLines = f.readlines()

        i = len(itpLines) - 1

        while i > 0:
            line = itpLines[i]
            parts = line.split(';')

            if len(parts) > 1:
                parts = parts[-1].split()
                if len(parts) > 1 and parts[0] == 'qtot':
                    fileCharge = float(parts[-1])
                    netCharge += fileCharge
                    #print(itpFile, fileCharge, netCharge)
                    break

            i -= 1

    return int(round(netCharge))

def itpChargeToMol2(itpFile, mol2File):

    with open(mol2File, 'r') as f:
        mol2Lines = f.readlines()

    molsys = MolSystem()
    molsys.readMol2(mol2Lines)

    with open(itpFile, 'r') as f:
        itpLines = f.readlines()

    inAtoms = False
    for line in itpLines:
        if line[0] == '[':
            if inAtoms:
                break
            if line.find('atoms') != -1:
                inAtoms = True
                continue
        line = line.strip()
        if len(line) <= 0 or line[0] == ';':
            continue

        if inAtoms:
            parts = line.split()
            if len(parts) < 7:
                continue

            atom = molsys.findAtomNumber(int(parts[0]))
            if not atom:
                print('Atom %s not found' % parts[0])
                return False

            atom.charge = float(parts[6])

    with open(mol2File, 'w') as f:
        molsys.writeMol2(f)

    return True

def mol2RenameToLig(mol2File):
    """ Read the specified mol2 file, rename the residues to LIG and
    save the file with the same name. """

    with open(mol2File) as f:
        contents = f.readlines()

    molsys = MolSystem()
    molsys.readMol2(contents)

    for molecule in molsys.molecules:
        for residue in molecule.residues:
            residue.setName('LIG')

    with open(mol2File, 'w') as f:
        molsys.writeMol2(f)

def getGroLigandName(groFile):
    """ Return the name of a molecule in a .gro file """

    with open(groFile) as f:
        lines = f.readlines()

    for line in lines[2:]:
        if len(line) > 10:
            nr = line[0:5].strip()
            name = line[5:10].strip()
            if nr.isdigit():
                return name

def findForceFieldsDir(forcefield):

    forcefieldDir = forcefield + '.ff'
    forcefieldsDir = os.environ.get('GMXLIB') or \
                     os.path.join(os.path.sep, 'usr', 'local', 'gromacs', 'share', 'gromacs' ,'top')
    if not os.path.exists(os.path.join(forcefieldsDir, forcefieldDir)):
        gmxPath = findPath('gmx'+gmxSuffix)
        print('gmxPath: ', gmxPath, os.path.dirname(gmxPath))
        forcefieldsDir = os.path.join(os.path.dirname(gmxPath), '..', 'share', 'gromacs', 'top')
        if not os.path.exists(os.path.join(forcefieldsDir, forcefieldDir)):
            progDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
            forcefieldsDir = os.path.join(progDir, 'forcefields')

    print('forcefieldsDir: ', forcefieldsDir)
    print('forcefieldDir: ', forcefieldDir)

    return forcefieldsDir

def setRtpLigandName(rtpFile, ligName):
    """ Set the name of a molecule in a .rtp file """

    with open(rtpFile) as f:
        lines = f.readlines()

    with open(rtpFile, 'w') as f:
        for line in lines:
            if line[0] == '[' and line.find('atoms') == -1 and \
               line.find('bonds') == -1 and line.find('impropers') == -1 and \
               line.find('bondedtypes') == -1:
                   line = '[ %s ]\n' % ligName
            f.write(line)

def runGamess(name, netCharge, suffix = '', gamessResultsFile = None, outputDir = None, onlyGenerateRunFiles = False, verbose = False):

    if not outputDir:
        outputDir = os.getcwd()

    currDir = os.getcwd()
    os.chdir(outputDir)

    inFileName = name + suffix + '.inp'
    outFileName = name + suffix + '.gamout'
    if gamessResultsFile:
        molFileName = None
    else:
        molFileName = name + '.mol2'

    if not (molFileName or gamessResultsFile):
        print('A mol2 file or a GAMESS results file must be used as input')
        os.chdir(currDir)
        return False

    #templateDir = os.path.join(os.path.dirname(__file__), '..', 'templates')
    templateDir = os.path.join(os.path.dirname(__file__), 'templates')
    templateFileName = os.path.join(templateDir, 'gamess' + suffix + '.inp')

    with open(templateFileName) as f:
        templateLines = f.readlines()

    babelConvert(molFileName or gamessResultsFile, inFileName)

    with open(inFileName) as f:
        lines = f.readlines()

    basisSetCompatible = True
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:
            if float(parts[1]) > 50.1:
                basisSetCompatible=False
                break

    with open(inFileName, 'w') as f:
        line = templateLines[0]
        line = line.replace('$END', 'ICHARG=%d $END' % netCharge)
        f.write(line)
        for line in templateLines[1:]:
            if not basisSetCompatible and '=CCT' in line:
                line = line.replace('=CCT', '=SPK-TZP')
            f.write(line)
        for line in lines[2:]:
            f.write(line)

    if onlyGenerateRunFiles:
        return True

    if isInPath('rungms'):
        binaryPath = 'rungms'
    else:
        progDir = os.path.dirname(__file__)
        #binaryPath = os.path.join(progDir, '..', 'external', 'gamess', 'rungms')
        binaryPath = os.path.join(progDir, 'external', 'gamess', 'rungms')
        if not isInPath(binaryPath):
            print('Cannot find rungms binary')
            os.chdir(currDir)
            return False

    scratchDir = None
    binaryPath = findPath(binaryPath)
    with open(binaryPath) as f:
        for line in f.readlines():
            if 'set USERSCR' in line:
                parts = line.strip().split('=')
                if len(parts) > 1:
                    scratchDir = os.path.expanduser(parts[-1].replace('$USER', os.environ['USER']))
                    break

    if scratchDir:
        scratchFiles = glob(os.path.join(scratchDir, name + suffix + '.*'))
        for scratchFile in scratchFiles:
            os.remove(scratchFile)

    gamessCommand = [binaryPath, inFileName]

    if verbose:
        print(' '.join(gamessCommand))

    try:
        result = getCommandOutput(gamessCommand)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(gamessCommand))
        print(e)
        os.chdir(currDir)
        return False

    with open(outFileName, 'w') as f:
        for line in result:
            f.write(line)

    shutil.copy(os.path.join(scratchDir, name + suffix + '.dat'), outputDir)

    os.chdir(currDir)

    return True

def runGaussian(name, netCharge, multiplier, suffix = '', gaussianResultsFile = None, outputDir = None, onlyGenerateRunFiles = False, verbose = False):

    if not outputDir:
        outputDir = os.getcwd()

    currDir = os.getcwd()
    os.chdir(outputDir)

    inFileName = name + suffix + '.com'
    outFileName = name + suffix + '.log'
    if gaussianResultsFile:
        molFileName = None
    else:
        molFileName = name + '.mol2'

    if not (molFileName or gaussianResultsFile):
        print('A mol2 file or a Gaussian results file must be used as input')
        os.chdir(currDir)
        return False

    #templateDir = os.path.join(os.path.dirname(__file__), '..', 'templates')
    templateDir = os.path.join(os.path.dirname(__file__), 'templates')
    templateFileName = os.path.join(templateDir, 'gaussian' + suffix + '.com')

    with open(templateFileName) as f:
        templateLines = f.readlines()

    babelConvert(molFileName or gaussianResultsFile, inFileName)

    with open(inFileName) as f:
        lines = f.readlines()

    with open(inFileName, 'w') as f:
        for line in templateLines:
            f.write(line)
        line = ('%d %d\n' % (netCharge, multiplier))
        f.write(line)
        for line in lines[6:]:
            f.write(line)

    if onlyGenerateRunFiles:
        return True

    if isInPath('g16'):
        binaryPath = 'g16'
    elif isInPath('g09'):
        binaryPath = 'g09'
    else:
        progDir = os.path.dirname(__file__)
        #binaryPath = os.path.join(progDir, '..', 'external', 'g16', 'rung16')
        binaryPath = os.path.join(progDir, 'external', 'g16', 'rung16')
        if not isInPath(binaryPath):
            print('Cannot find g16 binary')
            os.chdir(currDir)
            return False

    """
    scratchDir = None
    binaryPath = findPath(binaryPath)
    with open(binaryPath) as f:
        for line in f.readlines():
            if 'set USERSCR' in line:
                parts = line.strip().split('=')
                if len(parts) > 1:
                    scratchDir = os.path.expanduser(parts[-1].replace('$USER', os.environ['USER']))
                    break

    if scratchDir:
        scratchFiles = glob(os.path.join(scratchDir, name + suffix + '.*'))
        for scratchFile in scratchFiles:
            os.remove(scratchFile)
    """

    with open(inFile) as f:
        inFileContents = f.read()

    gaussianCommand = [binaryPath]

    if verbose:
        print(' '.join(gaussianCommand))

    try:
        result = getCommandOutput(gaussianCommand, inFileContents)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(gaussianCommand))
        print(e)
        os.chdir(currDir)
        return False

    with open(outFileName, 'w') as f:
        for line in result:
            f.write(line)

    os.chdir(currDir)

    return True

def runMakeresp(inFile, verbose = False):

    inFileBaseName = os.path.splitext(inFile)[0]

    if isInPath('makeresp'):
        binaryPath = 'makeresp'
    else:
        progDir = os.path.dirname(__file__)
        #binaryPath = os.path.join(progDir, '..', 'external', 'makeresp', 'makeresp')
        binaryPath = os.path.join(progDir, 'external', 'makeresp', 'makeresp')
        if not isInPath(binaryPath):
            print('Cannot find makeresp binary')
            return False

    makerespCommand = [binaryPath, inFile]

    if verbose:
        print(' '.join(makerespCommand))

    try:
        result = getCommandOutput(makerespCommand)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(makerespCommand))
        print(e)
        return False

    if os.path.isfile(inFileBaseName+'.res'):
        os.remove(inFileBaseName+'.res')
    if os.path.isfile(inFileBaseName+'_qout'):
        os.remove(inFileBaseName+'_qout')
    if os.path.isfile('punch'):
        os.remove('punch')

    respCommand = ['resp', '-i', inFileBaseName+'.in', '-e', inFileBaseName+'.esp', '-o', inFileBaseName+'.res', '-t', inFileBaseName+'_qout']

    if verbose:
        print(' '.join(respCommand))

    try:
        result = getCommandOutput(respCommand)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(respCommand))
        print(e)
        return False

    if os.path.isfile('punch'):
        os.remove('punch')

    return True

def runGmstoresp(inFile, verbose = False):

    inFileBaseName = os.path.splitext(inFile)[0]

    if isInPath('gmstoresp.sh'):
        binaryPath = 'gmstoresp.sh'
    else:
        progDir = os.path.dirname(__file__)
        #binaryPath = os.path.join(progDir, '..', 'external', 'gmstoresp', 'gmstoresp.sh')
        binaryPath = os.path.join(progDir, 'external', 'gmstoresp', 'gmstoresp.sh')
        if not isInPath(binaryPath):
            print('Cannot find gmstoresp.sh binary')
            return False

    with open(inFile) as f:
        inFileContents = f.read()

    gmstorespCommand = [binaryPath]

    if verbose:
        print(' '.join(gmstorespCommand))

    try:
        result = getCommandOutput(gmstorespCommand, inFileContents)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(gmstorespCommand))
        print(e)
        return False

    with open(inFileBaseName + '.ac', 'w') as f:
        f.write(result)

    return True

def generateCharges(inFile, method, netCharge, multiplier = 1.0, verbose = False):

    netCharge = int(round(netCharge,0))

    shutil.copy(inFile, inFile + '.bak')
    outFile =  inFile
    inFile = inFile + '.bak'

    if method == 'am1bcc-pol':
        molName = os.path.splitext(os.path.basename(inFile))[0]
        antechamberFileName = molName + '.ac'
        command = ['antechamber', '-fi', 'mol2', '-i', inFile, '-fo', 'ac', '-o', antechamberFileName, '-nc', '%d' % netCharge, '-c', 'mul', '-pf', 'y']

        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        #chargeSymPath = os.path.join(os.path.dirname(__file__), '..', 'external', 'chargesym.sh')
        chargeSymPath = os.path.join(os.path.dirname(__file__), 'external', 'chargesym.sh')

        command = [chargeSymPath]
        with open(antechamberFileName) as antechamberFile:
            contents = antechamberFile.read()
        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command, contents)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        with open(antechamberFileName, 'w') as antechamberFile:
            for line in result:
                antechamberFile.write(line)

        #parmFile = os.path.join(os.path.dirname(__file__), '..', 'ffOptimization', 'bcc-pol-parm.dat')
        parmFile = os.path.join(os.path.dirname(__file__), 'ffOptimization', 'bcc-pol-parm.dat')

        command = ['am1bcc', '-i', antechamberFileName, '-o', antechamberFileName, '-j', '4', '-p', parmFile]

        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        command = ['antechamber', '-fi', 'ac', '-i', antechamberFileName, '-fo', 'mol2', '-o', outFile, '-at', 'sybyl', '-pf', 'y']

        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

    elif method == 'am1bcc':
        command = ['antechamber', '-fi', 'mol2', '-i', inFile, '-fo', 'mol2', '-o', outFile, '-nc', '%d' % netCharge, '-c', 'bcc', '-pf', 'y']

        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        if abs(multiplier - 1.0) > 0.001:
            molSys = MolSystem()
            with open(outFile, 'r') as f:
                molSys.readMol2(f.readlines())
            molSys.multiplyCharges(multiplier)
            with open(outFile, 'w') as f:
                molSys.writeMol2(f)

    elif method in ['gasteiger', 'mmff94', 'eem', 'eqeq', 'qeq', 'qtpie', 'fromfile']:
        command = ['obabel', '-imol2', inFile, '-omol2', '-O', outFile, '--partialcharge', method]
        try:
            if verbose:
                print(' '.join(command))

            result = getCommandOutput(command)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(command))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Trying to continue anyhow.')

        if abs(multiplier - 1.0) > 0.001:
            molSys = MolSystem()
            with open(outFile, 'r') as f:
                molSys.readMol2(f.readlines())
            molSys.multiplyCharges(multiplier)
            with open(outFile, 'w') as f:
                molSys.writeMol2(f)

    d = os.path.dirname(inFile)
    rmFiles = glob(os.path.join(d, 'ANTECHAMBER*')) + glob(os.path.join(d, 'ATOMTYPE*')) + glob(os.path.join(d, 'sqm.*'))
    for f in rmFiles:
        os.remove(f)

def molecule2GamessDipole(molecule, outputDir, verbose = False):

    netCharge = int(getNetChargeOfMol2(molecule))
    if verbose:
        print('Net charge of molecule is %d' % netCharge)
        print('GAMESS gas optimization starting.')

    name = os.path.splitext(os.path.basename(molecule.split))[0]
    if not os.path.isfile(os.path.join(outputDir, name+'_gas_optimize.gamout')):
        if not runGamess(name, netCharge, '_gas_optimize', outputDir = outputDir, onlyGenerateRunFiles = False, verbose = verbose):
            print('GAMESS failed running gas optimization.')
            return None

    if verbose:
        print('GAMESS gas optimization finished.')
        print('GAMESS gas dipole and polarization calculations starting.')

    if not runGamess(name, netCharge, '_gas_dipole', name+'_gas_optimize.gamout', outputDir = outputDir, onlyGenerateRunFiles = False, verbose = verbose):
        print('GAMESS failed running dipole and polarization calculations.')
        return None

    if verbose:
        print('GAMESS gas dipole and polarization calculations finished.')

    return os.path.join(outputDir, name+'_gas_dipole.gamout')

def extractGamessDipoleData(gamessFile):

    dipole = alphaPol = center = None

    with open(gamessFile, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if '/D/' in line and 'APPLIED FIELD' in lines[i+7]:
                line = lines[i-1]
                parts = line.split()
                if len(parts) == 4:
                    center = array([[float(parts[0]), float(parts[1]), float(parts[2])]])
                    # Convert from Bohr to Angstroms
                    center *= 1.88971616463
                line = lines[i+1]
                parts = line.split()
                if len(parts) == 4:
                    dipole = array([[float(parts[0]), float(parts[1]), float(parts[2])]])
                    break

        for i in range(len(lines)-5, max(0, len(lines)-200), -1):
            line = lines[i]
            if 'DIPOLE BASED RESULTS' in line:
                alphaPol = empty((3,3))
                for j in range(3):
                    line = lines[i + 8 + j]
                    parts = line.split()
                    for k in range(3):
                        alphaPol[j,k] = float(parts[2 + k])

                break

    return(dipole, alphaPol, center)

def calculateMol2Dipole(mol2File, centerCoords = None):

    with open(mol2File, 'r') as f:
        mol2Lines = f.readlines()

    molsys = MolSystem()
    molsys.readMol2(mol2Lines)
    mol = molsys.molecules[0]
    dipole = array([mol.getDipole(centerCoords)])

    return dipole

def calculateGromacsDipole(trajectoryFile, runInputFile, indexFile, start = 0, verbose = False):

    outputDir = os.path.dirname(trajectoryFile)
    currDir = os.getcwd()
    if outputDir:
        os.chdir(outputDir)

    os.environ['GMX_MAXBACKUP'] = '0'

    outputFile = 'dipole.xvg'
    command = ['gmx'+gmxSuffix, 'g_dipoles', '-f', trajectoryFile, '-s', runInputFile, '-n', indexFile, '-b', '%d' % start, '-o', outputFile]

    try:
        result = getCommandOutput(command, 'LIG\n')

    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(command))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        os.chdir(currDir)
        return(None)

    with open(outputFile, 'r') as f:
        lines = f.readlines()

    dipoleSum = 0
    cnt = 0

    for line in lines:
        if line[0] == '#' or line[0] == '@':
            continue
        parts = line.split()
        if len(parts) == 5:
            dipoleSum += float(parts[4])
            cnt += 1

    dipole = dipoleSum/cnt

    os.chdir(currDir)

    return dipole

def replaceRtfCharges(rtfFile, mol2File):

    molSys = MolSystem()
    with open(mol2File, 'r') as f:
        molSys.readMol2(f.readlines())

    with open(rtfFile, 'r') as f:
        lines = f.readlines()

    with open(rtfFile, 'w') as f:
        for line in lines:
            parts = line.split()
            if len(parts) > 0:
                if parts[0] == "RESI":
                    resname = parts[1]
                    residues = []
                    for molecule in molSys.molecules:
                        residue = molecule.findResidue(resname)
                        if residue:
                            break
                if parts[0] == "ATOM":
                    atom = residue.findAtom(parts[1])
                    line = '%s %s %s %9f\n' % (parts[0], parts[1], parts[2], atom.getCharge())
            f.write(line)

def _findAtomInAtomList(atomList, nr):

    nAtoms = len(atomList)
    inr = int(nr)

    atom = atomList[inr-1]
    parts = atom.split()

    iatomnr = int(parts[0])

    if iatomnr == inr:
        return atom

    i = inr

    if iatomnr < inr:
        while iatomnr < inr and i < nAtoms:
            atom = atomList[i]
            iatomnr = int(parts[0])
            i += 1

        if iatomnr == inr:
            return atom

    i -= 1

    while iatomnr > inr and i >= 0:
        atom = atomList[i]
        iatomnr = int(parts[0])
        i -= 1

    if iatomnr == inr:
        return atom

    return None

def _setMassOfAtomInAtomList(atomList, atomNr, mass):

        atom = atomList[int(atomNr)-1]
        parts = atom.split(';')
        line = parts[0]
        if len(parts) > 1:
            comment = ';'.join(parts[1:])
        else:
            comment = None

        parts = line.split()
        if len(parts) >= 8:
            resnr = parts[2]
            resname = parts[3]
            atom = '%6s %8s %6s %4s %5s %5s %10s %12s' % (parts[0], parts[1],
                   resnr, resname, parts[4], parts[5], parts[6], mass)
            if comment:
                atom += ' ; %s' % comment

            atomList[int(atomNr)-1] = atom

def _findBondsToAtomInBondList(bondList, atomNr):

    bonds = []

    for i, bond in enumerate(bondList):
        parts = bond.split()
        if parts[0] == atomNr or parts[1] == atomNr:
            bonds.append(i)

    return bonds

def _findBondLengthsOfAtomTypes(bondTypeList, firstAtomType, lastAtomType):

    for bond in bondTypeList:
        bondParts = bond.split()
        if len(bondParts) < 5:
            continue
        if firstAtomType == bondParts[0] and lastAtomType == bondParts[1] or \
        firstAtomType == bondParts[1] and lastAtomType == bondParts[0]:
            return float(bondParts[3])

def _findBondLengthsOfAtoms(bondList, firstAtom, lastAtom):

    for bond in bondList:
        bondParts = bond.split()
        if len(bondParts) < 5:
            continue
        if firstAtom == bondParts[0] and lastAtom == bondParts[1] or \
        firstAtom == bondParts[1] and lastAtom == bondParts[0]:
            return float(bondParts[3])

def _optimizeMomentOfInertia(molecule, atom1, atom2, targetMomentOfInertia, accuracy = 0.00001):

    atom1Mass = atom1.getMass()
    atom2Mass = atom2.getMass()

    f = 1

    mom = molecule.getMomentOfInertia() / 100

    while abs(mom - targetMomentOfInertia) > accuracy and abs(f) > accuracy:

        #print 'diff', abs(mom - targetMomentOfInertia)
        #print 'f', f
        #print 'atom1Mass', atom1Mass
        #print 'atom2Mass', atom2Mass

        atom1Mass = atom1Mass + f
        atom2Mass = atom2Mass - f
        atom1.setMass(atom1Mass)
        atom2.setMass(atom2Mass)

        newMom = molecule.getMomentOfInertia() / 100

        if abs(targetMomentOfInertia - newMom) > abs(targetMomentOfInertia - mom):
            f = -f * 0.5

        mom = newMom

    return mom

def _optimizeMassAndDist(referenceAtomMass, centerOfMass, massCenterMass, momentOfInertia, accuracy = 0.00001):

    x = 0

    f = 1

    mom = pow(centerOfMass, 2) * x + pow(centerOfMass*x + centerOfMass*referenceAtomMass, 2) / (massCenterMass - x) + pow(centerOfMass, 2)*referenceAtomMass

    while abs(mom - momentOfInertia) > accuracy and abs(f) > accuracy:

        x += f

        #print 'diff', abs(mom - momentOfInertia)
        #print 'x', x, 'f', f

        newMom = pow(centerOfMass, 2) * x + pow(centerOfMass*x + centerOfMass*referenceAtomMass, 2) / (massCenterMass - x) + pow(centerOfMass, 2)*referenceAtomMass

        #print 'mom:', mom, 'newMom:', newMom

        if abs(momentOfInertia - newMom) > abs(momentOfInertia - mom):
            f = -f * 0.5

        mom = newMom

    massDiff = x

    massCenterPos = (centerOfMass * (referenceAtomMass + massDiff)) / (massCenterMass - massDiff)

    #print 'massDiff:', massDiff
    #print 'massCenterPos:', massCenterPos

    return (referenceAtomMass + massDiff, massCenterMass - massDiff, massCenterPos)

def _findIntersectionOnLine(coord1, coord2, eq, distToCoord2, accuracy = 0.00001):

    f = 1

    pos = 3

    dist = sqrt(sum(power(coord1 - coord2 + pos * eq, 2)))

    while abs(dist - distToCoord2) > accuracy and abs(f) > 0.001:

        #print 'pos', pos
        #print 'diff', abs(dist - distToCoord2)
        #print 'f', f

        pos += f

        newDist = sqrt(sum(power(coord1 - coord2 + pos * eq, 2)))

        if abs(distToCoord2 - newDist) > abs(distToCoord2 - dist):
            f = -f * 0.5

        dist = newDist

    return pos

def _readDataFields(lines, readableDataFields):

    currentDataField = None

    # Read the relevant data fields.
    for line in lines:
        modline = line.strip()
        if len(modline) < 1 or modline.startswith(';') or modline.startswith('#'):
            continue
        if modline.startswith('['):
            changedDataField = False
            for header in readableDataFields.keys():
                if header in line:
                    currentDataField = header
                    changedDataField = True
                    break
            else:
                currentDataField = None
            continue

        if currentDataField:
            readableDataFields[currentDataField].append(line)



def hydrogens2VirtualSites(outputDir, outputFileBaseName, verbose = False):
    """ Turns all hydrogens into virtual sites to allow longer timesteps. Experimental! """

    itpFile = os.path.join(outputDir, outputFileBaseName + '.itp')
    inGroFile = os.path.join(outputDir, outputFileBaseName + '.gro')
    outGroFile = os.path.join(outputDir, outputFileBaseName + '.gro')



def generateLinearVirtualSites(outputDir, outputFileBaseName, verbose = False):
    """ Turns collinear atoms into mass less virtual sites and creates new mass centers, if required. The total mass and
        moment of inertia is retained. """

    itpFile = os.path.join(outputDir, outputFileBaseName + '.itp')
    inGroFile = os.path.join(outputDir, outputFileBaseName + '.gro')
    outGroFile = os.path.join(outputDir, outputFileBaseName + '.gro')

    if not os.path.exists(inGroFile):
        inGroFile = os.path.join(outputDir, '..', outputFileBaseName + '.gro')

    # The names and atom types of the mass centers
    massCenterAtomType = 'MX'
    massCenterAtomName = 'MC'

    with open(itpFile) as f:
        lines = f.readlines()

    # Data fields that can be modified by this function
    dataFields = {'atomtypes' : [], 'bondtypes' : [], 'angletypes' : [], 'atoms' : [],
                 'bonds' : [], 'angles' : [], 'constraints' : [], 'virtual_sites2' : [],
                 'exclusions' : []}

    _readDataFields(lines, dataFields)

    anglesToRemove = []

    # Deduce what atoms are involved in the linear angles.
    affectedAtoms = []
    affectedAtomTypes = []
    if dataFields['angletypes']:
        for angle in dataFields['angletypes']:
            parts = angle.split()
            if len(parts) < 6:
                continue

            # This is only necessary for angles of type 1. Urey-Bradley works anyhow.
            if int(parts[3]) != 1:
                continue

            deg = float(parts[4])
            if abs(deg - 180.00) < 0.001:
                affectedAtomTypes.append(parts[0:3])
    else:
        for i, angle in enumerate(dataFields['angles']):
            parts = angle.split()
            if len(parts) < 6:
                continue
            deg = float(parts[4])
            if abs(deg - 180.00) < 0.001:
                affectedAtoms.append(parts[0:3])
                anglesToRemove.append(i)

    if not affectedAtomTypes and not affectedAtoms:
        return

    if verbose:
        print('Turning linear atoms into mass centers and virtual sites')

    if affectedAtomTypes:
        for i, angle in enumerate(dataFields['angles']):
            parts = angle.split()
            midAtom = _findAtomInAtomList(dataFields['atoms'], parts[1])
            if midAtom:
                midType = midAtom.split()[1]
                for affectedAtomTypeGroup in affectedAtomTypes:
                    if midType == affectedAtomTypeGroup[1]:
                        firstAtom = _findAtomInAtomList(dataFields['atoms'], parts[0])
                        lastAtom = _findAtomInAtomList(dataFields['atoms'], parts[2])
                        if firstAtom and lastAtom:
                            firstType = firstAtom.split()[1]
                            lastType = lastAtom.split()[1]

                            if firstType == affectedAtomTypeGroup[0] and lastType== affectedAtomTypeGroup[2] or \
                            firstType == affectedAtomTypeGroup[2] and lastType == affectedAtomTypeGroup[0]:
                                affectedAtoms.append([firstAtom.split()[0], midAtom.split()[0], lastAtom.split()[0]])
                                anglesToRemove.append(i)

    # Check if the atomtype used for mass centers is already specified.
    for atomtype in dataFields['atomtypes']:
        parts = atomtype.split()
        if parts[0] == massCenterAtomType:
            break
    else:
        massCenter = '%5s %5s          0.000   0.000    A              0.000      0.000\n' % (massCenterAtomType, massCenterAtomType)
        dataFields['atomtypes'].append(massCenter)

    nAtoms = len(dataFields['atoms'])
    molSys = MolSystem()
    molSys.readGroFile(inGroFile)
    mol = molSys.molecules[0]

    # Count the number of atoms with mass center atom names.
    # Also set the masses of the atoms in the molecular system.
    nMassCentraInFile = 0
    for i, atom in enumerate(dataFields['atoms']):
        parts = atom.split()
        name = ''.join(c for c in parts[4] if not c.isdigit())
        if name == massCenterAtomName:
            nMassCentraInFile += 1

        atomNr = int(parts[0])
        groAtom = molSys.findAtomNumber(atomNr)
        groAtom.setMass(float(parts[7]))

    targetMomentOfInertia = mol.getMomentOfInertia() / 100
    targetCenterOfMass = mol.getCenterOfMass() / 10

    if verbose:
        print('Old moment of inertia of molecule: %s' % targetMomentOfInertia)
        print('Old center of mass of molecule: %s' % targetCenterOfMass)
        print('Old mass of molecule: %s' % mol.getTotalMass())

    atomGroupsToRemove = []
    nAffectedAtoms = len(affectedAtoms)
    # If there are groups with two atoms in common in either end
    # make one group consisting of all atoms instead.
    for i in range(nAffectedAtoms):
        currAtomGroup = affectedAtoms[i]

        for j in range(i + 1, nAffectedAtoms):
            if j in atomGroupsToRemove:
                continue
            otherAtomGroup = affectedAtoms[j]
            if len(otherAtomGroup) > 3:
                continue


            if otherAtomGroup[0:2] == currAtomGroup[-2:]:
                currAtomGroup.append(otherAtomGroup[2])
                atomGroupsToRemove.append(j)
            elif otherAtomGroup[0:2] == currAtomGroup[1::-1]:
                currAtomGroup.insert(0, otherAtomGroup[2])
                atomGroupsToRemove.append(j)
            elif otherAtomGroup[1:3] == currAtomGroup[0:2]:
                currAtomGroup.insert(0, otherAtomGroup[0])
                atomGroupsToRemove.append(j)
            elif otherAtomGroup[1::-1] == currAtomGroup[0:2]:
                currAtomGroup.append(otherAtomGroup[2])
                atomGroupsToRemove.append(j)

    for i in reversed(atomGroupsToRemove):
        del(affectedAtoms[i])

    if verbose:
        print('Atoms involved in linear angles:', affectedAtoms)

    bondsToRemove = []
    for (groupNr, atomGroup) in enumerate(affectedAtoms):
        atoms = []
        for i in range(len(atomGroup)):
            try:
                atom = dataFields['atoms'][int(atomGroup[i])-1]
                atoms.append(atom)
            except Exception:
                break
        if len(atoms) < 3:
            continue

        bondsToFirstAtom = _findBondsToAtomInBondList(dataFields['bonds'], atomGroup[0])
        bondsToLastAtom = _findBondsToAtomInBondList(dataFields['bonds'], atomGroup[-1])

        nBondsToFirstAtom = len(bondsToFirstAtom)
        nBondsToLastAtom = len(bondsToLastAtom)
        nMassCentra = 1

        bondLengths = []

        if nBondsToFirstAtom > 1 and nBondsToLastAtom <= 1:
            virtualAtomsNumbers = atomGroup[1:]
            referenceAtomNumber = atomGroup[0]
            virtualAtoms = atoms[1:]
            anchorAtoms = [referenceAtomNumber]
            if affectedAtomTypes:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtomTypes(dataFields['bondtypes'], atoms[i].split()[1], atoms[i+1].split()[1]))
            else:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtoms(dataFields['bonds'], atoms[i].split()[0], atoms[i+1].split()[0]))
        elif nBondsToLastAtom > 1 and nBondsToFirstAtom <= 1:
            virtualAtomsNumbers = atomGroup[0:-1]
            virtualAtomsNumbers.reverse()
            referenceAtomNumber = atomGroup[-1]
            virtualAtoms = atoms[0:-1]
            virtualAtoms.reverse()
            anchorAtoms = [referenceAtomNumber]
            if affectedAtomTypes:
                for i in reversed(range(len(atoms) - 1)):
                    #print i+1, atoms[i+1].split()[1], atoms[i].split()[1]
                    bondLengths.append(_findBondLengthsOfAtomTypes(dataFields['bondtypes'], atoms[i+1].split()[1], atoms[i].split()[1]))
            else:
                for i in reversed(range(len(atoms) - 1)):
                    bondLengths.append(_findBondLengthsOfAtoms(dataFields['bonds'], atoms[i+1].split()[0], atoms[i].split()[0]))
        elif nBondsToFirstAtom > 1 and nBondsToLastAtom > 1:
            virtualAtomsNumbers = atomGroup[1:-1]
            referenceAtomNumber = atomGroup[0]
            virtualAtoms = atoms[1:-1]
            if len(virtualAtoms) <= 1:
                nMassCentra = 0
            anchorAtoms = [referenceAtomNumber, atomGroup[-1]]
            if affectedAtomTypes:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtomTypes(dataFields['bondtypes'], atoms[i].split()[1], atoms[i+1].split()[1]))
            else:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtoms(dataFields['bonds'], atoms[i].split()[0], atoms[i+1].split()[0]))
        else:
            # Use the first atom as the "anchor". Keep it as a real atom.
            virtualAtomsNumbers = atomGroup[1:]
            referenceAtomNumber = atomGroup[0]
            virtualAtoms = atoms[1:]
            anchorAtoms = [referenceAtomNumber]
            if affectedAtomTypes:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtomTypes(dataFields['bondtypes'], atoms[i].split()[1], atoms[i+1].split()[1]))
            else:
                for i in range(len(atoms) - 1):
                    bondLengths.append(_findBondLengthsOfAtoms(dataFields['bonds'], atoms[i].split()[0], atoms[i+1].split()[0]))

        # Remove masses of atoms (make them virtual sites)
        for atomNr in virtualAtomsNumbers:
            _setMassOfAtomInAtomList(dataFields['atoms'], atomNr, 0.00000)

        nAtomsInGroup = len(atomGroup)

        nVirtualAtoms = len(virtualAtoms)
        groReferenceAtom = molSys.findAtomNumber(int(referenceAtomNumber))

        referenceAtomPos = groReferenceAtom.getCoords()

        groVirtualAtoms = []
        for vAtomNr in virtualAtomsNumbers:
            groVirtualAtoms.append(molSys.findAtomNumber(int(vAtomNr)))


        # Get bond lengths and remove bonds that will not be used anymore
        iterList = virtualAtomsNumbers[:]
        if referenceAtomNumber != atomGroup[0]:
            iterList.reverse()

        nAtomsInGroup = len(atomGroup)
        nBonds = nAtomsInGroup - 1

        for i, vAtom in enumerate(iterList):
            vAtomBonds = _findBondsToAtomInBondList(dataFields['bonds'], vAtom)
            for bondNr in vAtomBonds:
                #bond = dataFields['bonds'][bondNr]

                if bondNr not in bondsToRemove:
                    bondsToRemove.append(bondNr)

            bondsToRemove.sort()

        # Add exclusions for atoms that are no longer connected through bonds.
        for i in range(nAtomsInGroup):
            excl = []
            excl.append('%7s' % atomGroup[i])
            for j in range(nAtomsInGroup):
                if j == i:
                    continue
                excl.append('%7s' % atomGroup[j])
            if len(excl) > 1:
                exclusion = ' '.join(excl) + '\n'
                dataFields['exclusions'].append(exclusion)

        #groupMomentOfInertia = 0
        #groupMassCenter = 0
        virtualAtomsMassCenter = 0
        accBondLength = 0
        totMass = 0

        for (i, virtualAtomNr) in enumerate(virtualAtomsNumbers):
            groAtom = molSys.findAtomNumber(int(virtualAtomNr))
            mass = groAtom.getMass()
            accBondLength += bondLengths[i]
            virtualAtomsMassCenter += mass * accBondLength
            totMass += mass

        virtualAtomsMassCenter /= totMass

        accBondLength = 0

        referenceAtomMass = groReferenceAtom.getMass()
        if len(anchorAtoms) == 2:
            groSecondAnchorAtom = molSys.findAtomNumber(int(anchorAtoms[1]))


        if nMassCentra:
            massCenterMass = 0
            for groVirtualAtom in groVirtualAtoms:
                mass = groVirtualAtom.getMass()
                massCenterMass += mass
                groVirtualAtom.setMass(0)
        else:
            secondAnchorAtomMass = groSecondAnchorAtom.getMass()
            for groVirtualAtom in groVirtualAtoms:
                mass = groVirtualAtom.getMass()
                referenceAtomMass += mass / 2
                secondAnchorAtomMass += mass / 2
                groVirtualAtom.setMass(0)
            groSecondAnchorAtom.setMass(secondAnchorAtomMass)
            groReferenceAtom.setMass(referenceAtomMass)

        # distMatrix is the difference in coordinates between the reference atom
        # and the last virtual atom. This is the internal coordinate axis along
        # which the mass center should be located.
        # distMatrix / dist is the "slope" of the line connecting the atoms.
        distMatrix = groReferenceAtom.getDistanceMatrixToAtom(groVirtualAtoms[-1])
        dist = groReferenceAtom.getDistanceToAtom(groVirtualAtoms[-1])

        if nMassCentra:
            nMassCentraInFile += 1
            nAtoms += 1
            massCenterAtomNumber = nAtoms
            atomName = '%s%d' % (massCenterAtomName, nMassCentraInFile)
            massCenterAtomName = atomName
            residue = groReferenceAtom.getResidue()
            resnr = residue.getNumber()
            resname = residue.getName()

            refDist = virtualAtomsMassCenter

            pos = referenceAtomPos + distMatrix * refDist * 10 / dist
            groReferenceAtom.setMass(referenceAtomMass)
            groMassCenterAtom = Atom(residue, massCenterAtomNumber, massCenterAtomName)
            groMassCenterAtom.setCoords(pos[0], pos[1], pos[2])
            groMassCenterAtom.setVels(None, None, None)
            groMassCenterAtom.setMass(massCenterMass)
            mol.atoms.append(groMassCenterAtom)

            # Calibrate the mass of the mass center atom and the reference atom to keep the moment of inertia correct.
            _optimizeMomentOfInertia(mol, groMassCenterAtom, groReferenceAtom, targetMomentOfInertia)

            # The masses of the reference atom and the mass center atom have probably changed after the calibration.
            massCenterMass = groMassCenterAtom.getMass()
            referenceAtomMass = groReferenceAtom.getMass()

            dataFields['atoms'].append('%6s %8s %6s %4s %5s %5s %10s %12s\n' % (nAtoms, massCenterAtomType,
                                      resnr, resname, atomName, nAtoms, 0.000, massCenterMass))
            dataFields['constraints'].append('%7s %7s %7s   %.6f\n' % (referenceAtomNumber, nAtoms, 1, refDist))

        else:
            dataFields['constraints'].append('%7s %7s %7s   %.6f\n' % (referenceAtomNumber, atomGroup[-1], 1, refDist))


        #referenceAtomMass = groReferenceAtom.getMass()
        _setMassOfAtomInAtomList(dataFields['atoms'], referenceAtomNumber, referenceAtomMass)
        if nMassCentra == 0:
            secondAnchorAtomMass = groSecondAnchorAtom.getMass()
            _setMassOfAtomInAtomList(dataFields['atoms'], atomGroup[-1], secondAnchorAtomMass)

        # If there is one mass center and anchor atoms in the respective end of the collinear group it is possible
        # that virtual atoms are placed on either side of the mass center (closest to either of the anchor atoms).
        if nMassCentra and len(anchorAtoms) == 2:
            altRefDist = groMassCenterAtom.getDistanceToAtom(groSecondAnchorAtom)

        for i, vAtomNr in enumerate(virtualAtomsNumbers):
            vAtomNr = int(vAtomNr)
            groVirtualAtom = groVirtualAtoms[i]
            distToReferenceAtom = groVirtualAtom.getDistanceToAtom(groReferenceAtom)
            distFactor = distToReferenceAtom / (10 * refDist)
            if nMassCentra:
                if len(anchorAtoms) == 2:
                    distToMassCenter = groVirtualAtom.getDistanceToAtom(groMassCenterAtom)
                    distToAnchorAtom = groVirtualAtom.getDistanceToAtom(groSecondAnchorAtom)
                    # If the virtual site is closer to the mass center than the reference atom and closer to the other
                    # anchor atom specify the virtual site relative to those closer atoms instead.
                    if distToMassCenter < distToReferenceAtom and distToAnchorAtom < distToReferenceAtom:
                        distFactor = distToMassCenter / (10 * altRefDist)
                        dataFields['virtual_sites2'].append('%7s %7s %7s %7s   %.6f\n' % (vAtomNr, massCenterAtomNumber, atomGroup[-1], 1, distFactor))
                    else:
                        dataFields['virtual_sites2'].append('%7s %7s %7s %7s   %.6f\n' % (vAtomNr, referenceAtomNumber, massCenterAtomNumber, 1, distFactor))
                else:
                    dataFields['virtual_sites2'].append('%7s %7s %7s %7s   %.6f\n' % (vAtomNr, referenceAtomNumber, massCenterAtomNumber, 1, distFactor))
            else:
                dataFields['virtual_sites2'].append('%7s %7s %7s %7s   %.6f\n' % (vAtomNr, referenceAtomNumber, atomGroup[-1], 1, distFactor))

    if verbose:
        print('New moment of inertia of molecule: %s' % (mol.getMomentOfInertia() / 100))
        print('New center of mass of molecule: %s' % (mol.getCenterOfMass() / 10))
        print('New mass of molecule: %s' % mol.getTotalMass())

    with open(outGroFile, 'w') as f:
        molSys.writeGro(f)

    for i in reversed(bondsToRemove):
        if verbose:
            print('Removing bond: %s' % dataFields['bonds'][i])
        del(dataFields['bonds'][i])

    for i in reversed(anglesToRemove):
        if verbose:
            print('Removing angle: %s' % dataFields['angles'][i])
        del(dataFields['angles'][i])


    # Write the updated itp file.
    with open(itpFile, 'w') as f:
        currentDataField = None
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith('['):
                # Some fields should be written before [ moleculetype ]
                if 'moleculetype' in line:
                    if dataFields.get('atomtypes'):
                        data = dataFields['atomtypes']
                        if data:
                            f.write('\n[ atomtypes ]\n')

                            for dataLine in data:
                                f.write(dataLine)
                        del(dataFields['atomtypes'])

                    if dataFields.get('bondtypes'):
                        data = dataFields['bondtypes']
                        if data:
                            f.write('\n[ bondtypes ]\n')

                            for dataLine in data:
                                f.write(dataLine)
                        del(dataFields['bondtypes'])

                    if dataFields.get('angletypes'):
                        data = dataFields['angletypes']
                        if data:
                            f.write('\n[ angletypes ]\n')

                            for dataLine in data:
                                f.write(dataLine)
                        del(dataFields['angletypes'])

                for header in dataFields.keys():
                    if header in line:
                        currentDataField = header
                        f.write('\n%s' % line)
                        while i < len(lines) and lines[i+1].startswith(';'):
                            i += 1
                            f.write(lines[i])
                        for dataLine in dataFields[header]:
                            f.write(dataLine)
                        del(dataFields[header])
                        break
                else:
                    f.write('\n')
                    currentDataField = None

            elif line.startswith('#'):
                if currentDataField:
                    f.write('\n')
                currentDataField = None

            if currentDataField == None:
                f.write('%s' % line)

            i += 1

        for header, data in dataFields.items():
            if data:
                f.write('\n[ %s ]\n' % header)

                for dataLine in data:
                    f.write(dataLine)

def copyItp(topology, outputDir, verbose = False):
    """ Copy user supplyied itp files to working directiry. """

    with open(topology) as topf:
        topologyLines = topf.readlines()

        for line in topologyLines:
            if '#include' in line and not '.ff/' in line:
                incF = line.split('"')[1]
                if verbose:
                    print('include file:', incF)
                if os.path.exists(incF):
                    shutil.copy(incF, os.path.join(outputDir, incF))

def modproteinItp(topology, outputDir, verbose = False):
    """ Modify user supplyied protein itp files in working directiry. """

    with open(topology) as topf:
        topologyLines = topf.readlines()

        itpFiles = []
        for line in topologyLines:
            if '#include' in line and not '.ff/' in line:
                incF = line.split('"')[1]
                if verbose:
                    print('include file:', incF)
                if os.path.exists(incF):
                    itpFiles.append(incF)

    currDir = os.getcwd()
    os.chdir(outputDir)

    if len(itpFiles) > 0:
        for itpFile in itpFiles:
            if not 'posre_' in itpFile:
                with open(itpFile) as f:
                    lines = f.readlines()

                with open(itpFile, 'w') as f:
                    for line in lines:
                        if line.strip() == '#ifdef POSRES':
                            #f.write('#ifdef POSRES || POSRES_PROT\n')
                            f.write('#ifdef POSRES_PROT\n')
                        else:
                            f.write(line)

    os.chdir(currDir)

def splitTopologyToItp(topology, verbose = False):
    """ Remove molecule specific information from a topology file to an itp
    file, which in turn is included by the topology file. Returns the name
    of the itp file."""

    with open(topology) as topf:
        topologyLines = topf.readlines()

    if verbose:
        print('topology: ', topology) 
        print('Splitting molecule information from topology to itp file')

    itpFile = os.path.splitext(topology)[0] + '.itp'
    if verbose:
        print('itpFile: ', itpFile) 
    if os.path.exists(itpFile):
        if verbose:
            print('Itp file already exists. Aborting.')
        return itpFile

    itpHeaders = ['atomtypes', 'moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals', 'cmap']

    toItp = False
    hasItp = False

    with open(topology, 'w') as topf:
        with open(itpFile, 'w') as itpf:
            for line in topologyLines:
                if line.startswith('['):
                    if any(itpHeader in line for itpHeader in itpHeaders):
                        toItp = True
                    else:
                        toItp = False
                if toItp:
                    if not hasItp:
                        hasItp = True
                        # Include the itp file in the topology file, but only once
                        topf.write('#include "%s"\n' % os.path.basename(itpFile))
                    itpf.write(line)
                else:
                    topf.write(line)

    if not hasItp:
        os.remove(itpFile)
        itpFile = None
        if verbose:
            print('No data to write to itp file.')

    return itpFile

def mergeTopologyFiles(ligandTop, proteinTop, outTop, verbose = False):
    """ Add information from the file proteinTop to the file ligandTop. """

    if verbose:
        print('ligandTop:', ligandTop)
        print('proteinTop:', proteinTop)
        print('OutTop:', outTop)

    with open(proteinTop) as f:
        proteinLines = f.readlines()

    with open(ligandTop) as f:
        ligandLines = f.readlines()

    proteinIncludeLines = []
    proteinSysName = ''
    proteinMoleculesLines = []

    outputDir = os.path.dirname(ligandTop)

    inSystem = 0
    inMolecules = 0

    for line in proteinLines:
        if inSystem:
            if inSystem == 1:
                inSystem += 1
                continue
            proteinSysName = line
            inSystem = 0
            continue

        if inMolecules:
            if inMolecules == 1:
                inMolecules += 1
                continue
            if line.startswith('['):
                inMolecules = 0
                continue
            proteinMoleculesLines.append(line)

        if '#include' in line and not '.ff/' in line:
            incF = line.split('"')[1]

            if line in ligandLines:
                continue

            #"if not os.path.isabs(incF):
            #    incF = os.path.join(outputDir, incF)

            line = '#include "%s"\n' % incF
            if line in ligandLines:
                continue

            proteinIncludeLines.append(line)

        if line.startswith('[ system ]') or line.startswith('[system]'):
            inSystem = 1
            continue

        if line.startswith('[ molecules ]') or line.startswith('[molecules]'):
            inMolecules = 1
            continue

    if verbose:
        if proteinIncludeLines:
            print('Copying the following include lines from %s to %s' % (proteinTop, outTop))
            for line in proteinIncludeLines:
                print(line,)
        if proteinMoleculesLines:
            print('Copying the following molecule lines from %s to %s' % (proteinTop, outTop))
            for line in proteinMoleculesLines:
                print(line,)

    with open(outTop, 'w') as f:
        inIncludes = 0
        inMolecules = 0
        for line in ligandLines:
            if inSystem:
                f.write(line.strip() + ' + '+ proteinSysName)
                inSystem = False
                continue

            if not inIncludes and '#include' in line:
                inIncludes = 1

            if line.startswith('[ system ]') or line.startswith('[system]'):
                inSystem = True

            if inIncludes == 1 and line.startswith('['):
                for includeLine in proteinIncludeLines:
                    f.write(includeLine)
                f.write('\n')
                inIncludes = 2

            f.write(line)

        for moleculeLine in proteinMoleculesLines:
            f.write(moleculeLine)

def mergeCoordinateFiles(ligandCoords, proteinCoords, outputCoords = None):
    """ Append one coordinate file to another (proteinCoords after ligandCoords).
    If outputCoords is not specified the results will be written to the
    ligandCoords file """

    if outputCoords == None:
        outputCoords = ligandCoords

    with open(proteinCoords) as f:
        proteinLines = f.readlines()

    with open(ligandCoords) as f:
        ligandLines = f.readlines()

    title = ligandLines[0].strip() + ' + ' + proteinLines[0].strip()
    nAtoms = int(ligandLines[1]) + int(proteinLines[1])

    with open(outputCoords, 'w') as f:
        f.write(title + '\n')
        f.write('%d\n' % nAtoms)

        for line in ligandLines[2:-1] + proteinLines[2:]:
            f.write(line)

def makeIndexRun(outputDir, outputFileBaseName, verbose = False):
    """ Run make_ndx to create an index file """

    indexFileName = os.path.join(outputDir, 'index.ndx')
    # Try to find a suitable coordinate file to use
    fileAlts = [outputFileBaseName + '_solvated_ionised.gro', outputFileBaseName + '_solvated.gro',
                outputFileBaseName + '_box.gro', outputFileBaseName + '.gro']

    for alt in fileAlts:
        coordsFile = os.path.join(outputDir, alt)
        if verbose:
            print('Searching for %s as .gro file for generating index.' % coordsFile)
        if os.path.exists(coordsFile):
            if verbose:
                print('Found gro file %s.' % coordsFile)
            break
    else:
        if verbose:
            print('Cannot find .gro file for generating index file. Searching in parent directory.')
        for alt in fileAlts:
            coordsFile = os.path.join(outputDir, '..', alt)
            if os.path.exists(coordsFile):
                break
        else:
            print('Cannot find .gro coordinate file for generating index file.')
            return

    if verbose:
        print('Making index file from %s' % coordsFile)
        print('indexFileName: ', indexFileName)

    # Make new groups:
    # 1) the whole system except the ligand
    command = '!2\n'
    # 2) The ligand without hydrogen
    command += '2 & ! a H*\n'
    # 3) One group with everything except water and ions (if present) for temp coupling
    # 4) One new group with water and ions (if present) for temp coupling
    if 'ionised' in coordsFile:
        command += '!"Water_and_ions"\n"Water_and_ions"\nq\n'
    elif 'solvated' in coordsFile:
        command += '!"Water"\n"Water"\nq\n'
    else:
        command += 'q\n'

    makeIndexCmd = ['gmx'+gmxSuffix, 'make_ndx', '-f', coordsFile, '-o', indexFileName]

    if verbose:
        print(' '.join(makeIndexCmd))

    try:
        result = getCommandOutput(makeIndexCmd, command)
        if verbose:
            print(result)

        if 'ionised' in coordsFile:
            groupA = '!Water_and_ions'
            groupB = 'Water_and_ions'
        else:
            groupA = '!Water'
            groupB = 'Water'
        groupC = 'LIG_&_!H*'
        groupD = '[ LIG ]'

        with open(indexFileName) as f:
            lines = f.readlines()

            # Rename the newly created temp coupling groups to get consistent names.
            if 'solvated' in coordsFile:

                foundA = False
                foundB = False
                for i in reversed(range(len(lines))):
                    line = lines[i]
                    if not foundA and groupA in line:
                        lines[i] = '[ TempCouplingA ]\n'
                        foundA = True
                        if foundB:
                            break
                    elif not foundB and groupB in line:
                        lines[i] = '[ TempCouplingB ]\n'
                        foundB = True
                        if foundA:
                            break

            # Rename the newly created ligand without hydrogen group to get consistent names.
            foundC = False
            for i in reversed(range(len(lines))):
                if not foundC and groupC in lines[i]:
                    lines[i] = '[ LIG_Heavy ]\n'
                    foundC = True
                    break

            # Rename the duplicated ligand group to unique names.
            countD = 0
            for i in range(len(lines)):
                if groupD in lines[i]:
                    if countD > 0:
                        lines[i] = '[ LIG_Dup{} ]\n'.format(countD)
                    countD += 1

            with open(indexFileName, 'w') as f:
                for line in lines:
                    f.write(line)

    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(makeIndexCmd))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')


def makeRestraintsRun(coordsFile, name = None, restrainHydrogens = False, verbose = False):
    """ run genrestr to create restraints for the ligand """

    outputDir = os.path.dirname(coordsFile)
    coordsBaseName = os.path.basename(coordsFile)
    if name:
        restrFile = os.path.join(outputDir, 'posre_' + name + '.itp')
    else:
        restrFile = os.path.join(outputDir, 'posre_' + os.path.splitext(coordsBaseName)[0] + '.itp')
    makeRestraintsCmd = ['gmx'+gmxSuffix, 'genrestr', '-f', coordsFile, '-o', restrFile]

    if verbose:
        print(' '.join(makeRestraintsCmd))

    try:
        # Restrain the ligand
        result = getCommandOutput(makeRestraintsCmd, '2\n')
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(makeRestraintsCmd))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    if not restrainHydrogens:
        try:
            fixRestraintsFile(restrFile, coordsFile)
        except IOError:
            print('Error!' + restrFile + ' not found. Trying again.')
            try:
                result = getCommandOutput(makeRestraintsCmd, '2\n')
                fixRestraintsFile(restrFile, coordsFile)
            except Exception:
                print('Still failing generating restraints. Trying to continue anyhow.')

def fixRestraintsFile(restrFile, groFile = None):
    """ Remove hydrogens from the restraints file """

    with open(restrFile) as rF:
        restrLines = rF.readlines()

    with open(groFile) as gF:
        groLines = gF.readlines()

    with open(restrFile, 'w') as rF:
        for line in restrLines:
            parts = line.split()
            if len(parts) == 5:
                if parts[0].isdigit():
                    nr = int(parts[0])

                    # Only look for a specific line and compare the
                    # atom number. Do not specifically search for the right
                    # atom in the .gro file.
                    if len(groLines) - 1 >= nr:
                        groLine = groLines[nr+1]
                        if len(groLine) > 20:
                            groNr = groLine[15:20].strip()
                            if groNr.isdigit():
                                name = groLine[10:15].strip()
                                groNr = int(groNr)
                                if groNr == nr and 'H' in name:
                                    continue
            rF.write(line)

def assignChargeGroups(mol2File):

    with open(mol2File) as f:
        contents = f.readlines()

    molsys = MolSystem()
    molsys.readMol2(contents)

    #FIXME: Finish this

def runPreMinimization(outputDir, outputFileBaseName, verbose = False):

    currDir = os.getcwd()
    os.chdir(outputDir)

    tprFileName = 'premin.tpr'
    topologyFileName = outputFileBaseName + '.top'
    groFileName = outputFileBaseName + '.gro'
    if not os.path.exists(groFileName):
        groFileName = os.path.join('..', outputFileBaseName + '.gro')
    #templateDir = os.path.join(os.path.dirname(__file__), '..', 'templates')
    templateDir = os.path.join(os.path.dirname(__file__), 'templates')
    templateFileName = os.path.join(templateDir, 'premin_template.mdp')

    gromppCommand = ['gmx'+gmxSuffix, 'grompp', '-f', templateFileName, '-p', topologyFileName,
                     '-c', groFileName, '-r', groFileName, '-o', tprFileName, '-maxwarn', '10']

    if verbose:
        print('Generating GROMACS input for early energy minimization of molecule')

    try:
        result = getCommandOutput(gromppCommand)
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(gromppCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')
        os.chdir(currDir)
        return

    mdrunCommand = ['gmx'+gmxSuffix, 'mdrun', '-deffnm', 'premin']
    try:
        result = getCommandOutput(mdrunCommand)
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(mdrunCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')
        os.chdir(currDir)
        return

    if os.path.exists('premin.gro'):
        shutil.move('premin.gro', outputFileBaseName + '.gro')

    os.chdir(currDir)

def solvateSystem(groFile, outputDir, outputFileBaseName, solvent, boxType = 'dodecahedron', 
                  bufferDist = 1.0, verbose = False):
    """ Solvate the system in a "box" of the specified solvent. The box shape is specified by boxType
    (default dodecahedron).
    The minimum distance from molecules in the system to the edges of the "box" is specified by bufferDist
    (default 1.0 nm). """

    topologyFileName = os.path.join(outputDir, outputFileBaseName + '.top')

    boxFileName = os.path.join(outputDir, outputFileBaseName + '_box.gro')
    solvatedFileName = os.path.join(outputDir, outputFileBaseName + '_solvated.gro')

    editConfCommand = ['gmx'+gmxSuffix, 'editconf', '-f', groFile, '-o', boxFileName,
    '-bt', boxType, '-d', '%.2f' % bufferDist]

    if verbose:
        print('Preparing solvent box')

    try:
        result = getCommandOutput(editConfCommand)
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(editConfCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    if solvent in ['tip4p', 'tip5p']:
        solventFile = solvent + '.gro'
    elif 'tip4p' in solvent or 'opc' in solvent:
        solventFile = 'tip4p.gro'
    else:
        solventFile = 'spc216.gro'

    genboxCommand = ['gmx'+gmxSuffix, 'solvate', '-cp', boxFileName, '-cs', solventFile,
    '-o', solvatedFileName, '-p', topologyFileName]

    if verbose:
        print(' '.join(genboxCommand))

    try:
        result = getCommandOutput(genboxCommand)
        if verbose:
            print(result)

        # If there was no output try the gromacs 4.x way of solvating
        if not os.path.exists(solvatedFileName) or 'This tool is not avaiable before Gromacs 5.0' in result:
            if verbose:
                print("No output from genbox, trying 'gmx genbox' instead.")
            genboxCommand = ['gmx'+gmxSuffix, 'genbox', '-cp', boxFileName, '-cs', solventFile,
                             '-o', solvatedFileName, '-p', topologyFileName]
            result = getCommandOutput(genboxCommand)
            if verbose:
                print(result)

    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(genboxCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    # For some reason a temp.top file is sometimes created where GROMMFAWrapper
    # is launched and the updated topology file is empty. Move the temp.top file
    # to replace the empty file.
    if os.path.exists('temp.top') and abs(os.path.getmtime('temp.top') - time.time()) < 5:
        shutil.move('temp.top', topologyFileName)

#def neutraliseSystem(outputDir, outputFileBaseName, systemNetCharge, pname = 'NA', nname = 'CL',
def neutraliseSystem(outputDir, outputFileBaseName, conc = 0.0, pname = 'NA', nname = 'CL',
                     verbose = False):
    """ Neutralise a solvated system by adding ions to counter the net charge """

    topologyFileName = os.path.join(outputDir, outputFileBaseName + '.top')
    #templateDir = os.path.join(os.path.dirname(__file__), '..', 'templates')
    templateDir = os.path.join(os.path.dirname(__file__), 'templates')
    templateFileName = os.path.join(templateDir, 'genion_template.mdp')

    #if verbose:
    #    print('Net charge: %d. Neutralising system by adding ions.' % systemNetCharge)

    fileNameAlts = [outputFileBaseName + '_solvated', outputFileBaseName]
    for fileName in fileNameAlts:
        solvatedFileName = os.path.join(outputDir, fileName + '.gro')
        if os.path.exists(solvatedFileName):
            break
    else:
        print('Cannot find input file for adding ions.')
        return

    ionisedFileName = os.path.join(outputDir, fileName + '_ionised.gro')
    tprFileName = os.path.join(outputDir, 'genion_input.tpr')

    gromppCommand = ['gmx'+gmxSuffix, 'grompp', '-f', templateFileName, '-c', solvatedFileName,
                     '-r', solvatedFileName, '-p', topologyFileName, '-o', tprFileName,
                     '-po', os.path.join(outputDir, 'mdout.mdp'), '-maxwarn', '10']

    if verbose:
        print('Generating ions')
        print(' '.join(gromppCommand))

    try:
        result = getCommandOutput(gromppCommand)
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(gromppCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    # For some reason a temp.top file is sometimes created where GROMMFAWrapper
    # is launched and the updated topology file is empty. Move the temp.top file
    # to replace the empty file.
    if os.path.exists('temp.top') and abs(os.path.getmtime('temp.top') - time.time()) < 5:
        shutil.move('temp.top', topologyFileName)

    genionCommand = ['gmx'+gmxSuffix, 'genion', '-s', tprFileName, '-p', topologyFileName, '-o',
                     ionisedFileName, '-conc', '%f' % conc, '-neutral', '-pname', pname, '-nname', nname]
    #                 ionisedFileName]
    #if systemNetCharge < 0:
    #    genionCommand += ['-np', '%d' % abs(round(systemNetCharge)), '-pname', pname]
    #else:
    #    genionCommand += ['-nn', '%d' % abs(round(systemNetCharge)), '-nname', nname]

    if verbose:
        print(' '.join(genionCommand))

    try:
        result = getCommandOutput(genionCommand, 'SOL')
        if verbose:
            print(result)
    except subprocess.CalledProcessError as e:
        print('Failed running', ' '.join(genionCommand))
        try:
            print(result, '\n')
        except Exception:
            pass
        print(e)
        print('Trying to continue anyhow.')

    # For some reason a temp.top file is sometimes created where GROMMFAWrapper
    # is launched and the updated topology file is empty. Move the temp.top file
    # to replace the empty file.
    if os.path.exists('temp.top') and abs(os.path.getmtime('temp.top') - time.time()) < 5:
        shutil.move('temp.top', topologyFileName)

    if os.path.exists(tprFileName): os.remove(tprFileName)
    mdpFileName = os.path.join(outputDir, 'mdout.mdp')
    if os.path.exists(mdpFileName): os.remove(mdpFileName)


def getForceFieldCalibration(calibrationFileName, forceFieldName, solvent, chargeModel):

    cal = dict()

    with open(calibrationFileName) as f:
        lines = f.readlines()

    i = 0
    nLines = len(lines)

    fName = None

    while i < nLines:
        line = lines[i]
        if len(line) < 1 or line[0] == '#' or line[0] == ';':
            i += 1
            continue
        if 'forcefield = ' in line.lower():
            parts = line.split('=')
            value = parts[1].strip().lower()
            if forceFieldName and value != forceFieldName:
                while i < nLines and lines[i] != '***':
                    i += 1
                fName = None
                continue
            fName = value
        elif 'solvent = ' in line.lower():
            parts = line.split('=')
            value = parts[1].strip().lower()
            if solvent and value != solvent:
                while i < nLines and lines[i] != '***':
                    i += 1
                fName = None
                continue
        elif 'charges = ' in line.lower():
            parts = line.split('=')
            value = parts[1].strip().lower()
            if chargeModel and value != chargeModel:
                while i < nLines and lines[i] != '***':
                    i += 1
                fName = None
                continue

        elif 'overallfactors = ' in line.lower():
            parts = line.split()
            if len(parts) >= 4:
                cal['overallfactors'] = [float(parts[2]), float(parts[3])]

        elif fName:
            parts = line.split()
            if len(parts) >= 3:
                cal[parts[0]] = [float(parts[1]), float(parts[2])]

        i += 1

    return cal


def calibrateVdW(outputDir, calibrationFileName, forceFieldName, solvent = None,
                 chargeModel = None, verbose = False):

    outputFileBaseName = os.path.basename(outputDir)
    ffDir = outputDir + '_%s' % forceFieldName
    itpFileName = os.path.join(ffDir, outputFileBaseName + '.itp')

    if forceFieldName:
        forceFieldName = forceFieldName.lower()
    if solvent:
        solvent = solvent.lower()
    if chargeModel:
        chargeModel = chargeModel.lower()

    atomCalibration = getForceFieldCalibration(calibrationFileName, forceFieldName, solvent, chargeModel)

    if not atomCalibration:
        return

    with open(itpFileName) as f:
        lines = f.readlines()

    inAtomTypes = False
    nLines = len(lines)
    for i in range(nLines):
        line = lines[i]
        if len(line) == 0 or line[0] == ';':
            continue

        if inAtomTypes:
            if line[0] == '[':
                break
            firstSplit = line.split(';')
            parts = firstSplit[0].strip().split()
            if len(parts) < 7:
                continue
            atomType = parts[0]
            sigma = float(parts[5])
            epsilon = float(parts[6])

            f = atomCalibration.get(atomType) or atomCalibration.get('overallfactors')
            if f and len(f) >= 2:
                #sigma *= (1 + 0.02 * (1 - f))
                #epsilon *= f
                sigma *= f[0]
                epsilon *= f[1]

                line = ' %-7s  %-11s %7s  %7s   %-4s  %.5e   %.5e' % (atomType, parts[1], parts[2],
                                                                parts[3], parts[4], sigma, epsilon)
                if len(firstSplit) > 1:
                    firstSplit[1] = firstSplit[1].strip()
                    firstSplit[0] = line + ' ;'
                    line = ' '.join(firstSplit)
                    line += '\n'
                lines[i] = line

        if '[ atomtypes ]' in line or '[atomtypes]' in line:
            inAtomTypes = True

    with open(itpFileName, 'w') as f:
        for line in lines:
            f.write(line)

def convertGmx2Amb(outputDir, outputFileBaseName, verbose = False):
    """ Run make_ndx to create an index file """
    import parmed as pmd

    topFile = os.path.join(outputDir, outputFileBaseName + '.top')
    # Try to find a suitable coordinate file to use
    fileAlts = [outputFileBaseName + '_solvated_ionised.gro', outputFileBaseName + '_solvated.gro',
                outputFileBaseName + '_box.gro', outputFileBaseName + '.gro']

    for alt in fileAlts:
        coordsFile = os.path.join(outputDir, alt)
        if verbose:
            print('Searching for %s as .gro file for generating index.' % coordsFile)
        if os.path.exists(coordsFile):
            if verbose:
                print('Found gro file %s.' % coordsFile)
            break
    else:
        if verbose:
            print('Cannot find .gro file for generating index file. Searching in parent directory.')
        for alt in fileAlts:
            coordsFile = os.path.join(outputDir, '..', alt)
            if os.path.exists(coordsFile):
                break
        else:
            print('Cannot find .gro coordinate file for generating index file.')
            return

    if verbose:
        print('Coordinate file from %s' % coordsFile)
        print('Topology file from %s' % topFile)

    gmx_top = pmd.load_file(topFile, xyz=coordsFile)
    # Write Amber coordinate and topology files
    coordsFile_save = os.path.splitext(coordsFile)[0] + '.inpcrd'
    topFile_save = os.path.splitext(topFile)[0] + '.prmtop'
    print('Coordinate file to %s' % coordsFile_save)
    print('Topology file to  %s' % topFile_save)
    gmx_top.save(coordsFile_save, format='rst7', overwrite=True)
    gmx_top.save(topFile_save, format='amber', overwrite=True)

    # Write Charmm coordinate and topology files
    #coordsFile_save = os.path.splitext(coordsFile)[0] + '.crd'
    #topFile_save = os.path.splitext(topFile)[0] + '.psf'
    #print('Coordinate file to %s' % coordsFile_save)
    #print('Topology file to  %s' % topFile_save)
    #gmx_top.save(coordsFile_save, format='charmmcrd', overwrite=True)
    #gmx_top.save(topFile_save, format='psf', overwrite=True)
