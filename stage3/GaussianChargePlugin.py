"""
    This file is part of STaGE3 - a wrapper for generating
    GROMACS topology files.

    Written by Michio Katouda

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

import os
import shutil
from glob import glob

from .ChargePlugin import ChargePlugin
from .molScripts.molecule import MolSystem
from .util import isInPath, getCommandOutput, runGaussian, subprocess

if not isInPath('antechamber'):
    raise ImportError('antechamber not available. Disabling charge plugin.')

class GaussianChargePlugin(ChargePlugin):

    def __init__(self):

        self.programName = 'g16'
        self.alternativeChargeMethods = {'hf': 'Hatree-Fock/6-31G(d) basis set followed by RESP'}
        self.order = 3

    def generateCharges(self, inFile, method, netCharge = None, multiplier = 1.0, verbose = False):

        if verbose:
            print('Generating %s charges.' % method)

        netCharge = int(round(netCharge,0))

        currDir = os.getcwd()
        molDir = os.path.dirname(inFile)
        os.chdir(molDir)

        molName = os.path.splitext(os.path.basename(inFile))[0]
        gaussianOutputFileName = molName + '_charge.log'
        if not os.path.isfile(gaussianOutputFileName):
            if verbose:
                print('Running GAUSSIAN')
            if not runGaussian(molName, netCharge, int(round(multiplier,0)), '_charge', onlyGenerateRunFiles = False, verbose = verbose):
                print('GAUSSIAN failed. Trying to continue anyhow.')

        tempMol2File = molName + '_charge_temp.mol2'

        antechamberCommand = ['antechamber', '-fi', 'gout', '-i', gaussianOutputFileName, 
                              '-fo', 'mol2', '-o', tempMol2File, '-c', 'resp', '-at', 'sybyl', '-pf', 'y']
        try:
            if verbose:
                print(' '.join(antechamberCommand))

            result = getCommandOutput(antechamberCommand)
            if verbose:
                print(result)

        except subprocess.CalledProcessError as e:
            print('Failed running', ' '.join(antechamberCommand))
            try:
                print(result, '\n')
            except Exception:
                pass
            print(e)
            print('Cannot continue without converting molecule from ac to mol2.')
            return

        molSys = MolSystem()
        with open(inFile) as f:
            molSys.readMol2(f.readlines())

        mol = molSys.molecules[0]
        atoms = sorted(mol.atoms, key = lambda a:a.number)

        tempMolSys = MolSystem()
        with open(tempMol2File) as f:
            tempMolSys.readMol2(f.readlines())

        tempMol = tempMolSys.molecules[0]
        tempAtoms = sorted(tempMol.atoms, key = lambda a:a.number)

        if len(atoms) != len(tempAtoms):
            print('Incompatible molecular systems after charge calculations')
            return

        for i, atom in enumerate(atoms):
            tempAtom = tempAtoms[i]
            atom.charge = tempAtom.charge

        if abs(multiplier - 1.0) > 0.001:
            molSys.multiplyCharges(multiplier)

        with open(inFile, 'w') as f:
            molSys.writeMol2(f)

        os.remove(tempMol2File)
        #if os.path.isfile('punch'):
        #    os.remove('punch')
        #if os.path.isfile('esout'):
        #    os.remove('esout')
        #if os.path.isfile('QOUT'):
        #    os.remove('QOUT')

        os.chdir(currDir)
