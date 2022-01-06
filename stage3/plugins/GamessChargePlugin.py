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
from glob import glob

from ChargePlugin import ChargePlugin
from molScripts.molecule import MolSystem
from util import isInPath, getCommandOutput, runGamess, runGmstoresp, subprocess

if not isInPath('rungms') or not isInPath('antechamber'):
    raise ImportError('rungms or antechamber not available. Disabling charge plugin.')

class GamessChargePlugin(ChargePlugin):

    def __init__(self):

        self.programName = 'rungms'
        self.alternativeChargeMethods = {'b3lyp/pcm': 'B3LYP with PCM and cc-pV(T+d)Z basis set followed by RESP'}
        self.order = 3

    def generateCharges(self, inFile, method, netCharge = None, multiplier = 1.0, verbose = False):

        if verbose:
            print('Generating %s charges.' % method)

        netCharge = int(round(netCharge,0))

        currDir = os.getcwd()
        molDir = os.path.dirname(inFile)
        os.chdir(molDir)

        molName = os.path.splitext(os.path.basename(inFile))[0]
        gamessOutputFileName = molName + '_charge.gamout'
        gmstorespOutputFileName = molName + '_charge.ac'
        if not os.path.isfile(gamessOutputFileName):
            if verbose:
                print('Running GAMESS')
            if not runGamess(molName, netCharge, '_charge', onlyGenerateRunFiles = False, verbose = verbose):
                print('GAMESS failed. Trying to continue anyhow.')
        if verbose:
            print('Running makeresp and resp')
        if not runGmstoresp(gamessOutputFileName, verbose):
            print('makeresp or resp failed. Trying to continue anyhow.')

        tempMol2File = molName + '_charge_temp.mol2'

        antechamberCommand = ['antechamber', '-fi', 'ac', '-i', gmstorespOutputFileName, '-fo', 'mol2', '-o', tempMol2File]

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

        os.chdir(currDir)
