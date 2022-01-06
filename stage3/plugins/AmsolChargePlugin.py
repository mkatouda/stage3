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
from util import isInPath, getCommandOutput

if not isInPath('amsol7.1'):
    raise ImportError('amsol7.1 not available. Disabling charge plugin.')

class AmsolChargePlugin(ChargePlugin):

    def __init__(self):

        self.programName = 'amsol7.1'
        self.alternativeChargeMethods = {'cm1a': 'CM1A',
                                         'cm3a': 'CM3A',
                                         'sm54a': 'AM1-SM5.4 is solvation model SM5 using class IV charges from AM1'}
        self.order = 2

    def generateCharges(self, inFile, method, netCharge = None, multiplier = 1.0, verbose = False):

        if verbose:
            print('Generating %s charges.' % method)

        netCharge = int(round(netCharge,0))

        shutil.copy(inFile, inFile + '.bak')

        method = method[:3].upper()

        molSys = MolSystem()
        with open(inFile, 'r') as f:
            molSys.readMol2(f.readlines())
        amsolFile = os.path.splitext(inFile)[0] + '.amin'
        with open(amsolFile, 'w') as f:
            if method == 'SM5':
                molSys.writeAmsol(f, 'SM5.4A AM1 CHARGE=%.1f CART SOLVNT=WATER TRUES HFCALC=OPT CYCLES=500 BFGS KICK=2' % float(netCharge))
                method = 'CM1'
            else:
                molSys.writeAmsol(f, 'AM1 %s CHARGE=%.1f CART' % (method, float(netCharge)))
        with open(amsolFile, 'r') as f:
            amsolInput = f.readlines()

        command = [self.programName]

        if verbose:
            print(' '.join(command))

        amsolOutput = getCommandOutput(command, ''.join(amsolInput))

        # If AMSOL failed it is worth trying reshuffling the first and the last atom. The problem seems to arise if there
        # is a 180 degrees angle starting from the first atom.
        if 'CALCULATION ABANDONED' in amsolOutput:
            if verbose:
                print('Amsol calculations failed. Fixing input file and trying again.')
            if 'STRAIGHT LINE' in amsolOutput:
                amsolFail = 1
                lastLine = amsolInput[-1]
                amsolInput[-1] = amsolInput[3]
                amsolInput[3] = lastLine
            amsolOutput = getCommandOutput(command, ''.join(amsolInput))
            if 'CALCULATION ABANDONED' not in amsolOutput:
                if verbose:
                    print('Second Amsol attempt succeeded.')
            else:
                print('Second Amsol attempt failed. Trying to continue anyhow.')
        else:
            amsolFail = None

        amsolFile = os.path.splitext(inFile)[0] + '.amout'
        with open(amsolFile, 'w') as f:
            f.write(amsolOutput)

        if verbose:
            print('Adding calculated charges to atoms')
        molSys.addAmsolChargesToAtoms(amsolOutput.split('\n'), method)

        # If the atoms were reshuffled above make sure the charges get added to the right atom.
        if amsolFail == 1:
            mol = molSys.molecules[0]
            atoms = sorted(mol.atoms, key = lambda a:a.number)
            lastAtomCharge = atoms[-1].charge
            atoms[-1].charge = atoms[0].charge
            atoms[0].charge = lastAtomCharge

        if abs(multiplier - 1.0) > 0.001:
            molSys.multiplyCharges(multiplier)

        if verbose:
            print('Writing', inFile)

        with open(inFile, 'w') as f:
            molSys.writeMol2(f)

        d = os.path.dirname(inFile)
        rmFiles = glob(os.path.join(d, 'fort.*'))

        if verbose and rmFiles:
            print('Cleaning up')

        for f in rmFiles:
            os.remove(f)
