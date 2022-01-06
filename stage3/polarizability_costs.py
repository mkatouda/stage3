#!/usr/bin/env python

"""
    A script to calculate polarization costs of a small molecule in a force
    field with fixed charges.

    Dipole moments in gas phase are obtained by generating an optimized structure
    using a B3LYP level of theory and a cc-pV(T+d)Z basis set. The dipole polarizability
    is computed by finite difference of energies and multipole moments
    each with different directions of small finite perturbing electrostatic
    potentials (0.001 au) or potential gradients (0.001 au).

    Calculations are performed according to eq 4 in
    Swope, W. C.; Horn, H. W.; Rice, J. E. Accounting for Polarization Cost When
    Using Fixed Charge Force Fields. II. Method and Application for Computing
    Effect of Polarization Cost on Free Energy of Hydration. J. Phys. Chem. B 2010,
    114: 8631-8645.


    Requires Open Babel and GAMESS.

                      VERSION 0.9

    Written by Magnus Lundborg
    Copyright (c) 2013-2014, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.

    The following applications must be available in the list of directories of
    executable programs (the PATH environment variable):

    rungms
    babel

    Requirements: Python 2.6 <= version < 3.0
                  Open Babel
                  GAMESS


    GAMESS
    Schmidt, M. W.; Baldridge, K. K.; Boatz, J. A.; Elbert, S. T.; Gordon,
    M. S.; Jensen, J. H.; Koseki, S.; Matsunaga, N.; Nguyen, K. A.; Su, S.;
    Windus, T. L.; Dupuis, M.; Montgomery, J. A. General atomic and
    molecular electronic structure system. J. Comput. Chem. 1993,
    14: 1347-1363.

    Open Babel
    O' Boyle, M. O., Banck, M., James, C. A., Morley, C., Vandermeersch, T.,
    Hutchison, G. R. Open Babel: An open checmical toolbox.
    J. Cheminform. 2011, 3: 33.

"""

import os
import argparse
import traceback
import shutil

from util import molecule2GamessDipole, extractGamessDipoleData, calculateMol2Dipole
from util import babelConvert, itpChargeToMol2
from numpy import transpose, power, sum
from math import sqrt

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A script to calculate '
    'polarization costs of a small molecule in a force field with fixed charges.\n\n'
    'Dipole moments in gas phase are obtained by generating an optimized structure using a B3LYP '
    'level of theory and a cc-pV(T+d)Z basis set. The dipole polarizability '
    'is thereafter computed by finite difference of energies and multipole moments '
    'each with different directions of small finite perturbing electrostatic '
    'potentials (0.001 au) or potential gradients (0.001 au).\n'
    'The gas phase optimized structure is used for calculating the dipole moment of the chosen '
    'charge model (by providing a Gromacs itp file).\n\n'
    'Calculations are performed according to eq 4 in \n'
    'Swope, W. C.; Horn, H. W.; Rice, J. E. Accounting for Polarization Cost When '
    'Using Fixed Charge Force Fields. II. Method and Application for Computing '
    'Effect of Polarization Cost on Free Energy of Hydration. J. Phys. Chem. B 2010, '
    '114: 8631-8645.')


    parser.add_argument('-g', '--gamess_dipole',
                        help = 'File name of GAMESS output with dipole and polarizability. '
                        'Generated if not specified')
    parser.add_argument('-m', '--molecule',
                        help = 'Molecule file (preferably .mol2, but can be in other formats)')
    parser.add_argument('-o', '--gamess_output_dir',
                        help = 'Directory to save GAMESS output files (if they are not provided).')
    parser.add_argument('-i', '--itp',
                        help = 'Gromacs itp file (containing partial charges in solvent). If an itp file is '
                        'not provided only GAMESS calculations will be performed.')
    parser.add_argument('-p', '--ph',
                        type = float,
                        help = 'Protonate the molecule according to this pH (float). '
                               'This does not always give correct results. It is safer '
                               'to provide correctly protonated input files.')
    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        help = 'Verbose output')

    args = parser.parse_args()

    if not args.gamess_dipole and not args.molecule:
        parser.error('If a GAMESS output file (--gamess_dipole) is not specified '
                     'a molecule file (--molecule) must be specified.')

    if (args.itp and ',' in args.itp) or (args.gamess_dipole and ',' in args.gamess_dipole) or \
        (args.molecule and ',' in args.molecule):
        parser.error('File names must not contain a comma (",")')

    itpFile = args.itp
    if itpFile and not os.path.isfile(itpFile):
        print('Cannot find Gromacs itp file')
        exit(1)


    if args.gamess_dipole:
        gamessFile = args.gamess_dipole
    else:
        if not os.path.isfile(args.molecule):
            print('Cannot find molecule file')
            exit(1)

        iFormat = os.path.splitext(args.molecule)[1].lower()
        if iFormat != 'mol2' or args.ph:
            mol2File = '.'.join(args.molecule.split('.')[:-1]) + '.mol2'
            inMolecules = babelConvert(inputFile = args.molecule, outputFile =
                                       mol2File, pH = args.ph,
                                       verbose = args.verbose)
        else:
            mol2File = args.molecule

        gamess_output_dir = args.gamess_output_dir or os.path.dirname(mol2File)

        # Perform gas optimization and polarization cost calculations using GAMESS.
        if not os.path.isdir(gamess_output_dir):
            if args.verbose:
                print('Cannot find directory to save GAMESS output. Creating it.')
            try:
                os.mkdir(gamess_output_dir)
            except OSError:
                print('Cannot create GAMESS output directory.')
                exit(1)
        if not os.path.isfile(os.path.join(gamess_output_dir, os.path.basename(mol2File))):
            if args.verbose:
                print('Copying mol2 file to GAMESS output dir.')
            shutil.copy(mol2File, gamess_output_dir)

        gamessFile = molecule2GamessDipole(mol2File, gamess_output_dir, verbose = args.verbose)
        if not gamessFile:
            print('Could not run GAMESS.')
            exit(1)

    if not itpFile:
        exit(0)

    # Generate a mol2 file from the optimized structure and apply the charges from the itp file to the
    # molecule.
    mol2File = os.path.splitext(gamessFile)[0]) + '.mol2'
    babelConvert(gamessFile, mol2File, verbose = args.verbose)
    itpChargeToMol2(itpFile, mol2File)

    gasDipole, alphaPol, dipoleCenter = extractGamessDipoleData(gamessFile)

    solventDipole = calculateMol2Dipole(mol2File, dipoleCenter)
    gasDipoleTot = power(gasDipole, 2)
    gasDipoleTot = sqrt(gasDipoleTot[0,0] + gasDipoleTot[0,1] + gasDipoleTot[0,2])
    solventDipoleTot = power(solventDipole, 2)
    solventDipoleTot = sqrt(solventDipoleTot[0,0] + solventDipoleTot[0,1] + solventDipoleTot[0,2])
    # 1D = 1e-21/299792458 C * m
    convConst = 1e-21/299792458
    gasDipoleConv = convConst * gasDipole
    solventDipoleConv = convConst *solventDipole

    gasDipoleTotConv = convConst * gasDipoleTot
    solventDipoleTotConv = convConst * solventDipoleTot
    #print('Alpha polarizability:', alphaPol)
    # Atomic unit of polarizability is 1.6487772754e-41 C^2 * m^2 * J^-1
    alphaPol *= 1.6487772754e-41
    polarizability = (alphaPol[0,0] + alphaPol[1,1] + alphaPol[2,2])/3
    if args.verbose:
        print('Gas dipole: %e D = %e C * m' % (gasDipoleTot, gasDipoleTotConv))
        print('Solvent dipole: %e D = %e C * m' % (solventDipoleTot, solventDipoleTotConv))
        print('Polarizability: %e C^2 * m^2 * J^-1' % polarizability)

    polarizationCost = 0.5 * transpose(solventDipoleConv-gasDipoleConv) * transpose(alphaPol ** -1) * \
    (solventDipoleConv-gasDipoleConv) * 6.02214129e23 / 1000

    polarizationCost = polarizationCost[0,0] + polarizationCost[1,1] + polarizationCost[2,2]

    print('The polarization cost (based on dipoles) is %f kJ/mol.' % polarizationCost)

    exit(0)