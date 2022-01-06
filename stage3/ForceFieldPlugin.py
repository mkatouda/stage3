"""
    This file is part of STaGE - a wrapper for generating
    GROMACS topology files.

    Written by Magnus Lundborg
    Copyright (c) 2013-2015, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

class ForceFieldPlugin():
    """ A plugin for generating a GROMACS topology for a force field. To generate a new
    plugin make a new file and class in the plugins directory inheriting this class. """

    def __init__(self):

        self.forceFieldName = None
        """ The name of the force field. """
        self.order = 9
        """ Used for listing and executing the plugins in a certain order.
        Plugins with low order are listed and executed before plugins with high order.
        If plugin A requires results from plugin B plugin A needs a higher order number than
        plugin B. """
        self.condition = None

    def generate(self, inputFile, output, keepMol2Charges = False, netCharge = None,
                  verbose = False):
        """ This method generates the general force field parameters. inputFile is a mol2 file
        describing the molecule. output is the basename of the output that will be generated, the
        force field parameters will be generated in a directory named output+'_'+self.forceFieldName.
        If keepMol2Charges is True the charges in the input file will be kept instead of generating
        new charges. netCharge is used to specify what the net charge of the molecule should be -
        this is only used when keepMol2Charges == False. If verbose == True more verbose output
        of the topology generation process will be given."""

        return None

    def fixAssignment(self, output, retainCharges = False, netCharge = 0, spreadCharges = False,
                      verbose = False):
        """ If any modifications need to be done to the generated topology it can be done in this
        method. output is the full path to the basename of the output. If retainCharges == True the
        charges of the atoms will not be modified. netCharge is used to specify what the net charge
        of the molecule should be - this is only used if retainCharges == False.
        If spreadCharges == True and the total charge of the atoms in the molecule does not match
        netCharge the difference will be compensated by spreading a charge over all atoms.
        If verbose == True more verbose output of the process will be given. The method does not
        return anything."""

        return None

    def convert2Gromacs(self, output, verbose = False):
        """ If the generated topology, so far, is not in GROMACS format this method should convert
        it to GROMACS format. output is the base name of the outputs (also used for finding
        the inputs), not including extensions
        and, e.g., force field directory names. If verbose == True more verbose output of the
        process will be given."""

        return None

    def genTop(self, output, solvent = None, verbose = False, forcefieldDir = None):
        """ Generates the actual topology (.top) file (which in turns includes the general force field
        parameters and the .itp files generated for the molecule and solvents etc). output is the base name
        of the output from STaGE. If a solvent is specified it will be added to the topology file. If
        verbose == True more verbose output of the process will be given. The file name of the
        topology is returned."""

        return None

    def coordsToTopology(self, output, coordsFile, verbose = False):
        """ Generates a topology using GROMACS tools (pdb2gmx) from a coordinate file (.pdb).
        This is used for e.g. generating a topology for a protein in a protein-ligand system.
        output is the base name of the output from STaGE. coordsFile is the pdb file to use
        as input. If verbose == Ture more verbose output of the process will be given.
        A tuple consisting of the name of the topology file and the generate .gro coordinate
        file is returned."""
        return None

    def finalClean(self, output):
        """ If there are files or directories that should be removed after generating the force field
        that should be done in this method. output is the base name of the files that should be removed."""

        return None
