"""
    This file is part of STaGE - a wrapper for generating
    GROMACS topology files.

    Written by Magnus Lundborg
    Copyright (c) 2013-2015, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

class ChargePlugin():
    """ A plugin for partial charge calculations. To generate a new plugin
    make a new file and class in the plugins directory inheriting this class. """

    def __init__(self):

        self.programName = None
        """ The name of the program used to generate the charges. """
        self.alternativeChargeMethods = {}
        """ A dictionary of charge methods implemented by this plugin, using the same infrastructure and
        a brief description of the charge method. Each alternative charge method should be handled in
        the generateCharges method """
        self.order = 9
        """ Used for listing the plugins in a certain order. Plugins with low order are listed
        before plugins with high order. """

    def generateCharges(self, inFile, method, netCharge, multiplier = 1.0, verbose = False):
        """ A method for generating partial charges. inFile is a mol2 file used as input and will
        also be populated with the calculated charges from this method. method
        is a string corresponding to one of the keys in self.alternativeChargeMethods.
        netCharge is the net charge of the molecule. multiplier is an optional factor to apply to
        the calculated charges. verbose is an option whether to generate verbose output or not.
        The method is not expected to return anything. """

        return None
