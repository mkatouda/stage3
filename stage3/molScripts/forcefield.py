"""
    This file is part of STaGE - a wrapper for generating
    GROMACS topology files.

    Written by Magnus Lundborg
    Copyright (c) 2013-2014, The GROMACS development team.
    Check out http://www.gromacs.org for more information.

    This code is released under GNU LGPL V3.0 and may be released under
    a later version of the LGPL.
"""

import math
import re

class Forcefield:
    def __init__(self, name=None, fieldType=None):

        self.name = name
        if fieldType is None:
            self.fieldType = None
        else:
            self.fieldType = fieldType.upper()
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.aliases = []

    def readGromacsItpFile(self, filename):

        with open(filename) as f:
            lines = f.readlines()

        self.readGromacsItp(lines)

    def readGromacsItp(self, lines):

        mode=None
        for line in lines:
            line=line.strip()
            if len(line) < 1 or line.startswith(';') or line.startswith('#'):
                continue

            if line.find('atomtypes') != -1 or \
            line.find('bondtypes') != -1 or \
            line.find('constrainttypes') != -1 or \
            line.find('angletypes') != -1 or \
            line.find('dihedraltypes') != -1:
                mode = line.strip('[] ')
                continue

            a = Atom('X', 0)
            a.kind = 'X'
            a.charge = 0
            self.atoms.append(a)

            if mode == 'atomtypes':
                parts = line.split()
                if len(parts) < 8:
                    continue
                name = parts[0]
                kind = parts[1]
                elementNr = int(parts[2])
                mass = float(parts[3])
                charge = float(parts[4])
                ptype = parts[5]
                sigma = float(parts[6])
                epsilon = float(parts[7])
                if len(parts) > 8 and parts[8][0] == ';':
                    parts[8] = parts[8].strip(';')
                    comment = ' '.join(parts[8:])
                else:
                    comment = None

                a = Atom(name, mass)
                a.kind = kind
                a.elementNr = elementNr
                a.charge = charge
                a.ptype = ptype
                a.sigma = sigma
                a.epsilon = epsilon
                a.comment = comment
                self.atoms.append(a)

            elif mode == 'bondtypes':
                parts = line.split()
                if len(parts) < 5:
                    continue
                fromAtom = self.findAtom(parts[0])
                toAtom = self.findAtom(parts[1])
                if not fromAtom or not toAtom:
                    print('Bond Atoms not found: %s %s' % (parts[0], parts[1]))
                    continue

                func = int(parts[2])
                if func in [1, 2, 6, 7]:
                    dist = float(parts[3])
                    force = float(parts[4])
                    if len(parts) > 5 and parts[5][0] == ';':
                        parts[5] = parts[5].strip(';')
                        comment = ' '.join(parts[5:])
                    else:
                        comment = None
                    b = Bond(fromAtom, toAtom)
                    b.func = func
                    b.dist = dist
                    b.force = force
                    b.comment = comment

                else:
                    print('Bond func not implemented')
                    continue

                self.bonds.append(b)

            elif mode == 'angletypes':
                parts = line.split()
                if len(parts) < 6:
                    continue
                atom1 = self.findAtom(parts[0])
                atom2 = self.findAtom(parts[1])
                atom3 = self.findAtom(parts[2])
                if not atom1 or not atom2 or not atom3:
                    print('Angle atoms not found: %s %s %s' % (parts[0], parts[1], parts[2]))
                    continue

                func = int(parts[3])
                if func in [1, 2, 6]:
                    angle = float(parts[4])
                    force = float(parts[5])
                    if len(parts) > 6 and parts[6][0] == ';':
                        parts[6] = parts[6].strip(';')
                        comment = ' '.join(parts[6:])
                    else:
                        comment = None
                    a = Angle(atom1, atom2, atom3)
                    a.func = func
                    a.angle = angle
                    a.force = force
                    a.comment = comment
                elif func == 5:
                    angle = float(parts[4])
                    force = float(parts[5])
                    r13 = float(parts[6])
                    kub = float(parts[7])
                    if len(parts) > 8 and parts[8][0] == ';':
                        parts[8] = parts[8].strip(';')
                        comment = ' '.join(parts[8:])
                    else:
                        comment = None
                    a = Angle(atom1, atom2, atom3)
                    a.func = func
                    a.angle = angle
                    a.force = force
                    a.r13 = r13
                    a.kub = kub
                    a.comment = comment

                else:
                    print('Angle func not implemented')
                    continue

                self.angles.append(a)

            elif mode == 'dihedraltypes':
                parts = line.split()
                if len(parts) < 7:
                    continue
                atom1 = self.findAtom(parts[0])
                atom2 = self.findAtom(parts[1])
                atom3 = self.findAtom(parts[2])
                atom4 = self.findAtom(parts[3])
                if not atom1 or not atom2 or not atom3 or not atom4:
                    if parts[0] != 'X' and parts[1] != 'X' and parts[2] != 'X' and parts[3] != 'X':
                        print('Dihedral atoms not found: %s %s %s %s' % (parts[0], parts[1], parts[2], parts[3]))
                        continue

                func = int(parts[4])
                if func in [1, 4, 9]:
                    phase = float(parts[5])
                    force = float(parts[6])
                    nMinima = int(parts[7])
                    if len(parts) > 8 and parts[8][0] == ';':
                        parts[8] = parts[8].strip(';')
                        comment = ' '.join(parts[8:])
                    else:
                        comment = None
                    d = Dihedral(atom1, atom2, atom3, atom4)
                    d.func = func
                    d.phase = phase
                    d.force = force
                    d.nMinima = nMinima
                    d.comment = comment

                # Improper dihedral
                elif func == 2:
                    phase = float(parts[5])
                    force = float(parts[6])
                    if len(parts) > 7 and parts[7][0] == ';':
                        parts[7] = parts[7].strip(';')
                        comment = ' '.join(parts[7:])
                    else:
                        comment = None
                    d = Dihedral(atom1, atom2, atom3, atom4)
                    d.func = func
                    d.phase = phase
                    d.force = force
                    d.comment = comment

                elif func == 3:
                    cs = []
                    for i in range(6):
                        cs.append(float(parts[5+i]))
                    d = Dihedral(atom1, atom2, atom3, atom4)
                    d.func = func
                    d.cs = cs
                    d.comment = comment


                else:
                    print('Dihedral func not implemented')
                    continue

                self.dihedrals.append(d)


    def findAtom(self, name):

        for atom in self.atoms:
            if atom.name == name:
                return atom

        for atom in self.atoms:
            if atom.kind == name:
                return atom

        return None

    def findAtomsOfKind(self, kind):

        atoms = []
        for atom in self.atoms:
            if atom.kind == kind:
                atoms.append(atom)

        return atoms

    def findBond(self, fromAtom, toAtom):

        for bond in self.bonds:
            if bond.fromAtom == fromAtom and bond.toAtom == toAtom:
                return bond
            # Allow reverse definition as well.
            if bond.fromAtom == toAtom and bond.toAtom == fromAtom:
                return bond

        # Check if the bond is using the same atom type/atom kind as the
        # requested atoms.
        if fromAtom.kind and toAtom.kind:
            for bond in self.bonds:
                if bond.fromAtom.kind == fromAtom.kind and bond.toAtom.kind == toAtom.kind:
                    return bond
                # Allow reverse definition as well.
                if bond.fromAtom.kind == toAtom.kind and bond.toAtom.kind == fromAtom.kind:
                    return bond

        return None

    def findAngle(self, a1, a2, a3):

        for angle in self.angles:
            if angle.atom1 == a1 and angle.atom2 == a2 and angle.atom3 == a3:
                return angle
            # Allow reverse definition as well.
            if angle.atom1 == a3 and angle.atom2 == a2 and angle.atom3 == a1:
                return angle

        # Check if the angle is using the same atom type/atom kind as the
        # requested atoms.
        if a1.kind and a2.kind and a3.kind:
            for angle in self.angles:
                if a1.kind == angle.atom1.kind and \
                a2.kind == angle.atom2.kind and \
                a3.kind == angle.atom3.kind:
                    return angle
                # Allow reverse definition as well.
                if a3.kind == angle.atom1.kind and \
                a2.kind == angle.atom2.kind and \
                a1.kind == angle.atom3.kind:
                    return angle

        return None

    def findDihedral(self, a1, a2, a3, a4, func = None):

        for dihedral in self.dihedrals:
            if dihedral.atom1 == a1 and dihedral.atom2 == a2 and dihedral.atom3 == a3 and dihedral.atom4 == a4:
                if func == None or dihedral.func == func:
                    return dihedral

            # Allow reverse definition as well.
            if dihedral.atom1 == a4 and dihedral.atom2 == a3 and dihedral.atom3 == a2 and dihedral.atom4 == a1:
                if func == None or dihedral.func == func:
                    return dihedral

        if a1.kind and a2.kind and a3.kind and a4.kind:
            for dihedral in self.dihedrals:
                if a1.kind == dihedral.atom1.kind and \
                a2.kind == dihedral.atom2.kind and \
                a3.kind == dihedral.atom3.kind and \
                a4.kind == dihedral.atom4.kind:
                    if func == None or dihedral.func == func:
                        return dihedral
                    # Allow reverse definition as well.
                if a4.kind == dihedral.atom1.kind and \
                a3.kind == dihedral.atom2.kind and \
                a2.kind == dihedral.atom3.kind and \
                a1.kind == dihedral.atom4.kind:
                    if func == None or dihedral.func == func:
                        return dihedral


        # Try with wildcards
        for dihedral in self.dihedrals:
            if (dihedral.atom1.name == 'X' or dihedral.atom1 == a1) and \
               (dihedral.atom2.name == 'X' or dihedral.atom2 == a2) and \
               (dihedral.atom3.name == 'X' or dihedral.atom3 == a3) and \
               (dihedral.atom4.name == 'X' or dihedral.atom4 == a4):
                if func == None or dihedral.func == func:
                    return dihedral
            if (dihedral.atom1.name == 'X' or dihedral.atom1 == a4) and \
            (dihedral.atom2.name == 'X' or dihedral.atom2 == a3) and \
            (dihedral.atom3.name == 'X' or dihedral.atom3 == a2) and \
            (dihedral.atom4.name == 'X' or dihedral.atom4 == a1):
                if func == None or dihedral.func == func:
                    return dihedral


        if a1.kind and a2.kind and a3.kind and a4.kind:
            for dihedral in self.dihedrals:
                if (dihedral.atom1.name == 'X' or dihedral.atom1.kind == a1.kind) and \
                (dihedral.atom2.name == 'X' or dihedral.atom2.kind == a2.kind) and \
                (dihedral.atom3.name == 'X' or dihedral.atom3.kind == a3.kind) and \
                (dihedral.atom4.name == 'X' or dihedral.atom4.kind == a4.kind):
                    if func == None or dihedral.func == func:
                        return dihedral
                if (dihedral.atom1.name == 'X' or dihedral.atom1.kind == a4.kind) and \
                (dihedral.atom2.name == 'X' or dihedral.atom2.kind == a3.kind) and \
                (dihedral.atom3.name == 'X' or dihedral.atom3.kind == a2.kind) and \
                (dihedral.atom4.name == 'X' or dihedral.atom4.kind == a1.kind):
                    if func == None or dihedral.func == func:
                        return dihedral

        return None

    def findDihedrals(self, a1, a2, a3, a4, func = None):

        results = []

        for dihedral in self.dihedrals:
            if dihedral.atom1 == a1 and dihedral.atom2 == a2 and dihedral.atom3 == a3 and dihedral.atom4 == a4:
                if func == None or dihedral.func == func:
                    results.append(dihedral)

            # Allow reverse definition as well.
            if dihedral.atom1 == a4 and dihedral.atom2 == a3 and dihedral.atom3 == a2 and dihedral.atom4 == a1:
                if func == None or dihedral.func == func:
                    results.append(dihedral)

        if a1.kind and a2.kind and a3.kind and a4.kind:
            for dihedral in self.dihedrals:
                if a1.kind == dihedral.atom1.kind and \
                a2.kind == dihedral.atom2.kind and \
                a3.kind == dihedral.atom3.kind and \
                a4.kind == dihedral.atom4.kind:
                    results.append(dihedral)

        # Try with wildcards
        for dihedral in self.dihedrals:
            if (dihedral.atom1.name == 'X' or dihedral.atom1 == a1) and \
            (dihedral.atom2.name == 'X' or dihedral.atom2 == a2) and \
            (dihedral.atom3.name == 'X' or dihedral.atom3 == a3) and \
            (dihedral.atom4.name == 'X' or dihedral.atom4 == a4):
                if func == None or dihedral.func == func:
                    results.append(dihedral)

        if a1.kind and a2.kind and a3.kind and a4.kind:
            for dihedral in self.dihedrals:
                if (dihedral.atom1.name == 'X' or dihedral.atom1.kind == a1.kind) and \
                (dihedral.atom2.name == 'X' or dihedral.atom2.kind == a2.kind) and \
                (dihedral.atom3.name == 'X' or dihedral.atom3.kind == a3.kind) and \
                (dihedral.atom4.name == 'X' or dihedral.atom4.kind == a4.kind):
                    if func == None or dihedral.func == func:
                        results.append(dihedral)

        return results

# Class for atoms, currently only containing data required for Q force field
class Atom:
    def __init__(self, name, mass, comment=None):
        self.name = name
        self.mass = mass
        self.kind = None
        self.elementNr = None
        self.charge = None
        self.ptype = None
        self.sigma = None
        self.epsilon = None
        self.r1 = None
        self.r2 = None
        self.r3 = None
        self.e1 = None
        self.e2 = None
        self.comment = comment

    # Not finished
    def getGromacsOutput(self, includeComment = True):
        if self.name == None or self.mass == None:
            return None
        str = '%-8s%11s' % (self.name, self.mass)
        if includeComment and self.comment != None:
            str = str + ' ; %s' % self.comment

    def setR1(self, r1):
        self.r1 = r1
    def setR2(self, r2):
        self.r2 = r2
    def setR3(self, r3):
        self.r3 = r3
    def setE1(self, e1):
        self.e1 = e1
    def setE2(self, e2):
        self.e2 = e2

    def setAllR(self, r1, r2, r3):
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
    def setAllE(self, e1, e2):
        self.e1 = e1
        self.e2 = e2
    def appendComment(self, comment):
        if self.comment != None:
            self.comment = self.comment+', '+comment
        else:
            self.comment = comment


class Bond:
    def __init__(self, fromAtom, toAtom, comment=None):

        self.fromAtom = fromAtom
        self.toAtom = toAtom
        self.force = None
        self.func = 1
        self.dist = None
        self.comment = comment

    # Not finished
    def getGromacsOutput(self, includeComment = True):

        if self.func in [1, 2, 6, 7]:
            str='%7s %7s %5d %10.4e %10.4e' % (self.fromAtom.name, self.toAtom.name, self.func, self.dist, self.force)
        else:
            print('Bond func not implemented')
            return ''

        if includeComment and self.comment != None:
            str=str + ' ; %s' % self.comment
        return str

class Angle:
    def __init__(self, atom1, atom2, atom3, comment = None):

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.func = 1
        self.force = None
        self.angle = None
        self.kub = None
        self.r13 = None
        self.s0 = None
        self.comment = comment

    # Not finished
    def getGromacsOutput(self, includeComment = True):

        if self.func in [1, 2, 6]:
            str = '%6s %6s %6s %5d %10.4e %10.4e' % (self.atom1.name, self.atom2.name,
            self.atom3.name, self.func, self.angle, self.force)
        elif self.func == 5:
            str = '%6s %6s %6s %7d %10f %10f %10f %10f' % (self.atom1.name, self.atom2.name,
            self.atom3.name, self.func, self.angle, self.force, self.r13, self.kub)
        else:
            print('Angle func not implemented')
            return ''

        if includeComment and self.comment != None:
            str=str + ' ; %s' % self.comment
        return str

class Dihedral:
    def __init__(self, atom1, atom2, atom3, atom4, comment = None):

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.func = 1
        self.force = None
        self.nMinima = None
        self.phase = None
        self.paths = None
        self.cs = []
        self.comment = comment

    def getGromacsOutput(self, includeComment = True):

        if self.func in [1, 4, 9]:
            str = '%6s %6s %6s %6s %5d %9.2f %9f %6d' % (self.atom1.name, self.atom2.name,
            self.atom3.name, self.atom4.name, self.func, self.phase, self.force, self.nMinima)

        # Improper dihedral
        elif self.func == 2:
            str = '%6s %6s %6s %6s %5d %9.2f %9f' % (self.atom1.name, self.atom2.name,
            self.atom3.name, self.atom4.name, self.func, self.phase, self.force)

        elif self.func == 3:
            str = '%6s %6s %6s %6s %5d ' % (self.atom1.name, self.atom2.name,
            self.atom3.name, self.atom4.name, self.func)
            for i in range(6):
                str += '%10.5f ' % self.cs[i]
        else:
            print('Dihedral func not implemented')
            return ''

        if includeComment and self.comment != None:
            str = str + ' ; %s' % self.comment
        return str

class Alias:
    def __init__(self, oldName, newName):
        self.oldName=oldName
        self.newName=newName
