from Bio.PDB import *
import numpy as np
import math

def alphaAngle(prevResidue, residue):
    vector1 = prevResidue['O3\''].get_vector()
    vector2 = residue['P'].get_vector()
    vector3 = residue['O5\''].get_vector()
    vector4 = residue['C5\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def betaAngle(residue):
    vector1 = residue['P'].get_vector()
    vector2 = residue['O5\''].get_vector()
    vector3 = residue['C5\''].get_vector()
    vector4 = residue['C4\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def gammaAngle(residue):
    vector1 = residue['O5\''].get_vector()
    vector2 = residue['C5\''].get_vector()
    vector3 = residue['C4\''].get_vector()
    vector4 = residue['C3\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def deltaAngle(residue):
    vector1 = residue['C5\''].get_vector()
    vector2 = residue['C4\''].get_vector()
    vector3 = residue['C3\''].get_vector()
    vector4 = residue['O3\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def epsilonAngle(residue, nextResidue):
    vector1 = residue['C4\''].get_vector()
    vector2 = residue['C3\''].get_vector()
    vector3 = residue['O3\''].get_vector()
    vector4 = nextResidue['P'].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def zetaAngle(residue, nextResidue):
    vector1 = residue['C3\''].get_vector()
    vector2 = residue['O3\''].get_vector()
    vector3 = nextResidue['P'].get_vector()
    vector4 = nextResidue['O5\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def chiAngle(residue):
    if(residue.get_resname() in ['A', 'G']):
        vector3 = residue['N9'].get_vector()
        vector4 = residue['C4'].get_vector()
    else:
        vector3 = residue['N1'].get_vector()
        vector4 = residue['C2'].get_vector()
    vector1 = residue['O4\''].get_vector()
    vector2 = residue['C1\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def v0Angle(residue):
    vector1 = residue['C4\''].get_vector()
    vector2 = residue['O4\''].get_vector()
    vector3 = residue['C1\''].get_vector()
    vector4 = residue['C2\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def v1Angle(residue):
    vector1 = residue['O4\''].get_vector()
    vector2 = residue['C1\''].get_vector()
    vector3 = residue['C2\''].get_vector()
    vector4 = residue['C3\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def v2Angle(residue):
    vector1 = residue['C1\''].get_vector()
    vector2 = residue['C2\''].get_vector()
    vector3 = residue['C3\''].get_vector()
    vector4 = residue['C4\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def v3Angle(residue):
    vector1 = residue['C2\''].get_vector()
    vector2 = residue['C3\''].get_vector()
    vector3 = residue['C4\''].get_vector()
    vector4 = residue['O4\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))

def v4Angle(residue):
    vector1 = residue['C3\''].get_vector()
    vector2 = residue['C4\''].get_vector()
    vector3 = residue['O4\''].get_vector()
    vector4 = residue['C1\''].get_vector()
    return math.degrees(calc_dihedral(vector1, vector2, vector3, vector4))