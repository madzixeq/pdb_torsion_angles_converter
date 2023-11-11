import numpy as np
import math
import json
from pdb_to_torsion_angles import getPDBFile

def getDataFromConfig():
    with open("config.json", "r") as read_file:
        config = json.load(read_file)
    return config

def getForthPoint(bondLength, bondAngle, torsionAngle, prevAtom):
    coords = list()
    coords.append(prevAtom.get_coord()[0] + bondLength * np.sin(np.radians(bondAngle)))
    coords.append(prevAtom.get_coord()[1] + bondLength * np.cos(np.radians(bondAngle)) * np.sin(np.radians(torsionAngle)))
    coords.append(prevAtom.get_coord()[2] + bondLength * np.cos(np.radians(bondAngle)) * np.cos(np.radians(torsionAngle)))
    return coords


if __name__ == "__main__":
    isFileCorrect, PDBFile = getPDBFile()
    residues = list()
    if isFileCorrect:
        for residue in PDBFile.get_residues():
            if(residue.get_resname() in ['A', 'C', 'U', 'G']):
                residues.append(residue)
    

