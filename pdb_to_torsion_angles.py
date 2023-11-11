import sys
import numpy as np
import pandas as pd
from Bio.PDB import *
from calculate_angles import *

def getPDBFile():
    isFileCorrect = True
    if "--pdb" in sys.argv and (sys.argv.index("--pdb")+1)<len(sys.argv):
        pdbFilepath = sys.argv[sys.argv.index("--pdb")+1]
        try:
            parser = PDBParser()
            PDBFile = parser.get_structure("Struct", pdbFilepath)
        except FileNotFoundError:
            print("Error! Wrong filepath.")
            isFileCorrect = False  
    else:
        print("Error! Please provide a file to run the algorithm.")
        isFileCorrect = False
        PDBFile = ""

    return isFileCorrect, PDBFile

def getTorsionAngleMatrix(PDBFile, includeSugarAngles):
    residues= [residue for residue in PDBFile.get_residues()]
    residuesWithoutHet = []
    for residue in residues:
        if(residue.get_resname() in ['A', 'C', 'U', 'G']):
            residuesWithoutHet.append(residue)

    angles = np.array([])
    anglesFromResidue = []
    for id, residue in enumerate(residuesWithoutHet):
        if id != 0:
            anglesFromResidue.append(alphaAngle(residuesWithoutHet[id-1],residue))
        else:
            anglesFromResidue.append(None)
        anglesFromResidue.append(betaAngle(residue))
        anglesFromResidue.append(gammaAngle(residue))
        anglesFromResidue.append(deltaAngle(residue))
        if id < (len(residuesWithoutHet) - 1):
            anglesFromResidue.append(epsilonAngle(residue, residuesWithoutHet[id+1]))
            anglesFromResidue.append(zetaAngle(residue, residuesWithoutHet[id+1]))
        else:
            anglesFromResidue.append(None)
            anglesFromResidue.append(None)
        anglesFromResidue.append(chiAngle(residue))
        if(includeSugarAngles):
            anglesFromResidue.append(v0Angle(residue))
            anglesFromResidue.append(v1Angle(residue))
            anglesFromResidue.append(v2Angle(residue))
            anglesFromResidue.append(v3Angle(residue))
            anglesFromResidue.append(v4Angle(residue))
        if id == 0:
            angles = np.array(anglesFromResidue)
        else:
            angles = np.vstack([angles,anglesFromResidue])
        anglesFromResidue=[]
    return angles

if __name__ == "__main__":
    isFileCorrect, PDBFile = getPDBFile()
    if isFileCorrect:
        angles = getTorsionAngleMatrix(PDBFile, True)
        DF = pd.DataFrame(angles)
        DF.to_csv("torsion_angles.csv", header=False, index=False)