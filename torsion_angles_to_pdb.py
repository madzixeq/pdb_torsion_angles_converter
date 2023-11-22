import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
import math
import json
import csv
from random import randint
from Bio.PDB import *
from pdb_to_torsion_angles import getPDBFile

def vectorLength(a):
    return math.sqrt(a[0]**2+a[1]**2+a[2]**2)

def getDataFromConfig():
    with open("config.json", "r") as read_file:
        config = json.load(read_file)
    return config

def getTorsionAngles():
    torsion_angles = list(csv.reader(open("torsion_angles.csv", "r"), delimiter=","))
    for i in range(len(torsion_angles)):
        for j in range(len(torsion_angles[0])):
            if torsion_angles[i][j]=='':
                torsion_angles[i][j]=None
            else:
                torsion_angles[i][j]=float(torsion_angles[i][j])
    return torsion_angles


def getForthPoint(distance_CD, angle_BCD, torsion_angle, A, B, C):
    #Równanie zostało wyznaczone na kartce i przepisane do funkcji
    #stąd specyficzne oznaczenia
    coords_A = np.array(A.get_coord())
    coords_B = np.array(B.get_coord())
    coords_C = np.array(C.get_coord())
    u1 = coords_B - coords_A
    u2 = coords_C - coords_B
    a = np.cross(u1, u2)
    b = vectorLength(a)*vectorLength(u2)*distance_CD*math.sin(np.deg2rad(angle_BCD))
    c1 = b*math.cos(np.deg2rad(torsion_angle))
    d1 = a[1]*u2[2]-a[2]*u2[1]
    e1 = a[2]*u2[0]-a[0]*u2[2]
    f1 = a[0]*u2[1] - a[1]*u2[0]
    c2 = b*vectorLength(u2)*math.sin(np.deg2rad(torsion_angle))
    d2 = u2[2]**2*a[0] + u2[1]**2*a[0]-u2[0]*a[1]*u2[1]-u2[0]*a[2]*u2[2]
    e2 = u2[0]**2*a[1] + u2[2]**2*a[1]-u2[1]*a[2]*u2[2]-u2[1]*a[0]*u2[0]
    f2 = u2[0]**2*a[2] + u2[1]**2*a[2]-u2[2]*a[0]*u2[0]-u2[2]*a[1]*u2[1]
    c3 = vectorLength(u2)*distance_CD*math.cos(np.deg2rad(angle_BCD))
    d3 = -u2[0]
    e3 = -u2[1]
    f3 = -u2[2]
    w = d1*e2*f3+e1*f2*d3+f1*d2*e3-f1*e2*d3-d1*f2*e3-e1*d2*f3
    wx = c1*e2*f3+e1*f2*c3+f1*c2*e3-f1*e2*c3-c1*f2*e3-e1*c2*f3
    wy = d1*c2*f3+c1*f2*d3+f1*d2*c3-f1*c2*d3-d1*f2*c3-c1*d2*f3
    wz = d1*e2*c3+e1*c2*d3+c1*d2*e3-c1*e2*d3-d1*c2*e3-e1*d2*c3
    if w!=0:
        u3 = np.array([wx/w,wy/w,wz/w])
        print(angle_BCD)
        print(np.rad2deg(np.arccos(np.dot(u2,u3)/(vectorLength(u2)*vectorLength(u3)))))
        return coords_C + u3
    else:
        print('wrong')
        return coords_C +[1,1,1]
    
#def getPointFromDistAnd2Angles()
#def getPointFromAngleDistAndPlane()
#def planeCoords()


if __name__ == "__main__":
    torsionAngles = getTorsionAngles()
    isFileCorrect, PDBFile = getPDBFile()
    residues = list()
    originalIndexes = list()
    if isFileCorrect:
        newPDBFile = PDBFile.copy()
        for id, residue in enumerate(PDBFile.get_residues()):
            if(residue.get_resname() in ['A', 'C', 'U', 'G']):
                residues.append(residue)
                originalIndexes.append(id)
    json = getDataFromConfig()
    for id, residue in enumerate(residues):
        if id!=0:
            residues[id]["C5'"].set_coord(getForthPoint(json["rybosine"]["distance"]["O5'C5'"], json["phosfor"]["angles"]["PO5'C5'"],torsionAngles[id][0],residues[id-1]["O3'"], residues[id]["P"], residues[id]["O5'"]))
        residues[id]["C4'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C5'C4'"], json["rybosine"]["angles"]["O5'C5'C4'"],torsionAngles[id][1],residues[id]["P"], residues[id]["O5'"], residues[id]["C5'"]))
        residues[id]["C3'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C4'C3'"], json["rybosine"]["angles"]["C5'C4'C3'"],torsionAngles[id][2],residues[id]["O5'"], residues[id]["C5'"], residues[id]["C4'"]))
        residues[id]["O3'"].set_coord(getForthPoint(json["rybosine"]["distance"]["O3'C3'"], json["rybosine"]["angles"]["C4'C3'O3'"],torsionAngles[id][3],residues[id]["C5'"], residues[id]["C4'"], residues[id]["C3'"]))
        if id < (len(residues) - 1):
            residues[id+1]["P"].set_coord(getForthPoint(json["phosfor"]["distance"]["O3'P"], json["phosfor"]["angles"]["C3'O3'P"],torsionAngles[id][4],residues[id]["C4'"], residues[id]["C3'"],residues[id]["O3'"]))
            residues[id+1]["O5'"].set_coord(getForthPoint(json["phosfor"]["distance"]["PO5'"], json["phosfor"]["angles"]["O3'PO5'"],torsionAngles[id][5],residues[id]["C3'"],residues[id]["O3'"], residues[id+1]["P"]))
        residues[id]["O4'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C4'O4'"], json["rybosine"]["angles"]["C3'C4'O4'"],torsionAngles[id][10],residue["C2'"],residue["C3'"],residue["C4'"]))
        residues[id]["C1'"].set_coord(getForthPoint(json["rybosine"]["distance"]["O4'C1'"], json["rybosine"]["angles"]["C4'O4'C1'"],torsionAngles[id][11],residue["C3'"],residue["C4'"],residue["O4'"]))
        residues[id]["C2'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C2'C1'"], json["rybosine"]["angles"]["C2'C1'O4'"],torsionAngles[id][1],residue["C4'"],residue["O4'"],residue["C1'"]))

        #Potem C4/C2 już normalnie
        #Zasada - wyznaczyć płaszczyznę i zrobić na niej punkty    

    newStructure = PDBFile.copy()
    tmp = 0
    for id, residue in enumerate(newStructure.get_residues()):
        if id in originalIndexes:
            for atom in residue.get_atoms():
                atom.set_coord(residues[tmp][atom.get_name()].get_coord())
            tmp+=1


    io=PDBIO()
    io.set_structure(newStructure)
    io.save("1ehz_from_torsion.pdb")  
            
    

