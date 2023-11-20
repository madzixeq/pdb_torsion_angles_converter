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
    diff=15
    coords_A = np.array(A.get_coord())
    coords_B = np.array(B.get_coord())
    coords_C = np.array(C.get_coord())
    u1 = coords_B - coords_A
    u2 = coords_C - coords_B
    a = np.cross(u1, u2)
    t = vectorLength(a)*vectorLength(u2)*math.sin(np.radians(angle_BCD))*distance_CD*math.cos(np.radians(torsion_angle))
    b = a[1]*u2[2] - a[2]*u2[1]
    c = a[2]*u2[0] - a[0]*u2[2]
    d = a[0]*u2[1] - a[1]*u2[0]
    s = vectorLength(u2)*distance_CD*math.cos(np.radians(angle_BCD))
    e = c/b - u2[1]/u2[0]
    f = u2[2]/u2[0] - d/b 
    g = -s/u2[0] + t/b 
    h = -((u2[1]/u2[0])*(f/e) + u2[2]*u2[0])
    j = s/u2[0] - (u2[1]/u2[0])*(g/e)
    l = h**2 + f**2/e**2 + 1
    m = 2*h*j +(2*f*g)/(e**2)
    n = j**2 + g**2/e**2 - distance_CD**2
    # print(l,m,n)
    delta = (m**2 - 4*l*n)
    if delta >= 0:
        u3z = (-m-math.sqrt(m**2 - 4*l*n))/(2*l)
        u3 = [h*u3z+j, (f/e)*u3z + g/e, u3z]
        u3 = np.array(u3)
        torsion_angle_check = np.rad2deg(np.arccos(np.dot(np.cross(u1,u2),np.cross(u2,u3))/(vectorLength(np.cross(u1,u2))*vectorLength(np.cross(u2,u3)))))
        if np.absolute(torsion_angle_check-torsion_angle)<=diff:
            coords_D = coords_C + u3
        else:
            u3z = (-m+math.sqrt(m**2 - 4*l*n))/(2*l)
            u3 = [h*u3z+j, (f/e)*u3z + g/e, u3z]
            u3 = np.array(u3)
            coords_D = coords_C + u3
        print("torsion", torsion_angle)
        if torsion_angle > 0:
            print(np.rad2deg(np.arccos(np.dot(np.cross(u1,u2),np.cross(u2,u3))/(vectorLength(np.cross(u1,u2))*vectorLength(np.cross(u2,u3))))))
        else:
            print(-np.rad2deg(np.arccos(np.dot(np.cross(u1,u2),np.cross(u2,u3))/(vectorLength(np.cross(u1,u2))*vectorLength(np.cross(u2,u3))))))
        return coords_D
    # print('wrong')
    return coords_C + [1,1,1]


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
            
    

