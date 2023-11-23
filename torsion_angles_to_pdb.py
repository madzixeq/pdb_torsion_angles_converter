import numpy as np
import pandas as pd
import math
import json
import csv
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

def planeCoords(a, b, c):
    u1 = b - a
    u2 = c - a
    unorm = np.cross(u1,u2)
    k = -(unorm[0]*a[0]+unorm[1]*a[1]+unorm[2]*a[2])
    return (unorm[0], unorm[1], unorm[2], k)


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
        return coords_C + u3
    else:
        print('wrong')
        return coords_C +[1,1,1]
    
def getPointFromDistAnd2Angles(distance_CD, angle_ACD, angle_BCD, A, B, C, plus):
    coords_A = np.array(A.get_coord())
    coords_B = np.array(B.get_coord())
    coords_C = np.array(C.get_coord())
    u1 = coords_C - coords_A
    u2 = coords_C - coords_B
    a = -vectorLength(u1)*distance_CD*math.cos(np.deg2rad(angle_ACD))
    b = -vectorLength(u2)*distance_CD*math.cos(np.deg2rad(angle_BCD))
    c = a/u1[0] - b/u2[0]
    d = u2[2]/u2[0] - u1[2]/u1[0]
    e = u1[1]/u1[0] - u2[1]/u2[0]
    f = -u1[2]/u1[0] - (u1[1]*d)/(u1[0]*e)
    g = a/u1[0] - (u1[1]*c)/(u1[0]*e)
    l = f**2 + d**2/e**2 + 1
    m = 2*f*g + 2*d*c/e**2
    n = g**2 + c**2/e**2 - distance_CD**2
    delta = m**2 - 4*l*n
    if delta>=0:
        if plus==True:
            u3z = (-m+math.sqrt(delta))/(2*l)
        else:
            u3z = (-m-math.sqrt(delta))/(2*l)
        u3x = f*u3z+g
        u3y = (d/e)*u3z + c/e
        u3 = np.array([u3x,u3y,u3z])
        return coords_C + u3
    else:
        print('wrong')
        return coords_C + [1,1,1]

def getPointFromAngleDistAndPlane(distance_BC,angle_ABC,A,B,p,plus):
    coords_A = np.array(A.get_coord())
    coords_B = np.array(B.get_coord())
    u1 = coords_B-coords_A
    a = -p[3]-p[0]*coords_B[0]-p[1]*coords_B[1]-p[2]*coords_B[2]
    b = -vectorLength(u1)*distance_BC*math.cos(np.deg2rad(angle_ABC))
    c = u1[1]/u1[0] - p[1]/p[0]
    d = b/u1[0] - a/p[0]
    e = p[2]/p[0] - u1[2]/u1[0]
    f = a/p[0] - (p[1]*d)/(p[0]*c)
    g = -((p[1]*e)/(p[0]*c)) - p[2]/p[0]
    l = g**2 + e**2/c**2 + 1
    m = 2*g*f + 2*e*d/c**2
    n = f**2 + d**2/c**2 - distance_BC**2
    delta = m**2 - 4*l*n
    if delta>=0:
        if plus==True:
            u2z = (-m+math.sqrt(delta))/(2*l)
        else:
            u2z = (-m-math.sqrt(delta))/(2*l)
        u2x = g*u2z+f
        u2y = (e/c)*u2z + d/c
        u2 = np.array([u2x,u2y,u2z])
        return coords_B + u2
    else:
        print('wrong')
        return coords_B + [1,1,1]

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
        residues[id]["C2'"].set_coord(getPointFromDistAnd2Angles(json["rybosine"]["distance"]["C3'C2'"],json["rybosine"]["angles"]["C2'C3'C4'"],json["rybosine"]["angles"]["O3'C3'C2'"],residues[id]["C4'"], residues[id]["O3'"], residues[id]["C3'"],True))
        residues[id]["O4'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C4'O4'"], json["rybosine"]["angles"]["C3'C4'O4'"],torsionAngles[id][10],residue["C2'"],residue["C3'"],residue["C4'"]))
        residues[id]["C1'"].set_coord(getForthPoint(json["rybosine"]["distance"]["O4'C1'"], json["rybosine"]["angles"]["C4'O4'C1'"],torsionAngles[id][11],residue["C3'"],residue["C4'"],residue["O4'"]))
        residues[id]["C2'"].set_coord(getForthPoint(json["rybosine"]["distance"]["C2'C1'"], json["rybosine"]["angles"]["C2'C1'O4'"],torsionAngles[id][7],residue["C4'"],residue["O4'"],residue["C1'"]))
        if residue.get_resname() in ['A','G']:
            residues[id]["N9"].set_coord(getPointFromDistAnd2Angles(json["rybosine"]["distance"]["C1'N"],json["rybosine"]["angles"]["C2'C1'N9"],json["rybosine"]["angles"]["O4'C1'N9"],residues[id]["C2'"], residues[id]["O4'"], residues[id]["C1'"],False))
            if residue.get_resname() == 'A':
                residues[id]["C4"].set_coord(getForthPoint(json["A"]["distance"]["C4N9"], json["A"]["angles"]["C1'N9C4"],torsionAngles[id][6],residue["O4'"],residue["C1'"],residue["N9"]))
                plane = planeCoords(residues[id]["C1'"].get_coord(),residues[id]["N9"].get_coord(),residues[id]["C4"].get_coord())
                residues[id]["N3"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["N3C4"],json["A"]["angles"]["N9C4N3"],residues[id]["N9"],residues[id]["C4"],plane,True))
                residues[id]["C2"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["C2N3"],json["A"]["angles"]["C4N3C2"],residues[id]["C4"],residues[id]["N3"],plane,False))
                residues[id]["N1"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["N1C2"],json["A"]["angles"]["N3C2N1"],residues[id]["N3"],residues[id]["C2"],plane,False))
                residues[id]["C6"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["C6N1"],json["A"]["angles"]["C2N1C6"],residues[id]["C2"],residues[id]["N1"],plane,False))
                residues[id]["N6"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["C6N6"],json["A"]["angles"]["N1C6N6"],residues[id]["N1"],residues[id]["C6"],plane,False))
                residues[id]["C5"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["C5C6"],json["A"]["angles"]["C5C6N1"],residues[id]["N1"],residues[id]["C6"],plane,True))
                residues[id]["N7"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["N7C5"],json["A"]["angles"]["N7C5C6"],residues[id]["C6"],residues[id]["C5"],plane,False))
                residues[id]["C8"].set_coord(getPointFromAngleDistAndPlane(json["A"]["distance"]["C8N7"],json["A"]["angles"]["C8N7C5"],residues[id]["C5"],residues[id]["N7"],plane,True))
            if residue.get_resname() == 'G':
                residues[id]["C4"].set_coord(getForthPoint(json["G"]["distance"]["C4N9"], json["G"]["angles"]["C1'N9C4"],torsionAngles[id][6],residue["O4'"],residue["C1'"],residue["N9"]))
                plane = planeCoords(residues[id]["C1'"].get_coord(),residues[id]["N9"].get_coord(),residues[id]["C4"].get_coord())
                residues[id]["N3"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["N3C4"],json["G"]["angles"]["N9C4N3"],residues[id]["N9"],residues[id]["C4"],plane,True))
                residues[id]["C2"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C2N3"],json["G"]["angles"]["C4N3C2"],residues[id]["C4"],residues[id]["N3"],plane,False))
                residues[id]["N2"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C2N2"],json["G"]["angles"]["N3C2N2"],residues[id]["N3"],residues[id]["C2"],plane,True))
                residues[id]["N1"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["N1C2"],json["G"]["angles"]["N3C2N1"],residues[id]["N3"],residues[id]["C2"],plane,False))
                residues[id]["C6"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C6N1"],json["G"]["angles"]["C2N1C6"],residues[id]["C2"],residues[id]["N1"],plane,False))
                residues[id]["O6"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C6O6"],json["G"]["angles"]["N1C6O6"],residues[id]["N1"],residues[id]["C6"],plane,False))
                residues[id]["C5"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C5C6"],json["G"]["angles"]["N1C6C5"],residues[id]["N1"],residues[id]["C6"],plane,True))
                residues[id]["N7"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["N7C5"],json["G"]["angles"]["C6C5N7"],residues[id]["C6"],residues[id]["C5"],plane,False))
                residues[id]["C8"].set_coord(getPointFromAngleDistAndPlane(json["G"]["distance"]["C8N7"],json["G"]["angles"]["C5N7C8"],residues[id]["C5"],residues[id]["N7"],plane,True))
        else:
            residues[id]["N1"].set_coord(getPointFromDistAnd2Angles(json["rybosine"]["distance"]["C1'N"],json["rybosine"]["angles"]["C2'C1'N1"],json["rybosine"]["angles"]["O4'C1'N1"],residues[id]["C2'"], residues[id]["O4'"], residues[id]["C1'"],False))
            if residue.get_resname() == 'C':
                residues[id]["C2"].set_coord(getForthPoint(json["C"]["distance"]["C2N1"], json["C"]["angles"]["C1'N1C2"],torsionAngles[id][6],residue["O4'"],residue["C1'"],residue["N1"]))
                plane = planeCoords(residues[id]["C1'"].get_coord(),residues[id]["N1"].get_coord(),residues[id]["C2"].get_coord())
                residues[id]["N3"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["N3C2"],json["C"]["angles"]["N1C2N3"],residues[id]["N1"],residues[id]["C2"],plane,True))
                residues[id]["O2"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["O2C2"],json["C"]["angles"]["N1C2O2"],residues[id]["N1"],residues[id]["C2"],plane,True))
                residues[id]["C4"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["C4N3"],json["C"]["angles"]["C2N3C4"],residues[id]["C2"],residues[id]["N3"],plane,True))
                residues[id]["N4"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["C4N4"],json["C"]["angles"]["N3C4N4"],residues[id]["N3"],residues[id]["C4"],plane,True))
                residues[id]["C5"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["C5C4"],json["C"]["angles"]["N3C4C5"],residues[id]["N3"],residues[id]["C4"],plane,True))
                residues[id]["C6"].set_coord(getPointFromAngleDistAndPlane(json["C"]["distance"]["C6C5"],json["C"]["angles"]["C4C5C6"],residues[id]["C4"],residues[id]["C5"],plane,True))
            if residue.get_resname() == 'U':
                residues[id]["C2"].set_coord(getForthPoint(json["U"]["distance"]["C2N1"], json["U"]["angles"]["C1'N1C2"],torsionAngles[id][6],residue["O4'"],residue["C1'"],residue["N1"]))
                plane = planeCoords(residues[id]["C1'"].get_coord(),residues[id]["N1"].get_coord(),residues[id]["C2"].get_coord())
                residues[id]["N3"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["N3C2"],json["U"]["angles"]["N1C2N3"],residues[id]["N1"],residues[id]["C2"],plane,True))
                residues[id]["O2"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["O2C2"],json["U"]["angles"]["N1C2O2"],residues[id]["N1"],residues[id]["C2"],plane,True))
                residues[id]["C4"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["C4N3"],json["U"]["angles"]["C2N3C4"],residues[id]["C2"],residues[id]["N3"],plane,True))
                residues[id]["O4"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["C4O4"],json["U"]["angles"]["N3C4O4"],residues[id]["N3"],residues[id]["C4"],plane,True))
                residues[id]["C5"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["C5C4"],json["U"]["angles"]["N3C4C5"],residues[id]["N3"],residues[id]["C4"],plane,True))
                residues[id]["C6"].set_coord(getPointFromAngleDistAndPlane(json["U"]["distance"]["C6C5"],json["U"]["angles"]["C4C5C6"],residues[id]["C4"],residues[id]["C5"],plane,True))

    newStructure = PDBFile.copy()
    tmp = 0
    for id, residue in enumerate(newStructure.get_residues()):
        if id in originalIndexes:
            for atom in residue.get_atoms():
                atom.set_coord(residues[tmp][atom.get_name()].get_coord())
            tmp+=1

    io=PDBIO()
    io.set_structure(newStructure)
    io.save("1ajf_from_torsion.pdb")  