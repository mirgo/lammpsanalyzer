import numpy as np
import tqdm
import re
import matplotlib.pyplot as plt

# Reformat lammpstrj script into numpy matrix of positions
# Matrix width 3 for x,y,z coordinates
# Matrix length number of atoms * timesteps
def readlammpstrjpositions(fileloc):
    #f = open(fileloc, 'r')
    b = 0
    index = 0
    addspace = False
    positions = []
    natoms = 0
    for line in tqdm.tqdm(fileloc):
        # Signal from atom count, append chunk of memory for next timestep
        if addspace:
            positions.extend(np.zeros((int(line), 3)))
            natoms = int(line)
            addspace = False
        # Line contains atom count
        if 'ITEM: NUMBER OF ATOMS' in line:
            addspace = True
        # Avoiding header space, take in coordinates
        if b >= 9:
            v = ([m.start() for m in re.finditer(' ', line)])[:4]
            positions[index][0] = float(line[v[0]:v[1]])
            positions[index][1] = float(line[v[1]:v[2]])
            positions[index][2] = float(line[v[2]:v[3]])
            index+=1
        b+=1
        if b == (9 + natoms):
            b=0
    #f.close()
    # Return matrix
    return np.array(positions), natoms

def rog(positions, natoms):
    rogs = []
    mat = positions
    # Chunking position block into each timestep
    for i in tqdm.tqdm(range(int(len(mat)/natoms))):
        allr = []
        beg = i*natoms
        # Calculate center of mass across each axis (x,y,z)
        comxyz = (np.average(mat[:][beg:(i+1)*natoms], axis=0))
        # Iterate through all atoms in timestep, taking square of distance from CoM
        for j in range(natoms):
            allr.append(((mat[beg+j][0] - comxyz[0])+(mat[beg+j][1] - comxyz[1])+(mat[beg+j][2] - comxyz[2]))**2)
        # Record the square of the sum of every atom's distance averaged
        rogs.append(np.sqrt(sum(allr)/natoms))
    return np.array(rogs)

# Distance formula between atoms that have x,y,z coordinates
def dist(f, s):
    return np.abs(np.sqrt(((f[0]-s[0])**2) + ((f[1]-s[1])**2) + ((f[2]-s[2])**2)))

# Average Displacement Per Atom Over Time
def avdis(positionarray, numberofatoms):
    atoms = numberofatoms
    pos = positionarray
    avgdp = []
    nf = int(len(pos)/atoms)
    # Chunk into timesteps
    for j in tqdm.tqdm(range(nf-1)):
        displacements = []
        # Get all self-displacement values with next timestep
        for i in range(atoms):
            displacements.append(dist(pos[(atoms*j) + i], pos[(atoms*j) + i+atoms]))  ## from position array, extract distance between same atom of original and next timestep
        # Record average of all those displacements
        avgdp.append(np.average(displacements))
    return avgdp
