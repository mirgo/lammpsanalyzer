# imports
import numpy as np
import matplotlib.pyplot as plt
import argparse, tqdm, sys, re, os

from supplementary import readlammpstrjpositions, rog, avdis

parser = argparse.ArgumentParser(description='Take in a LAMMPSTRJ file and output a variety of analyses.', epilog='Created July 31, 2021. Hope it helps! - Mihir')

parser.add_argument('infile', type=argparse.FileType('r'), help='.LAMMPSTRJ file to be analyzed. Specify with relative folder path.')
parser.add_argument('outfolder', help='Folder name to place output contents in. Example "output1" creates directory /output1.')
parser.add_argument('-q', action='store_true', help='quiet. Indicating flag will remove progress updates in the command line.')
parser.add_argument('-r', action='store_true', help='radius of gyration calculation. Will output chart over timesteps.')
parser.add_argument('-d', action='store_true', help='average self-displacement calculation. Will output chart over timesteps.')

parser.add_argument('-all', action='store_true', help='indicate ALL analyses to be outputted')

args = parser.parse_args()

# Create output folder, give option if about to overwrite
try:
    os.mkdir(args.outfolder)
except FileExistsError:
    prompt = False
    while not prompt:
        sys.stdout.write('Folder already exists, continue and overwrite contents? [y/n]: ')
        choice = input().lower()
        if choice == 'y' or choice == 'n':
            prompt = True
    if choice == 'n':
        sys.exit()

# Extract position matrix and number of atoms
pos, natoms = readlammpstrjpositions(args.infile)
print('Loaded Information!')

# Radius of Gyration calculation
if args.r or args.all:
    a = plt.plot(np.linspace(0, len(pos)/natoms, int(len(pos)/natoms)), rog(pos, natoms))
    plt.title('Radius of Gyration')
    plt.xlabel('Timestep')
    plt.ylabel('RoG')
    plt.savefig(args.outfolder + '/rog.png')
    print('Radius of Gyration Computed and Saved!')
    plt.clf()

if args.d or args.all:
    avgdp = avdis(pos, natoms)
    times = np.linspace(1, len(avgdp), len(avgdp))
    b = plt.plot(times, avgdp, linewidth=.5)
    plt.xlabel('Timestep')
    plt.ylabel('Average Displacement')
    plt.title('Average Displacement of Atoms per Timestep')
    plt.savefig(args.outfolder + '/avgdp.png')
    print('Average Displacement Computed and Saved!')
    plt.clf()
