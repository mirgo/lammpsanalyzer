# imports
import numpy as np
import matplotlib.pyplot as plt
import argparse, tqdm, sys, re, os, gif
import moviepy.editor as mp

from supplementary import readlammpstrjpositions, rog, avdis, animatedmovementcovmat, squaredisplacement

# Command Line Guide and Information
parser = argparse.ArgumentParser(description='Take in a LAMMPSTRJ file and output a variety of analyses.', epilog='Created July 31, 2021. Hope it helps! - Mihir')

parser.add_argument('infile', type=argparse.FileType('r'), help='.LAMMPSTRJ file to be analyzed. Specify with relative folder path.')
parser.add_argument('outfolder', help='Folder name to place output contents in. Example "output1" creates directory /output1.')

parser.add_argument('-gyration', action='store_true', help='radius of gyration calculation. Will output chart over timesteps.')
parser.add_argument('-displacement', action='store_true', help='average self-displacement calculation. Will output chart over timesteps.')
parser.add_argument('-covariance', type=int, metavar='formatinteger', help='distance covariance matrix. Will output either gif or mp4 illustrating movment covariance matrices over timesteps. For format enter 0 for "gif" or non-0 for "mp4". NOTE: "gif" option saves first 100 iterations, to keep it short. "MP4" takes all information.')
parser.add_argument('-rmsd', action='store_true', help='root mean square displacement calculation. Will output chart over timesteps.')

parser.add_argument('-all', action='store_true', help='indicate ALL analyses to be outputted.')

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
if args.gyration or args.all:
    plt.plot(np.linspace(0, len(pos)/natoms, int(len(pos)/natoms)), rog(pos, natoms))
    plt.title('Radius of Gyration')
    plt.xlabel('Timestep')
    plt.ylabel('RoG')
    plt.savefig(args.outfolder + '/rog.png')
    print('Radius of Gyration Computed and Saved!')
    plt.clf()

# Average Displacement Plot
if args.displacement or args.all:
    avgdp = avdis(pos, natoms)
    times = np.linspace(1, len(avgdp), len(avgdp))
    plt.plot(times, avgdp, linewidth=.5)
    plt.xlabel('Timestep')
    plt.ylabel('Average Displacement')
    plt.title('Average Displacement of Atoms per Timestep')
    plt.savefig(args.outfolder + '/avgdp.png')
    print('Average Displacement Computed and Saved!')
    plt.clf()

# Creating Movement Covariance GIF
if args.covariance == 0 or args.all:
    frames = animatedmovementcovmat(pos, natoms, 0)
    # Capping gif with 100 frames to 30 seconds
    durationofframe = 30
    gif.save(frames, args.outfolder + '/movementcm.gif', duration=durationofframe, unit='s', between='startend')
    print('Movement Covariance GIF Computed and Saved!')

# Creating Movement Covariance MP4
if args.covariance != 0 and args.covariance != None:
    frames = animatedmovementcovmat(pos, natoms, 1)
    # Making total time so that each frame has .2 seconds to display
    durationofframe = int(len(pos)/natoms) * .2
    # Create GIF to convert into MP4
    gif.save(frames, args.outfolder + '/movementcm.gif', duration=durationofframe, unit='s', between='startend')
    # Write MP4
    clip = mp.VideoFileClip(args.outfolder + '/movementcm.gif')
    clip.write_videofile(args.outfolder + '/movementcm.mp4')
    os.remove(args.outfolder + '/movementcm.gif')
    print('Movement Covariance MP4 Computed and Saved!')

# Creating Root Mean Square Displacement
if args.rmsd or args.all:
    plt.plot(np.linspace(1, len(pos)/natoms, int(len(pos)/natoms)-1), squaredisplacement(pos, natoms))
    plt.title('Root Mean Square Displacement over Time')
    plt.xlabel('Timestep')
    plt.ylabel('RMSD')
    plt.savefig(args.outfolder + '/rmsd.png')
    print('RMSD Computed and Saved!')
    plt.clf()
