# Molecule Simulation Analyzer: lammpsanalyzer

## Description

This repo contains files designed to parse through LAMMPS trajectory files, which contain
information on molecule positions over a series of timesteps.

**[Link to Short Video Description](https://youtu.be/H2jwzjH7r2Y)**

[LAMMPS](https://www.lammps.org/) is a molecule simulator, taking in structures of
molecules and allowing simulations under a variety of completely customizable conditions
(ex., implementing Lennard-Jones potential, and such).

I've designed code in this repo to output visuals and videos containing statistical analyses
of the trajectories of these molecules, to give insight to multiple features.

#### What is a lammpstrj file?

Take a look at `sampleinputs/sample24atoms.lammpstrj` to get a view! The other file might be too large to look at.

A lammpstrj file, also known as a LAMMPS trajectory file, describes the movement of a molecule over a set amount of timesteps, and classifies those movements with a series of x,y,z coordinates (per atom). You can probably imagine as timesteps and atoms increase, the complexity of the problem increases, as there are *n = timestep x atoms* iterations to go through.

## Navigation

`sampleinputs` contains example files to test. The structure of these files are consistent outputs by LAMMPS, and so a static comprehensive program in python can read it reliably.

There are two sample files, which have greatly different sizes, and so execution time will differ between them, which is why I appended a loading bar to the executable.

Input file | Number of Atoms | Timesteps | *n* Steps
--- | --- | --- | ---
`sample24atoms.lammpstrj` | 24 | 500 | 12000
`sample724atoms.lammpstrj` | 724 | 2000 | 1448000



## Installation Procedure

Assuming working on a linux machine. I used Ubuntu 20.02 OS.

1. Install python and git
```
sudo apt-get update
sudo install python3
sudo install git
```

2. Clone repository and navigate in
```
git clone []
cd /lammpsanalyzer
```

3. Install pip and prerequisite packages.
```
sudo install pip
pip install numpy matplotlib argparse tqdm gif pillow peakutils moviepy scipy
```

## Execution

Running `python3 analyzer.py -h` should give a help prompt, detailing all options
for execution in the command line.

`analyzer.py` is designed to be highly customizable and selective towards the outputs and analyses desired for a given molecule, because redundant calculations can become very costly with larger molecules. Thus, to output specific analyses, there are flag indicators that have to specified along with execution of the script for each analysis.

General execution follows the following format:
```
usage: analyzer.py [-h] [-gyration] [-displacement] [-covariance formatinteger] [-rmsd] [-all] infile outfolder
```

Where `infile` refers to the relative location for the input lammpstrj file, which should be (assuming you're in `/lammpsanalyzer`) in the `sampleinputs/` folder, corresponding either to `sampleinputs/sample24atoms.lammpstrj` or `sampleinputs/sample742atoms.lammpstrj`. **It is required.**

`outfolder` refers to the folder destination for the outputs desired. For each analysis indicated/flagged, there will be an output of some representation (i.e., a graph, gif, or mp4), that will have a name with the analysis (i.e., rsmd.png). **It is required.**

The following are optional arguments (at least 1 would be useful to be run):
Analysis | Flag | Additional Notes
--- | --- | ---
Radius of Gyration | -gyration |
Average Self-Displacement | -displacement |
Distance Covariance Matrix | -covariance (followed by 0 or 1) | Following with a 0 will create a GIF of the first 100 timesteps, to save memory. A 1 will create a MP4 of the entire simulation.
Root Mean Square Displacement | -rmsd |
Distance Covariance Rolling Average Peak | -peakaverage |
**Indicate all above analyses** | -all | Doesn't require any other flags to be indicated. Will not harm execution if so. Will produce 100 timestep GIF option for covariance.

Some examples with `sampleinputs/sample24atoms.lammpstrj` as infile and `output24` as outfolder:

1. To run all analyses:
```
python3 analyzer.py -all sampleinputs/sample24atoms.lammpstrj output24
```

2. To run just the radius of gyration and root mean square displacement:
```
python3 analyzer.py -gyration -rmsd sampleinputs/sample24atoms.lammpstrj output24
```

3. To run just distance covariance matrix as mp4:
```
python3 analyzer.py -covariance 1 sampleinputs/sample24atoms.lammpstrj output24
```

### Analyses Available

- Radius of Gyration

Describes the root mean square distance aggregate of all atoms from the center of mass for the molecule. Higher values indicate a more spread out polymer, while lower values indicate a more bundled-up, compact molecule.

- Average Self-Displacement

Takes the displacement of atoms through every timestep, averages per timestep, and outputs over time. Higher values indicate a molecule moving rapidly, while lower values descibe a less-fluctuating molecule (typically high suspension or compact conditions).

- Distance Covariance Matrix

Describes the association of movement for atoms between each other through every timestep. Higher values describe relational (similar) movement, while lower values describe limited association.

- Root Mean Square Displacement

Describes the movement of the molecule. Hails from Brownian random motion and diffusion.

- Distance Covariance Rolling Average Peak

An additional measure to the distance covariance matrix, but identifying and quantifying peaks - which would be bundles where a lot of atoms are moving in unison.

## References

- https://realpython.com/python-scipy-fft/
- https://en.wikipedia.org/wiki/Mean_squared_displacement
- https://docs.lammps.org/compute.html
- https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30
- https://datascienceplus.com/understanding-the-covariance-matrix/
