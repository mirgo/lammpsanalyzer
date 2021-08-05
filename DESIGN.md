# Design

This project was mainly an exploration into applying python as a tool for data analysis and visualization, exploring topics ranging from algorithm efficiency to command line interpretability to capabilities of visual generation.

I tried to condense my main design aspects into categories:
1. Overall Structure
2. Input and Output Handling
3. Command Line Interface
4. Analysis Creation
5. Algorithm Design

## Overall Structure

Two main files work together here, the main executable `analyzer.py` and the supplementary file `summplementary.py`.

`analyzer.py`, as designed to be the main script the user interacts with, has the most readable components of the program - it contains the argument information, and drop down if statements for each analysis option, clearly outlining the path taken if invoked at the command line.

`supplementary.py` contains the bulk of the data analysis and algorithmic operations. These are harder to immediately understand, but invoked as functions have clear outputs (like *movementcovariancematrix* or *readlammpstrjpositions*).

## Input and Output Handling

The rigid format of the input .lammpstrj files makes a rigid counterpart program for reading and extracting information very easy, hence a tight but efficient *readlammpstrjpositions* function. More on efficiency later.

To design an easy output generation method, I decided a folder with each analysis independently appended would be the best, to ensure maximum customizability. This is over, say a dashboard, where every analysis is created every time. The problem lies with increasing problem size and time required for each additional analysis.

The modules `matplotlib.pyplot`, `PIL`, `gif`, and `moviepy.editor` were extremely useful for generating output visuals. `matplotlib.pyplot` offered simple graph creation and handling, and the other modules used a frame decorator to sequentially create animations.

## Command Line Interface

`argparse` was a very neat module to create a rigid command line format. Through that I could create argument flags and descriptions, which I could follow with a if-then structure to `analyzer.py`.

With the help flag (-h or --help) a description of all variables and example usage message will appear in the command line.

Another very useful feature I implemented was the use of loading bars and print commands following the invocation of every if(flag) route, to indicate to the user through the command line of progress. It was very easy with the `tqdm` module to wrap around for loops and reflect on the progress/speed of an each analysis.

## Analysis Creation

All of the calculations for these analyses follow mathematical formulas described from following sources:
- [Radius of gyration](https://en.wikipedia.org/wiki/Radius_of_gyration)
- [Movement covariance matrix](https://en.wikipedia.org/wiki/Covariance_matrix)
- [Root mean sqaure distance](https://en.wikipedia.org/wiki/Mean_squared_displacement)
- [Fourier transform](https://numpy.org/doc/stable/reference/routines.fft.html#implementation-details)

And for efficient computing, most data transformations were done through numpy functions. 

## Algorithm Design

Creating a fast enough implementation was perphaps the most difficult part of the process. Through the use of numpy, allocating larger chunks of space at a time, and the minimal amount of runthroughs, I've created a fairly efficient result.

The following functions are in supplementary.py.
>*readlammpstrjpositions* utilizes `extend` over `append` to create larger chunks of list space to insert atom information into. Each line is iterated through once, and it retains flexibility for addressing any timestep amount.

>*rog* utilizes numpy `np.average` function to bring that calculation speed closer to C code execution. Numpy position matrix iterated through once.

>*movementcovariancematrix* utilizes all commands through numpy for covariance and linear algebra operations.
