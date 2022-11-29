# HIVtreeSimulations
## Introduction

This program was developed to study assess methods that infer HIV latent integration times.
This program simulates phylogenies and sequences for latent and non-latent viruses.
HIV\_treeSim simulates phylogenies based on models of within-host viral dynamics.
dnaSim simulated DNA sequences on trees produced by HIV\_treeSim.
No evolution occurs while sequences are latent.

Warning: The tree simulation program can use a lot RAM, depending on the simulation parameters. 
The memory usage scale approximately linearly with the mLBlood parameters in the control file. 
Users should not run large simulations on a personal computer. 
The memory usage is printed to the screen by default. 

## Download and Compile the Program
To clone the program from GitHub

```
git clone https://github.com/nage0178/HIVtreeSimulations.git
```

To compile the program,

```
cd ~/HIVtreeSimulations/treeSim
make
cd ~/HIVtreeSimulations/dnaSim
make
```

## Running the HIV_treeSim program

First you must run treeSim to generate the tree.
The user must specify a file with the sample times as either a command line argument, as in the first example, or in the control file.

```
./treeSim -i sampleTimes.csv
```

or

```
./treeSim -c control.ctl
```

## Options
Program options:
| Option                     | Description                                                                       |
| -------------------------- | --------------------------------------------------------------------------------- |
| **-c**                     | Name of control file. Does not have an effect if loading from a checkpoint.       |
| **-h**                     | Prints the program options and exits.		                                 |
| **-i**                     | Input file with sampling times. In csv file format without a header. First column is the time, second is number of active sequences to sample, third is the number of latent sequences to sample. Does not have an effect if loading from a checkpoint. |
| **-l**                     | Loads checkpoint files and restarts the run. Requires the prefix of the checkpoint files. If there is no prefix, use ""|
| **-o**                     | Output file prefix. |
| **-p**                     | Frequency of printing counts and memory usage. If 0, counts are not printed. Otherwise, every p events the counts are printed. |
| **-r**                     | Restart if population goes extinct. 0 if no, otherwise yes.  |
| **-s**                     | Starting seed. Does not have an effect if loading from a checkpoint. |
| **-t**                     | Time to checkpoint the program. The program cannot load a checkpoint and have another checkpoint later. |
| **-v**                     | Volume of the simulation in mL. Does not have an effect if loading from a checkpoint. |



## Control File

| Keyword            | Description                                                                |
| ------------------ | -------------------------------------------------------------------------- |
| **sample\_times**  | csv file with the sample times and number of viruses to sample |
| **seed**           | Random seed to the simulation |
| **volume**         | Simulation parameter: Blood volume to simulate in mL |
| **lamdba**         | Simulation parameter: Birth rate of uninfectied cell (cell per mL per day)  |
| **kappa**          | Simulation parameter: Transition rate from unifected to actively infected cells (mL per virion per day) |
| **d**              | Simulation parameter: Death rate of unifected cell (per day) |
| **delta**          | Simulation parameter: Death rate of actively infected cells (per day) |
| **pi**             | Simulation parameter: Viral birth rate (virions per cell per day) |
| **c**              | Simulation parameter: Viral clearance rate (per day) |
| **eta**            | Simulation parameter: Proportion of newly infected cells that are latent |
| **alpha**          | Simulation parameter: Rate of activation of replication-competent, latent cells (per day |
| **gamma**          | Simulation parameter: Proportion of viruses that are replication-incompetent |
| **sigma**          | Simulation parameter: Death rate of actively infected cells (per day) |
| **tau**            | Simulation parameter: Death rate of latent, replication-incompetent cells |
| **ARTstart**       | Simulation parameter: Time of at which ART is initiated. Kappa becomes zero at this time |


Options that can be specified in the control file or the command line (seed, volume, and sample\_times) can only be specified in one places.
If either is specified in both locations, the program will not run.
While there are default values for each of the simulation parameters, they are not "preferred" or "best" values.
The user should determine what parameters match their simulation needs based on the expected population sizes from the viral dynamic model.
The initial conditions and parameter values are printed to the log.txt file.
The user may wish to use the provided R script to see the expected population sizes given different parameter values.
For a full description of the model, please refer to the paper.

## Sample Times Files

The sample times file in in the format of a csv file without a header.
This file is required for the program to run.
The first column is the sample time in days.
The second column is the number of non-latent sequences sample.
The third column is the number of latent sequence to sample.
The sample times must be increasing in the file.
There can be no blanks in the file, but you can sample no sequences by having a zero in the sample number column.
The simulation will run until the last sample time provided in the sample times file.


## Running the dnaSim program

The dnaSim program is run after the tree simulation program.
It takes as input the output from the tree simulator.
To run the program,

```
./dnaSim -l LatentStates.txt -t SampledTree.txt
```
where the file names match the output from the tree simulator.
The default model is Jukes-Cantor model without among site rate variation, with a default mutation rate is 3.6 x 10 ^ -5.
By default, the program will simulate 1000 bases from a starting sequence drawn from the stationary distribution of the substitution model.

## Options
Program options:
| Option     | Description                                                                       |
| -------------------------- | ----------------------------------------------------------------- |
| **a**      |	Alpha for the +gamma model of among site rate variation |
| **b**      | 	number of bases. If the number of bases in the alignment does not match the input, the number of based in the alignment will be used. |
| **h**      | 	Prints options message |
| **f**      | 	Stationary frequencies followed by rate matrix parameters. The correct format is statFreq\_A:statFreq\_C:statFreq\_G:statFreq-rateParameter\_AC:rateParameter\_AG:rateParameter\_AT:rateParameter\_CG:rateParameter\_CT:rateParameter\_GT . Numbers must be in decimal format with a digit after the decimal. |
| **i**      | 	Input fasta file with ancestral DNA sequence. |
| **l**      | 	Latent history file from HIV\_treeSim |
| **o**      | 	Output file name |
| **r**      | 	Outgroup sequence for rooting |
| **s**      | 	Seed |
| **t**      | 	Tree file from HIV\_treeSim |
| **u**      | 	Substitution rate |

## Output Files
The tree simulator generates two files which are used in the tree simulation.
Additional intermediate files are generated if a checkpoint is used.
The file ending in "SampledTree.txt" has the tree topology and "LatentStates.txt" records which branches are latent.
Each row in the latent states file is in the format "node index : total time latent in the branch, sample time, and is latent".
The node index is an unique identification for each node that serves as the name in the tree file produced by the tree simulator.
The sample time is zero for internal nodes.
The is latent will be 0 if the node is a sampled latent tip and zero otherwise.
This file may be useful to check if there are internal branches that are latent.
The tree file is in newick format, but may not open in some tree viewers due to the lack of a semi-colon at the end.
The output files from the DNA simulator contain the same information and will open in many tree viewers.

The output files from the DNA simulator include a fasta file, a tree file in newick format, and possible a tree file with an outgroup if the outgroup option was used.
The sequences names are in the format "Node\_isLatent\_nodeIndex\_timeAtTransitionToLatency\_SampleTime".
Latent sequences have a one for the "isLatent" field.
Other sequences have a zero.
The node index is unique for each sequence.
The time at transition to latency will be equal to the sample time if the sequence is not latent.
The outgroup sequence name is "outgroup".
The branch length to the outgroup will be equal to the time of the last sample.
