# ReaxMD

ReaxMD is a program that analyzes trajectory files and displays bond formation/breaking events. Interatomic distances, dij, are compared to the bond length cutoff $D^{cutoff}_{ij}=\alpha*(r_i + r_j)$, where $\alpha$ was optimized to 1.30 for solvated and gas-phase dynamics and -a 1.0 is suggested for solid state. If Dij is lower than Bij, then there is a bond. If Dij increases above Bij, the bond is broken. If Dij decreases below Bij, a bond is formed.

It is also possible to analyze the Wannier center's (aka. "wanniers"). The distance between the wanniers and the atoms considers a tolerance of $d^\text{wannier}_i = \beta*(r_i + \gamma)$. $\beta$ and $\gamma$ were obtained empirically by analyzing trajectory files without reactions and minimizing the amount of wanniers movement.

Covalent radii were obtained from "B. Cordero et. al., Dalton Transactions, 2008, 21. https://doi.org/10.1039/B801115J" and included in the file "covalent_radii.ssv". An alternative covalent radii data file can be used with "-d path/to/file", but the ssv format (space separated values) must be kept as " Z    Element    r".

## INSTALATION
### Option 1: Download executable
- The ```reaxmd``` executable is provided. Copy it into your ```home/.local/bin``` folder and execute it from the command line as
```
reaxmd trajectory.xyz A B C
```

### Option 2: Download python script
- Download the file ```reaxmd.py``` and executing with Python3 as
```
python3 reaxmd.py trajectory.xyz A B C
```
- The python file can be compiled with PyInstaller (Recommended)

## USAGE
### Check version
```
reaxmd -v
```
Latest version: 260505

### Basic usage
```
reaxmd trajectory.xyz A B C
```
- ```trajectory.xyz```: trajectory file in .xyz format, or path to it.
- ```A B C```: cell dimensions in atomic units. ReaxMD was programmed for periodic boundary conditions. Single cells can be analyzed by specifying dimensions larger to the actual ones, preventing bonding between frontier atoms.
- It is suggested to include ```> trajectory.txt``` at the end to save the output into a text file.

### Batch analysis
```
reaxmd *.xyz A B C
```
- Several trajectory files can be provided to analyze with the same cell dimensions.
- To save the outputs in separated files, you can execute a command like (example for Bash):
```
for trajectory in *.xyz; do reaxmd $trajectory A B C > ${trajectory%.xyz}; done
```

### Atomic units
- Use the flag ```-au``` to indicate that the cell dimensions are given in atomic units:
```
reaxmd trajectory.xyz A B C -au
```

Some examples are provided in the folder "examples". All of those files were used for testing the program during the development. In every file, the cell dimensions are given in the second line.

### Example 1: Default usage
```
reaxmd bp3n/bp3n-tra3.xyz 12.7 12.7 12.7
```
Output
```
trajectory file: bp3n-tra3.xyz
atoms: 162
steps: 50
step 39: bond (34)O-H(73) is formed.
step 41: bond (2)O-H(73) is broken.
step 43: bond (35)O-H(138) is formed.
step 43: bond (39)O-H(91) is formed.
step 44: bond (34)O-H(138) is broken.
NET REACTION:
2 HO + C14H12O3  ->  C14H11O3 + H3O2
```
- The atomic IDs in parentesys corresponds with the order of the atoms in the trajectory file, starting from 1.

### Example 2: Save output
```
reaxmd bp4l-tra.xyz 12.7 12.7 12.7 -e C > bp4l-tra.log
```
then open the file bp4l-tra.log with any external file. In a bash command line you can use ```cat```:
```
cat bp4l/bp4l-tra.log
```
Output:
```
trajectory file: bp4l/bp4l-tra.xyz
atoms: 157
steps: 1900
step 1205: bond (62)C-O(16) is formed.
step 1471: bond (62)C-S(1) is broken.
NET REACTION:
CS  ->  CO
```
- There are several reactions in this trajectory file, but only the reactions involving C atoms are shown.
- In this case, the reaction "H2O(16) + Ar(62)-S(1)O3 -> Ar(62)-H2O(16)+ + S(1)O3-" is observed.

### Example 3: Elements
Use the flag ```-e E1-E2``` to specify elements:
- Give a single element ```E1``` to track reactions between ```E1``` with all other atoms. </br>
Command:
```
reaxmd bp3n/bp3n-tra.xyz 12.7 12.7 12.7 -e C | grep C
```
Output:
```
step 1236: bond (60)C-O(11) is formed.
step 3436: bond (62)C-O(41) is formed.
step 15422: bond (60)C-C(62) is broken.
```
- Give two elements ```E1-E2``` to track a particular kind of bond. The same element can also be provided as ```E1-E1```.</br>
Command:
```
reaxmd bp3n/bp3n-tra.xyz 12.7 12.7 12.7 -e C-O | grep C
```
Output:
```
step 1236: bond (60)C-O(11) is formed.
step 3436: bond (62)C-O(41) is formed.
```
- All of these cases can be used to 1. reduce the length of the output, and 2. reduce the execution time.

### Example 4: Wannier centers
```
reaxmd bp3n/bp3n-ion3.xyz 12.7 12.7 12.7 -o w > bp3n/bp3n-ion3.log
```
Output:
```
trajectory file: bp3n-ion3.xyz
atoms: 616
steps: 50
...
step 44: bond (34)O-H(138) is broken.
step 3: Wannier center 295 approached C(63).
step 3: Wannier center 522 approached C(63).
step 4: Wannier center 198 separated from H(160).
step 4: Wannier center 425 separated from H(160).
step 7: Wannier center 217 separated from C(57).
step 7: Wannier center 444 separated from C(57).
step 9: Wannier center 198 approached H(160).
step 9: Wannier center 295 separated from C(63).
step 9: Wannier center 385 approached C(63).
step 9: Wannier center 425 approached H(160).
step 9: Wannier center 522 separated from C(63).
step 9: Wannier center 612 approached C(63).
...
NET REACTION:
2 HO + C14H12O3  ->  C14H11O3 + H3O2
```
- Same reaction as in the previous example, but including the Wannier centers.
- Save the output and look for specific symbols and IDs using external (notes, kate, grep).

## Finally
This is the basic and intermediate usage. Don't hesitate to write me an email if you have further questions. 
-Luis Alvarez F.-
alvarez@theochem.uni-hannover.de

### Full help: 
```
reaxmd -h
```
Output
```
usage: reaxmd [-h] [-au] [-a ALPHA] [-d DELTA] [-b BETA] [-g GAMMA] [-e ELEMENT] [-o OPTIONS] [-r RADII] [-ro] [-s STRIDE] [-t TIMESTEP] [-v] [trajectory_files [trajectory_files ...]] A B C

positional arguments:
  trajectory_files      Path(s) to trajectory file(s).
  A                     A length of the cell in Angstrom.
  B                     B length of the cell in Angstrom.
  C                     C length of the cell in Angstrom.

optional arguments:
  -h, --help            show this help message and exit
  -au, --atomic_units   Specifies that the cell side lengths are in atomic units.
  -a ALPHA, --alpha ALPHA
                        Alpha value for bond lengths as alpha * (r_i + r_j) +
                        delta. Default value: 1.30
  -d DELTA, --delta DELTA
                        Delta value for bond lengths as alpha * (r_i + r_j) +
                        delta. Default value: 0.00
  -b BETA, --beta BETA  Beta value for Wannier centers proximities as beta *
                        r_i + gamma. Default value: 1.30
  -g GAMMA, --gamma GAMMA
                        Gamma value for Wannier centers proximities as beta *
                        r_i + gamma. Default value: 0.40
  -e ELEMENT, --element ELEMENT
                        "-e El" analyzes only the bonds involving element El.
                        "-e E1-E2" analyzes only bonds between elements E1 and
                        E2.
  -o OPTIONS, --options OPTIONS
                        "-o hmwf", where: h writes a trajectory file
                        "trajectory_file-H.xyz", where the atoms involved in
                        reactions are replaced by Si and the Wannier centers
                        (include "w" as option) by He. m writes a file
                        "trajectory-molecules.txt" with the reactive
                        molecules. w analyzes the Wannier centers. f filters
                        equilibrium reactions within 100 MD steps.
  -r RADII, --radii RADII
                        Full path to alternative covalent radii data file.
  -ro, --reaction_only  Print reaction only.
  -s STRIDE, --stride STRIDE
                        Jumps the specified amount of steps.
  -t TIMESTEP, --timestep TIMESTEP
                        Prints the time in ps next to the step number.
                        Timestep = (MD timestep) x (MD steps between
                        trajectory steps)
  -v, --version         Print version information
```

