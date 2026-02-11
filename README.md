# ReaxMD

mdanalyzer is a program that analyzes trajectory files and shows the bond formation and breaking events. Interatomic distances, Dij, are compared to the bonding distance Bij=(alpha)*(covalent_radii_1 + covalent_radii_2), where alpha=1.3. If Dij is lower than Bij, then there is a bond. If Dij increases above Bij, the bond is broken. If Dij decreases below Bij, a bond is formed.

It is also possible to analyze the Wannier center's (aka. "wanniers"). The distance between the wanniers and the atoms considers a margin of (beta)*(covalent_radii + gamma). Beta and gamma were obtained empirically by analyzing trajectory files without reactions and minimizing the amount of wanniers movement.

Covalent radii were obtained from "B. Cordero et. al., Dalton Transactions, 2008, 21. https://doi.org/10.1039/B801115J" and included in the file "covalent_radii.ssv". An alternative covalent radii data file can be used with "-d path/to/file", but the ssv format (space separated values) must be kept as " Z    Element    r".



#### INSTALATION
Although it can be executed directly with python ($python (or python3) mdanalyzer.py trajectoryfile.xyz A B C), an instalation with pyinstaller is recomended as follows:

mkdir ~/bin
pyinstaller --onefile mdanalyzer.py
mv dist/mdanalyzer ~/bin/



#### USAGE
$ mdanalyzer [-h] [-a ALPHA] [-b BETA] [-e ELEMENT] [-g GAMMA] [-o OPTIONS] [-r RADII] [-ro] [-t TIMESTEP] [--version] trajectory_file a b c

positional arguments:
  trajectory_file       name of the trajectory file.
  a                     A length of the cell in Angstrom.
  b                     B length of the cell in Angstrom.
  c                     C length of the cell in Angstrom.

optional arguments:
  -h, --help            show this help message and exit
  -a ALPHA, --alpha ALPHA
                        Alpha value for bond lengths as a*(r1+r2). Default value: 1.3
  -b BETA, --beta BETA  Beta value for Wannier centers proximities as b*(ri+g). Default value: 1.25
  -e ELEMENT, --element ELEMENT
                        "-e El" analyzes only the bonds involving element El. "-e E1-E2" analyzes only bonds between elements E1 and E2.
  -g GAMMA, --gamma GAMMA
                        Gamma value for Wannier centers proximities as b*(ri+g). Default value: 0.31
  -o OPTIONS, --options OPTIONS
                        "-o mtwf", where:
                        m writes a file "trajectory-molecules.txt" with the reactive molecules.
                        t writes a trajectory file "trajectory_file-H.xyz", where the atoms
                        involved in reactions are replaced by Si and the Wannier centers (if wt or tw ) by He.
                        w analyzes the Wannier centers.
                        f filters the Wannier centers close to
                        reactive atoms. Use with option "w"
  -r RADII, --radii RADII
                        Alternative full path to covalent radii data.
  -ro, --reaction_only  Print reaction only
  -t TIMESTEP, --timestep TIMESTEP
                        Prints the time in ps next to the step number. Timestep = (MD timestep) x (MD steps between trajectory steps)
  --version             Print version information



#### EXAMPLES
Some examples are provided in the folder "examples". All of those files were used for testing the program during the development. In every file, the cell dimensions are given in the second line.

## Example 1
> mdanalyzer bp3n/bp3n-tra3.xyz 12.7 12.7 12.7 > bp3n/bp3n-tra3.log
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

# The atomic IDs in parentesys corresponds with the order of the atoms in the trajectory file, starting from 1.


## Example 2
> mdanalyzer bp3n/bp3n-ion3.xyz 12.7 12.7 12.7 -o w > bp3n/bp3n-ion3.log
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
# This trajectory file has the same previous reaction. It is recomended to save the outputs for and then look for specific symbols and IDs.


## Example 3
> mdanalyzer bp4l/bp4l-tra.xyz 12.7 12.7 12.7 -e C > bp4l/bp4l-tra.log
trajectory file: bp4l/bp4l-tra.xyz
atoms: 157
steps: 1900
step 1205: bond (62)C-O(16) is formed.
step 1471: bond (62)C-S(1) is broken.
NET REACTION:
CS  ->  CO
# There are several reactions in this trajectory file, but only the reactions involving C atoms are shown. In this case, the reaction "H2O(16) + Ar(62)-S(1)O3 -> Ar(62)-H2O(16)+ + S(1)O3-" is observed.


## Example 4
bp3n-tra/bp3n-tra.xyz
    > mdanalyzer bp3n/bp3n-tra.xyz 12.7 12.7 12.7 > bp3n/full.log
    > mdanalyzer bp3n/bp3n-tra.xyz 12.7 12.7 12.7 -e C > bp3n/C.log
    > mdanalyzer bp3n/bp3n-tra.xyz 12.7 12.7 12.7 -e C-O > bp3n/C-O.log
# When looking for specific bonds or elements, the options "-e El" and "-e E1-E2" can be used to reduce the time. In the previous example, the time elapsed was 7.6 min for analyzing every bond, 1.3 min for bonds involving C, and 26 sec for C-O bonds only.
    > grep C bp3n/full.log
    step 1236: bond (11)O-C(60) is formed.
    step 3436: bond (41)O-C(62) is formed.
    step 15422: bond (60)C-C(62) is broken.
    > grep C bp3n/C.log
    step 1236: bond (60)C-O(11) is formed.
    step 3436: bond (62)C-O(41) is formed.
    step 15422: bond (60)C-C(62) is broken.
    > grep C bp3n/C-O.log
    step 1236: bond (60)C-O(11) is formed.
    step 3436: bond (62)C-O(41) is formed.


#### Finally
There are more options available, but this is the basic and intermedium usage. Don't hesitate to write me an email if you have further questions.

-Luis Alvarez F.-
alvarez@theochem.uni-hannover.de

