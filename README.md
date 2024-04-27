# Vibrational Analysis Using CNDO/2
Sam Wollenburg

Date: Apr 4, 2024

### Subdirectories

- basis/ : Contains the basis functions for each atom. Must be STO3G and follow the file naming convention, {Atomic symbol}_STO3G.txt.
- output/ : Has the same structure as the directory used for input. The output file name has the structure of {molecule name}.out
- sample_input/ : Contains sample input files to be run with each example. The input file name has the structure of {molecule name}.txt. It contains information on the number of atoms and the charge of the molecule. Following, each line is represented by an atomic number, and then x, y, and z coordinates. Coordinate values must be in a.u.**Contains additional example CO2 and NH3.**
- sample_output/ : Contains the expected output for the files in sample_input.

### Files

- scf_w_grad.cpp : Takes input molecule file and produces the overlap matrix, core hamiltonian matrix, CNDO/2 Fock Matrix, MO coefficients (C matrix) and the molecule's energies (total, nuclear repulsion, and electron). This also computes the overlap gradient matrix, the gamma gradient matrix, and the molecule's energy gradients (total, nuclear repulsion, and electron). Following, optimization of the molecule is performed using steepest gradient descent with backtracking (Armijo condition) to find the ideal bond length. To compile this file, use make target `scf_w_grad` as seen below. To input files directily use `./scf ./{input}/{molecule name}.txt` after compilation. The input folder must be at the same level as the `./sample_input/` subdirectory. All output files will go to the `./output/` subdirectory. The file names will be in the format `{molecule name}.out`. Input molecules can contain H, C, N, O and F atoms. The number of valence electrons will be split evenly between p and q with q having the extra electron if the total number of valence electrons are odd. The values of p and q can be set when calling the file. To do this use `./scf ./{input}/{molecule name}.txt {p} {q}` where p and q are integer values that sum to the total number of valence electrons in the system. 
- analysis.ipynb : This notebook was used to calculate the average bond distance for the output molecules. This notebook is included for completness but was primarily used for scratch work and calculations.
- Makefile : Contains targets to easily compile the code. Use target `scf_w_grad` to compile. To run the executable and all given example files use target `run_all`. To remove the executable, use target `clean`.
- README.md : This file.

### Analysis
Steepest descent with backtracking (Armijo condition) was used to optimize the location of the atoms in the molecule. These bond distances were then compared to Table 4.1b (Pople and Beveridge, 1970). The values provided in the book will be converted to atomic units (a.u.) to be consistent with the coordinates in the input files and output files.

| Molecule | Initial Bond Length (a.u.) | Optimized Bond Length (a.u.) | Reported CNDO Bond Length (a.u.)* | Obs. Bond Length (a.u.)* |
|---|---|---|---|---|
| H2 | 1.3984 | 1.377688 | 1.4097355875 | 1.4021766835 |
| HF | 1.72 | 1.811375 | 1.8897259886 | 1.7328787315 |
| HO | 1.83 | 1.887560 | 1.9388588643 | 1.8349239349 |
| CO2** | 2.1977513247 | 2.224459836 | 2.1958615987 | 2.1977513247 |
| NH3** | 1.8842912216117451 | 3.7184401648782806 | | 1.9275205084 |
|<td colspan=5> \*Values from Table 4.1b of Pople and Beveridge, 1970. Values converted from Angstrom to a.u. (1 Angstrom = 1.8897259886 a.u.). </td>|
|<td colspan=5> \*\*Bond lengths for molecules are average bond lengths. </td>|


As seen above, the HF and HO molecules were able to be optimized towards the reported CNDO bond length. On the other hand, the H2, CO2, and NH3 molecules performed worse and moved in the opposite direction than what was expected. I believe these results are caused by the steepest descent algorithm being used and the optimization landscape presented by the molecules tested. H2 and CO2 did not perform horribly, but NH3 had difficulty optimizing it's structure. This could be a result of the initial location of the atoms, or a result of the optimizer (or both!).

It should be noted that the CO2.txt molecule file was made by placing the carbon atom at (0, 0, 0) and then placing one atom in the -x direction at position (-2.1977513247, 0, 0) and one atom in the +x direction at position (2.1977513247, 0, 0). The NH3.txt file was made by using the first set of coordinates provided from: https://www.atomic-scale-physics.de/lattice/struk/NH3.html. As with the other molecules, the values of Angstrom were converted to a.u..