# Vibrational Analysis Using CNDO/2
Sam Wollenburg

Date: May 1, 2024

### Subdirectories

- basis/ : Contains the basis functions for each atom. Must be STO3G and follow the file naming convention, {Atomic symbol}_STO3G.txt.
- output/ : Has the same structure as the directory used for input. The output file name has the structure of {molecule name}.out
- input/ : Contains sample input files to be run with each example. The input file name has the structure of {molecule name}.txt. It contains information on the number of atoms and the charge of the molecule. Following, each line is represented by an atomic number, and then x, y, and z coordinates. Coordinate values must be in a.u.

### Files

- vibrational_analysis.cpp : Takes input molecule file and produces the overlap matrix, core hamiltonian matrix, CNDO/2 Fock Matrix, MO coefficients (C matrix) and the molecule's energies (total, nuclear repulsion, and electron). This also computes the overlap gradient matrix, the gamma gradient matrix, and the molecule's energy gradients (total, nuclear repulsion, and electron). Following, optimization of the molecule is performed using steepest gradient descent with backtracking (Armijo condition) to find the ideal bond length. Finally, vibrational analysis is performed and the frequencies and vibrational modes are reported. Vibrational analysis is performed using central difference theory. To compile this file, use make target `vibrational_analysis` as seen below. To input files directily use `./vibrational_analysis ./input/{molecule name}.txt` after compilation. All output files will go to the `./output/` subdirectory. The file names will be in the format `{molecule name}.out`. Input molecules can contain H, C, N, O and F atoms. The number of valence electrons will be split evenly between p and q with q having the extra electron if the total number of valence electrons are odd. The values of p and q can be set when calling the file. To do this use `./vibrational_analysis ./{input}/{molecule name}.txt {p} {q}` where p and q are integer values that sum to the total number of valence electrons in the system. 
- analysis.ipynb : This notebook was used to calculate the average bond distance for the output molecules. This notebook is included for completness but was primarily used for scratch work and calculations.
- Makefile : Contains targets to easily compile the code. Use target `vibrational_analysis` to compile. To run the executable and all given example files use target `run_all`. To remove the executable, use target `clean`.
- README.md : This file.

### Analysis
Vibrational frequency predictions were more accurate for bending modes, rather than for streching modes. This was seen across all molecules. For example, H2O bending mode was under predicted by less than 20%, while the streching modes were over predicted by more than 70% of the experimental value. This was also seen in CO2 where the bending modes were under predicted by ~22%, while the streching modes were over predicted by ~45%. NH3 shows similar results. H2 and O2 produced rotational modes that should not be possible as the CNDO/2 method is rotationally invariant. This is likely a result of the central difference approximation used to calculate the Hessian, or the optimizer being used during geometeric optimization. Further investigation is required.

The Hessian can be read as the column read by first recognizing that the matrix is constructed with the energy gradient with respect to x, y, and z for the first atom, then the next three numbers are the gradient with respect to x, y, z of the second atom, etc, for each column. These value in the first column is the gradient of the gradient with respect to the x of the first atom, the second column is the gradient of the gradient with respect to eh y of the first atom, etc for x, y, z, for each atom. This produces a symmetric matrix. 

After mass weighting and diagonalization, the normal modes (eigenvectors) can be read by grouping the first three numbers in the column as the gradient in the x, y, and z directions, respectively, for the first atom, the next set of three for the second atom, etc. in which the atoms will move during vibration. 

It should be noted that the CO2.txt molecule file was made by placing the carbon atom at (0, 0, 0) and then placing one atom in the -x direction at position (-2.1977513247, 0, 0) and one atom in the +x direction at position (2.1977513247, 0, 0). The NH3.txt file was made by using the first set of coordinates provided from: https://www.atomic-scale-physics.de/lattice/struk/NH3.html. The H2O.txt file was made be using the coordinates provided here https://github.com/susilehtola/erkale/blob/master/tests/xyz/H2O.xyz. As with the other molecules, the values of Angstrom were converted to a.u..