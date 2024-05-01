#include <armadillo>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <map>


const double EV_AU_CONV = 27.211;

std::map<int, std::string> periodicTable{{1,"H"}, {2,"He"}, {6,"C"}, {7,"N"}, {8, "O"}, {9, "F"}};
                                        // {3,"Li"}, {4,"Be"}, {5,"B"}, , {9, "F"}, {10,"Ne"}};

struct ContractedGaussian{
    int element, Z; // atomic num, num val e-
    double IA_2, beta; // 0.5*(Iu+Au), B
    std::vector<double> center; // (x, y, z)
    std::vector<int> L; // l, m, n
    std::vector<double> e; // e1, e2, e3
    std::vector<double> cc; // cc1, cc2, cc3
    std::vector<double> N;
};

struct Atom{
    int element, charge, L, Z;
    double x, y, z, ex, ey, ez, csx, csy, csz, cpx, cpy, cpz, mass;
    std::vector<std::vector<int>> orbitals; 
};

struct Molecule{
    std::vector<Atom> atoms;
    int numElements, charge;
    std::vector<ContractedGaussian> basisFunctions; // store
    int p, q, numE;
    double nuclearRepulsionEnergy, ElectronEnergy, TotalEnergy;
};

/*
 File Parsing
*/

void readBasisFile(Atom &atom, std::string basisSetType="STO3G") {
    /*
     Only works for STO3G currently
    */
    std::string filepath = "./basis/"+periodicTable[atom.element]+"_"+basisSetType+".txt";
    std::ifstream basisSet(filepath);
    std::string line;
    double num;

    if(basisSet.is_open()) {
        for(int i=0; i<3; ++i) {
            std::getline(basisSet, line);
            std::istringstream ss(line);
            if(atom.element==1) {
                if(i==0){
                    ss >> atom.ex >> atom.csx;
                } else if(i==1) {
                    ss >> atom.ey >> atom.csy;
                } else {
                    ss >> atom.ez >> atom.csz;
                }
                
            } else if(atom.element==6 || atom.element==7 || atom.element==8 || atom.element==9) {
                if(i==0){
                    ss >> atom.ex >> atom.csx >> atom.cpx;
                } else if(i==1) {
                    ss >> atom.ey >> atom.csy >> atom.cpy;
                } else {
                    ss >> atom.ez >> atom.csz >> atom.cpz;
                }
            } else {
                std::cout << "This atom is not supported." << std::endl;
            }
        }
    }
}

Molecule readMolecule(std::string filepath) {
    std::ifstream myfile(filepath);
    std::string line;
    Molecule molecule;
    int numElements, charge;
    
    if(myfile.is_open()) {
        std::getline(myfile, line);
        std::istringstream ss(line);
        ss >> numElements >> charge;
        molecule.numElements = numElements;
        molecule.charge = charge;

        for(int i=0; i<numElements; ++i) {
            Atom atom;
            std::getline(myfile, line);
            std::istringstream ss(line);
            ss >> atom.element >> atom.x >> atom.y >> atom.z;
            readBasisFile(atom);
            
            if(atom.element == 1) {
                atom.L = 0;
                atom.Z = 1;
                atom.mass = 1.0078;
            } else if(atom.element==6) {
                atom.L = 1;
                atom.Z = 4;
                atom.mass = 12.011;
            } else if(atom.element==7) {
                atom.L = 1;
                atom.Z = 5;
                atom.mass = 14.007;
            }else if(atom.element==8) {
                atom.L = 1;
                atom.Z = 6;
                atom.mass = 15.999;
            } else if(atom.element==9) {
                atom.L = 1;
                atom.Z = 7;
                atom.mass = 18.998;
            } else {
                std::cout << "Atom not supported." << std::endl;
            }
            
            molecule.atoms.push_back(atom);
        }
    }

    return molecule;
}

/*
 Helper function
*/
int numElectrons(Molecule molecule) {
    /*
     Returns number of valence electrons
    */
    int valE = 0;

    for(auto &atom:molecule.atoms) {
        valE += atom.Z;
    }
    // Adds charge
    valE -= molecule.charge;

    return valE;    
}

int numBasisFuncs(Molecule molecule) {
    int numH = 0;
    int numC = 0;
    int numO = 0;
    int numN = 0;
    int numF = 0;

    for(auto &atom:molecule.atoms) {
        if(atom.element==1) {
            ++numH;
        } else if(atom.element == 6) {
            ++numC;
        } else if(atom.element == 7) {
            ++numN;
        } else if(atom.element == 8) {
            ++numO;
        } else if(atom.element == 9) {
            ++numF;
        } else {
            std::cout << "Element number " << atom.element << " is not supported.";
        }
    }

    return 4*(numC + numN + numO + numF) + numH;
}

void printMol(Molecule molecule) {
    std::cout << "This molecule has " << molecule.atoms.size() << " atoms.";
    std::cout << "The atoms are:";
    for(int i=0; i<molecule.atoms.size(); ++i) {
        std::cout << " " << molecule.atoms[i].element;
    } 
    std::cout << std::endl;

    std::cout << "The atoms' angular momentums (L) are:";
    for(int i=0; i<molecule.atoms.size(); ++i) {
        std::cout << " " << molecule.atoms[i].L;
    } 
    std::cout << std::endl;

    std::cout << "This molecule has " << numBasisFuncs(molecule) << " basis functions." << std::endl;

    std::cout << "Each atoms' orbitals:" << std::endl;
    for(int i=0; i<molecule.atoms.size(); ++i) {
        for(auto &o:molecule.atoms[i].orbitals) {
            std::cout << "(";
            for(auto &oo:o) {
                std::cout << oo << " ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
        
    }

}

void printEXPandCC(Molecule molecule) {
    std::cout << "The Exponenets and Contraction Coefficients are: " << std::endl;
    std::cout << "(atom)" << std::endl;
    std::cout << "    (Exponent X)  (S Contr. Coeff. X)  (P Contr. Coeff. X)" << std::endl; 
    std::cout << "    (Exponent Y)  (S Contr. Coeff. Y)  (P Contr. Coeff. Y)" << std::endl; 
    std::cout << "    (Exponent Z)  (S Contr. Coeff. Z)  (P Contr. Coeff. Z)" << std::endl; 
    for(auto &atom:molecule.atoms) {
        std::cout << atom.element << std::endl;
        std::cout << "    " << atom.ex << "  " << atom.csx;
        if(atom.element > 4) {
            std::cout << "  " << atom.cpx;
        }
        std::cout << std::endl;

        std::cout << "    " << atom.ey << "  " << atom.csy;
        if(atom.element > 4) {
            std::cout << "  " << atom.cpy;
        }
        std::cout << std::endl;

        std::cout << "    " << atom.ez << "  " << atom.csz;
        if(atom.element > 4) {
            std::cout << "  " << atom.cpz;
        }
        std::cout << std::endl;
    }
    
}

void printBasisFuncs(Molecule molecule) {
    std::cout << "The basis functions for the molecule are: " << std::endl; 
    for(auto &basis:molecule.basisFunctions) {
        std::cout << "Element #: " << basis.element << std::endl;
        std::cout << "Center: (" << basis.center[0] << ", " << basis.center[1] << ", " << basis.center[2] << ")" << std::endl;
        std::cout << "L: (" << basis.L[0] << ", "  << basis.L[1] << ", " << basis.L[2] << ")" << std::endl;
        std::cout << "Exponents: " << basis.e[0] << ", " << basis.e[1] << ", " << basis.e[2] << std::endl;
        std::cout << "Cont. Coefs.: " << basis.cc[0] << ", " << basis.cc[1] << ", " << basis.cc[2] << std::endl;
        std::cout << "Norm. Coefs.: " << basis.N[0] << ", " << basis.N[1] << ", " << basis.N[2] << std::endl << std::endl;
    }
}

void assignShells(Molecule &molecule) {
    std::vector<std::vector<int>> s = {{0,0,0}};
    std::vector<std::vector<int>> sp = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};

    for(auto &atom:molecule.atoms) {
        if(atom.L==0) {
            atom.orbitals = s;
        } else if(atom.L==1) {
            atom.orbitals = sp;
        }
    } 
}

void createBasisFuncs(Molecule &molecule) {
    for(auto &atom:molecule.atoms) {
        int count = 0;
        for(auto &orbit:atom.orbitals) {
            ContractedGaussian cg;
            
            cg.element = atom.element;
            cg.center = {atom.x, atom.y, atom.z};
            cg.e = {atom.ex, atom.ey, atom.ez};
            cg.L = orbit;
            cg.Z = atom.Z;
            if(count == 0) {
                cg.cc = {atom.csx, atom.csy, atom.csz};
            } else {
                cg.cc = {atom.cpx, atom.cpy, atom.cpz};
            }
            
            molecule.basisFunctions.push_back(cg);
            ++count;
        }
    }
}

arma::Mat<double> createMatrix(std::vector<std::vector<double>> vecvec) {
    arma::Mat<double> mat(vecvec.size(), vecvec[0].size());
    int i = 0;
    for(auto &vec:vecvec) {
        int j = 0;
        for(auto &e:vec) {
            mat(i, j) = e;
            ++j;
        }
        ++i;
    }
    return mat;
}

/*
 Analytical Soln of Overlap
*/

int factorial(int x) {
    int value = 1;
    for(int i=1; i<=x; ++i) {
        value*=i;
    }
    return value;
}

int doubleFactorial(int x) {
    int value = 1;
    if(x%2==0) { // even
        for(auto i=2; i<=x; i+=2) {
            value*=i;
        }
    } else { // odd
        for(auto i=1; i<=x; i+=2) {
            value*=i;
        }
    }
    return value;
}

std::vector<double> weightedCenter(std::vector<double> center1, std::vector<double> center2, double e1, double e2) {
    // seperable
    std::vector<double> p;
    p.push_back((e1*center1[0] + e2*center2[0]) / (e1 + e2));
    p.push_back((e1*center1[1] + e2*center2[1]) / (e1 + e2));
    p.push_back((e1*center1[2] + e2*center2[2]) / (e1 + e2));
    return p;
}

double binomial(int m, int n) {
    return factorial(m)/(factorial(n)*factorial(m-n));
}

double analyticalEqPart(double alpha, double beta, double A, double B, int LA, int LB, double P) {
    double exp_part = std::exp(-((alpha*beta*std::pow((A-B), 2))/(alpha+beta)));
    double sqrt_part = std::sqrt(M_PI/(alpha + beta));
    double sum = 0;

    for(auto i=0; i<=LA; ++i) {
        for(auto j=0; j<=LB; ++j) {
            // only even (i+j) contributes non-0 result
            if((i+j)%2==0) {
                double numerator = doubleFactorial(i+j-1)*std::pow((P-A),(LA-i))*std::pow((P-B),(LB-j));
                double denominator = std::pow((2*(alpha+beta)), ((i+j)/2));
                sum += binomial(LA, i) * binomial(LB, j) * numerator/denominator;
            }
        }
    }
    return exp_part*sqrt_part*sum;
}

double analyticalEq(ContractedGaussian cg1,  ContractedGaussian cg2, int i, int j) {
    double product = 1;

    // weighted center
    Atom a1;
    // a1.x = 
    std::vector<double> P = weightedCenter(cg1.center, cg2.center, cg1.e[i], cg2.e[j]);

    for(auto idx=0; idx<3; ++idx) {
        if(idx==0) {
            product *= analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[0], cg2.center[0], cg1.L[0], cg2.L[0], P[0]);
        } else if (idx==1) {
            product *= analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[1], cg2.center[1], cg1.L[1], cg2.L[1], P[1]);
        } else {
            product *= analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[2], cg2.center[2], cg1.L[2], cg2.L[2], P[2]);
        }
    }
    return product;
}

/*
 Normalize and overlap
*/

void normalizeBasisFuncs(Molecule &molecule) {
    for(auto &basis:molecule.basisFunctions) {
        for(auto i=0; i<3; ++i) {
            // overlap
            double overlap, N;
            overlap = analyticalEq(basis, basis, i, i);
            
            // find and save N
            N =  1/std::sqrt(overlap);
            basis.N.push_back(1/std::sqrt(overlap)); 
        }
    }   
}

double calculateOverlap(ContractedGaussian cg1, ContractedGaussian cg2) {
    double sum=0;
    for(auto k=0; k<3; ++k) {
        for(auto l=0; l<3; ++l) {
            sum += cg1.cc[k]*cg2.cc[l]*cg1.N[k]*cg2.N[l]*analyticalEq(cg1, cg2, k, l);
        }
    }
    return sum;
}

std::vector<std::vector<double>> createOverlapMatrix(Molecule molecule) {
    std::vector<std::vector<double>> overlapMatrix;

    for(auto &basisY:molecule.basisFunctions) {
        std::vector<double> temp;
        for(auto &basisX:molecule.basisFunctions) {
            temp.push_back(calculateOverlap(basisY, basisX));
        }
        overlapMatrix.push_back(temp);
    }

    return overlapMatrix;
}

/*
 SCF
*/

void assignNumElectrons(Molecule &molecule) {
    molecule.numE = numElectrons(molecule);
    molecule.q = static_cast<int>(std::ceil(molecule.numE/2));
    molecule.p = molecule.numE - molecule.q;
}

void assignSemiEmpiricalParams(Molecule &molecule) {
    for(int i=0; i<molecule.basisFunctions.size(); ++i) {
        if(molecule.basisFunctions[i].element == 1) {
            molecule.basisFunctions[i].IA_2 = 7.176; // H(1s)
            molecule.basisFunctions[i].beta = -9;
        } else if(molecule.basisFunctions[i].element == 6) {
            molecule.basisFunctions[i].beta = -21;
            if(molecule.basisFunctions[i].L[0]+molecule.basisFunctions[i].L[1]+molecule.basisFunctions[i].L[2] == 0) {
                molecule.basisFunctions[i].IA_2 = 14.051; // C(2s)
            } else {
                molecule.basisFunctions[i].IA_2 = 5.572; // C(2p)
            }
        } else if(molecule.basisFunctions[i].element == 7) {
            molecule.basisFunctions[i].beta = -25;
            if(molecule.basisFunctions[i].L[0]+molecule.basisFunctions[i].L[1]+molecule.basisFunctions[i].L[2] == 0) {
                molecule.basisFunctions[i].IA_2 = 19.316; // N(2s)
            } else {
                molecule.basisFunctions[i].IA_2 = 7.275; // N(2p)
            }
        } else if(molecule.basisFunctions[i].element == 8) {
            molecule.basisFunctions[i].beta = -31;
            if(molecule.basisFunctions[i].L[0]+molecule.basisFunctions[i].L[1]+molecule.basisFunctions[i].L[2] == 0) {
                molecule.basisFunctions[i].IA_2 = 25.390; // O(2s)
            } else {
                molecule.basisFunctions[i].IA_2 = 9.111; // O(2p)
            }
        }  else if(molecule.basisFunctions[i].element == 9) {
            molecule.basisFunctions[i].beta = -39;
            if(molecule.basisFunctions[i].L[0]+molecule.basisFunctions[i].L[1]+molecule.basisFunctions[i].L[2] == 0) {
                molecule.basisFunctions[i].IA_2 = 32.272; // O(2s)
            } else {
                molecule.basisFunctions[i].IA_2 = 11.080; // O(2p)
            }
        } 
    }
}

double gammaDiag(ContractedGaussian cg1) {
    double sum = 0;
    for(int k=0; k<3; ++k) {
        for(int kk=0; kk<3; ++kk) {
            for(int l=0; l<3; ++l) {
                for(int ll=0; ll<3; ++ll) {
                    double sigmaA = 1/(cg1.e[k] + cg1.e[kk]);
                    double UA = std::pow((M_PI * sigmaA), 1.5);
                    double sigmaB = 1/(cg1.e[l] + cg1.e[ll]);
                    double UB = std::pow((M_PI * sigmaB), 1.5);
                    double V2 = 1/(sigmaA + sigmaB);
                    
                    double zerozero = UA*UB*std::sqrt(2*V2)*std::sqrt(2/M_PI);
                    sum += (cg1.cc[k]*cg1.N[k])*(cg1.cc[kk]*cg1.N[kk])*(cg1.cc[l]*cg1.N[l])*(cg1.cc[ll]*cg1.N[ll])*zerozero;
                }
            }
        }
    }
    return sum;
}

double gammaOffDiag(ContractedGaussian cg1, ContractedGaussian cg2) {
    double sum = 0;
    for(int k=0; k<3; ++k) {
        for(int kk=0; kk<3; ++kk) {
            for(int l=0; l<3; ++l) {
                for(int ll=0; ll<3; ++ll) {
                    double sigmaA = 1/(cg1.e[k] + cg1.e[kk]);
                    double UA = std::pow((M_PI * sigmaA), 1.5);
                    double sigmaB = 1/(cg2.e[l] + cg2.e[ll]);
                    double UB = std::pow((M_PI * sigmaB), 1.5);
                    double V2 = 1/(sigmaA + sigmaB);
                    double dist = std::sqrt(std::pow((cg1.center[0]-cg2.center[0]), 2)+
                                            std::pow((cg1.center[1]-cg2.center[1]), 2)+
                                            std::pow((cg1.center[2]-cg2.center[2]), 2));
                    double T = V2*std::pow(dist,2);

                    double zerozero = UA*UB*std::sqrt(1/std::pow(dist, 2))*std::erf(std::sqrt(T));
                    sum += (cg1.cc[k]*cg1.N[k])*(cg1.cc[kk]*cg1.N[kk])*(cg2.cc[l]*cg2.N[l])*(cg2.cc[ll]*cg2.N[ll]*zerozero);
                }
            }
        }
    }
    return sum;
}

std::vector<std::vector<double>> createGammaMatrix(Molecule molecule) {
    std::vector<std::vector<double>> gamma;
    int count1, count2;

    count1=0;
    for(auto &cg1:molecule.basisFunctions) {
        if(cg1.L!=std::vector{0,0,0}){continue;}
        count2=0;
        std::vector<double> temp;
        for(auto &cg2:molecule.basisFunctions) {
            if(cg2.L!=std::vector{0,0,0}){continue;}
            if(count1==count2) {
                temp.push_back(gammaDiag(cg1));
            } else {
                temp.push_back(gammaOffDiag(cg1, cg2));
            }
            ++count2;
        }
        gamma.push_back(temp);
        ++count1;
    }
    return gamma;
}

std::vector<std::vector<double>> createCoreHamiltonianMatrix(Molecule molecule, arma::Mat<double> gammaMat, arma::Mat<double> overlap) {
    std::vector<std::vector<double>> coreHam;
    int gcount1;

    gcount1 = -1;
    for(auto i=0; i<molecule.basisFunctions.size(); ++i) {
        if(molecule.basisFunctions[i].L == std::vector{0,0,0}){
            ++gcount1;
        }

        std::vector<double> temp;
        for(auto j=0; j<molecule.basisFunctions.size(); ++j) {
            if(i==j) {
                double sumZy = 0;
                int basis_idx = 0;
                for(auto k=0; k<gammaMat.n_cols; ++k) {
                    if(k != gcount1) {
                        sumZy += molecule.basisFunctions[basis_idx].Z*gammaMat(gcount1, k);
                    }
                    
                    if(molecule.basisFunctions[basis_idx].element==1) {
                        basis_idx += 1;
                    } else {
                        basis_idx += 4;  
                    }
                }
                temp.push_back(-molecule.basisFunctions[i].IA_2 - 
                               (molecule.basisFunctions[i].Z - 0.5) * gammaMat(gcount1, gcount1) - 
                               sumZy);
            } else {
                temp.push_back(0.5*(molecule.basisFunctions[i].beta+molecule.basisFunctions[j].beta)*overlap(i, j));
            }
        }

        coreHam.push_back(temp);
    }

    return coreHam;
}

arma::Mat<double> createDensityMat(arma::Mat<double> cMat, int pORq) {
    arma::Mat<double> pMat(cMat.n_rows, cMat.n_cols);
    
    for(int u=0; u<pMat.n_rows; ++u) {
        for(int v=0; v<pMat.n_cols; ++v) {
            pMat(u,v) = 0;
            for(int i=0; i<pORq; ++i) {
                pMat(u,v) += cMat(u,i)*cMat(v, i);
            }
        }
    }

    return pMat;
}

arma::Mat<double> createCNDO2FockMat(Molecule molecule, arma::Mat<double> ptotMat, arma::Mat<double> pMat, arma::Mat<double> gammaMat, arma::Mat<double> overlapMat) {
    arma::Mat<double> fockMat(pMat.n_rows, pMat.n_cols);

    int A = 0;
    for(int u=0; u<fockMat.n_rows; ++u) {
        int B = 0;
        for(int v=0; v<fockMat.n_cols; ++v) {

            if(u==v) { 
                double sumPBBZB = 0;
                for(int tempB=0; tempB<molecule.numElements; ++tempB) {
                    if(tempB!=A) {
                        sumPBBZB += (ptotMat(tempB, 0) - molecule.atoms[tempB].Z) * gammaMat(A, tempB);
                    }
                }

                fockMat(u,u) = -molecule.basisFunctions[u].IA_2+
                               ((ptotMat(A, 0) - molecule.basisFunctions[u].Z) - (pMat(u,u) - 0.5))*gammaMat(A, B)+
                               sumPBBZB;
                               
            } else {
                fockMat(u,v) = 0.5*(molecule.basisFunctions[u].beta + molecule.basisFunctions[v].beta)*overlapMat(u,v)-
                               pMat(u,v)*gammaMat(A, B);
            }

            // update B (v+1 might cause issues ***)
            if(molecule.basisFunctions[v+1].L == std::vector{0,0,0}) {
                ++B;
            }

        }
        // update A (v+1 might cause issues ***)
        if(molecule.basisFunctions[u+1].L == std::vector{0,0,0}) {
            ++A;
        }
    }

    return fockMat;
}

void scf(std::ofstream &out, Molecule &molecule, arma::Mat<double> coreHamiltonianMat, arma::Mat<double> gammaMat, arma::Mat<double> overlapMat, double cutoff=1e-6) {
    arma::Mat<double> currentMatAlpha = coreHamiltonianMat;
    arma::Mat<double> currentMatBeta = coreHamiltonianMat;
    arma::Mat<double> palphaOld = coreHamiltonianMat;
    arma::Mat<double> pbetaOld = coreHamiltonianMat;
    arma::Mat<double> palpha = coreHamiltonianMat;
    arma::Mat<double> pbeta = coreHamiltonianMat;
    arma::vec Ea, Eb;
    arma::mat Ca, Cb;

    double pDiff = 1;
    int count=0;
    while(pDiff>cutoff) {
        out << "Iteration: " << count << std::endl;
        
        // Print Fock Matrices 
        arma::Mat<double> Fa(currentMatAlpha);
        Fa.print(out, "Fa");
        arma::Mat<double> Fb(currentMatBeta);
        Fb.print(out, "Fb");

        out << "after solving eigen equation: " << count << std::endl;

        // Solve Eigenvectors and Eigenvalues
        arma::eig_sym(Ea, Ca, Fa);
        Ca.print(out, "Ca");
        arma::eig_sym(Eb, Cb, Fb);
        Cb.print(out, "Cb");

        out << "  p = " << molecule.p << "  q = " << molecule.q << std::endl;

        // Create Density Matrices, p
        palpha = createDensityMat(Ca, molecule.p);
        palpha.print(out, "Pa_new");
        pbeta = createDensityMat(Cb, molecule.q);
        pbeta.print(out, "Pb_new");

        // create ptot matrix
        arma::Mat<double> ptot(molecule.numElements, 1);
        int basis_idx = 0;
        int elementPos = 0;
        for(auto &bf:molecule.basisFunctions) {
            int numBasisElement = 0;  // number of basis funcs for element
            if(bf.L == std::vector{0,0,0}) {
                if(bf.element==1) {
                    numBasisElement = 1;
                } else {
                    numBasisElement = 4;  
                }

                double sum=0;
                for(auto k=0; k<numBasisElement; ++k) {
                    sum += palpha(basis_idx+k, basis_idx+k) + pbeta(basis_idx+k, basis_idx+k);
                }

                ptot(elementPos, 0) = sum;
                
                if(bf.element==1) {
                    basis_idx += 1;
                } else {
                    basis_idx += 4;  
                }
                ++elementPos;  
            } 
        }

        ptot.print(out, "P_t");

        // create CNDO/2 Fock Matrices
        arma::Mat<double> fockMatAlpha(palpha.n_rows, palpha.n_cols);
        currentMatAlpha = createCNDO2FockMat(molecule, ptot, palpha, gammaMat, overlapMat);

        arma::Mat<double> fockMatBeta(palpha.n_rows, palpha.n_cols);
        currentMatBeta = createCNDO2FockMat(molecule, ptot, pbeta, gammaMat, overlapMat);
        
        // Calculate difference between old and new density matrices, palpha and pbeta
        if(count!=0) {
            arma::Mat<double> palphaDiffMat(palphaOld);
            arma::Mat<double> pbetaDiffMat(pbetaOld);
            for(int i=0; i<palphaDiffMat.n_rows; ++i) {
                for(int j=0; j<palphaDiffMat.n_cols; ++j) {
                    palphaDiffMat(i,j) = std::abs(palphaOld(i,j)-palpha(i,j));
                    pbetaDiffMat(i,j) = std::abs(pbetaOld(i,j)-pbeta(i,j));
                }
            }
            
            // Find the maximum difference between the matrices to test for convergence
            pDiff = palphaDiffMat.max();
            if(pDiff<pbetaDiffMat.max()) {
                pDiff = pbetaDiffMat.max();
            }
        }
        palphaOld = palpha;
        pbetaOld = pbeta;
        
        ++count;
    }
    // Print Eigenvalues
    Ea.print(out, "Ea");
    Eb.print(out, "Eb");

    // Print Eigenvectors
    Ca.print(out, "Ca");
    Cb.print(out, "Cb");
    
    // Calculate Nuclear Repulsion Energy
    molecule.nuclearRepulsionEnergy = 0;
    for(int A=0; A<molecule.numElements; ++A) {
        for(int B=0; B<A; ++B) {
            double dist = std::sqrt(std::pow((molecule.atoms[A].x - molecule.atoms[B].x), 2)+
                                    std::pow((molecule.atoms[A].y - molecule.atoms[B].y), 2)+
                                    std::pow((molecule.atoms[A].z - molecule.atoms[B].z), 2));
                          
            molecule.nuclearRepulsionEnergy+=molecule.atoms[A].Z*molecule.atoms[B].Z/dist;
        }

    }
    molecule.nuclearRepulsionEnergy *= EV_AU_CONV;
    
    // Calculate Total Electron Energy
    double sumPHFalpha = 0;
    double sumPHFbeta = 0;
    for(int u=0; u<coreHamiltonianMat.n_cols; ++u) {
        for(int v=0; v<coreHamiltonianMat.n_rows; ++v) {
            sumPHFalpha += palpha(u,v)*(coreHamiltonianMat(u,v)+currentMatAlpha(u,v));
            sumPHFbeta += pbeta(u,v)*(coreHamiltonianMat(u,v)+currentMatBeta(u,v));
        }
        
    }
    molecule.ElectronEnergy = 0.5*sumPHFalpha + 0.5*sumPHFbeta;

    // Calculate Total Energy
    molecule.TotalEnergy = molecule.ElectronEnergy + molecule.nuclearRepulsionEnergy;
}

/*
 Gradient
*/

double analyticalEqXYZ(ContractedGaussian cg1,  ContractedGaussian cg2, int i, int j, int xyz) {
    // weighted center
    std::vector<double> P = weightedCenter(cg1.center, cg2.center, cg1.e[i], cg2.e[j]);

    double product = analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[xyz], cg2.center[xyz], cg1.L[xyz], cg2.L[xyz], P[xyz]);
    
    return product;
}

double analyticalEqGrad(ContractedGaussian cg1,  ContractedGaussian cg2, int i, int j, int xyz) {
    // weighted center
    std::vector<double> P = weightedCenter(cg1.center, cg2.center, cg1.e[i], cg2.e[j]);

    double product = -cg1.L[xyz]*(analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[xyz], cg2.center[xyz], cg1.L[xyz]-1, cg2.L[xyz], P[xyz]));
    product += 2*cg1.e[i]*(analyticalEqPart(cg1.e[i], cg2.e[j], cg1.center[xyz], cg2.center[xyz], cg1.L[xyz]+1, cg2.L[xyz], P[xyz]));
    
    return product;
}

double calculateOverlapGrad(ContractedGaussian cg1, ContractedGaussian cg2, int xyz) {
    
    // Determines remaining x, y, or z dimensions to find the real overlap instead of the gradient
    int xyz1, xyz2;
    if(xyz==0) {
        xyz1=1;
        xyz2=2;
    } else if(xyz==1) {
        xyz1=0;
        xyz2=2;
    } else {
        xyz1=0;
        xyz2=1;
    }
    
    double sum=0;
    for(auto k=0; k<3; ++k) {
        for(auto l=0; l<3; ++l) {
            sum += cg1.cc[k]*cg1.N[k]*cg2.cc[l]*cg2.N[l]*analyticalEqGrad(cg1, cg2, k, l, xyz)*\
                                                         analyticalEqXYZ(cg1, cg2, k, l, xyz1)*\
                                                         analyticalEqXYZ(cg1, cg2, k, l, xyz2);
        }
    }
    return sum;
}

arma::Mat<double> createOverlapDerivMat(Molecule molecule) {
    int numBasis = numBasisFuncs(molecule);
    arma::Mat<double> overlapDeriv(3, numBasis*numBasis);

    for(auto xyz=0; xyz<3; ++xyz) {

        int elementCount_i = -1; // Tracks element id. Immediately gets set to zero upon entering loop as all elements begin with basisFunction with L = (0, 0, 0).
        for(auto i=0; i<molecule.basisFunctions.size(); ++i) {
            if(molecule.basisFunctions[i].L == std::vector{0,0,0}) {
                ++elementCount_i;
            }

            int elementCount_j = -1; // Tracks element id
            for(auto j=0; j<molecule.basisFunctions.size(); ++j) {
                if(molecule.basisFunctions[j].L == std::vector{0,0,0}) {
                    ++elementCount_j;
                }
                
                if(elementCount_i == elementCount_j) {
                    overlapDeriv(xyz, numBasis*i+j) = 0;
                } else {
                    overlapDeriv(xyz, numBasis*i+j) = calculateOverlapGrad(molecule.basisFunctions[i], molecule.basisFunctions[j], xyz);
                }
            }
        }
    }
    
    return overlapDeriv;
}

double gammaOffDiagDeriv(ContractedGaussian cg1, ContractedGaussian cg2, int xyz) {
    double distXYZ = cg1.center[xyz] - cg2.center[xyz];
    double dist = std::sqrt(std::pow((cg1.center[0]-cg2.center[0]), 2)+ \
                            std::pow((cg1.center[1]-cg2.center[1]), 2)+ \
                            std::pow((cg1.center[2]-cg2.center[2]), 2));
    double sum = 0;
    for(int k=0; k<3; ++k) {
        for(int kk=0; kk<3; ++kk) {
            for(int l=0; l<3; ++l) {
                for(int ll=0; ll<3; ++ll) {
                    double sigmaA = 1/(cg1.e[k] + cg1.e[kk]);
                    double UA = std::pow((M_PI * sigmaA), 1.5);
                    double sigmaB = 1/(cg2.e[l] + cg2.e[ll]);
                    double UB = std::pow((M_PI * sigmaB), 1.5);
                    double V2 = 1/(sigmaA + sigmaB);

                    double T = V2*std::pow(dist,2);

                    double zerozero = (UA*UB*distXYZ/std::pow(dist,2)) * \
                                      ((-std::erf(std::sqrt(T))/dist) + \
                                      2*std::sqrt(V2)*std::exp(-T)/std::sqrt(M_PI));
                    sum += (cg1.cc[k]*cg1.N[k])*(cg1.cc[kk]*cg1.N[kk])*(cg2.cc[l]*cg2.N[l])*(cg2.cc[ll]*cg2.N[ll]*zerozero);
                }
            }
        }
    }
    return sum;
}

arma::Mat<double> createGammaMatrixGrad(Molecule molecule) {
    int numAtoms = molecule.atoms.size();
    arma::Mat<double> gammaGrad(3, numAtoms*numAtoms);

    int count1, count2;
    for(int xyz=0; xyz<3; ++xyz) {
        count1=0;   // keeps track of element id
        for(int i=0; i<molecule.basisFunctions.size(); ++i) {
            if(molecule.basisFunctions[i].L!=std::vector{0,0,0}){continue;} // do I need these? ***
            count2=0;   // keeps track of element id
            for(int j=0; j<molecule.basisFunctions.size(); ++j) {
                if(molecule.basisFunctions[j].L!=std::vector{0,0,0}){continue;} // do I need these? ***
                if(count1==count2) {
                    gammaGrad(xyz, numAtoms*count1+count2) = 0;
                } else {
                    gammaGrad(xyz, numAtoms*count1+count2) = gammaOffDiagDeriv(molecule.basisFunctions[i], molecule.basisFunctions[j], xyz);
                }
                ++count2;
            }
            ++count1;
        }
    }

    return gammaGrad;
}

arma::Mat<double> createNuclearRepulsionGradMat(Molecule molecule) {
    int numAtoms = molecule.atoms.size();
    arma::Mat<double> VGrad(3, numAtoms);
    double distXYZ, dist, i, j;

    for(int xyz=0; xyz<3; ++xyz) {
        
        for(i=0; i<numAtoms; ++i) {
            VGrad(xyz, i) = 0;
            for(j=0; j<numAtoms; ++j) {
                if(i==j){continue;}
                double dist = std::sqrt(std::pow((molecule.atoms[i].x-molecule.atoms[j].x), 2)+
                                        std::pow((molecule.atoms[i].y-molecule.atoms[j].y), 2)+
                                        std::pow((molecule.atoms[i].z-molecule.atoms[j].z), 2));

                // select dimension
                if(xyz==0) {
                    distXYZ = molecule.atoms[i].x - molecule.atoms[j].x;
                } else if(xyz==1) {
                    distXYZ = molecule.atoms[i].y - molecule.atoms[j].y;
                } else if(xyz==2) {
                    distXYZ = molecule.atoms[i].z - molecule.atoms[j].z;
                }
                VGrad(xyz, i) += molecule.atoms[i].Z*molecule.atoms[j].Z*distXYZ/std::pow(dist, 3);
            }
            VGrad(xyz, i) *= -1;
        }
    }
    VGrad *= EV_AU_CONV;
    
    return VGrad;
}

arma::Mat<double> createXMat(Molecule molecule, arma::Mat<double> ptotal) {
    arma::Mat<double> x(ptotal.n_rows, ptotal.n_cols);

    for(int i=0; i<molecule.basisFunctions.size(); ++i) {
        for(int j=0; j<molecule.basisFunctions.size(); ++j) {
            x(i,j) = (molecule.basisFunctions[i].beta + molecule.basisFunctions[j].beta)*ptotal(i,j);
        }
    }

    return x;
}

arma::Mat<double> createYMat(Molecule molecule, arma::Mat<double> ptot, arma::Mat<double> ptotal, arma::Mat<double> palpha, arma::Mat<double> pbeta) {
    int numBasisF = molecule.basisFunctions.size();
    int numAtoms = molecule.atoms.size();
    arma::Mat<double> y(numAtoms, numAtoms);

    
    int basis_i, basis_j; // idx for basis functions
    basis_i = 0;
    for(int i=0; i<numAtoms; ++i) {
        basis_j = 0;
        for(int j=0; j<numAtoms; ++j) {
            double sumsum = 0;
            for(int u=basis_i; u<basis_i+molecule.atoms[i].orbitals.size(); ++u) {
                for(int v=basis_j; v<basis_j+molecule.atoms[j].orbitals.size(); ++v) {
                    // std::cout << molecule.atoms[i].orbitals.size() << " " << molecule.atoms[j].orbitals.size() << "    " << u << " " << v << std::endl;
                    sumsum += palpha(u,v)*palpha(u,v) + pbeta(u,v)*pbeta(u,v);
                }
            }
            // std::cout << std::endl;

            y(i, j) = ptot(i,0)*ptot(j,0) - molecule.atoms[j].Z*ptot(i,0) - molecule.atoms[i].Z*ptot(j,0) - sumsum;
            basis_j += molecule.atoms[j].orbitals.size();
        }
        basis_i += molecule.atoms[i].orbitals.size();
    }  

    return y;
} 

arma::Mat<double> createElectronGradMat(Molecule molecule, arma::Mat<double> x, arma::Mat<double> overlapGrad, arma::Mat<double> y, arma::Mat<double> gammaGrad) {
    int numBasisF = molecule.basisFunctions.size();
    int numAtoms = molecule.atoms.size();
    arma::Mat<double> electronGrad(3, numAtoms);

    // i and j are the element IDs, while basis_i and basis_j are the corresponding basisFunction IDs
    for(int xyz=0; xyz<3; ++xyz) {
        int basis_i = 0; 
        for(int i=0; i<numAtoms; ++i) { // loop through elements
            int basis_j = 0; 
            double sum = 0;
            double sumsum = 0;
            for(int j=0; j<numAtoms; ++j) { // loop through elements
                if(i!=j){ // ignores when u is of A, but B is A (B == A when i == j)
                    
                    for(int u=basis_i; u<basis_i+molecule.atoms[i].orbitals.size(); ++u) { // loop through u of A
                        for(int v=basis_j; v<basis_j+molecule.atoms[j].orbitals.size(); ++v) { // loop through v of B != A
                            sumsum += x(u,v)*overlapGrad(xyz, numBasisF*u+v);
                        }
                    }

                    sum += y(i,j)*gammaGrad(xyz, numAtoms*i+j);
                
                } 
                
                basis_j += molecule.atoms[j].orbitals.size();
            }
            basis_i += molecule.atoms[i].orbitals.size();
            electronGrad(xyz, i) = sumsum + sum;
        }
    }

    return electronGrad;
}

void scfWithGrad(std::ofstream &out, Molecule &molecule, arma::Mat<double> coreHamiltonianMat, \
                 arma::Mat<double> gammaMat, arma::Mat<double> overlapMat, 
                 arma::Mat<double> gammaMatDeriv, arma::Mat<double> overlapMatDeriv, double cutoff=1e-6) {
    arma::Mat<double> currentMatAlpha = coreHamiltonianMat;
    arma::Mat<double> currentMatBeta = coreHamiltonianMat;
    arma::Mat<double> palphaOld = coreHamiltonianMat;
    arma::Mat<double> pbetaOld = coreHamiltonianMat;
    arma::Mat<double> palpha = coreHamiltonianMat;
    arma::Mat<double> pbeta = coreHamiltonianMat;
    arma::Mat<double> ptot(molecule.numElements, 1);
    arma::Mat<double> ptotal; // not to be confused with ptot; ptotal = palpha + pbeta
    arma::vec Ea, Eb;
    arma::mat Ca, Cb;

    double pDiff = 1;
    int count = 0;
    while(pDiff>cutoff) {
        // Fock Matrices 
        arma::Mat<double> Fa(currentMatAlpha);
        arma::Mat<double> Fb(currentMatBeta);

        // Solve Eigenvectors and Eigenvalues
        arma::eig_sym(Ea, Ca, Fa);
        arma::eig_sym(Eb, Cb, Fb);

        // Create Density Matrices, p
        palpha = createDensityMat(Ca, molecule.p);
        pbeta = createDensityMat(Cb, molecule.q);

        // create ptot matrix
        int basis_idx = 0;
        int elementPos = 0;
        for(auto &bf:molecule.basisFunctions) {
            int numBasisElement = 0;  // number of basis funcs for element
            if(bf.L == std::vector{0,0,0}) {
                if(bf.element==1) {
                    numBasisElement = 1;
                } else {
                    numBasisElement = 4;  
                }

                double sum=0;
                for(auto k=0; k<numBasisElement; ++k) {
                    sum += palpha(basis_idx+k, basis_idx+k) + pbeta(basis_idx+k, basis_idx+k);
                }

                ptot(elementPos, 0) = sum;
                
                if(bf.element==1) {
                    basis_idx += 1;
                } else {
                    basis_idx += 4;  
                }
                ++elementPos;  
            } 
        }

        // create CNDO/2 Fock Matrices
        arma::Mat<double> fockMatAlpha(palpha.n_rows, palpha.n_cols);
        currentMatAlpha = createCNDO2FockMat(molecule, ptot, palpha, gammaMat, overlapMat);

        arma::Mat<double> fockMatBeta(palpha.n_rows, palpha.n_cols);
        currentMatBeta = createCNDO2FockMat(molecule, ptot, pbeta, gammaMat, overlapMat);
        
        // Calculate difference between old and new density matrices, palpha and pbeta
        if(count!=0) {
            arma::Mat<double> palphaDiffMat(palphaOld);
            arma::Mat<double> pbetaDiffMat(pbetaOld);
            for(int i=0; i<palphaDiffMat.n_rows; ++i) {
                for(int j=0; j<palphaDiffMat.n_cols; ++j) {
                    palphaDiffMat(i,j) = std::abs(palphaOld(i,j)-palpha(i,j));
                    pbetaDiffMat(i,j) = std::abs(pbetaOld(i,j)-pbeta(i,j));
                }
            }
            
            // Find the maximum difference between the matrices to test for convergence
            pDiff = palphaDiffMat.max();
            if(pDiff<pbetaDiffMat.max()) {
                pDiff = pbetaDiffMat.max();
            }
        }
        palphaOld = palpha;
        pbetaOld = pbeta;
        ptotal = palpha + pbeta;

        ++count;
    }
    
    // Calculate Nuclear Repulsion Energy
    molecule.nuclearRepulsionEnergy = 0;
    for(int A=0; A<molecule.numElements; ++A) {
        for(int B=0; B<A; ++B) {
            double dist = std::sqrt(std::pow((molecule.atoms[A].x - molecule.atoms[B].x), 2)+
                                    std::pow((molecule.atoms[A].y - molecule.atoms[B].y), 2)+
                                    std::pow((molecule.atoms[A].z - molecule.atoms[B].z), 2));
                          
            molecule.nuclearRepulsionEnergy+=molecule.atoms[A].Z*molecule.atoms[B].Z/dist;
        }

    }
    molecule.nuclearRepulsionEnergy *= EV_AU_CONV;

    // Calculate Nuclear Repulsion Grad
    arma::Mat<double> nuclearRepulsionGrad = createNuclearRepulsionGradMat(molecule);
    nuclearRepulsionGrad.print(out, "Nuclear Repulsion Grad Matrix (V):");
    
    // Calculate Total Electron Energy
    double sumPHFalpha = 0;
    double sumPHFbeta = 0;
    for(int u=0; u<coreHamiltonianMat.n_cols; ++u) {
        for(int v=0; v<coreHamiltonianMat.n_rows; ++v) {
            sumPHFalpha += palpha(u,v)*(coreHamiltonianMat(u,v)+currentMatAlpha(u,v));
            sumPHFbeta += pbeta(u,v)*(coreHamiltonianMat(u,v)+currentMatBeta(u,v));
        }
        
    }
    molecule.ElectronEnergy = 0.5*sumPHFalpha + 0.5*sumPHFbeta;

    // Calculate Total Electron Grad
    arma::Mat<double> x = createXMat(molecule, ptotal);
    arma::Mat<double> y = createYMat(molecule, ptot, ptotal, palpha, pbeta);
    arma::Mat<double> totalElectronGrad = createElectronGradMat(molecule, x, overlapMatDeriv, y, gammaMatDeriv);
    totalElectronGrad.print(out, "Electron Grad Matrix:");

    // Calculate Total Energy
    molecule.TotalEnergy = molecule.ElectronEnergy + molecule.nuclearRepulsionEnergy;

    // Calculate Total Grad
    arma::Mat<double> totalGrad = totalElectronGrad + nuclearRepulsionGrad;
    totalGrad.print(out, "Total Grad Matrix:");
}

arma::Mat<double> scfWithGradOut(std::ofstream &out, Molecule &molecule, arma::Mat<double> coreHamiltonianMat, \
                                 arma::Mat<double> gammaMat, arma::Mat<double> overlapMat, 
                                 arma::Mat<double> gammaMatDeriv, arma::Mat<double> overlapMatDeriv, double cutoff=1e-6) {
    arma::Mat<double> currentMatAlpha = coreHamiltonianMat;
    arma::Mat<double> currentMatBeta = coreHamiltonianMat;
    arma::Mat<double> palphaOld = coreHamiltonianMat;
    arma::Mat<double> pbetaOld = coreHamiltonianMat;
    arma::Mat<double> palpha = coreHamiltonianMat;
    arma::Mat<double> pbeta = coreHamiltonianMat;
    arma::Mat<double> ptot(molecule.numElements, 1);
    arma::Mat<double> ptotal; // not to be confused with ptot; ptotal = palpha + pbeta
    arma::vec Ea, Eb;
    arma::mat Ca, Cb;

    double pDiff = 1;
    int count = 0;
    while(pDiff>cutoff) {
        // std::cout << pDiff << std::endl;
        // Fock Matrices 
        arma::Mat<double> Fa(currentMatAlpha);
        arma::Mat<double> Fb(currentMatBeta);

        // Solve Eigenvectors and Eigenvalues
        arma::eig_sym(Ea, Ca, Fa);
        arma::eig_sym(Eb, Cb, Fb);

        // Create Density Matrices, p
        palpha = createDensityMat(Ca, molecule.p);
        pbeta = createDensityMat(Cb, molecule.q);

        // create ptot matrix
        int basis_idx = 0;
        int elementPos = 0;
        for(auto &bf:molecule.basisFunctions) {
            int numBasisElement = 0;  // number of basis funcs for element
            if(bf.L == std::vector{0,0,0}) {
                if(bf.element==1) {
                    numBasisElement = 1;
                } else {
                    numBasisElement = 4;  
                }

                double sum=0;
                for(auto k=0; k<numBasisElement; ++k) {
                    sum += palpha(basis_idx+k, basis_idx+k) + pbeta(basis_idx+k, basis_idx+k);
                }

                ptot(elementPos, 0) = sum;
                
                if(bf.element==1) {
                    basis_idx += 1;
                } else {
                    basis_idx += 4;  
                }
                ++elementPos;  
            } 
        }

        // create CNDO/2 Fock Matrices
        arma::Mat<double> fockMatAlpha(palpha.n_rows, palpha.n_cols);
        currentMatAlpha = createCNDO2FockMat(molecule, ptot, palpha, gammaMat, overlapMat);

        arma::Mat<double> fockMatBeta(palpha.n_rows, palpha.n_cols);
        currentMatBeta = createCNDO2FockMat(molecule, ptot, pbeta, gammaMat, overlapMat);
        
        // Calculate difference between old and new density matrices, palpha and pbeta
        if(count!=0) {
            arma::Mat<double> palphaDiffMat(palphaOld);
            arma::Mat<double> pbetaDiffMat(pbetaOld);
            for(int i=0; i<palphaDiffMat.n_rows; ++i) {
                for(int j=0; j<palphaDiffMat.n_cols; ++j) {
                    palphaDiffMat(i,j) = std::abs(palphaOld(i,j)-palpha(i,j));
                    pbetaDiffMat(i,j) = std::abs(pbetaOld(i,j)-pbeta(i,j));
                }
            }
            
            // Find the maximum difference between the matrices to test for convergence
            pDiff = palphaDiffMat.max();
            if(pDiff<pbetaDiffMat.max()) {
                pDiff = pbetaDiffMat.max();
            }
        }
        palphaOld = palpha;
        pbetaOld = pbeta;
        ptotal = palpha + pbeta;

        ++count;
    }
    
    // Calculate Nuclear Repulsion Energy
    molecule.nuclearRepulsionEnergy = 0;
    for(int A=0; A<molecule.numElements; ++A) {
        for(int B=0; B<A; ++B) {
            double dist = std::sqrt(std::pow((molecule.atoms[A].x - molecule.atoms[B].x), 2)+
                                    std::pow((molecule.atoms[A].y - molecule.atoms[B].y), 2)+
                                    std::pow((molecule.atoms[A].z - molecule.atoms[B].z), 2));
                          
            molecule.nuclearRepulsionEnergy+=molecule.atoms[A].Z*molecule.atoms[B].Z/dist;
        }

    }
    molecule.nuclearRepulsionEnergy *= EV_AU_CONV;

    // Calculate Nuclear Repulsion Grad
    arma::Mat<double> nuclearRepulsionGrad = createNuclearRepulsionGradMat(molecule);
    // nuclearRepulsionGrad.print(out, "Nuclear Repulsion Grad Matrix (V):");
    
    // Calculate Total Electron Energy
    double sumPHFalpha = 0;
    double sumPHFbeta = 0;
    for(int u=0; u<coreHamiltonianMat.n_cols; ++u) {
        for(int v=0; v<coreHamiltonianMat.n_rows; ++v) {
            sumPHFalpha += palpha(u,v)*(coreHamiltonianMat(u,v)+currentMatAlpha(u,v));
            sumPHFbeta += pbeta(u,v)*(coreHamiltonianMat(u,v)+currentMatBeta(u,v));
        }
        
    }
    molecule.ElectronEnergy = 0.5*sumPHFalpha + 0.5*sumPHFbeta;

    // Calculate Total Electron Grad
    arma::Mat<double> x = createXMat(molecule, ptotal);
    arma::Mat<double> y = createYMat(molecule, ptot, ptotal, palpha, pbeta);
    arma::Mat<double> totalElectronGrad = createElectronGradMat(molecule, x, overlapMatDeriv, y, gammaMatDeriv);
    // totalElectronGrad.print(out, "Electron Grad Matrix:");

    // Calculate Total Energy
    molecule.TotalEnergy = molecule.ElectronEnergy + molecule.nuclearRepulsionEnergy;

    // Calculate Total Grad
    arma::Mat<double> totalGrad = totalElectronGrad + nuclearRepulsionGrad;
    // totalGrad.print(out, "Total Grad Matrix:");

    return totalGrad;
}

/*
 Geometric Optimization
*/
double l2_norm(arma::Mat<double> totalGrad) {
    /*
    Calculate the Euclidean norm of the force vector

    Parameters
    ----------
    atom : Atom
        The atom to calculate the L2 norm of.
    
    Returns
    -------
    l2_norm : double 
        The norm of the force vector

    */
    double l2_norm = 0;
    for(int i=0; i<totalGrad.n_cols; ++i) {
        l2_norm += std::sqrt(std::pow(totalGrad(0, i), 2) + std::pow(totalGrad(1, i), 2) + std::pow(totalGrad(2, i), 2));
    }
    return l2_norm;
}

double backtrackingLineSearch(std::ofstream &out, Molecule molecule, arma::Mat<double> coreHamiltonianMat, \
                              arma::Mat<double> gammaMat, arma::Mat<double> overlapMat, \
                              arma::Mat<double> gammaMatDeriv, arma::Mat<double> overlapMatDeriv) {
    double alpha=0.00001; // max step size
    double beta=0.1; // proportion to backtrack step
    double gamma=0.0001; // Armijo variable for sufficient decrease

    // Calculate force grad
    arma::Mat<double> totalGrad = scfWithGradOut(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, 1e-6);

    // Create force gradient and energy
    arma::Mat<double> forceGrad = totalGrad*-1;

    // Calculate norm of the force
    double forceNorm = l2_norm(forceGrad);

    // Enter loop until an appropriate step is found
    // Potential for infinite loop, possible correction
    //     -Set limit for alpha and once below, set alpha to larger than original.
    //     -Possibly use recursive calls to enlarge alpha. 
    for(auto i=0; i<1000; ++i) {
        // Create copy of molecule to adjust
        Molecule testMolecule = molecule;

        for(auto i=0; i<testMolecule.atoms.size(); ++i) {
            molecule.atoms[i].x += alpha*forceGrad(0,i);
            molecule.atoms[i].y += alpha*forceGrad(1,i);
            molecule.atoms[i].z += alpha*forceGrad(2,i);
        }

        // update set of basisFuncs
        int elementID = -1;
        for(auto i=0; i<testMolecule.basisFunctions.size(); ++i) {
            if(molecule.basisFunctions[i].L == std::vector{0,0,0}){++elementID;}
            molecule.basisFunctions[i].center[0] += alpha*forceGrad(0,elementID);
            molecule.basisFunctions[i].center[1] += alpha*forceGrad(1,elementID);
            molecule.basisFunctions[i].center[2] += alpha*forceGrad(2,elementID);
        }

        // Calculate energy after adjustment
        scfWithGradOut(out, testMolecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, 1e-6);

        // Check Armijo condition (sufficient decrease)
        if(testMolecule.TotalEnergy <= molecule.TotalEnergy + gamma*alpha*forceNorm) {
            return alpha;
        
        // If the condition is not satisfied, try a smaller step
        } else {
            alpha *= beta;
        }

    }

    return alpha;
}

void steepestDescentBacktracking(std::ofstream &out, Molecule &molecule, arma::Mat<double> coreHamiltonianMat, \
                                 arma::Mat<double> gammaMat, arma::Mat<double> overlapMat, \
                                 arma::Mat<double> gammaMatDeriv, arma::Mat<double> overlapMatDeriv, \
                                 int maxIter=1000,  double convergLimit=0.01) {
    arma::Mat<double> totalGrad;
    for(auto i=0; i<maxIter; ++i){
        // Calculate force grad and energy
        totalGrad = scfWithGradOut(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, 1e-6);

        // Create unit gradient
        arma::Mat<double> forceGrad = totalGrad*-1;

        // Caclulate norm Force
        double normForce = l2_norm(forceGrad);
        out << "Step = " << i << "; E = " << molecule.TotalEnergy << "; F_norm = " << normForce << std::endl;

        // Use backtracking line search to find stepsize
        double stepsize = backtrackingLineSearch(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv);
        

        // cut off should be for normForce or should cut off be for each atom?
        if(normForce > convergLimit) {
            
            // Update via steepest descent
            for(auto ii=0; ii<molecule.atoms.size(); ++ii) {
                molecule.atoms[ii].x += stepsize*forceGrad(0, ii);
                molecule.atoms[ii].y += stepsize*forceGrad(1, ii);
                molecule.atoms[ii].z += stepsize*forceGrad(2, ii);
            }

            // update set of test basisFuncs
            int elementID = -1;
            for(auto i=0; i<molecule.basisFunctions.size(); ++i) {
                if(molecule.basisFunctions[i].L == std::vector{0,0,0}){++elementID;}
                molecule.basisFunctions[i].center[0] += stepsize*forceGrad(0,elementID);
                molecule.basisFunctions[i].center[1] += stepsize*forceGrad(1,elementID);
                molecule.basisFunctions[i].center[2] += stepsize*forceGrad(2,elementID);
            }

        } else {
            out << "Converged in " << i+1 << " steps. To a cluster energy of " << molecule.TotalEnergy << std::endl;
            break;  
        }
    }
    
    out << "----Optimized Positions----" << std::endl;
    for(auto i=0; i<molecule.atoms.size(); i++) {
        out << molecule.atoms[i].x << ", " << molecule.atoms[i].y << ", " << molecule.atoms[i].z << std::endl;
    }
    totalGrad.print(out, "Optimized Grad");
}

/*
 Vibrational Analysis
*/
void makeMassWeightedHessian(Molecule molecule, arma::Mat<double> &Hessian) {
    int N = molecule.atoms.size();
    int A = 3*N;
    int x, y;
    double massX, massY, mass;
    for(int atomX=0; atomX<N; ++atomX) {
        for(int dirX=0; dirX<3; ++dirX) {
            for(int atomY=0; atomY<N; ++atomY) {
                for(int dirY=0; dirY<3; ++dirY) {
                    x = atomX*3 + dirX;
                    y = atomY*3 + dirY; 
                    massX = molecule.atoms[atomX].mass;
                    massY = molecule.atoms[atomY].mass;
                    mass = std::sqrt(massX*massY);
                    Hessian(x,y) = Hessian(x,y)/mass;
                }
            }
        }
    }
}

void vibrationalAnalysis(std::ofstream &out, Molecule molecule, double h=1e-4, double cutoff=1e-6) {
    
    // make copy so I don't mess up original molecule object (just in case)
    // Molecule copymol = molecule;
    int A = molecule.atoms.size()*3;
    arma::Mat<double> Hessian(A, A);
    int numBasisF = numBasisFuncs(molecule);

    int basisStart = 0;
    for(int atom=0; atom<molecule.atoms.size(); ++atom) {
        for(int xyz=0; xyz<3; ++xyz) {
            // std::cout << "TESTxyz" << std::endl;
            // // print statment to check for atom position. (Paste where needed.)
            // if(xyz==0) {
            //     std::cout << molecule.atoms[atom].x << "  ";
            // } else if(xyz==1) {
            //     std::cout << molecule.atoms[atom].y << "  ";
            // } else if(xyz==2) {
            //     std::cout << molecule.atoms[atom].z << "  ";
            // }
            // for(int bf=basisStart; bf<basisStart+molecule.atoms[atom].orbitals.size(); ++bf) {
            //     std::cout << molecule.basisFunctions[bf].center[xyz] << std::endl;
            // }

            // update atom coordinate PLUS
            if(xyz==0) {
                molecule.atoms[atom].x += h;
            } else if(xyz==1) {
                molecule.atoms[atom].y += h;
            } else if(xyz==2) {
                molecule.atoms[atom].z += h;
            }
            //update basis functions coordinates for atom
            for(int bf=basisStart; bf<basisStart+molecule.atoms[atom].orbitals.size(); ++bf) {
                molecule.basisFunctions[bf].center[xyz] += h;
            }
            
            // Overlap Matrix (create vector of vectors, feed into matrix)
            arma::Mat<double> overlapMat(numBasisF, numBasisF);
            overlapMat = createMatrix(createOverlapMatrix(molecule));

            // Overlap Matrix Gradient
            arma::Mat<double> overlapMatDeriv(3, numBasisF*numBasisF);
            overlapMatDeriv = createOverlapDerivMat(molecule);
            // overlapMatDeriv.print(out, "Overlap Grad Matrix:");

            // Create gamma
            arma::Mat<double> gammaMat = createMatrix(createGammaMatrix(molecule));
            gammaMat *= EV_AU_CONV;

            // Gamma Matrix Gradient
            arma::Mat<double> gammaMatDeriv(3, numBasisF*numBasisF);
            gammaMatDeriv = createGammaMatrixGrad(molecule);
            gammaMatDeriv *= EV_AU_CONV;
            // gammaMatDeriv.print(out, "Gamma Grad Matrix:");

            // Create Core Hamiltonian Matrix
            arma::Mat<double> coreHamiltonianMat(numBasisF, numBasisF);
            coreHamiltonianMat = createMatrix(createCoreHamiltonianMatrix(molecule, gammaMat, overlapMat));

            // calculate PLUS grad
            arma::Mat<double> totalGrad_plus = scfWithGradOut(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, cutoff);
            // update atom coordinate to MINUS
            if(xyz==0) {
                molecule.atoms[atom].x -= 2*h;
            } else if(xyz==1) {
                molecule.atoms[atom].y -= 2*h;
            } else if(xyz==2) {
                molecule.atoms[atom].z -= 2*h;
            }
            //update basis functions coordinates for atom
            for(int bf=basisStart; bf<basisStart+molecule.atoms[atom].orbitals.size(); ++bf) {
                molecule.basisFunctions[bf].center[xyz] -= 2*h;
            }

            // Overlap Matrix (create vector of vectors, feed into matrix)
            overlapMat = createMatrix(createOverlapMatrix(molecule));

            // Overlap Matrix Gradient
            overlapMatDeriv = createOverlapDerivMat(molecule);
            // overlapMatDeriv.print(out, "Overlap Grad Matrix:");

            // Create gamma
            gammaMat = createMatrix(createGammaMatrix(molecule));
            gammaMat *= EV_AU_CONV;

            // Gamma Matrix Gradient
            gammaMatDeriv = createGammaMatrixGrad(molecule);
            gammaMatDeriv *= EV_AU_CONV;
            // gammaMatDeriv.print(out, "Gamma Grad Matrix:");

            // Create Core Hamiltonian Matrix
            coreHamiltonianMat = createMatrix(createCoreHamiltonianMatrix(molecule, gammaMat, overlapMat));

            // calculate PLUS grad
            arma::Mat<double> totalGrad_minus = scfWithGradOut(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, cutoff);
            // RESET COORDS
            // update atom coordinate to ORIGIN
            if(xyz==0) {
                molecule.atoms[atom].x += h;
            } else if(xyz==1) {
                molecule.atoms[atom].y += h;
            } else if(xyz==2) {
                molecule.atoms[atom].z += h;
            }
            //update basis functions coordinates for atom
            for(int bf=basisStart; bf<basisStart+molecule.atoms[atom].orbitals.size(); ++bf) {
                molecule.basisFunctions[bf].center[xyz] += h;
            }

            // FIND DIFFERENCE
            arma::Mat<double> doublederiv = totalGrad_plus - totalGrad_minus;
            doublederiv = doublederiv/(2*h);
            doublederiv.reshape(A,1);
            
            for(int i=0; i<A; ++i) {
                Hessian(i,atom*3+xyz) = doublederiv(i,0);
            }
        }
        basisStart += molecule.atoms[atom].orbitals.size();
    }
    // average error ***
    for(int i=0; i<A; ++i) {
        for(int j=0; j<i; ++j) {
            double avg = (Hessian(i,j) + Hessian(j,i)) / 2;
            Hessian(i,j) = avg;
            Hessian(j,i) = avg;
        }
    }

    /*
        CONVERT UNITS
    */
    // Hessian = Hessian/EV_AU_CONV; // Conversion from EV to AU
    // Hessian = Hessian*627.5; // Conversion to kcal/mol/au^2
    Hessian = Hessian*23.06035; // conversion from eV to kcal/mol
    Hessian = Hessian/(0.529177249*0.529177249); // Conversion to kcal/mol/A^2; 1 a.u. = 0.529177249 Angstrom; kcal/mol/A^2
    Hessian = ((1e8 * 4184)/(1e-10 * 6.022e23)) * Hessian; // (J/m to millidynes)(Kcal to J)/(Ang to m)(Mol to molecule) = millidynes/A

    // Hessian = Hessian*1.602176e-19; // (from eV to J)
    // Hessian = Hessian/(0.529177e-10*0.529177e-10); // (/au^2 to /m^2)
    // Hessian = Hessian/1.66054e-27; // (/amu to /kg)

    // Mass weighted Hessian
    Hessian.print(out, "Hessian");
    makeMassWeightedHessian(molecule, Hessian); // millidynes/A/amu
    Hessian.print(out, "Mass Weighted");

    

    // // Diagnolize(molecule, Hessian)
    // arma::Mat<double> diagMat = arma::diagmat(Hessian);
    // diagMat.print("Diagonalized Hessian");  

    arma::vec E;
    arma::Mat<double> V;
    arma::eig_sym(E, V, Hessian);
    E.print(out, "Force constants");
    V.print(out, "Modes");
    // arma::eig_sym(E, V, Hessian);
    // E.print("H - Frequencies");
    // V.print("H - Modes");


    // Calculate frequency
    for(auto i=0; i<A; ++i) {
        
        if(E(i) > 1e-10) {
            // E(i) = 4.13567e16/2.9979e10 * std::sqrt(E(i)); 
            E(i) = 1/(2*M_PI*2.9979e10) * std::sqrt(1e5*6.022e23*E(i)); // (1e2)*(1e3)*(6.022e23); from millidynes/A (1 mDy/A = 1e2 J/m^2; 1 amu = 1e-3 / 6.022e23 kg) = N/m/kg
            // E(i) = 1302.79 * std::sqrt(E(i)); //
        } else {
            E(i) = 0;
        }
        
         
    }

    E.print(out, "Frequencies (cm^-1)");




}
/*
 Solving
*/

int main(int argc, char** argv) {
    Molecule molecule;
    // User input
    std::string filepath = argv[1];
    int p, q;
    if(argc==4) {
        p = std::stoi(argv[2]);
        q = std::stoi(argv[3]);
    }
    std::string location, filename;
    std::stringstream path(filepath);
    std::string outpath;
    std::vector<std::string> pathList, fileList;
    
    while(std::getline(path, location, '/')) {
        pathList.push_back(location);
    }

    // change to output 
    pathList[1] = "output";

    // change .txt to .out
    std::stringstream path2(pathList[2]);
    while(std::getline(path2, location, '.')) {
        fileList.push_back(location);
    }
    fileList[1] = "out";
    filename = fileList[0] + "." + fileList[1];
    pathList[2] = filename;

    for(auto i=0; i<pathList.size()-1; i++) {
        outpath += pathList[i];
        outpath += "/";
    }
    outpath += pathList.back();

    // create output file stream
    std::ofstream out(outpath);

    // Read file, create molecule
    molecule = readMolecule(filepath);

    // Some prep
    assignShells(molecule);
    createBasisFuncs(molecule);
    normalizeBasisFuncs(molecule);
    int numBasisF = numBasisFuncs(molecule);

    // Assign semi-empirical parameters to basisfunctions
    assignSemiEmpiricalParams(molecule);

    // Overlap Matrix (create vector of vectors, feed into matrix)
    arma::Mat<double> overlapMat(numBasisF, numBasisF);
    overlapMat = createMatrix(createOverlapMatrix(molecule));
    
    // Overlap Matrix Gradient
    arma::Mat<double> overlapMatDeriv(3, numBasisF*numBasisF);
    overlapMatDeriv = createOverlapDerivMat(molecule);
    overlapMatDeriv.print(out, "Overlap Grad Matrix:");

    // Create gamma
    arma::Mat<double> gammaMat = createMatrix(createGammaMatrix(molecule));
    gammaMat *= EV_AU_CONV;

    // Gamma Matrix Gradient
    arma::Mat<double> gammaMatDeriv(3, numBasisF*numBasisF);
    gammaMatDeriv = createGammaMatrixGrad(molecule);
    gammaMatDeriv *= EV_AU_CONV;
    gammaMatDeriv.print(out, "Gamma Grad Matrix:");

    // Assign or calculate p and q
    if(argc==4) {
        molecule.numE = p + q;
        molecule.p = p;
        molecule.q = q;
    } else {
        assignNumElectrons(molecule);
    }

    // Create Core Hamiltonian Matrix
    arma::Mat<double> coreHamiltonianMat(numBasisF, numBasisF);
    coreHamiltonianMat = createMatrix(createCoreHamiltonianMatrix(molecule, gammaMat, overlapMat));

    // begin algorithm
    scfWithGrad(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, 1e-6);

    // print energies
    out << "Nuclear Repulsion Energy is " << molecule.nuclearRepulsionEnergy << " eV." << std::endl;
    out << "Electron Energy is " << molecule.ElectronEnergy << " eV." << std::endl;
    out << "The molecule in file ." << filepath << " has energy " << std::setprecision(10) << molecule.TotalEnergy << " eV." << std::endl;
    out << std::endl;

    // Optimization with Steepest Descent with Backtracking
    out << "Beginning Steepest Descent with Backtracking (Armijo condition) to optimize location." << std::endl;
    steepestDescentBacktracking(out, molecule, coreHamiltonianMat, gammaMat, overlapMat, gammaMatDeriv, overlapMatDeriv, 1000, 0.001);
    out << std::endl;

    // Begin vibrational analysis
    out << "Beginning Vibrational Analysis." << std::endl;
    vibrationalAnalysis(out, molecule, 1e-5, 1e-6);
    out << std::endl;

    return 0;
}