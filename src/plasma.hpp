#ifndef PLASMA_HPP
#define PLASMA_HPP

#include <gsl/integration>

/*
Dimensions:
    NPTSF = number of points in finite-differenced electron energy
            distribution function
    NSPECIES = number of different chemical species that the code can handle
    NLEV = maximum number of states or levels allowed for each species
*/
#define NPTSF 10000
#define EMAX 50000
#define EMIN 0
// measured in eV
// Consider gridding this differently



/*
Represents the electron plasma in the vicinity of the ionised molecule
f[i] - sclar electron distribution function
*/
class Plasma{
public:
    Plasma(double vmax);
    void simulate(int nts);
    void simulate(double T);
private:
    double f[NPTSF];
    double dt;
    void ts();
};

#endif
