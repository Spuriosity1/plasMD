#ifndef MOLECULE_CPP_H
#define MOLECULE_CPP_H

#define GTO_EXPANSION_ORDER 10
#define BOUND_MAX_N 4
/*
Num (nl) states -- sum (2n+1)
*/
// Number of available configurations increases with N^2
// This class contains all of the electronic information involved in the molecular solution

#include "nuclei.hpp"

typedef double data_t;
typedef double occ_t;


typedef struct{
    double state[2*(2*BOUND_MAX_N+1)];
} orbital_t;
// Structure for storing population associated with a particular nucleus

class Molecule : public Nuclei {
public:
    Molecule();
protected:
    std::vector<orbital_t> orbitals;

private:
    orbital_t *allocate_orb(int Z);
};


/*
STATE STORAGE STRUCTURE
Would like this to be VERY FAST - store by an index

*/

#endif /* end of include guard: MOLECULE_CPP_H */
