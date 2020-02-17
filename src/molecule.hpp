#ifndef MOLECULE_CPP_H
#define MOLECULE_CPP_H

#define GTO_EXPANSION_ORDER 10
#define BOUND_MAX_N 6

// This class contains all of the electronic information involved in the molecular solution

#include "nuclei.hpp"

typedef double data_t;
typedef double occ_t;


typedef struct{
    occ_t* bound[BOUND_MAX_N];
} orbital_t;
// Structure for storing population associated with a particular nucleus

class Molecule : public Nuclei {
public:
    Molecule();
protected:
    std::vector<orbital_t> orbitals;

private:
    orbital_t allocate(int num);
};


#endif /* end of include guard: MOLECULE_CPP_H */
