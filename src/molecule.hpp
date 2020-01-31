#ifndef MOLECULE_CPP_H
#define MOLECULE_CPP_H

// #include "adhocmath.h"

#define GTO_EXPANSION_ORDER 10

// This class contains all of the electronic information involved in the molecular solution

#include "nuclei.hpp"

// Structure for storing population associated with a particular nucleus
typedef struct{
    uint Z;
} electron_t;


class Molecule : public Nuclei {
public:
    Molecule();
protected:
    std::vector<electron_t> orbitals;
};


#endif /* end of include guard: MOLECULE_CPP_H */
