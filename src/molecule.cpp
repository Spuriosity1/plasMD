#include "molecule.hpp"

Molecule::Molecule(){
    orbitals.resize(nuclei.size());
    for (size_t i = 0; i < nuclei.size(); i++) {
        orbitals[i]=allocate(nuclei[i].Z);
    }
}

void Molecule::comp_charges(){
    for (size_t i = 0; i < orbitals.size(); i++) {
        orbitals[i].q = nuclei[i].Z - orb_charge(orbitals[i]);
    }
}

double orb_charge(orbital_t *orb){
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        for (size_t l = 0; l < n; l++) {
            q -= orb->bound[n-1][l];
        }
    }
}

orbital_t *Molecule::allocate_orb(int Z){
    // Initialises to a ground state containing num electrons
    orbital_t *orb = new orbital_t;
    orb.bound.resize(BOUND_MAX_N);
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        orb.bound[n-1] = new occ_t[n];
    }
    orb.q=Z;
    // Degeneracy given by 2(2l+1)
    // Allocate in order 1s, 2s, 2p, 3s, 3p, 4s, 3d,
    // TODO: allocate in appropriate order
    return &orb;
}

void Molecule::deallocate_orb(orbital_t* orb){
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        delete[] orb->bound[n-1];
    }
}

void Molecule::timestep(){
    for (size_t a = 0; a < orbitals.size(); a++) {
        orb_update(orbitals[i]);
    }
}

void orb_update(orbital_t orb){
    double delta = 0;
    for (size_t n = 0; n < BOUND_MAX_N; n++) {
        for (size_t l = 0; l < n; l++) {
            orb.bound[n-1][l] += calc_orb_delta(n,l);
        }
    }
}

double calc_orb_delta(orbital_t orb, unsigned n, unsigned l){

}

//
// void print_bound(orbital_t orb){
//     printf("\n--BOUND--\n");
//     printf("n | l\n  | 1  2  3  4  5  6  7  8  9  10\n");
//     for (size_t n = 0; n < BOUND_MAX_N; n++) {
//         printf("%2d|",n);
//         for (size_t l = 0; l < n; l++) {
//             printf(" %2d",orb.bound[n-1][l]);
//         }
//         printf("\n");
//     }
//     printf("\n--FREE--\n");
// }
