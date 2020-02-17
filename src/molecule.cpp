#include "molecule.hpp"

Molecule::Molecule(){
    orbitals.resize(nuclei.size());
    for (size_t i = 0; i < nuclei.size(); i++) {
        orbitals[i]=allocate(nuclei[i].Z);
    }
}



orbital_t allocate(int Z){
    // Initialises to a ground state containing num electrons
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        bound[n-1] = (occ_t*) malloc(sizeof(occ_t)*n);
    }
    // Degeneracy given by 2(2l+1)
    // Allocate in order 1s, 2s, 2p, 3s, 3p, 4s, 3d,
    // TODO: allocate in appropriate order
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
