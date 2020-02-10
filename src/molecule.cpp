#include "molecule.hpp"

Molecule::Molecule(){
    orbitals.resize(nuclei.size());
    // Initialises to a ground state containing num electrons
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        bound[n-1] = (occ_t*) malloc(sizeof(occ_t)*n);
        for (occ_t l = 0; l < n && num > 0; l++) {
            if (num < 2*(2*l+1)){
                bound[n-1][l] = num;
                num = 0;
            } else {
                bound[n-1][l] = 2*(2*l+1);
                num -= 2*(2*l+1);
            }
        }
    }
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
