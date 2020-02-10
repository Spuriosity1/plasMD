#include "molecule.hpp"

Molecule::Molecule(){
    orbitals.resize(nuclei.size());
    // Initialises to a ground state containing num electrons
    for (size_t n = 1; n <= BOUND_MAX_N; n++) {
        bound[n-1] = (occ_t*) malloc(sizeof(occ_t)*n);
        for (occ_t l = 0; l < n && num > 0; l++) {
            num -= 2*(2*l+1);
            if (num < 0){
                bound[n-1][l] = num;
                break;
            } else {
                bound[n-1][l] = 2*(2*l+1);
            }
        }
    }
}

void print(){
    printf("n | l\n  | 1  2  3  4  5  6  7  8  9  10\n");
    for (size_t n = 0; n < BOUND_MAX_N; n++) {
        printf("%2d|",n);
        for (size_t l = 0; l < n; l++) {
            printf(" %2d",l);
        }
        printf("\n");
    }
}
