#include <iostream>
#include "molecule.hpp"





#include "molecule.hpp"
#include <iostream>
// CALL SEQUENCE
// plasMD input_file.dat
// input_file has the 3D coordinates of all atoms, with some syntax specified
// in moleculeio.h
int main(int argc, char const *argv[]) {
    if (argc < 2){
        std::cerr << "Please specify an input file." <<std::endl;
        std::cerr << FORMAT_SPECIFICATION <<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string path(argv[1]);

    Molecule m;
    m.import_from(path);
    m.print();




    return 0;
}
