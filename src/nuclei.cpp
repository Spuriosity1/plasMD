#include "nuclei.hpp"


Nuclei::Nuclei(const std::string &path){
    import_from(path);
}

// Prints a list of all coordinates of nuclear positions stored in nuclei[]
void Nuclei::print_nuclei(){
    printf("LIST OF ALL NUCLEAR POSITIONS:\n");
    for (uint i=0; i<nuclei.size(); i++) {
        vector3_t r = nuclei[i].r;
        printf("Z=%2d, (%3.5lf,%3.5lf,%3.5lf)\n",nuclei[i].Z,r.v[0],r.v[1],r.v[2]);
    }
}


bool zcomp(nucleus_t n1, nucleus_t n2){
    return n1.Z < n2.Z;
}

// Saves new nuclear positions to memory
void Nuclei::save_as(const std::string &filename){
    std::ofstream outdata;
    outdata.open(filename);

    // Prepend some date-time info
    std::time_t time = std::time(nullptr);
    outdata<<"# Saved "<< std::ctime(&time) <<std::endl;

    // Sort by Z to be sure nothing suspicious has happened

    std::sort(nuclei.begin(), nuclei.end(), zcomp);

    uint Z = 0;
    for (size_t i = 0; i < nuclei.size(); i++) {
        if (nuclei[i].Z != Z) {
            // New Z value, make a note of it
            Z = nuclei[i].Z;
            outdata<<"Z:"<<Z<<std::endl;
        }
        vector3_t r = nuclei[i].r;
        outdata<<r.v[0]<<" "<<r.v[1]<<" "<<r.v[2]<<std::endl;
    }
    outdata.close();
}


// Do this the old fashioned way
// This is quite possibly some of the worst code ever written
void Nuclei::import_from(const std::string &path){
    printf("Importing...\n");
    std::ifstream indata;
    indata.open(path);
    if (!indata.is_open()){
        fprintf(stderr, "Failed to load file `%s` !\n",path.c_str());
        return;
    }

    std::string line;
    uint lineno = 0;

    nelectrons = 0;
    // Structure:
    // nuclei[Z][idx][x,y,z]
    uint Zval = 0;
    nuclei.resize(1); // Put in a dummy element so the damn thing doesn't
    // segfault immediately

    nucleus_t nucleus;
    while (getline(indata, line)) {
        lineno++;
        if (line[0] == '#'){
            // Ignore comments
            continue;
        } else if (line[0] == 'Z'){
            // New Z value
            Zval = stoi(line.substr(2));
        } else {
            nelectrons += Zval;
            // Assume we are reading 3 coordinates
            nucleus.r = read_vector3(line);
            nucleus.Z = Zval;

            nuclei.push_back(nucleus);
        }
    }
}


// Reads 3 Euclidean coords like 1 2 -0.222e-4 from a std::string line,
// returns as vector3_t
vector3_t Nuclei::read_vector3(std::string line){
    std::stringstream lineStream(line);
    std::string cell;

    uint pos=0;
    vector3_t R;

    while (getline(lineStream, cell, ' ')) {
        if (pos < 3){
            R.v[pos] = stod(cell);
            ++pos;
        } else {
            fprintf(stderr,"Found extra data!\n");
            // Ignore any superfluous data
        }
    }

    return R;
}
