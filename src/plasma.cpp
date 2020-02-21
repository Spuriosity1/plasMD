#include "plasma.hpp"

using namespace Plasma;

Plasma(){
    // Error bounds should go here

}

// dt of f
void timestep(){
    for (size_t i = 0; i < NPTSF; i++) {
        f[i] += QB(v(i));
    }
}

double QB(double e){
    double integral = 0;
    for (size_t i = 0; i < NPTSF; i++) {
        f[i]
    }
    return pow(e,-0.5)*integral;
}

double v(int i){
    // Returns velocity indexed by integer i
    // Here, not very interesting
    return 1.*(EMAX-EMIN)*i/NPTSF;
}
