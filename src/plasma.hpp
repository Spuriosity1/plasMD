

/*
Dimensions:                                                         ABLX0062
    NPTSF = number of points in finite-differenced electron energy     ABLX0063
            distribution function                                      ABLX0064
    NSPECIES = number of different chemical species that the code can handleABLX0069
    NLEV = maximum number of states or levels allowed for each species ABLX0070
*/

class Plasma{
public:
    Plasma(int nptsf);
private:
    double f[];
}
