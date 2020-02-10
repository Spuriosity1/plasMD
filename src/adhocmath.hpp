#ifndef ADHOCMATH_CPP_H
#define ADHOCMATH_CPP_H


// template recursion -- this looks a little cursed to be honest
template <unsigned int exponent>
inline double intpow(double base) {
    return intpow<exponent-1>(base) * base;
}

template <>
inline double intpow<0>(double base) {
    return 1;
}


// // small brain
// inline double i_fact(uint n){
//     double retval = 1;
//     for (size_t i = n; i >1; i--) {
//         retval *= i;
//     }
//     return retval;
// }
//
// // big brain
// inline double r_fact(uint n){
//     if (n <= 2){
//         return n;
//     } else {
//         return n*r_fact(n-1);
//     }
// }

// galaxy brain
const static double a_fact[20] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600, 6227020800, 87178291200,
    1307674368000, 20922789888000, 355687428096000, 6402373705728000,
    121645100408832000};

// returns integral of x^(n) exp(-ax^2)
inline double _gint(uint n, double a){
    if (n%2 ==1) return 0;
    n /= 2;
    return pow(a,-n-0.5)*a_fact[2*n]*SQRT_PI/(intpow(4,n) * a_fact[n]);
}

#endif
