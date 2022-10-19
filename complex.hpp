#ifndef COMPLEX_HPP_INCLUDE
#define COMPLEX_HPP_INCLUDE

#include <iostream>
#include <complex>
#include <cmath>

typedef std::complex<double> complex_d_t;

// Output stream of a complex number.
std::ostream &operator<<(std::ostream &output, const complex_d_t &z)
{
    if (z.imag() >= 0)
        output << z.real() << " + " << std::abs(z.imag()) << "j";
    else
        output << z.real() << " - " << std::abs(z.imag()) << "j";
    return output;
}

#endif
