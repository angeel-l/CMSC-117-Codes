/* Sample implementation of Cardano's formula. */
#include <iomanip>
#include <cstdio>
#include "cardano.hpp"

int main()
{
    // coefficients
    double a = 1, b = 2, c = 3, d = 4;
    // pre-allocate array of solutions
    complex_d_t x[3];
    // apply Cardano's formula
    cardano(a, b, c, d, x);
    freopen("cardano.txt", "w+", stdout);
    std::cout << "The solutions of the cubic equation $f(x) = "
              << a << "x^3 + " << b << "x^2 + " << c << "x + " << d << " = 0$ "
              << "are:" << std::endl;
    // print solutions and function values
    std::cout << std::scientific << std::setprecision(8);
    for (int k = 0; k < 3; k++)
    {
        std::cout << "\\[x = " << x[k] << "\t";
        std::cout << "\\hspace{1cm} f(x) = " << evalCubicPoly(a, b, c, d, x[k]) << "\\]" << std::endl;
    }
    fclose(stdout);
    return 0;
}