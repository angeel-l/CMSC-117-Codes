/* Sample implementation of bisection method. */
#include "rootscalar.hpp"

double f(double x)
{
    // return std::pow(x, 3) - 6 * std::sqrt(2);
    return 0.25 * std::pow(std::cos(2. * x), 2) - std::pow(x, 2);
}

double df(double x)
{
    return -std::cos(2. * x) * std::sin(2. * x) - 2 * x;
}
double dplusf(double x, double h)
{
    return (f(x + h) - f(x)) / h;
}
double dminusf(double x, double h)
{
    return (f(x) - f(x - h)) / h;
}
double d0f(double x, double h)
{
    return (f(x + h) - f(x - h)) / (2 * h);
}

double g(double x)
{
    return 0.5 * std::cos(2 * x);
}

int main()
{
    // print filename
    std::cout << "File: " << __FILE__ << std::endl;
    // parameter object
    root::scalar::param parameter;
    parameter.tol = 1e-15;
    parameter.maxit = 100;
    parameter.wa = 0.9;
    parameter.wf = 0.1;
    // approximate a root by bisection method
    std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
    std::cout << "METHOD \t\tAPPROXIMATE ROOT \tFUNCTION VALUE \t\t\tERROR \t\t\t\t\tNITERS" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
    root::scalar::RootScalarResult _result;
    // _result = bisection(f, 2., 2.5, parameter);
    // _result.print();
    _result = chord(f, 0., 0.5, 1, parameter);
    _result.print_edit();
    _result = regulafalsi(f, 0.5, 0., parameter);
    _result.print_edit();
    _result = secant(f, 0.5, 0., parameter);
    _result.print_edit();
    _result = newton(f, df, 0.5, parameter);
    _result.print_edit();
    _result = innewton(f, dplusf, 0.5, "Inexact Newton", parameter);
    _result.print_edit();
    _result = steffensen(f, 0.5, parameter);
    _result.print_edit();
    _result = fixpoint(g, 0.5, parameter);
    _result.print();
    // _result = aitken(g, 0.5, parameter);
    // _result.print();
    return 0;
}