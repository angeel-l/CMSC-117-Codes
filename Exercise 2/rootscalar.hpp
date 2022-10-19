/* C++ header file for root-finding algorithms. */
#ifndef ROOTSCALAR_HPP_INCLUDE
#define ROOTSCALAR_HPP_INCLUDE
// standard library includes
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>
#include <vector>

// local include
#include "timer.hpp"
#include "../complex.hpp"

namespace root
{
    namespace scalar
    {

        // type for a function with input and output having the type double
        using UniVarFunction = double(double);
        // type for the dervative of the function and output having the type double
        using DerivFunction = double(double, double);

        // struct for parameters in scalar-root finding algorithms
        struct param
        {
            double tol = std::numeric_limits<double>::epsilon();
            double h = std::sqrt(std::numeric_limits<double>::epsilon());
            int maxit = 1000;
            double wa = 0.5;
            double wf = 0.5;
        };

        // struct for solution to scalar-root finding problems
        struct RootScalarResult
        {
            int numit;
            int maxit;
            double x;
            double funval;
            double error;
            double tol;
            double elapsed_time;
            std::string method_name;
            std::string termination_flag;

            // default constructor
            RootScalarResult() {}

            // user-defined constructor
            RootScalarResult(const int &numit, const int &maxit, const double &x,
                             const double &funval, const double &error, const double &tol,
                             const double &elapsed_time, const std::string &method_name,
                             const std::string &termination_flag)
            {
                this->numit = numit;
                this->maxit = maxit;
                this->x = x;
                this->funval = funval;
                this->error = error;
                this->tol = tol;
                this->elapsed_time = elapsed_time;
                this->method_name = method_name;
                this->termination_flag = termination_flag;
            }

            void print()
            {
                std::cout << "ROOT FINDER: "
                          << method_name << std::endl;
                std::cout << std::setprecision(16);
                std::cout << std::fixed;
                std::cout << "APPROXIMATE ROOT / LAST ITERATE: "
                          << x << std::endl;
                std::cout << "TERMINATION: "
                          << termination_flag << std::endl;
                std::cout << std::scientific;
                std::cout << "FUNCTION VALUE: "
                          << funval << std::endl;
                std::cout << "ERROR: "
                          << error << std::endl;
                std::cout << "TOLERANCE: "
                          << tol << std::endl;
                std::cout << "NUM ITERATIONS: "
                          << numit << std::endl;
                std::cout << "MAX ITERATIONS: "
                          << maxit << std::endl;
                std::cout << "ELAPSED TIME: "
                          << elapsed_time
                          << " seconds" << std::endl;
                std::cout << std::defaultfloat;
            }
            void print_edit()
            {
                std::cout << method_name << "\t"
                          << std::setprecision(16)
                          << std::fixed
                          << x << "\t"
                          << std::scientific
                          << funval << "\t"
                          << error << "\t"
                          << numit << std::endl;
                std::cout << std::defaultfloat;
            }
        };

        // bisection method for approximating a solution of scalar equation f(x) = 0.
        RootScalarResult bisection(UniVarFunction &f, double a, double b, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = b - a;
            double fa = f(a);
            double fb = f(b);
            int k = 0;
            double c;
            if (fa == 0)
            {
                c = a;
                err = 0;
            }
            if (fb == 0)
            {
                c = b;
                err = 0;
            }
            if (fa * fb > 0)
            {
                std::cerr << "Method Fails!" << std::endl;
                return RootScalarResult();
            }
            // main loop
            double fc;
            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                c = (a + b) / 2.0;
                fc = f(c);
                if (fc * fa > 0)
                {
                    a = c;
                    fa = fc;
                }
                else
                {
                    b = c;
                }
                err = std::abs(b - a);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }
            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, c, fc, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Bisection\t", term_flag);
        }

        RootScalarResult chord(UniVarFunction &f, double a, double b, double x, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double q = (f(b) - f(a)) / (b - a);
            double fx = f(x);
            double x0;
            int k = 0;

            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                x = x - fx / q;
                fx = f(x);
                err = parameter.wa * std::abs(x - x0) + parameter.wf * std::abs(fx);
                k++;
            }

            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }

            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Chord\t\t", term_flag);
        }

        RootScalarResult secant(UniVarFunction &f, double x0, double x1, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double f0 = f(x0);
            double f1 = f(x1);
            double q, xtemp;
            int k = 1;
            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                q = (f1 - f0) / (x1 - x0);
                xtemp = x1;
                x1 = x1 - f1 / q;
                x0 = xtemp;
                f0 = f1;
                f1 = f(x1);
                err = parameter.wa * std::abs(x1 - x0) + parameter.wf * std::abs(f1);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }
            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x1, f1, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Secant\t\t", term_flag);
        }

        RootScalarResult regulafalsi(UniVarFunction &f, double x0, double x1, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            std::vector<double> x{x0, x1}, fx{f(x0), f(x1)};
            double fc, ft, q, xr, xc, xt;
            int k = 1, kt;

            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                xc = x.at(k);
                fc = fx.at(k);
                kt = k - 1;
                xt = x.at(kt);
                ft = fx.at(kt);
                while ((ft * fc >= 0) && (kt > 1))
                {
                    kt--;
                    xt = x.at(kt);
                    ft = fx.at(kt);
                }
                q = (fc - ft) / (xc - xt);
                xr = xc - fc / q;
                x.push_back(xr);
                fx.push_back(f(xr));
                err = parameter.wa * std::abs(xr - xc) + parameter.wf * std::abs(f(xr));
                k++;
            }

            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }

            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x.at(k), fx.at(k), err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "RegulaFalsi\t", term_flag);
        }

        RootScalarResult newton(UniVarFunction &f, UniVarFunction &df, double x, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double fx = f(x);
            double x0;
            int k = 0;

            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                x = x - fx / df(x);
                fx = f(x);
                err = parameter.wa * std::abs(x - x0) + parameter.wf * std::abs(fx);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }

            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Newton\t\t", term_flag);
        }

        RootScalarResult innewton(UniVarFunction &f, DerivFunction &df, double x, std::string method, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double fx = f(x);
            double x0;
            int k = 0;

            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                x = x - fx / df(x, parameter.h);
                fx = f(x);
                err = parameter.wa * std::abs(x - x0) + parameter.wf * std::abs(fx);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }

            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), method, term_flag);
        }

        RootScalarResult steffensen(UniVarFunction &f, double x, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double fx = f(x);
            double x0, q;
            int k = 0;
            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                q = (f(x + fx) - fx) / fx;
                x = x - fx / q;
                fx = f(x);
                err = parameter.wa * std::abs(x - x0) + parameter.wf * std::abs(fx);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }
            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, fx, err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Steffensen\t", term_flag);
        }

        RootScalarResult fixpoint(UniVarFunction &g, double x, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double x0;
            int k = 0;
            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                x = g(x);
                err = std::abs(x - x0);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }
            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, g(x), err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "FixPoint\t", term_flag);
        }

        RootScalarResult aitken(UniVarFunction &g, double x, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err = parameter.tol + 1;
            double x0, x1, x2, x3;
            int k = 0;
            while ((err > parameter.tol) && (k < parameter.maxit))
            {
                x0 = x;
                x1 = g(x);
                x2 = g(x1);
                x3 = x2 - x1;
                x = x2 - (x3 * x3) / (x2 - (2 * x1) + x);
                err = std::abs(x - x0);
                k++;
            }
            if ((err > parameter.tol) && (k == parameter.maxit))
            {
                term_flag = "Fail";
            }
            stopwatch.stop();
            return RootScalarResult(k, parameter.maxit, x, g(x), err, parameter.tol,
                                    stopwatch.get_elapsed_time(), "Aitken\t", term_flag);
        }

    } // end of namespace root::scalar
    namespace poly
    {
        using PolynomialFunc = complex_d_t(complex_d_t);

        struct param
        {
            double tol = std::numeric_limits<double>::epsilon();
            double reftol = std::sqrt(std::numeric_limits<double>::epsilon());
            int maxit = 100;
            int refmax = 1000;
            bool ref;
        };
        struct RootPolyResult
        {
            std::vector<int> numit;
            std::vector<int> numref;
            complex_d_t z;
            complex_d_t pz;
            std::vector<complex_d_t> pzm;
            std::vector<complex_d_t> qz;
            std::vector<complex_d_t> zm;

            RootPolyResult() {}

            RootPolyResult(const complex_d_t &pz,
                           std::vector<complex_d_t> &qz)
            {
                this->pz = pz;
                this->qz = qz;
            }
            RootPolyResult(const std::vector<int> &numit, const std::vector<int> &numref, const complex_d_t &z, const std::vector<complex_d_t> &pzm,
                           std::vector<complex_d_t> &zm)
            {
                this->numit = numit;
                this->numref = numref;
                this->z = z;
                this->pzm = pzm;
                this->zm = zm;
            }

            void print_horner()
            {
                std::cout << pz << "\t"
                          << "[";
                for (int i = 0; i < qz.size(); i++)
                {
                    std::cout << qz.at(i);
                    if (i <= qz.size() - 2)
                    {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
            }

            void print()
            {
                std::cout << std::setprecision(16)
                          << std::fixed
                          << std::scientific;
                for (int i = 0; i < zm.size(); i++)
                {
                    std::cout << zm.at(i).real() << "\t" << zm.at(i).imag()
                              << "\t" << std::abs(pzm.at(i))
                              << "\t" << numit.at(i) << "\n";
                }
            }

            void print_edit()
            {
                std::cout << std::setprecision(16)
                          << std::fixed
                          << std::scientific;
                for (int i = 0; i < zm.size(); i++)
                {
                    std::cout << zm.at(i).real() << "\t" << zm.at(i).imag()
                              << "\t" << std::abs(pzm.at(i))
                              << "\t" << numit.at(i) << "\t\t"
                              << numref.at(i) << "\n";
                }
            }
            ~RootPolyResult() {}
        };

        RootPolyResult horner(std::vector<complex_d_t> &pn, complex_d_t z)
        {
            timer stopwatch;
            stopwatch.start();

            std::vector<complex_d_t> bn = pn;
            for (int k = pn.size() - 2; k >= 0; k--)
            {
                bn.at(k) = pn.at(k) + bn.at(k + 1) * z;
            }
            std::vector<complex_d_t> qn = bn;
            qn.erase(qn.begin());
            stopwatch.stop();
            return RootPolyResult(bn.at(0), qn);
        }
        RootPolyResult newtonhorner(PolynomialFunc &px, std::vector<complex_d_t> pn, complex_d_t z, param &parameter)
        {
            timer stopwatch;
            stopwatch.start();
            std::string term_flag = "Success";
            double err_ref = parameter.tol * parameter.reftol;
            double err;
            RootPolyResult rpr;
            complex_d_t z0, z_ref, z_ref2;
            std::vector<complex_d_t> pnm = pn;
            complex_d_t pz, qz;
            std::vector<complex_d_t> qnm, qn, pzm, zm;
            std::vector<int> numit, numref;
            int k, k_ref, n = pn.size() - 1;
            for (int m = 0; m <= n - 1; m++)
            {
                k = 0;
                z = complex_d_t(z.real(), z.imag());
                err = parameter.tol + 1;
                if (m == n - 1)
                {
                    k++;
                    z = -pn.at(0) / pn.at(1);
                }
                else
                {
                    // pnm.resize(pn.size() - m - 1);
                    while ((err > parameter.tol) && (k < parameter.maxit))
                    {
                        k++;
                        z0 = z;
                        rpr = horner(pnm, z);
                        pz = rpr.pz;
                        qnm = rpr.qz;
                        rpr = horner(qnm, z);
                        qz = rpr.pz;
                        if (std::abs(qz) > std::numeric_limits<float>::epsilon())
                        {
                            z = z0 - pz / qz;
                            err = std::max(std::abs(z - z0), std::abs(pz));
                        }
                        else
                        {
                            term_flag = "Need to divide by small number.";
                            err = 0;
                            z = z0;
                        }
                    }
                }
                if (parameter.ref == true)
                {
                    k_ref = 0;
                    z_ref = z;
                    err = parameter.tol + 1;
                    while ((err > err_ref) && (k_ref < parameter.refmax))
                    {
                        k_ref++;
                        rpr = horner(pn, z_ref);
                        pz = rpr.pz;
                        qn = rpr.qz;
                        rpr = horner(qn, z_ref);
                        qz = rpr.pz;
                        if (std::abs(qz) > std::numeric_limits<double>::epsilon())
                        {
                            z_ref2 = z_ref - pz / qz;
                            err = std::max(std::abs(z_ref - z_ref2), std::abs(pz));
                            z_ref = z_ref2;
                        }
                        else
                        {
                            term_flag = "Need to divide by small number.";
                            err = 0;
                        }
                    }
                    z = z_ref;
                    numref.push_back(k_ref);
                }
                numit.push_back(k);
                zm.push_back(z);
                pzm.push_back(pz);
                rpr = horner(pnm, z);
                pnm = rpr.qz;
            }
            return RootPolyResult(numit, numref, z, pzm, zm);
        }

    } // end of namespace root::poly
}

#endif