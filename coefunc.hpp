#ifndef COEF_HPP
#define COEF_HPP

double ss2(const double s, const double v, const double r) {return (0.5 * s * s * v);}
double vv2(const double s, const double v, const double r) {return (0.5 * sig1 * sig1 * v);}
double rr2(const double s, const double v, const double r) {return (0.5 * sig2 * sig2);}
double sv2(const double s, const double v, const double r) {return (rho12 * sig1 * s * v);}
double sr2(const double s, const double v, const double r) {return (rho13 * sig2 * s * pow(v, 0.5));}
double vr2(const double s, const double v, const double r) {return (rho23 * sig1 * sig2 * pow(v, 0.5));}
double s1(const double s, const double v, const double r) {return ((r - fairfee) * s);}
double v1(const double s, const double v, const double r) {return (kappa * (Gamma - v));}
//r1 is time-dependent
//double r1(const double s, const double v, const double r, const double t) {return (a * (c1 - c2 * exp(-c3 * (t))) - r);}
double r1(const double s, const double v, const double r, const double t) {return (rkappa * (rGamma + (sig2 * sig2 / (2.0 * rkappa * rkappa)) * (1.0 - exp(- 2.0 * rkappa * (Tto - t))) - r));}
double zero0(const double s, const double v, const double r) {return (-r);}

#endif