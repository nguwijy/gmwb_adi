#ifndef PDE_CPP
#define PDE_CPP

#include "parabolicpde.hpp"

ParabolicPDE::ParabolicPDE() {
	saxis = Srange();
	vaxis = Vrange();
	raxis = Rrange();
	taxis = Trange();
	baxis = Trange();

	//ss, vv, rr, sv, sr, vr, s, v, zero = double (*) (const double, const double, const double) ();
	//r = double (*) (const double, const double, const double, const double) ();
	//ic = double (*) (const double, const double) ();
	//extremev = bool (*) (const double) ();
}

ParabolicPDE::ParabolicPDE(const Srange& SRANGE, const Vrange& VRANGE, const Rrange& RRANGE, const Trange& TRANGE, const Trange& BRANGE,
			double (*SS) (const double, const double, const double), double (*VV) (const double, const double, const double), double (*RR) (const double, const double, const double),
			double (*SV) (const double, const double, const double), double (*SR) (const double, const double, const double), double (*VR) (const double, const double, const double), double (*S) (const double, const double, const double),
			double (*V) (const double, const double, const double), double (*R) (const double, const double, const double, const double), double (*ZERO) (const double, const double, const double), double (*IC) (const double, const double),
			bool (*EXTREMEV) (const double)) {
	saxis = SRANGE;
	vaxis = VRANGE;
	raxis = RRANGE;
	taxis = TRANGE;
	baxis = BRANGE;

	ss = SS;
	vv = VV;
	rr = RR;
	sv = SV;
	sr = SR;
	vr = VR;
	s = S;
	v = V;
	r = R;
	zero = ZERO;
	ic = IC;

	extremev = EXTREMEV;
}

ParabolicPDE::~ParabolicPDE() {;}

void ParabolicPDE::setS(Srange& input) {
	saxis = input;
}

void ParabolicPDE::setV(Vrange& input) {
	vaxis = input;
}

void ParabolicPDE::setR(Rrange& input) {
	raxis = input;
}

void ParabolicPDE::setT(Trange& input) {
	taxis = input;
}

void ParabolicPDE::setB(Trange& input) {
	baxis = input;
}


Srange ParabolicPDE::getS() {
	return saxis;
}

Vrange ParabolicPDE::getV() {
	return vaxis;
}

Rrange ParabolicPDE::getR() {
	return raxis;
}

Trange ParabolicPDE::getT() {
	return taxis;
}

Trange ParabolicPDE::getB() {
	return baxis;
}

double ParabolicPDE::SS(const double svar, const double vvar, const double rvar) {
	return (*ss) (svar, vvar, rvar);
}

double ParabolicPDE::VV(const double svar, const double vvar, const double rvar) {
	return (*vv) (svar, vvar, rvar);
}

double ParabolicPDE::RR(const double svar, const double vvar, const double rvar) {
	return (*rr) (svar, vvar, rvar);
}

double ParabolicPDE::SV(const double svar, const double vvar, const double rvar) {
	return (*sv) (svar, vvar, rvar);
}

double ParabolicPDE::SR(const double svar, const double vvar, const double rvar) {
	return (*sr) (svar, vvar, rvar);
}

double ParabolicPDE::VR(const double svar, const double vvar, const double rvar) {
	return (*vr) (svar, vvar, rvar);
}

double ParabolicPDE::S(const double svar, const double vvar, const double rvar) {
	return (*s) (svar, vvar, rvar);
}

double ParabolicPDE::V(const double svar, const double vvar, const double rvar) {
	return (*v) (svar, vvar, rvar);
}

double ParabolicPDE::R(const double svar, const double vvar, const double rvar, const double t) {
	return (*r) (svar, vvar, rvar, t);
}

double ParabolicPDE::ZERO(const double svar, const double vvar, const double rvar) {
	return (*zero) (svar, vvar, rvar);
}

double ParabolicPDE::IC(const double svar, const double B) {
	return (*ic) (svar, B);
}

bool ParabolicPDE::EXTREMEV(const double vvar) {
	return (*extremev) (vvar);
}

#endif