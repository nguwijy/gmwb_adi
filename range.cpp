#ifndef Range_CPP
#define Range_CPP

#include "range.hpp"

Range::Range() {
	double def = 0.0;
	lo = def;
	hi = def;
}

Range::Range(const double low, const double high) {
	lo = low;
	hi = high;
}

Range::Range(const Range& ran2) {
	lo = ran2.lo;
	hi = ran2.hi;
}

Range::~Range() {

}

void Range::low(const double& t1) {
	lo = t1;
}

void Range::high(const double& t1) {
	hi = t1;
}

double Range::low() const {
	return lo;
}

double Range::high() const {
	return hi;
}

Range& Range::operator = (const Range& ran2) {
	lo = ran2.lo;
	hi = ran2.hi;

	return *this;
}

Srange::Srange() {
	double def = 0.0;
	lo = def;
	hi = def;
}

Srange::Srange(const double low, const double high) {
	lo = low;
	hi = high;
}

std::vector<double> Srange::mesh(const int nSteps, const double sImp) const {
	std::vector<double> meshes(nSteps+1);
	double sleft = 0.8 * sImp;
	double sright = 1.2 * sImp;
	double d1 = sImp / 20.0;
	double smin = asinh(-sleft / d1);
	double sint = (sright - sleft) / d1;
	double smax = sint + asinh((hi - sright) / d1);
	double del = (smax - smin) / nSteps;
	for (int i = 0; i < (nSteps + 1); i++) {
		meshes[i] = smin + i * del;
		if (meshes[i] < 0) {meshes[i] = sleft + d1 * sinh(meshes[i]);}
		else if (meshes[i] < sint) {meshes[i] = sleft + d1 * meshes[i];}
		else {meshes[i] = sright + d1 * sinh(meshes[i] - sint);}
	}
	return meshes;
}

Vrange::Vrange() {
	double def = 0.0;
	lo = def;
	hi = def;
}

Vrange::Vrange(const double low, const double high) {
	lo = low;
	hi = high;
}

std::vector<double> Vrange::mesh(const int nSteps) const {
	std::vector<double> meshes(nSteps+1);
	double d2 = hi / 500.0;
	double del = asinh(hi / d2) / nSteps;
	for (int i = 0; i < (nSteps + 1); i++) {
		meshes[i] = i * del;
		meshes[i] = d2 * sinh(meshes[i]);
	}
	return meshes;
}

Rrange::Rrange() {
	double def = 0.0;
	lo = def;
	hi = def;
}

Rrange::Rrange(const double low, const double high) {
	lo = low;
	hi = high;
}

std::vector<double> Rrange::mesh(const int nSteps, const double rImp) const {
	std::vector<double> meshes(nSteps+1);
	double d3 = hi / 400.0;
	double del = (asinh((hi - rImp) / d3) - asinh((lo - rImp) / d3)) / nSteps;
	for (int i = 0; i < (nSteps + 1); i++) {
		meshes[i] = asinh((lo-rImp) / d3) + i * del;
		meshes[i] = rImp + d3 * sinh(meshes[i]);
	}
	return meshes;
}

Trange::Trange() {
	double def = 0.0;
	lo = def;
	hi = def;
}

Trange::Trange(const double low, const double high) {
	lo = low;
	hi = high;
}

std::vector<double> Trange::mesh(const int nSteps) const {
	std::vector<double> meshes(nSteps+1);
	double del = (hi - lo) / nSteps;
	for (int i = 0; i < (nSteps + 1); i++) {
		meshes[i] = lo + i * del;
	}
	return meshes;
}

#endif