#ifndef LU_CPP
#define LU_CPP


#include "lu.hpp"

LU::LU () {
	//empty ok
}

LU::LU (const LU& s2) {
	//empty ok
}

LU::LU(const std::vector<double>& a, const std::vector<double>& b, 
	const std::vector<double>& c, const std::vector<double>& d,
	const std::vector<double>& f) {
	aa = a;
	bb = b;
	cc = c;
	dd = d;
	ff = f;
}

LU::~LU () {
	//empty ok
}

/*
LU& LU::operator = (const LU& i2) {
	//empty ok
}*/

// result 
std::vector<double> LU::result() const
{ // Code to actually create the solution to the tridiagonal system

	size_t N =  aa.size();

	std::vector<double> alpha(N);
	std::vector<double> beta(N);
	std::vector<double> phi(N);

	alpha[0] = 0;
	alpha[1] = bb[1];

	beta[0] = cc[0];

	phi[0] = dd[0] / beta[0];

	beta[1] = cc[1] - alpha[1]*phi[0];
	phi[1] = dd[1] / beta[1];

	for (long x = 2; x < N; x++)
	{
		alpha[x] = bb[x] - aa[x]*phi[x-2];
		beta[x] = cc[x] - alpha[x]*phi[x-1];
		phi[x] = dd[x] / beta[x];
	}

	std::vector<double> z(N);
	z[0] = ff[0] / beta[0];
	z[1] = (ff[1] - z[0]*alpha[1]) / beta[1];

	for (long x = 2; x < N; x++)
	{
		z[x] = (ff[x] - z[x-1]*alpha[x] - z[x-2]*aa[x]) / beta[x];
	}


	std::vector<double> u(N); 
	u[N - 1] = z[N - 1];

	for (long x = N-2; x >= 0; x--)
	{
		u[x] = z[x] - phi[x]*u[x+1];
	}
	
	return u;

}

#endif