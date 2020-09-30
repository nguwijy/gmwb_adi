#include "stdc++.h"
#include "characteristics.hpp"
#include "range.cpp"
#include "coefunc.hpp"
#include "terbounfunc.hpp"
#include "parabolicpde.cpp"
#include "lu.cpp"
#include "luforA2.cpp"
#include "spline.hpp"		//from http://kluge.in-chemnitz.de/opensource/spline/
#include "parabolicfdm.cpp"
#include "findfair.cpp"

int main() {
	clock_t t1, t2;

	t1 = clock();

	Srange srange(Sfrom, Sto);
	Vrange vrange(Vfrom, Vto);
	Rrange rrange(Rfrom, Rto);
	Trange trange(Tfrom, Tto);
	Trange brange(Bfrom, Bto);		//Trange used because we do equivalent size mesh

	double (*ss) (const double,const double, const double) = &ss2;
	double (*vv) (const double,const double, const double) = &vv2;
	double (*rr) (const double,const double, const double) = &rr2;
	double (*sv) (const double,const double, const double) = &sv2;
	double (*sr) (const double,const double, const double) = &sr2;
	double (*vr) (const double,const double, const double) = &vr2;
	double (*s) (const double,const double, const double) = &s1;
	double (*v) (const double,const double, const double) = &v1;
	double (*r) (const double,const double, const double, const double) = &r1;
	double (*zero) (const double,const double, const double) = &zero0;

	double (*ic) (const double, const double) = &Ic;
	bool (*extremev) (const double) = &ExtremeV;

	ParabolicPDE pde(srange, vrange, rrange, trange, brange, ss, vv, rr, sv, sr, vr, s, v, r, zero, ic, extremev);

	//choose between Douglas, CS, MCS, HV
	ParabolicFDM fdm(pde, Snum, Vnum, Rnum, Tnum, Bnum, theta, Douglas);

	long maxiter = 50;
	double tol = 1e-2;

	double starget = P0;
	double vtarget = Gamma;
	double rtarget = c1;

	

	FindFair fairVal(fdm, starget, vtarget, rtarget, Bnum, maxiter, tol);
	fairVal.start();

	double sWeight = fairVal.getsWeight();
	double vWeight = fairVal.getvWeight();
	double rWeight = fairVal.getrWeight();
	long sIndex = fairVal.getsIndex();
	long vIndex = fairVal.getvIndex();
	long rIndex = fairVal.getrIndex();
	double interestedV = fairVal.getInterestedV();

	fdm.start();

	std::vector<double> sarr = fdm.getSARR();
	std::vector<double> varr = fdm.getVARR();
	std::vector<double> rarr = fdm.getRARR();
	std::vector<double> barr = fdm.getBARR();
	std::vector<double> final = fdm.getres();

	std::ofstream output_file("./5Y_douglas_Anew.txt");

	t2 = clock();
	float diff ((float)t2 - (float)t1);
	diff = diff / CLOCKS_PER_SEC;

    std::cout << "The calculated fair fee is " << fairfee << "\n";

	/* output_file << diff << "\t" << (1 - sWeight) * sarr[sIndex] + sWeight * sarr[sIndex + 1] << "\t" << (1 - vWeight) * varr[vIndex] + vWeight * varr[vIndex + 1] << "\t" << (1 - rWeight) * rarr[rIndex] + rWeight * rarr[rIndex + 1] << "\t" << interestedV << "\t" << fairfee << "\n"; */
	output_file << std::setw(20) << "b: " << std::setw(20) << "r: " << std::setw(20) << "v: " << std::setw(20) << "s: " << std::setw(20) << "value: " << "\n";

	for (long bb = 0; bb< (Bnum + 1); bb++) {
		for (long i = 0; i < (Rnum + 1); i++) {
			for (long j = 0; j < (Vnum + 1); j++) {
				for (long k = 0; k < (Snum + 1); k++) {
					output_file << std::setw(20) << barr[bb] << std::setw(20) << rarr[i] << std::setw(20) << varr[j] << std::setw(20) << sarr[k] << std::setw(20) << final[bb * (Rnum + 1) * (Vnum + 1) * (Snum + 1) + i * (Vnum + 1) * (Snum + 1) + j * (Snum + 1) + k] << "\n";
				}
			}
		}
	}

	return 0;
}
