#ifndef FDM_CPP
#define FDM_CPP

#include "parabolicfdm.hpp"

void ParabolicFDM::init() {
	SARR = std::vector<double> (snum + 1);
	VARR = std::vector<double> (vnum + 1);
	RARR = std::vector<double> (rnum + 1);
	TARR = std::vector<double> (tnum + 1);
	BARR = std::vector<double> (bnum + 1);
	tmp = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	res = std::vector<double> ((bnum + 1) * (snum + 1) * (vnum + 1) * (rnum + 1));


	SARR = pde.getS().mesh(snum, P0);
	VARR = pde.getV().mesh(vnum);
	RARR = pde.getR().mesh(rnum, Rint);
	BARR = pde.getB().mesh(bnum);
	TARR = pde.getT().mesh(tnum);

	h = (pde.getT().high() - pde.getT().low()) / tnum;

}

ParabolicFDM::ParabolicFDM() {
	pde = ParabolicPDE();

	current = pde.getT().low();
}

ParabolicFDM::ParabolicFDM(const ParabolicPDE& context, unsigned long Sintervals, unsigned long Vintervals, 
				unsigned long Rintervals, unsigned long Tintervals, unsigned long Bintervals, const double& theta, schemeType type) {
	pde = context;

	typ = type;

	current = pde.getT().low();

	snum = Sintervals;
	vnum = Vintervals;
	rnum = Rintervals;
	tnum = Tintervals;
	bnum = Bintervals;
	the = theta;

	init();
}

void ParabolicFDM::start() {
	for (long bb = 0; bb < (bnum + 1); bb++) {
		for (long i = 0; i < (rnum + 1); i++) {
			for (long j = 0; j < (vnum + 1); j++) {
				for (long k = 0; k < (snum + 1); k++) {
					res[bb * (rnum + 1) * (vnum + 1) * (snum + 1) + i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = pde.IC(SARR[k], BARR[bb]);		//****** need to change 0 to B
				}
			}
		}
	}

	A1next = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A1last = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A1now = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A2next = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A2last = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A2lastlast = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A2now = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3next_now = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3last_now = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3now_now = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3next_next = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3last_next = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	A3now_next = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));

	snext = std::vector<double> (snum - 1);
	slast = std::vector<double> (snum - 1);
	snow = std::vector<double> (snum - 1);
	vnext = std::vector<double> (vnum - 1);
	vlast = std::vector<double> (vnum - 1);
	vnow = std::vector<double> (vnum - 1);
	rnext = std::vector<double> (rnum - 1);
	rlast = std::vector<double> (rnum - 1);
	rnow = std::vector<double> (rnum - 1);

	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 1; k < snum; k++) {
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k + 1] - SARR[k]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k])) + pde.S(SARR[k], VARR[j], RARR[i]) * (SARR[k] - SARR[k - 1]) / ((SARR[k + 1] - SARR[k]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k])));
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k] - SARR[k - 1]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k])) - pde.S(SARR[k], VARR[j], RARR[i]) * (SARR[k + 1] - SARR[k ]) / ((SARR[k] - SARR[k - 1]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k])));
				A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k + 1] - SARR[k]) * (SARR[k] - SARR[k - 1])) + pde.S(SARR[k], VARR[j], RARR[i]) * (SARR[k + 1] - SARR[k] - SARR[k] + SARR[k - 1]) / ((SARR[k + 1] - SARR[k]) * (SARR[k] - SARR[k - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			}
			//first derivative added later
			long k = 0;
			A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])));
			A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])));
			A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.SS(SARR[k], VARR[j], RARR[i]) * 2.0 / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);

			//not used
			k = snum;
			A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
		}
	}

	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			for (long j = 1; j < vnum; j++) {
				if (pde.EXTREMEV(VARR[j])) {
					A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])));
					A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])) - pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j - 1] - VARR[j - 2] + VARR[j] - VARR[j - 1])/((VARR[j - 1] - VARR[j - 2]) * (VARR[j] - VARR[j - 1])));
					A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j] - VARR[j-1]) / ((VARR[j - 1] - VARR[j - 2]) * (VARR[j - 1] - VARR[j - 2] + VARR[j] - VARR[j - 1])));
					A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1])) + pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j - 1] - VARR[j - 2] + 2.0 * (VARR[j] - VARR[j - 1]))/((VARR[j] - VARR[j - 1]) * (VARR[j - 1] - VARR[j - 2] + VARR[j] - VARR[j - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				} else {
					A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])) + pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j] - VARR[j - 1])/((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])));
					A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])) - pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j + 1] - VARR[j])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j])));
					A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
					A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1])) + pde.V(SARR[k], VARR[j], RARR[i]) * (VARR[j + 1] - VARR[j] - VARR[j] + VARR[j - 1])/((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				}
			}
			//not used
			long j = 0;
			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;

			//first derivative added later
			j = vnum;
			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])));
			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])));
			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * 0;
			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.VV(SARR[k], VARR[j], RARR[i]) * 2.0 / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
		}
	}

	for (long k = 1; k < snum; k++) {
		snext[k - 1] = (SARR[k] - SARR[k - 1]) / ((SARR[k + 1] - SARR[k]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k]));
		slast[k - 1] = (SARR[k] - SARR[k + 1]) / ((SARR[k] - SARR[k - 1]) * (SARR[k] - SARR[k - 1] + SARR[k + 1] - SARR[k]));
		snow[k - 1] = (SARR[k + 1] - SARR[k] - SARR[k] + SARR[k - 1]) / ((SARR[k] - SARR[k - 1]) * (SARR[k + 1] - SARR[k]));
	}

	for (long j = 1; j < vnum; j++) {
		vnext[j - 1] = (VARR[j] - VARR[j - 1]) / ((VARR[j + 1] - VARR[j]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j]));
		vlast[j - 1] = (VARR[j] - VARR[j + 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j + 1] - VARR[j]));
		vnow[j - 1] = (VARR[j + 1] - VARR[j] - VARR[j] + VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j + 1] - VARR[j]));
	}

	for (long i = 1; i < rnum; i++) {
		rnext[i - 1] = (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i]));
		rlast[i - 1] = (RARR[i] - RARR[i + 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i]));
		rnow[i - 1] = (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i + 1] - RARR[i]));
	}

	while (!finished()) {
		for (long bb = 0; bb < (bnum + 1); bb++) {
			for (long i = 0; i < (rnum + 1); i++) {
				for (long j = 0; j < (vnum + 1); j++) {
					for (long k = 0; k < (snum + 1); k++) {
						tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = res[bb * (rnum + 1) * (vnum + 1) * (snum + 1) + i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
					}
				}
			}

			if (typ == Douglas) {
				douglas();
			}

			if (typ == CS) {
				cs();
			}

			if (typ == MCS) {
				mcs();
			}

			if (typ == HV) {
				hv();
			}

			for (long i = 0; i < (rnum + 1); i++) {
				for (long j = 0; j < (vnum + 1); j++) {
					for (long k = 0; k < (snum + 1); k++) {
						res[bb * (rnum + 1) * (vnum + 1) * (snum + 1) + i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
					}
				}
			}
		}

		current += h;

		if (((trunc(current + 1e-5) + 1e-5) >= current) && ((trunc(current + 1e-5) - 1e-5) <= current) && (!finished())) {
			opti();
		}
	}
}

void ParabolicFDM::reset() {
	SARR = std::vector<double> (snum + 1);
	VARR = std::vector<double> (vnum + 1);
	RARR = std::vector<double> (rnum + 1);
	TARR = std::vector<double> (tnum + 1);
	BARR = std::vector<double> (bnum + 1);
	tmp = std::vector<double> ((snum + 1) * (vnum + 1) * (rnum + 1));
	res = std::vector<double> ((bnum + 1) * (snum + 1) * (vnum + 1) * (rnum + 1));


	SARR = pde.getS().mesh(snum, P0);
	VARR = pde.getV().mesh(vnum);
	RARR = pde.getR().mesh(rnum, Rint);
	BARR = pde.getB().mesh(bnum);
	TARR = pde.getT().mesh(tnum);

	h = (pde.getT().high() - pde.getT().low()) / tnum;
	
	current = pde.getT().low();
}

void ParabolicFDM::douglas() {
	//setting up
	double t = current;

	std::vector<double> Y0 = tmp;
	std::vector<double> Y1((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3((snum + 1) * (vnum + 1) * (rnum + 1));

	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			for (long i = 1; i < rnum; i++) {
				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			}
			//first derivative added later
			long i = 0;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);

			//first derivative added later
			i = rnum;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));;
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
		}
	}

	//Douglas first step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]);

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
					A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
					A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
					A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)];
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + h * (pde.V(SARR[k], VARR[j], RARR[i])*(
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
						RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];
				} else if (j < vnum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k];
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
            		
            		sr = false;
            		vr = false;
				} else {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				}

				if (k == snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]);
				}

				if (sr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]);
				}

				if (vr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]);
				}
			}
		}
	}

	//Douglas second step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//Douglas third step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//Douglas forth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	tmp = Y3;
}

void ParabolicFDM::cs() {
	//setting up
	double t = current;

	std::vector<double> Y0 = tmp;
	std::vector<double> Y1((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y0squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y1squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3squiggle((snum + 1) * (vnum + 1) * (rnum + 1));

	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			for (long i = 1; i < rnum; i++) {
				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			}
			//first derivative added later
			long i = 0;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);

			//first derivative added later
			i = rnum;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));;
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
		}
	}

	//CS first step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]);

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
					A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
					A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
					A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)];
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + h * (pde.V(SARR[k], VARR[j], RARR[i])*(
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
						RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];
				} else if (j < vnum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k];
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
            		
            		sr = false;
            		vr = false;
				} else {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				}

				if (k == snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]);
				}

				if (sr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]);
				}

				if (vr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]);
				}
			}
		}
	}

	//CS second step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//CS third step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//CS forth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	//CS fifth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];

				bool sv = true;
				bool sr = true;
				bool vr = true;

				if ((k == 0) || (k == snum)) {
					sv = false;
					sr = false;
				} 

				if ((j == 0) || (j == vnum)) {
					sv = false;
					vr = false;
				}

				if ((i == 0) || (i == rnum)) {
            		sr = false;
            		vr = false;
				}

				if (sv) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)]) +
						snext[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)]) +
						snext[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)]) +
						slast[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)]) +
						slast[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)]) +
						slast[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)]) +
						snow[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)]) +
						snow[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)]) +
						snow[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]));
				}

				if (sr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						snext[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						snext[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						slast[k - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						slast[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						slast[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						snow[k - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) +
						snow[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) +
						snow[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]));
				}

				if (vr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vnext[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vnext[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vlast[j - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vlast[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vlast[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vnow[j - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) +
						vnow[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) +
						vnow[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]));
				}
			}
		}
	}

	//CS sixth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//CS seventh step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//CS eighth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	tmp = Y3squiggle;
}

void ParabolicFDM::mcs() {
	//setting up
	double t = current;

	std::vector<double> Y0 = tmp;
	std::vector<double> Y1((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y0hat((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y0squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y1squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3squiggle((snum + 1) * (vnum + 1) * (rnum + 1));

	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			for (long i = 1; i < rnum; i++) {
				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			}
			//first derivative added later
			long i = 0;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);

			//first derivative added later
			i = rnum;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));;
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
		}
	}

	//MCS first step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]);

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
					A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
					A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
					A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)];
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + h * (pde.V(SARR[k], VARR[j], RARR[i])*(
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
						RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];
				} else if (j < vnum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k];
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
            		
            		sr = false;
            		vr = false;
				} else {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				}

				if (k == snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]);
				}

				if (sr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]);
				}

				if (vr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]);
				}
			}
		}
	}

	//MCS second step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//MCS third step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//MCS forth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	//MCS fifth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];

				bool sv = true;
				bool sr = true;
				bool vr = true;

				if ((k == 0) || (k == snum)) {
					sv = false;
					sr = false;
				} 

				if ((j == 0) || (j == vnum)) {
					sv = false;
					vr = false;
				}

				if ((i == 0) || (i == rnum)) {
            		sr = false;
            		vr = false;
				}

				if (sv) {
					Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + theta * pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)]) +
						snext[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)]) +
						snext[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)]) +
						slast[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)]) +
						slast[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)]) +
						slast[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)]) +
						snow[k - 1] * vnext[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)]) +
						snow[k - 1] * vlast[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)]) +
						snow[k - 1] * vnow[j - 1] * (Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)] - tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]));
				}

				if (sr) {
					Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + theta * pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						snext[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						snext[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						slast[k - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						slast[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						slast[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]) +
						snow[k - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) +
						snow[k - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) +
						snow[k - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]));
				}

				if (vr) {
					Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + theta * pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vnext[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vnext[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
						vlast[j - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vlast[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vlast[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
						vnow[j - 1] * rnext[i - 1] * (Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) +
						vnow[j - 1] * rlast[i - 1] * (Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) +
						vnow[j - 1] * rnow[i - 1] * (Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] - tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]));
				}
			}
		}
	}

	//MCS sixth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0hat[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];

				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) + (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)])));

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) * 
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]));
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) * h * (pde.V(SARR[k], VARR[j], RARR[i]) * (
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k])) -
						RARR[i] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k])) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]));

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]));
				} else if (j < vnum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]));
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])));

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));
            		
            		sr = false;
            		vr = false;
				} else {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));
				}

				if (k == snum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) * pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						((snext[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]) - 
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]));
				}

				if (sr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) * pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						((snext[k - 1] * rnext[i- 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) - 
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]));
				}

				if (vr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (0.5 - theta) * pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						((vnext[j - 1] * rnext[i- 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) - 
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]));
				}
			}
		}
	}

	//MCS seventh step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//MCS eighth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//MCS ninth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	tmp = Y3squiggle;
}

void ParabolicFDM::hv() {
	//setting up
	double t = current;

	std::vector<double> Y0 = tmp;
	std::vector<double> Y1((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y0squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y1squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y2squiggle((snum + 1) * (vnum + 1) * (rnum + 1));
	std::vector<double> Y3squiggle((snum + 1) * (vnum + 1) * (rnum + 1));

	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			for (long i = 1; i < rnum; i++) {
				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], current) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i] - RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])) - pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i + 1] - RARR[i])));
				A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.R(SARR[k], VARR[j], RARR[i], (current + h)) * (RARR[i + 1] - RARR[i] - RARR[i] + RARR[i - 1]) / ((RARR[i + 1] - RARR[i]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			}
			//first derivative added later
			long i = 0;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])));
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);

			//first derivative added later
			i = rnum;
			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
			A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));
			A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])));;
			A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = h * (- pde.RR(SARR[k], VARR[j], RARR[i]) * 2.0 / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) + pde.ZERO(SARR[k], VARR[j], RARR[i]) / 3.0);
		}
	}

	//HV first step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]);

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
					A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
					A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
					A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)];
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + h * (pde.V(SARR[k], VARR[j], RARR[i])*(
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
						RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k];
				} else if (j < vnum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
            			A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k];
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
            			((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
            		
            		sr = false;
            		vr = false;
				} else {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				}

				if (k == snum) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]);
				}

				if (sr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]);
				}

				if (vr) {
					Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]);
				}
			}
		}
	}

	//HV second step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//HV third step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//HV forth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	//HV fifth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			for (long k = 0; k < (snum + 1); k++) {
				Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];

				bool sv = true;
				bool sr = true;
				bool vr = true;

				//s = 0, du / ds = 0
				if (k == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
						((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) + (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)])));

					sv = false;
					sr = false;
				} else if (k < snum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * 
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]) +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]));
				}

				//v = 0, all v - related derivatives vanish
				if (j == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * h * (pde.V(SARR[k], VARR[j], RARR[i]) * (
						(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) + 
						((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) - 
						(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k])) -
						RARR[i] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) / 3.0);

					sv = false;
					vr = false;
				} else if (j == vnum) {
					//v = max, du / dv = 0
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k])) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]));

            		sv = false;
            		vr = false;
				} else if ((j < vnum) && (pde.EXTREMEV(VARR[j]))) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]));
				} else if (j < vnum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * (Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* (Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] - tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]));
				}

				//r = min, du / dr = 0
				if (i == 0) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])));

            		sr = false;
            		vr = false;
				} else if (i == rnum) {
					//r = max, du / dr = 0
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           				((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1])) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));
            		
            		sr = false;
            		vr = false;
				} else {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 *
           				((A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]) - 
           				(A3now_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));
				}

				if (k == snum) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));

					sv = false;
					sr = false;
				}

				if (sv) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.SV(SARR[k], VARR[j], RARR[i]) * h *
						((snext[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * Y3[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]) - 
						(snext[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 1)] +
						snext[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 1)] +
						slast[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k - 1)] +
						slast[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k - 1)] +
						snow[k - 1] * vnext[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vlast[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + (k + 0)] +
						snow[k - 1] * vnow[j - 1] * tmp[i * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + (k + 0)]));
				}

				if (sr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.SR(SARR[k], VARR[j], RARR[i]) * h *
						((snext[k - 1] * rnext[i- 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]) - 
						(snext[k - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						snext[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						slast[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						slast[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)] +
						snow[k - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)] +
						snow[k - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 0)]));
				}

				if (vr) {
					Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 0.5 * pde.VR(SARR[k], VARR[j], RARR[i]) * h *
						((vnext[j - 1] * rnext[i- 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * Y3[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]) - 
						(vnext[j - 1] * rnext[i- 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vnext[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
						vlast[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vlast[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
						vnow[j - 1] * rnext[i - 1] * tmp[(i + 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rlast[i - 1] * tmp[(i - 1) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k] +
						vnow[j - 1] * rnow[i - 1] * tmp[(i + 0) * (vnum + 1) * (snum + 1) + (j + 0) * (snum + 1) + k]));
				}
			}
		}
	}

	//MCS sixth step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (vnum + 1); j++) {
			std::vector<double> temp(snum + 1);
			std::vector<double> tempnext(snum + 1);
			std::vector<double> templast(snum + 1);
			std::vector<double> tempnow(snum + 1);
			std::vector<double> templastlast(snum + 1, 0);

			for (long k = 0; k < (snum + 1); k++) {
				if ((k != 0) && (k != snum)) {
					temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
						(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
						A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
						A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k - 1)]);
				}

				tempnext[k] = - theta * A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[k] = - theta * A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[k] = 1.0 - theta * A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//s = max
			long k = snum;
			temp[k] = SARR[k] * exp(- fairfee * (current + h));
			tempnext[k] = 0;
			templast[k] = 0;
			tempnow[k] = 1;

			//s = min, du / ds = 0
			k = 0;
			temp[k] = Y0squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
				(A1now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
				A1next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)] +
				A1last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
				((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * ((SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]+(SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + (k + 1)]));
			tempnext[k] = tempnext[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]+SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k] - SARR[k + 1] + SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k]));
			tempnow[k] = tempnow[k] + templast[k] * ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]) / (SARR[k + 1] - SARR[k])) * (SARR[k + 1] - SARR[k]) / ((SARR[k + 1] - SARR[k]) * (SARR[k + 1] - SARR[k] + SARR[k + 1] - SARR[k]));
      		templast[k] = 0;

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long k = 0; k < (snum + 1); k++) {
      			Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[k];
      		}
		}
	}

	//MCS seventh step
	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(vnum + 1);
			std::vector<double> tempnext(vnum + 1);
			std::vector<double> templast(vnum + 1);
			std::vector<double> tempnow(vnum + 1);
			std::vector<double> templastlast(vnum + 1);

			for (long j = 1; j < (vnum + 1); j++) {
				if (pde.EXTREMEV(VARR[j])) {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);
				} else {
					temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]* Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]);
				}

				tempnext[j] = - theta * A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[j] = - theta * A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templastlast[j] = - theta * A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[j] = 1.0 - theta * A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//v = 0, all v - derivatives vanish
			long j = 0;
			tempnow[j] = 1.0 - theta * h * ((- 2.0 * (VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]) - RARR[i] / 3.0);
    		tempnext[j] = - theta * h * ((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);
    		templast[j] = 0;
    		templastlast[j] = 0;
    		temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * h * (pde.V(SARR[k], VARR[j], RARR[i])*(
				(- 2.0 *(VARR[j + 1] - VARR[j]) - (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + 
				((VARR[j + 1] - VARR[j]) + (VARR[j + 2] - VARR[j + 1])) / ((VARR[j + 1] - VARR[j]) * (VARR[j + 2] - VARR[j + 1])) * Y3[i * (vnum + 1) * (snum + 1) + (j + 1) * (snum + 1) + k] - 
				(VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * Y3[i * (vnum + 1) * (snum + 1) + (j + 2) * (snum + 1) + k]) -
				RARR[i] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] / 3.0);
    		//negative * negative = positive
    		double A2nextnext = theta * h * (VARR[j + 1] - VARR[j]) / ((VARR[j + 2] - VARR[j + 1]) * (VARR[j + 1] - VARR[j] + VARR[j + 2] - VARR[j + 1])) * pde.V(SARR[k], VARR[j], RARR[i]);

			//v = max, du / dv = 0
			j = vnum;
			templast[j] = templast[j] - tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1])/((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]));
			tempnow[j] = tempnow[j] + tempnext[j] * (VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1]));
			tempnext[j] = 0;
			temp[j] = Y1squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
            			(A2now[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
            			A2next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] *
            			(VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / (VARR[j] - VARR[j - 1] + 2.0 * (VARR[j] - VARR[j - 1])) * ((VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (VARR[j] - VARR[j - 1]) / ((VARR[j] - VARR[j - 1]) * (VARR[j] - VARR[j - 1] + VARR[j] - VARR[j - 1])) * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k]) +
            			A2last[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + (j - 1) * (snum + 1) + k] +
            			A2lastlast[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + (j - 2) * (snum + 1) + k]);

      		LUforA2 restemp(templastlast, templast, tempnow, tempnext, temp, A2nextnext);
      		temp = restemp.result();

      		for (long j = 0; j < (vnum + 1); j++) {
      			if (k == snum) {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[j];
      			}
      		}
		}
	}

	//MCS eighth step
	for (long j = 0; j < (vnum + 1); j++) {
		for (long k = 0; k < (snum + 1); k++) {
			std::vector<double> temp(rnum + 1);
			std::vector<double> tempnext(rnum + 1);
			std::vector<double> templast(rnum + 1);
			std::vector<double> tempnow(rnum + 1);
			std::vector<double> templastlast(rnum + 1, 0);

			for (long i = 0; i < (rnum + 1); i++) {
				if ((i != 0) && (i != rnum)) {
					temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           				(A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           				A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);
				}

				tempnext[i] = - theta * A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				templast[i] = - theta * A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
				tempnow[i] = 1.0 - theta * A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
			}

			//r = min, du / dr = 0
			long i = 0;
			tempnext[i] = tempnext[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i]));
    		tempnow[i] = tempnow[i] + templast[i] * ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]));
    		templast[i] = 0;
    		temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta *
           		(A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i]) / (RARR[i + 1] - RARR[i])) * ((RARR[i + 1] - RARR[i] - RARR[i + 1] + RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] + (RARR[i + 1] - RARR[i]) / ((RARR[i + 1] - RARR[i]) * (RARR[i + 1] - RARR[i] + RARR[i + 1] - RARR[i])) * Y3[(i + 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]));

			//r = max, du / dr = 0
			i = rnum;
			templast[i] = templast[i] + tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]) / ((RARR[i]-RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]));
			tempnow[i] = tempnow[i] - tempnext[i] * (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1]));
			tempnext[i] = 0;
			temp[i] = Y2squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - theta * 
           		(A3now_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] +
           		A3next_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * 
           		((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1]) / (RARR[i] - RARR[i - 1]) * ((RARR[i] - RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1] + RARR[i] - RARR[i - 1])) * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] - (RARR[i] - RARR[i - 1] - RARR[i] + RARR[i - 1]) / ((RARR[i] - RARR[i - 1]) * (RARR[i] - RARR[i - 1])) * Y3[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k])) +
           		A3last_next[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] * Y3[(i - 1) * (vnum + 1) * (snum + 1) + j * (snum + 1) + k]);

      		LU restemp(templastlast, templast, tempnow, tempnext, temp);
      		temp = restemp.result();

      		for (long i = 0; i < (rnum + 1); i++) {
      			if (k == snum) {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = SARR[k] * exp(- fairfee * (current + h));
      			} else {
      				Y3squiggle[i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = temp[i];
      			}
      		}
		}
	}

	tmp = Y3squiggle;
}

void ParabolicFDM::opti() {
	std::vector<double> opttemp = res;

	for (long bb = 1; bb < (bnum + 1); bb++) {
		for (long i = 0; i < (rnum + 1); i++) {
			for (long j = 0; j < (vnum + 1); j++) {
				std::vector<double> bestvalue((snum + 1), 0);
				//for (long bbsub = 0; bbsub < (bb + 1); bbsub++) {
					//double withdrawn = BARR[bb] - BARR[bbsub];
					long bbsub = bb - 1;
					double withdrawn = G;
					tk::spline sp;
					std::vector<double> sp_x = SARR;
					std::vector<double> sp_y(snum + 1);
					for (long k = 0; k < (snum + 1); k++) {
						sp_x[k] = (sp_x[k] > withdrawn)? (sp_x[k] - withdrawn): SARR[0]; 
						sp_y[k] = opttemp[bbsub * (rnum + 1) * (vnum + 1) * (snum + 1) + i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k];
					}
					sp.set_points(SARR, sp_y);
					for (long k = 0; k < (snum + 1); k++) {
						double valuenow = sp(sp_x[k]) + withdrawn - withdrawalPen * ((withdrawn > G)? (withdrawn - G): 0);
						bestvalue[k] = (bestvalue[k] > valuenow)? bestvalue[k]: valuenow;
					}
				//}
				for (long k = 0; k < (snum + 1); k++) {
					res[bb * (rnum + 1) * (vnum + 1) * (snum + 1) + i * (vnum + 1) * (snum + 1) + j * (snum + 1) + k] = bestvalue[k];
				}
			}
		}
	}
}

bool ParabolicFDM::finished() {
	if (current >= (pde.getT().high() - 1e-5)) {return true;}
	else {return false;}
}

void ParabolicFDM::setS(unsigned long& input) {
	snum = input;
}

void ParabolicFDM::setV(unsigned long& input) {
	vnum = input;
}

void ParabolicFDM::setR(unsigned long& input) {
	rnum = input;
}

void ParabolicFDM::setT(unsigned long& input) {
	tnum = input;
}

void ParabolicFDM::setB(unsigned long& input) {
	bnum = input;
}

void ParabolicFDM::setPDE(ParabolicPDE& input) {
	pde = input;
}

std::vector<double> ParabolicFDM::getSARR() {
	return SARR;
}

std::vector<double> ParabolicFDM::getVARR() {
	return VARR;
}

std::vector<double> ParabolicFDM::getRARR() {
	return RARR;
}

std::vector<double> ParabolicFDM::getBARR() {
	return BARR;
}

std::vector<double> ParabolicFDM::getres() {
	return res;
}

#endif