#ifndef FFAIR_CPP
#define FFAIR_CPP

#include "findfair.hpp"

void FindFair::init() {
	sarr = scheme.getSARR();
	varr = scheme.getVARR();
	rarr = scheme.getRARR();

	snum = sarr.size() - 1;
	vnum = varr.size() - 1;
	rnum = rarr.size() - 1;

	fairfee = 0.0;
	lastlastfair = fairfee;

	scheme.start();
	arr = scheme.getres();

	for (long ss = 0; ss < (snum + 1); ss++) {
		if (sarr[ss] < sTarget) {
			sIndex = ss;
		} else {
			break;
		}
	}
	sWeight = (sTarget - sarr[sIndex]) / (sarr[sIndex + 1] - sarr[sIndex]);

	for (long vv = 0; vv < (vnum + 1); vv++) {
		if (varr[vv] < vTarget) {
			vIndex = vv;
		} else {
			break;
		}
	}
	vWeight = (vTarget - varr[vIndex]) / (varr[vIndex + 1] - varr[vIndex]);

	for (long rr = 0; rr < (rnum + 1); rr++) {
		if (rarr[rr] < rTarget) {
			rIndex = rr;
		} else {
			break;
		}
	}
	rWeight = (rTarget - rarr[rIndex]) / (rarr[rIndex + 1] - rarr[rIndex]);

	lastlastinterestedV = (1 - sWeight) * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
						+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex])
						+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
						+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex]))
							+ sWeight * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)])
							+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)]));

	fairfee = 0.02;
	lastfair = fairfee;

	scheme.reset();
	scheme.start();
	arr = scheme.getres();

	for (long ss = 0; ss < (snum + 1); ss++) {
		if (sarr[ss] < sTarget) {
			sIndex = ss;
		} else {
			break;
		}
	}
	sWeight = (sTarget - sarr[sIndex]) / (sarr[sIndex + 1] - sarr[sIndex]);

	for (long vv = 0; vv < (vnum + 1); vv++) {
		if (varr[vv] < vTarget) {
			vIndex = vv;
		} else {
			break;
		}
	}
	vWeight = (vTarget - varr[vIndex]) / (varr[vIndex + 1] - varr[vIndex]);

	for (long rr = 0; rr < (rnum + 1); rr++) {
		if (rarr[rr] < rTarget) {
			rIndex = rr;
		} else {
			break;
		}
	}
	rWeight = (rTarget - rarr[rIndex]) / (rarr[rIndex + 1] - rarr[rIndex]);

	lastinterestedV = (1 - sWeight) * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
						+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex])
						+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
						+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex]))
							+ sWeight * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)])
							+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)]));

	interestedV = lastinterestedV;

	counter = 0;
}

FindFair::FindFair(const ParabolicFDM& FDM, double STARGET, double VTARGET, double RTARGET, long BINDEX,
			long MAXITER, double TOL) {
	scheme = FDM;

	sTarget = STARGET;
	vTarget = VTARGET;
	rTarget = RTARGET;
	bIndex = BINDEX;

	maxIter = MAXITER;
	Tol = TOL;

	init();
}

void FindFair::reset() {
	init();
}

void FindFair::start() {
	while ((counter < maxIter) && (((interestedV - sTarget) > Tol) || ((sTarget - interestedV) > Tol))) {
		fairfee = (lastlastfair * (lastinterestedV - sTarget) - lastfair * (lastlastinterestedV - sTarget)) / ((lastinterestedV - sTarget) - (lastlastinterestedV - sTarget));
		scheme.reset();
		scheme.start();
		arr = scheme.getres();
		interestedV = (1 - sWeight) * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex])
							+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
							+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + sIndex]))
								+ sWeight * ((1 - rWeight) * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
								+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)])
								+ rWeight * ((1 - vWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
								+ vWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + (vIndex + 1) * (snum + 1) + (sIndex + 1)]));

		lastlastfair = lastfair;
		lastlastinterestedV = lastinterestedV;

		lastfair = fairfee;
		lastinterestedV = interestedV;
		
		counter ++;
	}
}

void FindFair::setFDM(ParabolicFDM& input) {
	scheme = input;
}

double FindFair::getsWeight() {
	return sWeight;
}

double FindFair::getvWeight() {
	return vWeight;
}

double FindFair::getrWeight() {
	return rWeight;
}

long FindFair::getsIndex() {
	return sIndex;
}

long FindFair::getvIndex() {
	return vIndex;
}

long FindFair::getrIndex() {
	return rIndex;
}

double FindFair::getInterestedV() {
	return interestedV;
}

#endif