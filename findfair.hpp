#ifndef FFAIR_HPP
#define FFAIR_HPP

class FindFair {
	private:
		ParabolicFDM scheme;

		std::vector<double> sarr;
		std::vector<double> varr;
		std::vector<double> rarr;

		std::vector<double> arr;

		double sTarget;
		double vTarget;
		double rTarget;

		double sWeight;
		double vWeight;
		double rWeight;

		long snum;
		long vnum;
		long rnum;

		long sIndex;
		long vIndex;
		long rIndex;
		long bIndex;

		long maxIter;
		double Tol;

		double interestedV;
		double lastinterestedV;
		double lastlastinterestedV;

		double lastfair;
		double lastlastfair;

		long counter;

		void init();

	public:
		FindFair(const ParabolicFDM& FDM, double STARGET, double VTARGET, double RTARGET, long BINDEX,
			long MAXITER, double TOL);
		void reset();
		void start();

		//setter
		void setFDM(ParabolicFDM& input);

		//getter
		double getsWeight();
		double getvWeight();
		double getrWeight();

		long getsIndex();
		long getvIndex();
		long getrIndex();

		double getInterestedV();
};

#endif