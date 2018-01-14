#ifndef FDM_HPP
#define FDM_HPP

enum schemeType {Douglas, CS, MCS, HV};			//types of scheme

class FDM {};

class ParabolicFDM : public FDM {
	private:
		ParabolicPDE pde;

		schemeType typ;

		unsigned long snum;
		unsigned long vnum;
		unsigned long rnum;
		unsigned long tnum;
		unsigned long bnum;
		double the;

		double h;					//del T

		double current;				//time to maturity, start from maturity, current = 0

		std::vector<double> SARR;
		std::vector<double> VARR;
		std::vector<double> RARR;
		std::vector<double> TARR;
		std::vector<double> BARR;

		std::vector<double> tmp;	//value in current time as specified by B
		std::vector<double> res;	//value in current time in overall

		std::vector<double> A1next;
		std::vector<double> A1last;
		std::vector<double> A1now;
		std::vector<double> A2next;
		std::vector<double> A2last;
		std::vector<double> A2lastlast;
		std::vector<double> A2now;
		std::vector<double> A3next_now;
		std::vector<double> A3last_now;
		std::vector<double> A3now_now;
		std::vector<double> A3next_next;
		std::vector<double> A3last_next;
		std::vector<double> A3now_next;

		std::vector<double> snext;
		std::vector<double> slast;
		std::vector<double> snow;
		std::vector<double> vnext;
		std::vector<double> vlast;
		std::vector<double> vnow;
		std::vector<double> rnext;
		std::vector<double> rlast;
		std::vector<double> rnow;

		void init();				//initialize all matrices

	public:
		ParabolicFDM();
		ParabolicFDM(const ParabolicPDE& context, unsigned long Sintervals, unsigned long Vintervals, 
				unsigned long Rintervals, unsigned long Tintervals, unsigned long Bintervals, const double& theta, schemeType type);

		void start();
		void reset();
		void douglas();
		void cs();
		void mcs();
		void hv();
		void opti();
		bool finished();
		
		//setter
		void setS(unsigned long& input);
		void setV(unsigned long& input);
		void setR(unsigned long& input);
		void setT(unsigned long& input);
		void setB(unsigned long& input);
		void setPDE(ParabolicPDE& input);

		//getter
		std::vector<double> getSARR();
		std::vector<double> getVARR();
		std::vector<double> getRARR();
		std::vector<double> getBARR();
		std::vector<double> getres();
};

#endif