#ifndef PDE_HPP
#define PDE_HPP

class PDE {};

class ParabolicPDE : public PDE {
	private:
		Srange saxis;
		Vrange vaxis;
		Rrange raxis;
		Trange taxis;
		Trange baxis;
		double (*ss) (const double, const double, const double);
		double (*vv) (const double, const double, const double);
		double (*rr) (const double, const double, const double);
		double (*sv) (const double, const double, const double);
		double (*sr) (const double, const double, const double);
		double (*vr) (const double, const double, const double);
		double (*s) (const double, const double, const double);
		double (*v) (const double, const double, const double);
		double (*zero) (const double, const double, const double);
		double (*r) (const double, const double, const double, const double);
		double (*ic) (const double, const double);
		bool (*extremev) (const double);

	public:
		//constructor
		ParabolicPDE();
		ParabolicPDE(const Srange& SRANGE, const Vrange& VRANGE, const Rrange& RRANGE, const Trange& TRANGE, const Trange& BRANGE, 
			double (*SS) (const double, const double, const double), double (*VV) (const double, const double, const double), double (*RR) (const double, const double, const double),
			double (*SV) (const double, const double, const double), double (*SR) (const double, const double, const double), double (*VR) (const double, const double, const double), double (*S) (const double, const double, const double),
			double (*V) (const double, const double, const double), double (*R) (const double, const double, const double, const double), double (*ZERO) (const double, const double, const double), double (*IC) (const double, const double),
			bool (*EXTREMEV) (const double));

		//destructor
		virtual ~ParabolicPDE();

		//setter
		void setS(Srange& input);
		void setV(Vrange& input);
		void setR(Rrange& input);
		void setT(Trange& input);
		void setB(Trange& input);

		//getter
		Srange getS();
		Vrange getV();
		Rrange getR();
		Trange getT();
		Trange getB();
		double SS(const double svar, const double vvar, const double rvar);
		double VV(const double svar, const double vvar, const double rvar);
		double RR(const double svar, const double vvar, const double rvar);
		double SV(const double svar, const double vvar, const double rvar);
		double SR(const double svar, const double vvar, const double rvar);
		double VR(const double svar, const double vvar, const double rvar);
		double S(const double svar, const double vvar, const double rvar);
		double V(const double svar, const double vvar, const double rvar);
		double R(const double s, const double v, const double r, const double t);
		double ZERO(const double svar, const double vvar, const double rvar);
		double IC(const double svar, const double B);
		bool EXTREMEV(const double vvar);

};

#endif