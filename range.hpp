#ifndef Range_HPP
#define Range_HPP

class Range {
	protected:
		double lo;
		double hi;

	public:
		//constructor
		Range();
		Range(const double low, const double high);
		Range(const Range& ran2);

		//destructor
		virtual ~Range();

		//setter
		void low(const double& t1);
		void high(const double& t1);

		//getter
		double low() const;
		double high() const;

		//assignment operator
		Range& operator = (const Range& ran2);

		//mesher
		std::vector<double> mesh(const int nSteps) const;
};

class Srange : public Range {
	private:

	public:
		//need to define constructor because default constructor not working properly
		Srange();
		Srange(const double low, const double high);
		//mesher
		std::vector<double> mesh(const int nSteps, const double sImp) const;
};

class Vrange : public Range {
	private:

	public:
		//need to define constructor because default constructor not working properly
		Vrange();
		Vrange(const double low, const double high);
		//mesher
		std::vector<double> mesh(const int nSteps) const;
};

class Rrange : public Range {
	private:

	public:
		//need to define constructor because default constructor not working properly
		Rrange();
		Rrange(const double low, const double high);
		//mesher
		std::vector<double> mesh(const int nSteps, const double rImp) const;
};

class Trange : public Range {
	private:

	public:
		//need to define constructor because default constructor not working properly
		Trange();
		Trange(const double low, const double high);
		//mesher
		std::vector<double> mesh(const int nSteps) const;
};

#endif