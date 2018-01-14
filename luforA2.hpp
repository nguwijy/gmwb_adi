#ifndef LUFORA2_HPP
#define LUFORA2_HPP

class LUforA2
{ // The Balayage method from Godunov

private:
	
		// The vectors
		std::vector<double> aa, bb, cc, dd, ff;
		double ee;

public:
		// Constructors and destructor
		LUforA2 ();							
		LUforA2 (const LU& s2);	

		// Create members to initialise input for AU = F, A = (a,b,c)
		LUforA2 (const std::vector<double>& a, const std::vector<double>& b, 
			const std::vector<double>& c, const std::vector<double>& d,
			const std::vector<double>& f, const double& e);
		virtual ~LUforA2();


		// Operator overloading
		//LU& operator = (const LU& i2);

		// Result; this is a vector in range [1, J-1]
		std::vector<double> result() const;
};
#endif