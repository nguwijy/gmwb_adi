#ifndef LU_HPP
#define LU_HPP

class LU
{ // The Balayage method from Godunov

private:
	
		// The vectors
		std::vector<double> aa, bb, cc, dd, ff;

public:
		// Constructors and destructor
		LU();							
		LU (const LU& s2);	

		// Create members to initialise input for AU = F, A = (a,b,c)
		LU(const std::vector<double>& a, const std::vector<double>& b, 
			const std::vector<double>& c, const std::vector<double>& d,
			const std::vector<double>& f);
		virtual ~LU();


		// Operator overloading
		//LU& operator = (const LU& i2);

		// Result; this is a vector in range [1, J-1]
		std::vector<double> result() const;
};
#endif