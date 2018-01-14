#ifndef TBF_HPP
#define TBF_HPP

double Ic(const double s, const double B) {
	double minLeft = ((G > B)? B: G);
	double cf = ((1 - withdrawalPen) * B + withdrawalPen * minLeft);
	return ((cf > s)? cf: s);}

bool ExtremeV(const double v) {return ((v > Gamma)? true: false);}

#endif