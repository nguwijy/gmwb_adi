#ifndef CHAR_HPP
#define CHAR_HPP

const double kappa = 1.0;		//v - mean reversion
const double Gamma = .2 * .2;	//v - mean
const double sig1 = .2;			//v - vol
const double a = 1.0;			//r - mean reversion
/*const double c1 = .05;			//r - mean, for b(t)
const double c2 = 0.0;			//r - mean, for b(t)
const double c3 = 0.0;*/			//r - mean, for b(t)
const double c1 = .05;			//r - mean, for b(t)
const double c2 = 0.01;			//r - mean, for b(t)
const double c3 = 1.0;			//r - mean, for b(t)
const double rkappa = 1.0;
const double rGamma = .05;
const double sig2 = .2;			//r - vol
const double rho12 = - .5;		//relationship sv
const double rho13 = - .5;		//relationship sr
const double rho23 = 0.0;		//relationship vr
/*const double rho12 = -.5;		//relationship sv
const double rho13 = -。5;		//relationship sr
const double rho23 = 0;*/		//relationship vr
double fairfee = .01;		//fair fee charged for insurance company
const double managefee = .01;	//management fee charged by mutual fund
//const double fairfee = 0.0117 + managefee;

const double Vfrom = 0.0;		//volatility domain
const double Vto = 4.0;			//volatility domain
//min(max(100*gamma,1),5)
const double Rto = 0.5;			//interest rate domain
//changed from -1.0 to -0.1
const double Rfrom = - Rto;		//interest rate domain
const double Tfrom = 0.0;		//time domain
double Tto = 5.0;			//time domain
//const double Tto = 5.0;			//time domain
unsigned long Snum = 20;	//number of stock price
unsigned long Vnum = 4;	//number of volatility
unsigned long Rnum = 10;	//number of interest rate
unsigned long Tnum = 10;	//number of time
unsigned long Bnum = 5;	//number of guaranteed account

const double withdrawalPen = .1;//withdrawal penalty
const double P0 = 100.0;		//initial premium
const double WF = 1.0;			//withdrawal frequence
const double G = P0 / (Tto * WF);	//guaranteed minimum withdrawal

const double Sfrom = 0.0;		//stock price domain
//changed from 100 * Tto * P0 to 14 * P0
const double Sto = 10.0 * P0;	//stock price domain
const double Bfrom = 0.0;		//guaranteed account domain
const double Bto = P0;			//guaranteed account domain

const double Rint = 0.05;		//the point in R that we are particularly interested in

const double theta = 2.0/3.0;

#endif