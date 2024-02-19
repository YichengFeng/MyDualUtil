#ifndef MyToolkit_H
#define MyToolkit_H

#include <iostream>
#include <vector>
#include <cmath>

#include <TGraphErrors.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>


namespace MyToolkit
{ // namespace


//------------------------------------------------------------------------------
const TGraphErrors operator+(const double c, const TGraphErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		y.push_back(y1[i]+c);
		ye.push_back(y1e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors operator+(const TGraphErrors &g1, const double c) {

	return c+g1;
}


//------------------------------------------------------------------------------
const TGraphErrors operator*(const double c, const TGraphErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		y.push_back(y1[i]*c);
		ye.push_back(fabs(y1e*c));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors operator*(const TGraphErrors &g1, const double c) {

	return c*g1;
}


//------------------------------------------------------------------------------
const TGraphErrors operator-(const double c, const TGraphErrors &g1) {

	return c+(-1.0*g1);
}


//------------------------------------------------------------------------------
const TGraphErrors operator-(const TGraphErrors &g1, const double c) {

	return g1+(-1.0*c);
}


//------------------------------------------------------------------------------
const TGraphErrors operator/(const TGraphErrors &g1, const double c) {

	return g1*(1.0/c);
}


//------------------------------------------------------------------------------
const TGraphErrors operator/(const double c, const TGraphErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		double tmpy = c/y1[i];
		double tmpye = c*y1e/y1[i]/y1[i];
		if(c==0) {
			tmpy  = 0;
			tmpye = 0;
		}
		y.push_back(tmpy);
		ye.push_back(tmpye);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}

//------------------------------------------------------------------------------
const TGraphErrors operator+(const TGraphErrors &g1, const TGraphErrors &g2) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	double *x2 = g2.GetX();
	double *y2 = g2.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		double y2e = g2.GetErrorY(i);
		y.push_back(y1[i] + y2[i]);
		ye.push_back(sqrt(y1e*y1e + y2e*y2e));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors operator-(const TGraphErrors &g1, const TGraphErrors &g2) {

	return g1+(-1.0*g2);
}


//------------------------------------------------------------------------------
const TGraphErrors operator*(const TGraphErrors &g1, const TGraphErrors &g2) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	double *x2 = g2.GetX();
	double *y2 = g2.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		double y2e = g2.GetErrorY(i);
		y.push_back(y1[i]*y2[i]);
		ye.push_back(fabs(y1[i]*y2[i])*sqrt(y1e*y1e/y1[i]/y1[i] + y2e*y2e/y2[i]/y2[i]));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors operator/(const TGraphErrors &g1, const TGraphErrors &g2) {

	return (g1 * (1.0/g2));
}


//------------------------------------------------------------------------------
const TGraphErrors operator-(const TGraphErrors &g, const TF1 &f) {

	const int n= g.GetN();
	double *gx = g.GetX();
	double *gy = g.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(gx[i]);
		xe.push_back(0);
		double gye = g.GetErrorY(i);
		y.push_back(gy[i] - f.Eval(gx[i]));
		ye.push_back(gye);
	}

	return TGraphErrors(x.size(), &x[0], &y[0], &xe[0], &ye[0]);
}


//------------------------------------------------------------------------------
const TGraphErrors GraphAbs(const TGraphErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		y.push_back(fabs(y1[i]));
		ye.push_back(fabs(y1e));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphPow(const TGraphErrors &g1, double c) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		y.push_back(pow(y1[i],c));
		ye.push_back(y1e*c*pow(y1[i],c-1));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphSqrt(const TGraphErrors &g1) {

	return GraphPow(g1, 0.5);
}


//------------------------------------------------------------------------------
const TGraphErrors GraphFrom(const TH1D *h) {

	if(h==nullptr) return TGraphErrors();

	const int n = h->GetSize()-2;
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;
	for(int i=0; i<n; i++) {
		x.push_back(h->GetBinCenter(i+1));
		xe.push_back(0);
		y.push_back(h->GetBinContent(i+1));
		ye.push_back(h->GetBinError(i+1));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphFrom(const TH1F *h) {

	if(h==nullptr) return TGraphErrors();

	const int n = h->GetSize()-2;
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;
	for(int i=0; i<n; i++) {
		x.push_back(h->GetBinCenter(i+1));
		xe.push_back(0);
		y.push_back(h->GetBinContent(i+1));
		ye.push_back(h->GetBinError(i+1));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphFrom(const TGraphErrors *g1) {

	if(g1==nullptr) return TGraphErrors();

	TGraphErrors g(g1->GetN(), g1->GetX(), g1->GetY(), g1->GetEX(), g1->GetEY());

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphCent(const TGraphErrors &g1) {

	double cent[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
	//double cente[9]= { 5,  5,  5,  5,  5,  5,  5, 2.5, 2.5};
	double cente[9]= {0};

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;

	int di = 0;
	for(int i=0; i<n; i++) {
		if(n == 10 && i == 0) {
			di = 1;
			continue;
		}
		double y1e = g1.GetErrorY(i);
		x.push_back(cent[i-di]);
		//xe.push_back(cente[i-di]);
		xe.push_back(g1.GetErrorX(i));
		y.push_back(y1[i]);
		ye.push_back(y1e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphCent8(const TGraphErrors &g1) {

	double cent[8] = {75, 65, 55, 45, 35, 25, 15, 5};
	//double cente[8]= { 5,  5,  5,  5,  5,  5,  5, 5};
	double cente[8]= {0};

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;

	int di = 0;
	for(int i=0; i<8; i++) {
		double y1e = g1.GetErrorY(i);
		x.push_back(cent[i-di]);
		xe.push_back(cente[i-di]);
		y.push_back(y1[i]);
		ye.push_back(y1e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphCent16(const TGraphErrors &g1) {

	double cent[16] = {77.5, 72.5, 67.5, 62.5, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.5, 2.5};
	//double cente[16]= {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
	double cente[16]= {0};

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;

	int di = 0;
	for(int i=0; i<n; i++) {
		if(n == 17 && i == 0) {
			di = 1;
			continue;
		}
		double y1e = g1.GetErrorY(i);
		x.push_back(cent[i-di]);
		xe.push_back(cente[i-di]);
		y.push_back(y1[i]);
		ye.push_back(y1e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphAsymmErrors GraphCent(const TGraphAsymmErrors &g1) {

	double cent[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
	//double cente[9]= { 5,  5,  5,  5,  5,  5,  5, 2.5, 2.5};
	double cente[9]= {0};

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> xle;
	std::vector<double> xhe;
	std::vector<double> y;
	std::vector<double> yle;
	std::vector<double> yhe;

	int di = 0;
	for(int i=0; i<n; i++) {
		if(n == 10 && i == 0) {
			di = 1;
			continue;
		}
		x.push_back(cent[i-di]);
		xle.push_back(g1.GetErrorXlow(i));
		xhe.push_back(g1.GetErrorXhigh(i));
		y.push_back(y1[i]);
		yle.push_back(g1.GetErrorYlow(i));
		yhe.push_back(g1.GetErrorYhigh(i));
	}

	TGraphAsymmErrors g(x.size(), &x[0], &y[0], &xle[0],&xhe[0], &yle[0],&yhe[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraph GraphEdge(const TGraphErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> y;

	for(int i=0; i<n; i++) {
		double y1e = g1.GetErrorY(i);
		x.push_back(x1[i]);
		y.push_back(y1[i]+y1e);
	}
	for(int i=n-1; i>=0; i--) {
		double y1e = g1.GetErrorY(i);
		x.push_back(x1[i]);
		y.push_back(y1[i]-y1e);
	}
	x.push_back(x[0]);
	y.push_back(y[0]);

	TGraph g(x.size(), &x[0], &y[0]);
	g.SetLineColor(g1.GetLineColor());
	g.SetMarkerColor(g1.GetMarkerColor());

	return g;
}


//------------------------------------------------------------------------------
const TGraph GraphEdge(const TGraphAsymmErrors &g1) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	std::vector<double> x;
	std::vector<double> y;

	for(int i=0; i<n; i++) {
		double y1e = g1.GetErrorYhigh(i);
		x.push_back(x1[i]);
		y.push_back(y1[i]+y1e);
	}
	for(int i=n-1; i>=0; i--) {
		double y1e = g1.GetErrorYlow(i);
		x.push_back(x1[i]);
		y.push_back(y1[i]-y1e);
	}
	x.push_back(x[0]);
	y.push_back(y[0]);

	TGraph g(x.size(), &x[0], &y[0]);
	g.SetLineColor(g1.GetLineColor());
	g.SetMarkerColor(g1.GetMarkerColor());

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphOne(double con, double err) {

	double x  = 0;
	double xe = 0;
	double y  = con;
	double ye = err;

	TGraphErrors g(1, &x, &y, &xe, &ye);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphOne(double vx, double ex, double con, double err) {

	double x  = vx;
	double xe = ex;
	double y  = con;
	double ye = err;

	TGraphErrors g(1, &x, &y, &xe, &ye);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors RemoveZero(const TGraphErrors &g) {

	const int n = g.GetN();
	double *x = g.GetX();
	double *y = g.GetY();
	double xe[n];
	double ye[n];

	int on = 0;
	double ox[n];
	double oy[n];
	double oxe[n];
	double oye[n];

	for(int i=0; i<n; i++) {
		xe[i] = g.GetErrorX(i);
		ye[i] = g.GetErrorY(i);

		if(y[i]!=0 && (ye[i]!=0 && !isnan(ye[i]) && !isinf(ye[i]))) {
			ox[on] = x[i];
			oy[on] = y[i];
			oxe[on] = xe[i];
			oye[on] = ye[i];
			on ++;
		}
	}

	TGraphErrors og(on, ox, oy, oxe, oye);

	return og;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphMap(const TGraphErrors &g1, const TGraphErrors &g2) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	double *x2 = g2.GetX();
	double *y2 = g2.GetY();
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;

	for(int i=0; i<n; i++) {
		double y1e = g1.GetErrorY(i);
		double y2e = g2.GetErrorY(i);
		x.push_back(y1[i]);
		xe.push_back(y1e);
		y.push_back(y2[i]);
		ye.push_back(y2e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphMerge(const TGraphErrors &g1, const TGraphErrors &g2) {

	const int n1= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	const int n2= g2.GetN();
	double *x2 = g2.GetX();
	double *y2 = g2.GetY();
	std::vector<double> x;
	std::vector<double> xe;
	std::vector<double> y;
	std::vector<double> ye;

	for(int i=0; i<n1; i++) {
		double x1e = g1.GetErrorX(i);
		double y1e = g1.GetErrorY(i);
		x.push_back(x1[i]);
		xe.push_back(x1e);
		y.push_back(y1[i]);
		ye.push_back(y1e);
	}
	for(int i=0; i<n2; i++) {
		double x2e = g2.GetErrorX(i);
		double y2e = g2.GetErrorY(i);
		x.push_back(x2[i]);
		xe.push_back(x2e);
		y.push_back(y2[i]);
		ye.push_back(y2e);
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}


//------------------------------------------------------------------------------
const TGraphErrors GraphErrToVal(const TGraphErrors &g1) {
	TGraphErrors g(g1.GetN(), g1.GetX(), g1.GetEY(), g1.GetEX(), g1.GetEY());
	for(int i=0; i<g.GetN(); i++) g.GetEY()[i] = 0;
	return g;
}


////------------------------------------------------------------------------------
//const TGraphErrors GraphMoveSigma(const TGraphErrors &g1, double c) {
//	TGraphErrors ge = GraphErrToVal(g1);
//	return (g1 + c*ge);
//}


//------------------------------------------------------------------------------
double CalcChiFromReso(double res, double (*CalcResoFromChi)(double)) {
	// Calculates chi from the event plane resolution
	double chi   = 2.0;
	double delta = 1.0;

	for (int i = 0; i < 15; i++) {
		chi   = (CalcResoFromChi(chi) < res) ? chi + delta : chi - delta;
		delta = delta / 2.;
	}

	return chi;
}

//------------------------------------------------------------------------------
double CalcResoPsiK2(double chi) {
	// Calculates the event plane resolution as a function of chi
	//  for the case k=2.

	double con = sqrt(M_PI/2)/2; //0.626657;
	double arg = chi*chi/4.;
	double halfpi = M_PI/2; //1.570796;
	double besselOneHalf = sqrt(arg/halfpi) * sinh(arg)/arg;
	double besselThreeHalfs = sqrt(arg/halfpi) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
	double res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

	// Approximations.
	//res = 0.25*chi*chi-0.01141*chi*chi*chi-0.034726*chi*chi*chi*chi+0.006815*chi*chi*chi*chi*chi;

	return res;
}

//------------------------------------------------------------------------------
double CalcResoPsiK1(double chi) {
	// Calculates the event plane resolution as a function of chi

	double con = sqrt(M_PI/2)/2; //0.626657;
	double arg = chi*chi/4.;

	double res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

	return res;
}

//------------------------------------------------------------------------------
double CalcDerivative(double x, double (*f)(double)) {
	// very simple derivative calculator

	double xl;
	double xh;

	if(x!=0) {
		xl = x*0.99;
		xh = x*1.01;
	} else {
		xl = x-1e-4;
		xh = x+1e-4;
	}

	return (f(xh)-f(xl))/(xh-xl);
}

//------------------------------------------------------------------------------
const TGraphErrors ResHalfToFull(const TGraphErrors &g) {
	// input half resolution, output full resolution

	const int n = g.GetN();

	double *x = g.GetX();
	double *y = g.GetY();
	double xe[n];
	double ye[n];
	double oy[n];
	double oye[n];

	for(int i=0; i<n; i++) {
		xe[i] = g.GetErrorX(i);
		ye[i] = g.GetErrorY(i);
		double TmpRes = y[i];
		double TmpResErr = ye[i];
		double TmpChi = CalcChiFromReso(TmpRes, CalcResoPsiK1);
		double TmpChiErr = TmpResErr/fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
		TmpChi    *= sqrt(2.0);
		TmpChiErr *= sqrt(2.0);
		double TmpFullRes = CalcResoPsiK1(TmpChi);
		double TmpFullResErr = TmpChiErr*fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
		oy[i]  = TmpFullRes;
		oye[i] = TmpFullResErr;
		if(ye[i]==0) oye[i]=0;
		if(y[i] ==0) oy[i] =0;
	}

	TGraphErrors og(n, x, oy, xe, oye);

	return og;
}

//------------------------------------------------------------------------------
const TGraphErrors ResHalfToFullK2(const TGraphErrors &g) {
	// input half resolution, output full resolution

	const int n = g.GetN();

	double *x = g.GetX();
	double *y = g.GetY();
	double xe[n];
	double ye[n];
	double oy[n];
	double oye[n];

	for(int i=0; i<n; i++) {
		xe[i] = g.GetErrorX(i);
		ye[i] = g.GetErrorY(i);
		double TmpRes = y[i];
		double TmpResErr = ye[i];
		double TmpChi = CalcChiFromReso(TmpRes, CalcResoPsiK1);
		double TmpChiErr = TmpResErr/fabs(CalcDerivative(TmpChi, CalcResoPsiK1));
		TmpChi    *= sqrt(2.0);
		TmpChiErr *= sqrt(2.0);
		double TmpFullRes = CalcResoPsiK2(TmpChi);
		double TmpFullResErr = TmpChiErr*fabs(CalcDerivative(TmpChi, CalcResoPsiK2));
		oy[i]  = TmpFullRes;
		oye[i] = TmpFullResErr;
		if(ye[i]==0) oye[i]=0;
		if(y[i] ==0) oy[i] =0;
	}

	TGraphErrors og(n, x, oy, xe, oye);

	return og;
}

//------------------------------------------------------------------------------
const TGraphErrors GraphBarlowDiff(const TGraphErrors &g1, const TGraphErrors &g2) {

	const int n= g1.GetN();
	double *x1 = g1.GetX();
	double *y1 = g1.GetY();
	double *x2 = g2.GetX();
	double *y2 = g2.GetY();
	vector<double> x;
	vector<double> xe;
	vector<double> y;
	vector<double> ye;

	for(int i=0; i<n; i++) {
		x.push_back(x1[i]);
		xe.push_back(0);
		double y1e = g1.GetErrorY(i);
		double y2e = g2.GetErrorY(i);
		y.push_back(y1[i] - y2[i]);
		ye.push_back(sqrt(fabs(y1e*y1e - y2e*y2e)));
	}

	TGraphErrors g(x.size(), &x[0], &y[0], &xe[0], &ye[0]);

	return g;
}

//------------------------------------------------------------------------------
} // namespace

#endif
