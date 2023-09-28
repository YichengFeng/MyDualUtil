/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 * Note: dual number to calculate precise values of first-order derivative
 *       of multi-variable functions with common math formula
 *       this class is wrapped with ROOT for I/O and plotting 
 **************************************************************************/

#ifndef MyDualGraph_H
#define MyDualGraph_H

#include <iostream>
#include <cmath>
#include <vector>

#include "MyDualNumber.h"
#include "MyDualMultiv.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1.h"


class MyDualPoint
{
public:
	MyDualNumber Px;
	MyDualMultiv Py;

	MyDualPoint() {
	}
	MyDualPoint(MyDualNumber px) {
		Px = px;
	}
	MyDualPoint(MyDualNumber px, MyDualMultiv py) {
		Px = px;
		Py = py;
	}
	~MyDualPoint() {
	}
};


class MyDualGraph
{
private:
	std::vector<MyDualPoint> Points;
	bool IsUpdated;

public:
	TGraphErrors Graph;

	static int IdxNow; // the latest index
	static const int IdxMax = 1000000; // maximal index

	bool GetIsUpdated() const {
		return IsUpdated;
	}

	bool CheckSize(int n) const {
		if(n!=(int)Points.size()) {
			std::cout << "MyDualGraph: warning: size not match!" << std::endl;
			return false;
		}
		return true;
	}
	bool CheckIndex(int i) const {
		int n = (int)Points.size();
		if(i<0 || i>=n) {
			std::cout << "MyDualGraph:: warning: index out of range! " << i << "/" << n << std::endl;
			return false;
		}
		return true;
	}

	// normally SetPoints() is not used.
	void SetPoints(const std::vector<MyDualPoint> &points) {
		Points = points;
		IsUpdated = false;
	}

	// normally SetDual() is not used.
	void SetDual(int indx, const TH1 &h, bool xreset, bool xerror, bool yreset) {
		int n = h.GetXaxis()->GetNbins();
		if(!CheckSize(n)) return;
		for(int i=0; i<n; i++) {
			if(xreset) Points[i].Px = MyDualNumber(h.GetXaxis()->GetBinCenter(i+1), xerror?(h.GetXaxis()->GetBinWidth(i+1)*0.5):0);
			if(yreset) Points[i].Py.SetValu(h.GetBinContent(i+1));
			Points[i].Py.SetDual(indx, h.GetBinError(i+1));
		}
		IsUpdated = false;
	}
	void SetDual(int indx, const TH1 &h) {
		SetDual(indx, h, false, false, false);
	}
	void SetDual(int indx, const TGraphErrors &g, bool xreset, bool xerror, bool yreset) {
		int n = g.GetN();
		if(!CheckSize(n)) return;
		for(int i=0; i<n; i++) {
			if(xreset) Points[i].Px = MyDualNumber(g.GetX()[i], xerror?(g.GetEX()[i]):0);
			if(yreset) Points[i].Py.SetValu(g.GetY()[i]);
			Points[i].Py.SetDual(indx, g.GetEY()[i]);
		}
		IsUpdated = false;
	}
	void SetDual(int indx, const TGraphErrors &g) {
		SetDual(indx, g, false, false, false);
	}

	void Reset(int n=0) {
		Points.clear();
		for(int i=0; i<n; i++) {
			MyDualPoint dp;
			Points.push_back(dp);
		}
		IsUpdated = false;
	}

	static int AutoNewIdx() { // IdxNow range 1 ~ IdxMax-1
		if(IdxNow>=IdxMax-1) {
			std::cout << "MyDualGraph::IdxNow out of range! " << IdxNow << "/" << IdxMax-1 << std::endl;
		} else {
			IdxNow ++;
		}
		return IdxNow;
	}

	MyDualGraph() {
		Reset(0);
	}
	MyDualGraph(int n) {
		Reset(n);
	}
	MyDualGraph(const TH1 &h, bool xerror=false) {
		Reset(h.GetXaxis()->GetNbins());
		AutoNewIdx();
		SetDual(IdxNow, h, true, xerror, true);
	}
	MyDualGraph(const TGraphErrors &g, bool xerror=false) {
		Reset(g.GetN());
		AutoNewIdx();
		SetDual(IdxNow, g, true, xerror, true);
	}
	MyDualGraph(int indx, const TH1 &h, bool xerror=false) {
		Reset(h.GetXaxis()->GetNbins());
		SetDual(indx, h, true, xerror, true);
	}
	MyDualGraph(int indx, const TGraphErrors &g, bool xerror=false) {
		Reset(g.GetN());
		SetDual(indx, g, true, xerror, true);
	}
	MyDualGraph(const std::vector<MyDualPoint> &points) {
		Points = points;
		IsUpdated = false;
	}
	MyDualGraph(const MyDualGraph &dg) {
		//Reset(0);
		//Points = dg.Points;
		*this = dg;
	}
	~MyDualGraph() {
		Reset(0);
	}

	int GetIdxOne() const {
		if(Points.size()<=0) {
			std::cout << "MyDualGraph::GetIdxOne() empty Points!" << std::endl;
			return 0; // 0 is not used as an index
		}
		std::vector<int> list = Points[0].Py.GetList();
		if(list.size()!=1) {
			std::cout << "MyDualGraph::GetIdxOne() not the single index!" << std::endl;
			return 0;
		} else {
			return list[0];
		}
	}

	int GetN() const {
		return (int)Points.size();
	}

	MyDualPoint GetPoint(int i) const {
		if(!CheckIndex(i)) return MyDualPoint();
		return Points[i];
	}

	void SetPoint(int i, const MyDualPoint &dp) {
		if(!CheckIndex(i)) return;
		Points[i] = dp;
	}

	void Calc() {
		int n = (int)Points.size();
		double vx[n];
		double ex[n];
		double vy[n];
		double ey[n];
		for(int i=0; i<n; i++) {
			vx[i] = Points[i].Px.GetValu();
			ex[i] = Points[i].Px.GetDual();
			vy[i] = Points[i].Py.GetValu();
			ey[i] = Points[i].Py.GetUnce();
		}
		Graph = TGraphErrors(n, vx, vy, ex, ey);
		IsUpdated = true;
	}

	void ShiftX() {
		if(!IsUpdated) Calc();
		for(int i=0; i<Graph.GetN(); i++) {
			Graph.GetX()[i] = Points[i].Px.GetValu();
		}
	}

	void ShiftX(double dx) {
		if(!IsUpdated) Calc();
		for(int i=0; i<Graph.GetN(); i++) {
			Graph.GetX()[i] += dx;
		}
	}

	TString StrLatex(int i=0, string symboltype="R") const {
		TString strlatex;
		if(CheckIndex(i)) strlatex = Points[i].Py.StrLatex(symboltype);
		return strlatex;
	}

	// merge
	// must use different bins for merge, so they should be independent
	// this code will take care of that: #bin * -IdxMax
	const MyDualGraph AveBin(std::vector<int> range) const {
		std::sort( range.begin(), range.end() );
		range.erase( unique( range.begin(), range.end() ), range.end() );
		int n = (int)range.size();
		if(n<=0) {
			std::cout << "MyDualGraph::AveBin() invalid range!" << std::endl;
			return MyDualGraph();
		}
		if(!CheckIndex(range[0]) || !CheckIndex(range[n-1])) {
			std::cout << "MyDualGraph::AveBin() invalid range!" << std::endl;
			return MyDualGraph();
		}
		MyDualNumber avex(0,0);
		MyDualMultiv avey;
		for(int k=0; k<n; k++) {
			int i = range[k];
			avex = avex + Points[i].Px;
			MyDualMultiv dm;
			dm.SetValu(Points[i].Py.GetValu());
			std::vector<int> list = Points[i].Py.GetList();
			for(int j=0; j<(int)list.size(); j++) {
				int idx = list[j];
				dm.SetDual(idx-IdxMax*i, Points[i].Py.GetDual(idx));
			}
			avey = k==0?dm:AvePlus(avey, dm);
		}
		avex = (avex / (1.0*n));
		std::vector<MyDualPoint> points;
		points.push_back(MyDualPoint(avex, avey));
		return MyDualGraph(points);
	}

	const MyDualGraph AveBin(int bl, int bh) const {
		if(!CheckIndex(bl) || !CheckIndex(bh) || bl>bh) {
			std::cout << "MyDualGraph::AveBin() invalid range!" << std::endl;
			return MyDualGraph();
		}
		std::vector<int> range;
		for(int i=bl; i<=bh; i++) {
			range.push_back(i);
		}
		return AveBin(range);
	}

	const MyDualGraph MergeUnce() {
		if(!IsUpdated) Calc();
		return MyDualGraph(Graph, true);
	}

	void Write(TFile *f, TString name) const {
		f->WriteObjectAny(this, "MyDualGraph", name);
	}
};


int MyDualGraph::IdxNow = 0;


// math
const MyDualGraph operator+(const MyDualGraph &dg1, double c) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, (dp.Py + c)));
	}
	return dg2;
}

const MyDualGraph operator+(double c, const MyDualGraph &dg1) {
	return (dg1 + c);
}

const MyDualGraph operator*(const MyDualGraph &dg1, double c) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, (dp.Py * c)));
	}
	return dg2;
}

const MyDualGraph operator*(double c, const MyDualGraph &dg1) {
	return (dg1 * c);
}

const MyDualGraph operator-(const MyDualGraph &dg1, double c) {
	return (dg1 + (-1.0*c));
}

const MyDualGraph operator-(double c, const MyDualGraph &dg1) {
	return (-1.0*dg1 + c);
}

const MyDualGraph operator/(const MyDualGraph &dg1, double c) {
	return (dg1 * (1.0/c));
}

const MyDualGraph operator/(double c, const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, (c / dp.Py)));
	}
	return dg2;
}

const MyDualGraph operator+(const MyDualGraph &dg1) {
	return dg1;
}

const MyDualGraph operator-(const MyDualGraph &dg1) {
	return (0.0 - dg1);
}

const MyDualGraph operator+(const MyDualGraph &dg1, const MyDualGraph &dg2) {
	int n = dg1.GetN();
	if(!dg2.CheckSize(n)) return MyDualGraph();
	MyDualGraph dg3(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp1 = dg1.GetPoint(i);
		MyDualPoint dp2 = dg2.GetPoint(i);
		dg3.SetPoint(i, MyDualPoint(dp1.Px, (dp1.Py + dp2.Py)));
	}
	return dg3;
}

const MyDualGraph operator-(const MyDualGraph &dg1, const MyDualGraph &dg2) {
	return (dg1 + (-1.0*dg2));
}

const MyDualGraph operator*(const MyDualGraph &dg1, const MyDualGraph &dg2) {
	int n = dg1.GetN();
	if(!dg2.CheckSize(n)) return MyDualGraph();
	MyDualGraph dg3(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp1 = dg1.GetPoint(i);
		MyDualPoint dp2 = dg2.GetPoint(i);
		dg3.SetPoint(i, MyDualPoint(dp1.Px, (dp1.Py * dp2.Py)));
	}
	return dg3;
}

const MyDualGraph operator/(const MyDualGraph &dg1, const MyDualGraph &dg2) {
	return (dg1 * (1.0/dg2));
}

const MyDualGraph pow(const MyDualGraph &dg1, double c) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, pow(dp.Py, c)));
	}
	return dg2;
}

const MyDualGraph sqrt(const MyDualGraph &dg1) {
	return pow(dg1, 0.5);
}

const MyDualGraph abs(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, abs(dp.Py)));
	}
	return dg2;
}

const MyDualGraph exp(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, exp(dp.Py)));
	}
	return dg2;
}

const MyDualGraph log(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, log(dp.Py)));
	}
	return dg2;
}

const MyDualGraph sin(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, sin(dp.Py)));
	}
	return dg2;
}

const MyDualGraph cos(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, cos(dp.Py)));
	}
	return dg2;
}

const MyDualGraph tan(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, tan(dp.Py)));
	}
	return dg2;
}

const MyDualGraph atan(const MyDualGraph &dg1) {
	int n = dg1.GetN();
	MyDualGraph dg2(n);
	for(int i=0; i<n; i++) {
		MyDualPoint dp = dg1.GetPoint(i);
		dg2.SetPoint(i, MyDualPoint(dp.Px, atan(dp.Py)));
	}
	return dg2;
}


#endif
