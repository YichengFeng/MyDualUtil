#ifndef MySystGraph_H
#define MySystGraph_H

#include <iostream>
#include <cmath>

#include "MyDualNumber.h"
#include "MySystGraph.h"
#include "MyDualGraph.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"


class MyPackGraph
{
public:
	MyDualGraph g; // graph
	TString s;
	double w; // weight
	int m; // mode: 1: s^2; 2: s^2-e^2

	MyPackGraph() {
		w = 0;
		m = 1;
	}
	MyPackGraph(const MyDualGraph &gg, TString ss="", double ww=1.0, int mm=1) {
		g = gg;
		s = ss;
		w = ww;
		if(mm==1 || mm==2) {
			m = mm;
		} else {
			m = 1;
		}
	}
	MyPackGraph(const MyPackGraph &pg) {
		g = pg.g;
		s = pg.s;
		w = pg.w;
		m = pg.m;
	}
};


class MySystGraph
{
private:
	int Mode;

	// data are saved as MyDualGraph
	MyDualGraph Def;
	std::map<int,MyPackGraph> Var;

	// systematic box width
	double Width;

	// is graphs updated
	bool IsUpdated;

public:
	// graphs for plot
	TGraphErrors gStat;
	TGraphErrors gSyst;
	TGraphAsymmErrors gAsym;

	bool GetIsUpdated() const {
		return IsUpdated;
	}

	void Reset() {
		Var.clear();
		Width = 0.8;
		IsUpdated = false;
	}

	MySystGraph() {
		Reset();
		SetMode(0);
	}
	MySystGraph(const MyDualGraph &def, int mode=0) {
		Reset();
		SetMode(mode);
		Def = def;
	}
	MySystGraph(const MySystGraph &sg1) {
		//Reset();
		//Mode = sg1.Mode;
		//Def = sg1.Def;
		//Var = sg1.Var;
		//Width = sg1.Width;
		*this = sg1;
	}
	~MySystGraph() {
		Reset();
	}

	void SetMode(int mode) { // 0: as mode; 1: s^2; 2: s^2-e^2
		if(mode>=0 && mode<3) {
			Mode = mode;
		} else {
			std::cout << "MySystGraph::SetMode() invalid mode " << mode << ". Reset to default 0." << std::endl;
			Mode = 0;
		}
		IsUpdated = false;
	}
	int GetMode() const {
		return Mode;
	}

	void SetWidth(double width) {
		Width = fabs(width);
	}
	double GetWidth() const {
		return Width;
	}

	void SetDef(const MyDualGraph &def) {
		Def = def;
		IsUpdated = false;
	}
	const MyDualGraph GetDef() const {
		return Def;
	}

	std::vector<int> GetList() const {
		std::vector<int> list;
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			list.push_back(it->first);
		}
		return list;
	}

	void AddVar(int idx, const MyPackGraph &pg) {
		if(Var.count(idx)==0) {
			Var.insert({idx, pg});
		} else {
			std::cout << "MySystGraph::AddVar() index " << idx << " already exists!" << std::endl;
		}
		IsUpdated = false;
	}
	void ChgVar(int idx, const MyPackGraph &pg) {
		if(Var.count(idx)==0) {
			std::cout << "MySystGraph::ChgVar() index " << idx << " does not exist!" << std::endl;
		} else {
			Var[idx] = pg;
		}
		IsUpdated = false;
	}
	void SetVar(int idx, const MyPackGraph &pg) {
		if(Var.count(idx)==0) {
			Var.insert({idx, pg});
		} else {
			Var[idx] = pg;
		}
		IsUpdated = false;
	}
	const MyPackGraph GetVar(int idx) const {
		if(Var.count(idx)==0) {
			std::cout << "MySystGraph::GetVar() index " << idx << " does not exist!" << std::endl;
		}
		return Var.at(idx);
	}

	void SetVarAll(const std::map<int,MyPackGraph> &var) {
		Var = var;
		IsUpdated = false;
	}
	const std::map<int,MyPackGraph> GetVarAll() const {
		return Var;
	}

	void Calc() {
		int n = Def.GetN();
		double sx[n];
		double sy[n];
		double syl[n];
		double syh[n];
		for(int i=0; i<n; i++) {
			sx[i] = Width;
			sy[i] = 0;
			syl[i] = 0;
			syh[i] = 0;
		}
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			const int &idx = it->first;
			MyPackGraph &pg = it->second;
			if(!pg.g.GetIsUpdated()) pg.g.Calc();
			if(!pg.g.CheckSize(n)) return;
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = mm==2?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
			}
		}
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
		}
		if(!Def.GetIsUpdated()) Def.Calc();
		gStat = Def.Graph;
		gSyst = TGraphErrors(n, gStat.GetX(), gStat.GetY(), sx, sy);
		gAsym = TGraphAsymmErrors(n, gStat.GetX(), gStat.GetY(), sx, sx, syl, syh);
		IsUpdated = true;
	}

	TString StrLatex(int i=0) {
		if(!Def.CheckIndex(i)) return TString();
		if(!IsUpdated) Calc();

		double e_val, e_err;
		int e_n10, e_np;
		MyDualNumber e_dn(gStat.GetY()[i], gStat.GetEY()[i]);
		e_dn.StrLatex(e_val, e_err, e_n10, e_np);

		double s_val, s_err;
		int s_n10, s_np;
		MyDualNumber s_dn(gSyst.GetY()[i], gSyst.GetEY()[i]);
		s_dn.StrLatex(s_val, s_err, s_n10, s_np);

		double val, err, sys;
		int n10, np;
		//if(s_n10 < e_n10) {
		if(gSyst.GetEY()[i] < gStat.GetEY()[i]) {
			val = s_val;
			err = e_err * pow(10.0, 1.0*(e_n10-s_n10));
			sys = s_err;
			n10 = s_n10;
			np = s_np;
		} else {
			val = e_val;
			err = e_err;
			sys = s_err * pow(10.0, 1.0*(s_n10-e_n10));
			n10 = e_n10;
			np = e_np;
		}

		std::stringstream strval;
		std::stringstream strerr;
		std::stringstream strsys;
		strval << std::fixed << std::setprecision(np) << val;
		strerr << std::fixed << std::setprecision(np) << err;
		strsys << std::fixed << std::setprecision(np) << sys;
		std::stringstream strlatex;
		if(n10==0) {
			strlatex << strval.str() << "\\pm" << strerr.str() << "\\pm" << strsys.str();
		} else {
			strlatex << "(" << strval.str() << "\\pm" << strerr.str() << "\\pm" << strsys.str() << ")\\times10^{" << n10 << "}";
		}
		TString tstrlatex = strlatex.str();

		return tstrlatex;
	}

	TString StrLatexAsym(int i=0) {
		if(!Def.CheckIndex(i)) return TString();
		if(!IsUpdated) Calc();

		double e_val, e_err;
		int e_n10, e_np;
		MyDualNumber e_dn(gStat.GetY()[i], gStat.GetEY()[i]);
		e_dn.StrLatex(e_val, e_err, e_n10, e_np);

		double sl_val, sl_err;
		int sl_n10, sl_np;
		MyDualNumber sl_dn(gAsym.GetY()[i], gAsym.GetEYlow()[i]);
		sl_dn.StrLatex(sl_val, sl_err, sl_n10, sl_np);

		double sh_val, sh_err;
		int sh_n10, sh_np;
		MyDualNumber sh_dn(gAsym.GetY()[i], gAsym.GetEYhigh()[i]);
		sh_dn.StrLatex(sh_val, sh_err, sh_n10, sh_np);

		double val, err, syl, syh;
		int n10, np;
		//if(sl_n10 < e_n10 && sl_n10 < sh_n10) {
		if(gAsym.GetEYlow()[i] < gStat.GetEY()[i] && gAsym.GetEYlow()[i] < gAsym.GetEYhigh()[i]) {
			val = sl_val;
			err = e_err * pow(10.0, 1.0*(e_n10-sl_n10));
			syl = sl_err;
			syh = sh_err * pow(10.0, 1.0*(sh_n10-sl_n10));
			n10 = sl_n10;
			np = sl_np;
		//} else if(sh_n10 < e_n10 && sh_n10 < sl_n10) {
		} else if(gAsym.GetEYhigh()[i] < gStat.GetEY()[i] && gAsym.GetEYhigh()[i] < gAsym.GetEYlow()[i]) {
			val = sh_val;
			err = e_err * pow(10.0, 1.0*(e_n10-sh_n10));
			syl = sl_err * pow(10.0, 1.0*(sl_n10-sh_n10));
			syh = sh_err;
			n10 = sh_n10;
			np = sh_np;
		} else {
			val = e_val;
			err = e_err;
			syl = sl_err * pow(10.0, 1.0*(sl_n10-e_n10));
			syh = sh_err * pow(10.0, 1.0*(sh_n10-e_n10));
			n10 = e_n10;
			np = e_np;
		}

		std::stringstream strval;
		std::stringstream strerr;
		std::stringstream strsyl;
		std::stringstream strsyh;
		strval << std::fixed << std::setprecision(np) << val;
		strerr << std::fixed << std::setprecision(np) << err;
		strsyl << std::fixed << std::setprecision(np) << syl;
		strsyh << std::fixed << std::setprecision(np) << syh;
		std::stringstream strlatex;
		if(n10==0) {
			strlatex << strval.str() << "\\pm" << strerr.str() << "_{-" << strsyl.str() << "}^{+" << strsyh.str() << "}";
		} else {
			strlatex << "(" << strval.str() << "\\pm" << strerr.str() << "_{-" << strsyl.str() << "}^{+" << strsyh.str() << "})\\times10^{" << n10 << "}";
		}
		TString tstrlatex = strlatex.str();

		return tstrlatex;
	}

	void ShiftX() {
		Def.ShiftX();
		if(!IsUpdated) Calc();
		for(int i=0; i<gStat.GetN(); i++) {
			gStat.GetX()[i] = Def.Graph.GetX()[i];
			gSyst.GetX()[i] = Def.Graph.GetX()[i];
			gAsym.GetX()[i] = Def.Graph.GetX()[i];
		}
	}

	void ShiftX(double dx) {
		Def.ShiftX();
		if(!IsUpdated) Calc();
		for(int i=0; i<gStat.GetN(); i++) {
			gStat.GetX()[i] = Def.Graph.GetX()[i] + dx;
			gSyst.GetX()[i] = Def.Graph.GetX()[i] + dx;
			gAsym.GetX()[i] = Def.Graph.GetX()[i] + dx;
		}
	}

	void MakePlot(TString name, TH2D* hFrame, TString StrPath="systplot") {
		TCanvas *cTmp = new TCanvas(name, name);
		hFrame->Draw();
		double xl = hFrame->GetXaxis()->GetXmin();
		double xh = hFrame->GetXaxis()->GetXmax();
		double yl = hFrame->GetYaxis()->GetXmin();
		double yh = hFrame->GetYaxis()->GetXmax();
		TLine *lineTmp;
		if(yl<0 && yh>0) { lineTmp = new TLine(xl,0, xh,0); lineTmp->SetLineColor(kGray); lineTmp->Draw("same"); }
		if(yl<1 && yh>1) { lineTmp = new TLine(xl,1, xh,1); lineTmp->SetLineColor(kGray); lineTmp->Draw("same"); }
		TGraphErrors gTmp = gStat;
		gTmp.SetMarkerStyle(20);
		gTmp.Draw("PL same");
		TLegend *lTmp = new TLegend(0.60,0.60,0.95,0.88);
		lTmp->SetFillStyle(0);
		lTmp->SetBorderSize(0);
		lTmp->AddEntry(&gTmp, "default", "lp");
		int i = 0;
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			const int &idx = it->first;
			MyPackGraph &pg = it->second;
			pg.g.Graph.SetMarkerStyle(24+i);
			pg.g.Graph.SetMarkerColor(i==3?kOrange+1:i+2);
			pg.g.Graph.SetLineColor(i==3?kOrange+1:i+2);
			pg.g.Graph.Draw("PX0 same");
			lTmp->AddEntry(&(pg.g.Graph), pg.s, "lp");
			i++;
		}
		lTmp->Draw("same");
		gPad->SetTicks();
		cTmp->SaveAs(StrPath+"/"+name+".pdf");
		delete cTmp;
	}

	void MakePlot(TString name="", TString StrPath="systplot") {
		TString tmpname = name;
		if(tmpname=="") tmpname = gStat.GetName();
		int n = gStat.GetN();
		double xl = TMath::MinElement(n, gStat.GetX());
		double xh = TMath::MaxElement(n, gStat.GetX());
		double yl = TMath::MinElement(n, gStat.GetY());
		double yh = TMath::MaxElement(n, gStat.GetY());
		double xw = xh - xl;
		double yw = yh - yl;
		xl = xl - 0.05*xw;
		xh = xh + 0.05*xw;
		yl = yl - 0.05*yw;
		yh = yh + 0.15*yw;
		xw = xh - xl;
		yw = yh - yl;
		TH2D *hTmpFrame = new TH2D(name+"Frame", name+"Frame", 100,xl,xh, 100,yl,yh);
		hTmpFrame->GetXaxis()->SetTitle(gStat.GetXaxis()->GetTitle());
		hTmpFrame->GetYaxis()->SetTitle(gStat.GetYaxis()->GetTitle());
		MakePlot(tmpname, hTmpFrame, StrPath);
		delete hTmpFrame;
	}

	const MySystGraph AveBin(std::vector<int> range) const {
		MySystGraph sg(*this);
		sg.SetDef(Def.AveBin(range));
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			int idx = it->first;
			MyPackGraph pg = it->second;
			pg.g = pg.g.AveBin(range);
			sg.SetVar(idx, pg);
		}
		return sg;
	}

	const MySystGraph AveBin(int bl, int bh) const {
		MySystGraph sg(*this);
		sg.SetDef(Def.AveBin(bl, bh));
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			int idx = it->first;
			MyPackGraph pg = it->second;
			pg.g = pg.g.AveBin(bl, bh);
			sg.SetVar(idx, pg);
		}
		return sg;
	}

	void Write(TFile *f, TString name) const {
		f->WriteObjectAny(this, "MySystGraph", name);
	}
};


const MySystGraph operator+(const MySystGraph &sg1, double c) {
	MySystGraph sg2(sg1);
	sg2.SetDef(sg1.GetDef() + c);
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = (pg.g + c);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph operator+(double c, const MySystGraph &sg1) {
	return (sg1 + c);
}

const MySystGraph operator*(const MySystGraph &sg1, double c) {
	MySystGraph sg2(sg1);
	sg2.SetDef(sg1.GetDef() * c);
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = (pg.g * c);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph operator*(double c, const MySystGraph &sg1) {
	return (sg1 * c);
}

const MySystGraph operator-(const MySystGraph &sg1, double c) {
	return (sg1 + (-c));
}

const MySystGraph operator-(double c, const MySystGraph &sg1) {
	return (-1.0*sg1 + c);
}

const MySystGraph operator/(const MySystGraph &sg1, double c) {
	return (sg1 * (1.0/c));
}

const MySystGraph operator/(double c, const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(c / sg1.GetDef());
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = (c / pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph operator+(const MySystGraph &sg1) {
	return sg1;
}

const MySystGraph operator-(const MySystGraph &sg1) {
	return (0.0 - sg1);
}

const MySystGraph operator+(const MySystGraph &sg1, const MySystGraph &sg2) {
	if(sg1.GetMode() != sg2.GetMode()) std::cout << "MySystGraph: Mode not match!" << std::endl;
	MySystGraph sg3(sg1);
	sg3.SetDef(sg1.GetDef() + sg2.GetDef());
	std::vector<int> list1 = sg1.GetList();
	std::vector<int> list2 = sg2.GetList();
	int i = 0;
	int j = 0;
	while(i<(int)list1.size() && j<(int)list2.size()) {
		if(list1[i] < list2[j]) {
			int idx = list1[i];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = (pg.g + sg2.GetDef());
			sg3.SetVar(idx, pg);
			i++;
		} else if(list1[i] > list2[j]) {
			int idx = list2[j];
			MyPackGraph pg = sg2.GetVar(idx);
			pg.g = (sg1.GetDef() + pg.g);
			sg3.SetVar(idx, pg);
			j++;
		} else {
			int idx = list1[i]; // = list2[j];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = (pg.g + sg2.GetVar(idx).g);
			sg3.SetVar(idx, pg);
			i++;
			j++;
		}
	}
	while(i<(int)list1.size()) {
		int idx = list1[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = (pg.g + sg2.GetDef());
		sg3.SetVar(idx, pg);
		i++;
	}
	while(j<(int)list2.size()) {
		int idx = list2[j];
		MyPackGraph pg = sg2.GetVar(idx);
		pg.g = (sg1.GetDef() + pg.g);
		sg3.SetVar(idx, pg);
		j++;
	}
	return sg3;
}

const MySystGraph operator-(const MySystGraph &sg1, const MySystGraph &sg2) {
	return (sg1 + (-1.0*sg2));
}

const MySystGraph operator*(const MySystGraph &sg1, const MySystGraph &sg2) {
	if(sg1.GetMode() != sg2.GetMode()) std::cout << "MySystGraph: Mode not match!" << std::endl;
	MySystGraph sg3(sg1);
	sg3.SetDef(sg1.GetDef() * sg2.GetDef());
	std::vector<int> list1 = sg1.GetList();
	std::vector<int> list2 = sg2.GetList();
	int i = 0;
	int j = 0;
	while(i<(int)list1.size() && j<(int)list2.size()) {
		if(list1[i] < list2[j]) {
			int idx = list1[i];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = (pg.g * sg2.GetDef());
			sg3.SetVar(idx, pg);
			i++;
		} else if(list1[i] > list2[j]) {
			int idx = list2[j];
			MyPackGraph pg = sg2.GetVar(idx);
			pg.g = (sg1.GetDef() * pg.g);
			sg3.SetVar(idx, pg);
			j++;
		} else {
			int idx = list1[i]; // = list2[j];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = (pg.g * sg2.GetVar(idx).g);
			sg3.SetVar(idx, pg);
			i++;
			j++;
		}
	}
	while(i<(int)list1.size()) {
		int idx = list1[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = (pg.g * sg2.GetDef());
		sg3.SetVar(idx, pg);
		i++;
	}
	while(j<(int)list2.size()) {
		int idx = list2[j];
		MyPackGraph pg = sg2.GetVar(idx);
		pg.g = (sg1.GetDef() * pg.g);
		sg3.SetVar(idx, pg);
		j++;
	}
	return sg3;
}

const MySystGraph operator/(const MySystGraph &sg1, const MySystGraph &sg2) {
	return (sg1 * (1.0/sg2));
}

const MySystGraph pow(const MySystGraph &sg1, double c) {
	MySystGraph sg2(sg1);
	sg2.SetDef(pow(sg1.GetDef(), c));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = pow(pg.g, c);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph sqrt(const MySystGraph &sg1) {
	return pow(sg1, 0.5);
}

const MySystGraph abs(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(abs(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = abs(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph exp(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(exp(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = exp(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph log(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(log(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = log(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph sin(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(sin(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = sin(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph cos(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(cos(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = cos(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph tan(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(tan(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = tan(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}

const MySystGraph atan(const MySystGraph &sg1) {
	MySystGraph sg2(sg1);
	sg2.SetDef(atan(sg1.GetDef()));
	std::vector<int> list = sg1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int idx = list[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = atan(pg.g);
		sg2.SetVar(idx, pg);
	}
	return sg2;
}


#endif
