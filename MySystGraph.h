/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 * Note: dual number to calculate precise values of first-order derivative
 *       of multi-variable functions with common math formula
 *       this class is wrapped with ROOT for I/O and plotting
 *       this class can also calculate systematic uncertainties
 *       with Barlow's prescription: https://arxiv.org/abs/hep-ex/0207026 
 **************************************************************************/

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
	int m; // mode: 1: s^2; 2: s^2-(e^2-e^2); 3: s^2-(e^2+e^2)

	MyPackGraph() {
		w = 0;
		m = 1;
	}
	MyPackGraph(const MyDualGraph &gg, TString ss="", double ww=1.0, int mm=1) {
		g = gg;
		s = ss;
		w = ww;
		if(mm==1 || mm==2 || mm==3) {
			m = mm;
		} else {
			m = 1;
		}
	}
	MyPackGraph(const TGraphErrors &gg, TString ss="", double ww=1.0, int mm=1) {
		g = MyDualGraph(gg);
		s = ss;
		w = ww;
		if(mm==1 || mm==2 || mm==3) {
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

	vector< vector<double> > syVar;  // syst^2 for each variation
	vector< vector<double> > sylVar; // syst^2 for each variation
	vector< vector<double> > syhVar; // syst^2 for each variation

	// systematic box width
	double Width;

	// is graphs updated
	bool IsUpdated;
	bool IsHoldForDraw;

public:
	// graphs for plot
	TGraphErrors gStat;
	TGraphErrors gSyst;
	TGraphErrors gComb;
	TGraphAsymmErrors gAsym;

	TGraphErrors GetStat() { CalcIfNotUpdated(); return gStat;}
	TGraphErrors GetSyst() { CalcIfNotUpdated(); return gSyst;}
	TGraphErrors GetComb() { CalcIfNotUpdated(); return gComb;}
	TGraphAsymmErrors GetAsym() { CalcIfNotUpdated(); return gAsym; }

	static int VidNow; // the latest index
	static const int VidMax = 1000000; // maximal index

	bool GetIsUpdated() const {
		return IsUpdated;
	}

	void Reset() {
		Var.clear();
		Width = 0.8;
		IsUpdated = false;
		IsHoldForDraw = false;
	}

	MySystGraph() {
		Reset();
		SetMode(0);
	}
	MySystGraph(const MyDualGraph &def, int mode=0) {
		Reset();
		SetMode(mode);
		Def = def;
		if(Def.GetN()>=2) {
			//Width = 0.1*fabs(def.Graph.GetX()[0] - def.Graph.GetX()[1]);
			Width = 0.1*fabs(def.GetPoint(0).Px.GetValu() - def.GetPoint(1).Px.GetValu());
		}
	}
	MySystGraph(const TGraphErrors &def, int mode=0) {
		Reset();
		SetMode(mode);
		Def = MyDualGraph(def);
		if(Def.GetN()>=2) {
			//Width = 0.1*fabs(def.Graph.GetX()[0] - def.Graph.GetX()[1]);
			Width = 0.1*fabs(Def.GetPoint(0).Px.GetValu() - Def.GetPoint(1).Px.GetValu());
		}
	}
	MySystGraph(const MySystGraph &sg1) {
		//Reset();
		//Mode = sg1.Mode;
		//Def = sg1.Def;
		//Var = sg1.Var;
		//Width = sg1.Width;
		*this = sg1;
		if(!IsUpdated) Calc();
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
		for(int i=0; i<gSyst.GetN(); i++) {
			gSyst.GetEX()[i] = Width;
			gComb.GetEX()[i] = Width;
		}
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

	static int AutoNewVid() { // VidNow range 1 ~ VidMax-1
		if(VidNow>=VidMax-1) {
			std::cout << "MySystGraph::VidNow out of range! " << VidNow << "/" << VidMax-1 << std::endl;
		} else {
			VidNow ++;
		}
		return VidNow;
	}

	bool CheckVar(int idx) const {
		return (Var.count(idx)!=0);
	}
	void AddVar(int idx, const MyPackGraph &pg) {
		if(Var.count(idx)==0) {
			Var.insert({idx, pg});
		} else {
			std::cout << "MySystGraph::AddVar() index " << idx << " already exists!" << std::endl;
		}
		IsUpdated = false;
	}
	void AddVar(const MyPackGraph &pg) {
		AutoNewVid();
		AddVar(VidNow, pg);
	}
	void AddUnc(int idxl, int idxh, const TGraphErrors &g, TString ss="input") { // only for DrawSyst()
		const int n = Def.GetN();
		Def.Calc();
		if(g.GetN() != n) std::cout << "MySystGraph::AddUnc() size does not match!" << std::endl;
		double tmpyl[n];
		double tmpyh[n];
		for(int i=0; i<n; i++) {
			tmpyl[i] = Def.Graph.GetY()[i] - g.GetEY()[i];
			tmpyh[i] = Def.Graph.GetY()[i] + g.GetEY()[i];
		}
		TGraphErrors gl(n, Def.Graph.GetX(), tmpyl, Def.Graph.GetEX(), Def.Graph.GetEY());
		TGraphErrors gh(n, Def.Graph.GetX(), tmpyh, Def.Graph.GetEX(), Def.Graph.GetEY());
		MyPackGraph pgl(gl, ss+" d-s", 0.5, 1);
		MyPackGraph pgh(gh, ss+" d+s", 0.5, 1);
		AddVar(idxl, pgl);
		AddVar(idxh, pgh);
	} 
	void AddUnc(const TGraphErrors &g, TString ss="input") { // only for DrawSyst()
		int idxl = AutoNewVid();
		int idxh = AutoNewVid();
		AddUnc(idxl, idxh, g, ss);
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
			//Var.insert({idx, pg});
			AddVar(idx, pg);
		} else {
			//Var[idx] = pg;
			ChgVar(idx, pg);
		}
		IsUpdated = false;
	}
	void EraseVar(int idx) {
		if(Var.count(idx)==0) {
			//std::cout << "MySystGraph::ChgVar() index " << idx << " does not exist!" << std::endl;
		} else {
			Var.erase(idx);
		}
		IsUpdated = false;
	}
	void EraseVar(std::vector<int> vidx) {
		for(int i=0; i<vidx.size(); i++) {
			EraseVar(vidx[i]);
		}
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

	void Calc(std::vector<int> vidlist) {
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
		syVar.clear();
		sylVar.clear();
		syhVar.clear();
//		for(auto it=Var.begin(); it!=Var.end(); it++) {
//			const int &idx = it->first;
//			MyPackGraph &pg = it->second;
		for(int j=0; j<vidlist.size(); j++) {
			int idx = vidlist[j];
			if(!CheckVar(idx)) continue;
			MyPackGraph pg = Var[idx];
			if(!pg.g.GetIsUpdated()) pg.g.Calc();
			if(!pg.g.CheckSize(n)) return;
			vector<double> sylVarOne;
			vector<double> syhVarOne;
			vector<double> syVarOne;
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || mm==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				syVarOne.push_back(w*ss);
				if(d<0) {
					sylVarOne.push_back(w*ss);
					syhVarOne.push_back(0);
				} else {
					sylVarOne.push_back(0);
					syhVarOne.push_back(w*ss);
				}
			}
			syVar.push_back(syVarOne);
			sylVar.push_back(sylVarOne);
			syhVar.push_back(syhVarOne);
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

		double syc[n];
		for(int i=0; i<n; i++) {
			syc[i] = sqrt(pow(gStat.GetEY()[i],2) + pow(gSyst.GetEY()[i],2));
		}
		gComb = TGraphErrors(n, gStat.GetX(), gStat.GetY(), sx, syc);

		IsUpdated = true;
	}

	void Calc() {
		Calc(GetList());
	}

	void CalcIfNotUpdated() {
		if(!IsUpdated) Calc();
	}

	TString StrLatex(int i=0, std::string symboltype="R") {
		if(!Def.CheckIndex(i)) return TString();
		if(!IsUpdated) Calc();
		MyValuErrs ve(gStat.GetY()[i], gStat.GetEY()[i], gSyst.GetEY()[i]);
		return TString(ve.PrintErr2(symboltype));
	}

	MyValuErrs ValuStatSyst(int i=0) {
		if(!Def.CheckIndex(i)) return MyValuErrs();
		if(!IsUpdated) Calc();
		MyValuErrs ve(gStat.GetY()[i], gStat.GetEY()[i], gSyst.GetEY()[i]);
		return ve;
	}

	TString StrLatexAsym(int i=0, std::string symboltype="R") {
		if(!Def.CheckIndex(i)) return TString();
		if(!IsUpdated) Calc();
		MyValuErrs ve(gStat.GetY()[i], gStat.GetEY()[i], gAsym.GetEYlow()[i], gAsym.GetEYhigh()[i]);
		return TString(ve.PrintErr3(symboltype));
	}

	MyValuErrs ValuStatAsym(int i=0) {
		if(!Def.CheckIndex(i)) return MyValuErrs();
		if(!IsUpdated) Calc();
		MyValuErrs ve(gStat.GetY()[i], gStat.GetEY()[i], gAsym.GetEYlow()[i], gAsym.GetEYhigh()[i]);
		return ve;
	}

	void ShiftX() {
		Def.ShiftX();
		if(!IsUpdated) Calc();
		for(int i=0; i<gStat.GetN(); i++) {
			gStat.GetX()[i] = Def.Graph.GetX()[i];
			gSyst.GetX()[i] = Def.Graph.GetX()[i];
			gComb.GetX()[i] = Def.Graph.GetX()[i];
			gAsym.GetX()[i] = Def.Graph.GetX()[i];
		}
	}

	void ShiftX(double dx) {
		Def.ShiftX();
		if(!IsUpdated) Calc();
		for(int i=0; i<gStat.GetN(); i++) {
			gStat.GetX()[i] = Def.Graph.GetX()[i] + dx;
			gSyst.GetX()[i] = Def.Graph.GetX()[i] + dx;
			gComb.GetX()[i] = Def.Graph.GetX()[i] + dx;
			gAsym.GetX()[i] = Def.Graph.GetX()[i] + dx;
		}
	}

	void SetStyleColor(Int_t style, Int_t color) {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gSyst.SetFillStyle(0);
		gComb.SetFillStyle(0);
		gAsym.SetFillStyle(0);
		gSyst.SetFillColor(color);
		gComb.SetFillColor(color);
		gAsym.SetFillColor(color);
		gSyst.SetMarkerStyle(style);
		gComb.SetMarkerStyle(style);
		gAsym.SetMarkerStyle(style);
		gStat.SetMarkerStyle(style);
		gSyst.SetMarkerColor(color);
		gComb.SetMarkerColor(color);
		gAsym.SetMarkerColor(color);
		gStat.SetMarkerColor(color);
		gSyst.SetLineColor(color);
		gComb.SetLineColor(color);
		gAsym.SetLineColor(color);
		gStat.SetLineColor(color);
	}

	void SetXaxisRange(double xl, double xh) {
		if(!IsUpdated) Calc();
		gSyst.GetXaxis()->SetLimits(xl, xh);
		gComb.GetXaxis()->SetLimits(xl, xh);
		gAsym.GetXaxis()->SetLimits(xl, xh);
		gStat.GetXaxis()->SetLimits(xl, xh);
	}

	void SetXaxisLimit(double xl, double xh) {
		if(!IsUpdated) Calc();
		vector<double> vx, vy, ex, ey, exl, exh, eyl, eyh;

		vx.clear();
		vy.clear();
		ex.clear();
		ey.clear();
		for(int i=0; i<gSyst.GetN(); i++) {
			if(gSyst.GetX()[i]<xl || gSyst.GetX()[i]>=xh) continue;
			vx.push_back(gSyst.GetX()[i]);
			vy.push_back(gSyst.GetY()[i]);
			ex.push_back(gSyst.GetEX()[i]);
			ey.push_back(gSyst.GetEY()[i]);
		}
		gSyst = TGraphErrors(vx.size(), &vx[0], &vy[0], &ex[0], &ey[0]);

		vx.clear();
		vy.clear();
		ex.clear();
		ey.clear();
		for(int i=0; i<gComb.GetN(); i++) {
			if(gComb.GetX()[i]<xl || gComb.GetX()[i]>=xh) continue;
			vx.push_back(gComb.GetX()[i]);
			vy.push_back(gComb.GetY()[i]);
			ex.push_back(gComb.GetEX()[i]);
			ey.push_back(gComb.GetEY()[i]);
		}
		gComb = TGraphErrors(vx.size(), &vx[0], &vy[0], &ex[0], &ey[0]);

		vx.clear();
		vy.clear();
		ex.clear();
		ey.clear();
		for(int i=0; i<gStat.GetN(); i++) {
			if(gStat.GetX()[i]<xl || gStat.GetX()[i]>=xh) continue;
			vx.push_back(gStat.GetX()[i]);
			vy.push_back(gStat.GetY()[i]);
			ex.push_back(gStat.GetEX()[i]);
			ey.push_back(gStat.GetEY()[i]);
		}
		gStat = TGraphErrors(vx.size(), &vx[0], &vy[0], &ex[0], &ey[0]);

		vx.clear();
		vy.clear();
		exl.clear();
		exh.clear();
		eyl.clear();
		eyh.clear();
		for(int i=0; i<gAsym.GetN(); i++) {
			if(gAsym.GetX()[i]<xl || gAsym.GetX()[i]>=xh) continue;
			vx.push_back(gAsym.GetX()[i]);
			vy.push_back(gAsym.GetY()[i]);
			exl.push_back(gAsym.GetEXlow ()[i]);
			exh.push_back(gAsym.GetEXhigh()[i]);
			eyl.push_back(gAsym.GetEYlow ()[i]);
			eyh.push_back(gAsym.GetEYhigh()[i]);
		}
		gAsym = TGraphAsymmErrors(vx.size(), &vx[0], &vy[0], &exl[0],&exh[0], &eyl[0],&eyh[0]);

		IsUpdated = false;
		IsHoldForDraw = true;
	}

	void DrawStat(TString opt = "PL same") {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		//gStat.Draw("PL"+opt+" same");
		gStat.Draw(opt);
	}

	void DrawComb(TString opt = "PL same") {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gComb.Draw(opt);
	}

	void DrawSyst(TString opt = "") {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gSyst.Draw("P5"+opt+" same");
		gStat.Draw("PL"+opt+" same");
	}

	void DrawSyst(TString optstat, TString optsyst) {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gSyst.Draw(optsyst);
		gStat.Draw(optstat);
	}

	void DrawAsym(TString opt = "") {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gAsym.Draw("P5"+opt+" same");
		gStat.Draw("PL"+opt+" same");
	}

	void DrawAsym(TString optstat, TString optsyst) {
		if(!IsUpdated && !IsHoldForDraw) Calc();
		gAsym.Draw(optsyst);
		gStat.Draw(optstat);
	}

	void MakePlot(TString name, TH2D* hFrame, TString StrPath="systplot", vector<TString> FileType = vector<TString>{"pdf"}) {
		ShiftX();

		int vcol[10] = {2, 3, 4, kOrange, 6, 7, 8, 9, 16, 28};
		int vmrk[10] = {20, 24, 21, 25, 22, 26, 23, 32, 29, 30};
		// support 100 variations

		int nvar = Var.size();
		int ncln = nvar / 20 + 1;
		if(nvar==0) ncln = 0;
		int wpxl = 700 + ncln*175;
		int hpxl = 500;
		TCanvas *cTmp = new TCanvas("SystPlot"+name, "SystPlot"+name, wpxl,hpxl);
		gPad->SetLeftMargin(0.18*700/wpxl);
		gPad->SetRightMargin((0.04*700+ncln*175)/wpxl);
		gPad->SetTopMargin(0.10);
		gPad->SetBottomMargin(0.16);

		double xl = hFrame->GetXaxis()->GetXmin();
		double xh = hFrame->GetXaxis()->GetXmax();
		double yl = hFrame->GetYaxis()->GetXmin();
		double yh = hFrame->GetYaxis()->GetXmax();
		TH2D *hFrameTmp = new TH2D("SystPlotFrame"+name, "SystPlotFrame"+name, 100,xl,xh, 100,yl,yh);
		hFrameTmp->GetXaxis()->SetTitle(hFrame->GetXaxis()->GetTitle());
		hFrameTmp->GetYaxis()->SetTitle(hFrame->GetYaxis()->GetTitle());
		hFrameTmp->GetYaxis()->SetTitleOffset(1.0/wpxl*700);
		hFrameTmp->Draw();

		TLine *lineTmp;
		if(yl<0 && yh>0) { lineTmp = new TLine(xl,0, xh,0); lineTmp->SetLineColor(kGray); lineTmp->Draw("same"); }
		if(yl<1 && yh>1) { lineTmp = new TLine(xl,1, xh,1); lineTmp->SetLineColor(kGray); lineTmp->Draw("same"); }

		TGraphErrors gTmp = gStat;
		gTmp.SetLineColor(kBlack);
		gTmp.SetMarkerColor(kBlack);
		gTmp.SetMarkerStyle(20);
		gTmp.Draw("PL same");

		double bw = fabs(xh-xl) / (gTmp.GetN()+1.0);
		for(int i=1; i<gTmp.GetN(); i++) {
			double tmpbw = fabs(gTmp.GetX()[i]-gTmp.GetX()[i-1]);
			if(i==1) {
				bw = tmpbw;
			} else {
				bw = bw<tmpbw?bw:tmpbw;
			}
		}

		TLegend *lTmp = new TLegend(700.0/wpxl,0.01, 1-0.01*hpxl/wpxl,0.99);
		lTmp->SetFillStyle(0);
		lTmp->SetBorderSize(1);
		lTmp->SetNColumns(ncln);
		lTmp->SetTextSize(0.04);
		lTmp->AddEntry(&gTmp, "default", "lp");

		int nsx = (int)Var.size() + 2;
		int isx = 1;
		int i = 0;
		double sw = 0.06*nsx;
		if(sw<0.6) sw = 0.6;
		if(sw>1.0) sw = 1.0;
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			const int &idx = it->first;
			MyPackGraph &pg = it->second;
			TGraphErrors &gVar = pg.g.Graph;
			for(int j=0; j<gTmp.GetN(); j++) gVar.GetX()[j] += sw*bw*isx/nsx;
			isx++;
			int col = vcol[i%10];
			int imrk = i / 10;
			while(imrk>=10) imrk /= 10;
			int mrk = vmrk[imrk];
			gVar.SetMarkerStyle(mrk);
			gVar.SetMarkerColor(col);
			gVar.SetLineColor(col);
			gVar.Draw("P same");
			lTmp->AddEntry(&(pg.g.Graph), pg.s, "lp");
			i++;
		}
		lTmp->Draw("same");
		gPad->SetTicks();
		for(int i=0; i<FileType.size(); i++) cTmp->SaveAs(StrPath+"/Syst"+name+"."+FileType[i]);
		delete cTmp;
	}

	void MakePlot(TString name="", TString StrPath="systplot", vector<TString> FileType = vector<TString>{"pdf"}) {
		TString tmpname = name;
		if(tmpname=="") tmpname = gStat.GetName();
		int n = gStat.GetN();
		double xl = TMath::MinElement(n, gStat.GetX());
		double xh = TMath::MaxElement(n, gStat.GetX());
		double yl = TMath::MinElement(n, gStat.GetY());
		double yh = TMath::MaxElement(n, gStat.GetY());
		double xw = xh - xl;
		double yw = yh - yl;
		if(n==1) {
			xw = 10;
			yw = 40*gStat.GetEY()[0];
		}
		xl = xl - 0.05*xw;
		xh = xh + 0.05*xw;
		yl = yl - 0.05*yw;
		yh = yh + 0.15*yw;
		xw = xh - xl;
		yw = yh - yl;
		TH2D *hTmpFrame = new TH2D("FramePlot"+name, "FramePlot"+name, 100,xl,xh, 100,yl,yh);
		hTmpFrame->GetXaxis()->SetTitle(gStat.GetXaxis()->GetTitle());
		hTmpFrame->GetYaxis()->SetTitle(gStat.GetYaxis()->GetTitle());
		MakePlot(tmpname, hTmpFrame, StrPath, FileType);
		delete hTmpFrame;
	}

	void MakeHist(TString name, TH2D* hFrame, bool ispercent=false, TString StrPath="systplot", vector<TString> FileType = vector<TString>{"pdf"}) {
		ShiftX();

		if(ispercent) name = name+"Percent";

		int vcol[10] = {2, 3, 4, kOrange, 6, 7, 8, 9, 16, 28};
		int vfls[10] = {1001, 3001, 3002, 3004, 3244, 3003, 3344, 3345, 3444, 3354};
		// support 100 variations

		int nvar = Var.size();
		int ncln = nvar / 20 + 1;
		int wpxl = 700 + ncln*175;
		int hpxl = 500;
		TCanvas *cTmp = new TCanvas("SystHist"+name, "SystHist"+name, wpxl,hpxl);
		gPad->SetLeftMargin(0.18*700/wpxl);
		gPad->SetRightMargin((0.04*700+ncln*175)/wpxl);
		gPad->SetTopMargin(0.10);
		gPad->SetBottomMargin(0.16);

		double xl = hFrame->GetXaxis()->GetXmin();
		double xh = hFrame->GetXaxis()->GetXmax();
		TString xtitle = TString(hFrame->GetXaxis()->GetTitle());
		TString ytitle = TString(hFrame->GetYaxis()->GetTitle()) + Form(": syst^{2}");
		if(ispercent) ytitle = ytitle + Form(" (%%)");

		TGraphErrors gTmp = gStat;
		gTmp.SetLineColor(kBlack);
		gTmp.SetMarkerColor(kBlack);
		gTmp.SetMarkerStyle(20);
		gTmp.Draw("PL same");

		int nBin = gTmp.GetN();
		if(nBin<=0) return;
		double bw = fabs(xh-xl) / (nBin+1.0);
		for(int i=1; i<nBin; i++) {
			double tmpbw = fabs(gTmp.GetX()[i]-gTmp.GetX()[i-1]);
			if(i==1) {
				bw = tmpbw;
			} else {
				bw = bw<tmpbw?bw:tmpbw;
			}
		}
		bool isgoodxedge = true;
		vector<double> xcenter;
		for(int i=0; i<nBin; i++) xcenter.push_back(gTmp.GetX()[i]);
		std::sort(xcenter.begin(), xcenter.end());
		double xedge[nBin+1];
		xedge[0] = xl;
		for(int i=0; i<nBin; i++) {
			xedge[i+1] = xedge[i] + (xcenter[i] - xedge[i])*2;
			if(xedge[i+1] <= xedge[i]) isgoodxedge = false;
		}
		if(!isgoodxedge) {
			xedge[0] = xcenter[0] - 0.5*bw;
			for(int i=0; i<nBin; i++) xedge[i+1] = xcenter[i] + 0.5*bw;
		}

		TLegend *lTmp = new TLegend(700.0/wpxl,0.01, 1-0.01*hpxl/wpxl,0.99);
		lTmp->SetFillStyle(0);
		lTmp->SetBorderSize(1);
		lTmp->SetNColumns(ncln);
		lTmp->SetTextSize(0.04);

		std::vector<int> list = GetList();
		TBox syBox[nBin][syVar.size()];
		double syMax = 0;
		for(int j=0; j<nBin; j++) {
			double syBin = 0;
			//for(int i=0; i<syVar.size(); i++) {
			for(int i=syVar.size()-1; i>=0; i--) {
				double tmpvalue = syVar[i][j];
				double ttlvalue = pow(gSyst.GetEY()[j],2);
				if(ispercent && ttlvalue!=0) tmpvalue = 100 * tmpvalue / ttlvalue;
				syBox[j][i] = TBox(xedge[j],syBin, xedge[j+1],syBin+tmpvalue);
				syBin += tmpvalue;
				int col = vcol[i%10];
				int ifls = i / 10;
				while(ifls>=10) ifls /= 10;
				int fls = vfls[ifls];
				syBox[j][i].SetLineColor(col);
				syBox[j][i].SetFillStyle(fls);
				syBox[j][i].SetFillColor(col);
			}
			if(syBin>syMax) syMax = syBin;
		}

		TH2D *hFrameTmp = new TH2D("SystHistFrame"+name, "SystHistFrame"+name, 100,xl,xh, 100,0,syMax);
		hFrameTmp->GetXaxis()->SetTitle(xtitle);
		hFrameTmp->GetYaxis()->SetTitle(ytitle);
		hFrameTmp->GetYaxis()->SetTitleOffset(1.0/wpxl*700);
		hFrameTmp->Draw();

		for(int i=0; i<syVar.size(); i++) lTmp->AddEntry(&(syBox[0][i]), Var[list[i]].s, "f");
		for(int j=0; j<nBin; j++) {
			for(int i=0; i<syVar.size(); i++) {
				syBox[j][i].Draw();
			}
		}

		lTmp->Draw("same");

		gPad->SetTicks();
		for(int i=0; i<FileType.size(); i++) cTmp->SaveAs(StrPath+"/Syst"+name+"Hist."+FileType[i]);
		delete cTmp;
	}

	void MakeHist(TString name="", bool ispercent=false, TString StrPath="systplot", vector<TString> FileType = vector<TString>{"pdf"}) {
		TString tmpname = name;
		if(tmpname=="") tmpname = gStat.GetName();
		int n = gStat.GetN();
		double xl = TMath::MinElement(n, gStat.GetX());
		double xh = TMath::MaxElement(n, gStat.GetX());
		double yl = TMath::MinElement(n, gStat.GetY());
		double yh = TMath::MaxElement(n, gStat.GetY());
		double xw = xh - xl;
		double yw = yh - yl;
		if(n==1) {
			xw = 10;
			yw = 40*gStat.GetEY()[0];
		}
		xl = xl - 0.05*xw;
		xh = xh + 0.05*xw;
		yl = yl - 0.05*yw;
		yh = yh + 0.15*yw;
		xw = xh - xl;
		yw = yh - yl;
		TH2D *hTmpFrame = new TH2D("FramePlot"+name, "FramePlot"+name, 100,xl,xh, 100,yl,yh);
		hTmpFrame->GetXaxis()->SetTitle(gStat.GetXaxis()->GetTitle());
		hTmpFrame->GetYaxis()->SetTitle(gStat.GetYaxis()->GetTitle());
		MakeHist(tmpname, hTmpFrame, StrPath, ispercent, FileType);
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

	const MySystGraph SelectBin(std::vector<int> range) const {
		MySystGraph sg(*this);
		sg.SetDef(Def.SelectBin(range));
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			int idx = it->first;
			MyPackGraph pg = it->second;
			pg.g = pg.g.SelectBin(range);
			sg.SetVar(idx, pg);
		}
		return sg;
	}

	const MySystGraph SelectBin(int bl, int bh) const {
		MySystGraph sg(*this);
		sg.SetDef(Def.SelectBin(bl, bh));
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			int idx = it->first;
			MyPackGraph pg = it->second;
			pg.g = pg.g.SelectBin(bl, bh);
			sg.SetVar(idx, pg);
		}
		return sg;
	}

	const MySystGraph SelectBin(int bl) const {
		return SelectBin(bl, bl);
	}

	const MySystGraph MergeUnce() {
		MySystGraph sg(*this);
		sg.SetDef(Def.MergeUnce());
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			int idx = it->first;
			MyPackGraph pg = it->second;
			pg.g = pg.g.MergeUnce();
			sg.SetVar(idx, pg);
		}
		return sg;
	}

	static TString Sn(int n) {
		TString s;
		for(int i=0; i<n; i++) s = s + ' ';
		return s;
	}

	static TString Sn(TString ss, int n=16) {
		TString s;
		for(int i=0; i<n; i++) {
			if(i<ss.Length()) {
				s = s + ss[i];
			} else {
				s = s + ' ';
			}
		}
		return s;
	}

	static TString Sn(double a, int n=16) {
		TString ss = Form("%f", a);
		return Sn(ss, n);
	}

	void PrintAsym(int width = 176) {
		if(!IsUpdated) Calc();
		TString header = Sn("x-value") + Sn("y-value") + Sn("y-error");
		TString weight = Sn("weight") + Sn("") + Sn("");
		int n = Def.GetN();
		vector<TString> number(n);
		for(int i=0; i<n; i++) {
			number[i] = Sn(Def.GetPoint(i).Px.GetValu()) + Sn(Def.GetPoint(i).Py.GetValu()) + Sn(Def.GetPoint(i).Py.GetUnce());
		}
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
			header = header + Sn(pg.s);
			weight = weight + Sn(pg.w);
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || mm==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				int sign = d<0?-1:1;
				number[i] = number[i] + Sn(sign*sqrt(fabs(ss)));
			}
		}
		header = header + Sn("total", 32);
		weight = weight + Sn("", 32);
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
			number[i] = number[i] + Sn(Form("-%f+%f", syl[i], syh[i]), 32);
		}
		int nLine = header.Length()/width + (header.Length()%width==0?0:1);
		for(int iLine=0; iLine<nLine; iLine++) {
			TString tmph;
			for(int j=0; j<width && iLine*width+j<header.Length(); j++) tmph = tmph + header[iLine*width+j];
			std::cout << tmph << std::endl;
			TString tmpw;
			for(int j=0; j<width && iLine*width+j<weight.Length(); j++) tmpw = tmpw + weight[iLine*width+j];
			std::cout << tmpw << std::endl;
			for(int i=0; i<n; i++) {
				TString tmpn;
				for(int j=0; j<width && iLine*width+j<number[i].Length(); j++) tmpn = tmpn + number[i][iLine*width+j];
				std::cout << tmpn << std::endl;
			}
		}
	}

	void PrintSyst(int width = 176) {
		if(!IsUpdated) Calc();
		TString header = Sn("x-value") + Sn("y-value") + Sn("y-error");
		TString weight = Sn("weight") + Sn("") + Sn("");
		int n = Def.GetN();
		vector<TString> number(n);
		for(int i=0; i<n; i++) {
			number[i] = Sn(Def.GetPoint(i).Px.GetValu()) + Sn(Def.GetPoint(i).Py.GetValu()) + Sn(Def.GetPoint(i).Py.GetUnce());
		}
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
			header = header + Sn(pg.s);
			weight = weight + Sn(pg.w);
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || mm==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				int sign = d<0?-1:1;
				number[i] = number[i] + Sn(sqrt(fabs(ss)));
			}
		}
		header = header + Sn("total");
		weight = weight + Sn("");
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
			number[i] = number[i] + Sn(sy[i]);
		}
		int nLine = header.Length()/width + (header.Length()%width==0?0:1);
		for(int iLine=0; iLine<nLine; iLine++) {
			TString tmph;
			for(int j=0; j<width && iLine*width+j<header.Length(); j++) tmph = tmph + header[iLine*width+j];
			std::cout << tmph << std::endl;
			TString tmpw;
			for(int j=0; j<width && iLine*width+j<weight.Length(); j++) tmpw = tmpw + weight[iLine*width+j];
			std::cout << tmpw << std::endl;
			for(int i=0; i<n; i++) {
				TString tmpn;
				for(int j=0; j<width && iLine*width+j<number[i].Length(); j++) tmpn = tmpn + number[i][iLine*width+j];
				std::cout << tmpn << std::endl;
			}
		}
	}

	void Print(int width=176) {
		PrintAsym(width);
	}

	void Table(bool isHeader=true) {
		if(!IsUpdated) Calc();
		TString tn = " & ";
		TString header = Sn("x-value") + tn + Sn("y-value") + tn + Sn("y-error");
		TString weight = Sn("weight") + tn + Sn("") + tn + Sn("");
		int n = Def.GetN();
		vector<TString> number(n);
		for(int i=0; i<n; i++) {
			double intx;
			number[i] = (std::modf(Def.GetPoint(i).Px.GetValu()*10,&intx)==0.0?Sn(Form("%.2g",intx/10)):Sn(Def.GetPoint(i).Px.GetValu())) + tn + Sn(Def.GetPoint(i).Py.GetValu()) + tn + Sn(Form("%.2g",Def.GetPoint(i).Py.GetUnce()));
		}
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
			header = header + tn + Sn(pg.s);
			weight = weight + tn + Sn(Form("%.2g",pg.w));
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || mm==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				int sign = d<0?-1:1;
				//number[i] = number[i] + tn + Sn(Form("%.2g",sign*sqrt(fabs(ss))));
				number[i] = number[i] + tn + Sn(Form("%.2g",d));
			}
		}
		header = header + tn + Sn("total", 32) + tn + Sn("sym-total", 16);
		weight = weight + tn + Sn("",32) + tn + Sn("",16);
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
			number[i] = number[i] + tn + Sn(Form("-%.2g+%.2g", syl[i], syh[i]), 32) + tn + Sn(Form("%.2g",sy[i]));
		}

		if(isHeader) std::cout << header << " \\\\ \\hline" << std::endl;
		if(isHeader) std::cout << weight << " \\\\ \\hline" << std::endl;
		for(int i=0; i<n; i++) {
			std::cout << number[i] << " \\\\ " << std::endl;
		}
		std::cout << "\\hline" << std::endl;
	}

	void TableErr(bool isHeader=true, bool isInline = false) {
		if(!IsUpdated) Calc();
		TString tn = " & ";
		TString header = Sn("x-value") + tn + Sn("y-value") + tn + Sn("y-error");
		TString weight = Sn("weight") + tn + Sn("") + tn + Sn("");
		int n = Def.GetN();
		vector<TString> number(n);
		vector<TString> diferr(n);
		vector<TString> numerr(n);
		for(int i=0; i<n; i++) {
			double intx;
			number[i] = (std::modf(Def.GetPoint(i).Px.GetValu()*10,&intx)==0.0?Sn(Form("%.2g",intx/10)):Sn(Def.GetPoint(i).Px.GetValu())) + tn + Sn(Def.GetPoint(i).Py.GetValu()) + tn + Sn(Form("%.2g",Def.GetPoint(i).Py.GetUnce()));
			diferr[i] = Sn("") + tn + Sn("") + tn + Sn("");
			numerr[i] = number[i];
		}
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
			header = header + tn + Sn(pg.s, isInline?24:16);
			weight = weight + tn + Sn(Form("%.2g",pg.w), isInline?24:16);
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || mm==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				int sign = d<0?-1:1;
				//number[i] = number[i] + tn + Sn(Form("%.2g",sign*sqrt(fabs(ss))));
				number[i] = number[i] + tn + Sn(Form("%.2g",d));
				diferr[i] = diferr[i] + tn + Sn(Form("$\\pm$%.2g",sqrt(ee)));
				//numerr[i] = numerr[i] + tn + Sn(Form("%.2g$\\pm$%.2g", sign*sqrt(fabs(ss)), sqrt(ee)), 24);
				numerr[i] = numerr[i] + tn + Sn(Form("%.2g$\\pm$%.2g", d, sqrt(ee)), 24);
			}
		}
		header = header + tn + Sn("total");
		weight = weight + tn + Sn("");
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
			number[i] = number[i] + tn + Sn(Form("$\\pm$%.2g",sy[i]));
			diferr[i] = diferr[i] + tn + Sn("");
			numerr[i] = numerr[i] + tn + Sn(Form("$\\pm$%.2g",sy[i]));
		}

		if(isHeader) std::cout << header << " \\\\ \\hline" << std::endl;
		if(isHeader) std::cout << weight << " \\\\ \\hline" << std::endl;
		for(int i=0; i<n; i++) {
			if(isInline) {
				std::cout << numerr[i] << " \\\\ \\hline" << std::endl;
			} else {
				std::cout << number[i] << " \\\\ " << std::endl;
				std::cout << diferr[i] << " \\\\ \\hline" << std::endl;
			}
		}
		std::cout << "\\hline" << std::endl;
	}
/*
	void TableOpt() {
		if(!IsUpdated) Calc();
		TString tn = " & ";
		int n = Def.GetN();
		int m = Var.size();
		std::vector<TString> mt;
		std::vector<double> mw;
		std::vector<double> ms[n];
		std::vector<double> me[n];
		mt.push_back("x-value"); mw.push_back(-999);
		mt.push_back("y-value"); mw.push_back(-999);
		mt.push_back("y-error"); mw.push_back(-999);
		for(int i=0; i<n; i++) {
			ms[i].push_back(Def.GetPoint(i).Px.GetValu()); me[i].push_back(-999);
			ms[i].push_back(Def.GetPoint(i).Py.GetValu()); me[i].push_back(-999);
			ms[i].push_back(Def.GetPoint(i).Py.GetUnce()); me[i].push_back(-999);
		}
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
			mt.push_back(pg.s);
			mw.push_back(pg.w);
			for(int i=0; i<n; i++) {
				int mm = Mode==0?pg.m:Mode;
				double w = pg.w;
				double d = pg.g.GetPoint(i).Py.GetValu() - Def.GetPoint(i).Py.GetValu();
				double dd = d*d;
				double ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) - pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				if(mm==3) ee = fabs(pow(pg.g.GetPoint(i).Py.GetUnce(),2.0) + pow(Def.GetPoint(i).Py.GetUnce(),2.0));
				double ss = (mm==2 || m==3)?(dd-ee):dd;
				if(ss<0) ss = 0;
				sy[i] += w*ss;
				if(d<0) syl[i] += w*ss;
				if(d>0) syh[i] += w*ss;
				int sign = d<0?-1:1;
				//ms[i].push_back(sign*sqrt(fabs(ss)));
				ms[i].push_back(d);
				me[i].push_back(sqrt(ee));
			}
		}
		mt.push_back("total");
		mw.push_back(-999);
		for(int i=0; i<n; i++) {
			sy[i] = sqrt(sy[i]);
			syl[i] = sqrt(syl[i]);
			syh[i] = sqrt(syh[i]);
			ms[i].push_back(sy[i]);
			me[i].push_back(-999);
		}

		bool isVert = true;
		if(isVert) {
			for(int j=0; j<mt.size(); j++) {
				TString tmpend = (j==(mt.size()-1))?" \\\\ \\hline":tn;
				std::cout << Sn(mt[j]) << tmpend;
			}
			std::cout << std::endl;
			for(int j=0; j<mw.size(); j++) {
				TString tmpend = (j==(mw.size()-1))?" \\\\ \\hline":tn;
				std::cout << (j==0?Sn("weight"):(mw[j]<0?Sn(""):Sn(Form("%.2g",mw[j])))) << tmpend;
			}
			std::cout << std::endl;
			for(int i=0; i<n; i++) {
				for(int j=0; j<ms[i].size(); j++) {
				TString tmpend = (j==(mw.size()-1))?" \\\\ \\hline":tn;
					
				}
			}
		}
	}
*/

	void Write(TFile *f, TString name) const {
		f->WriteObjectAny(this, "MySystGraph", name);
	}

	void WriteGraph(TString name="", bool isstat=true, bool issyst=true, bool isasym=false) {
		//if(!IsUpdated) Calc();
		ShiftX();
		if(isstat) gStat.Write(name+"Stat");
		if(issyst) gSyst.Write(name+"Syst");
		if(isasym) gAsym.Write(name+"Asym");
	}

	void WriteDefVar(TString name="") {
		if(!Def.GetIsUpdated()) Def.Calc();
		Def.Graph.Write(name+"Def");
		for(auto it=Var.begin(); it!=Var.end(); it++) {
			TString StrSyst = Form("Syst%d", it->first);
			if(!it->second.g.GetIsUpdated()) it->second.g.Calc();
			it->second.g.Graph.Write(name+StrSyst);
		}
	}
};


int MySystGraph::VidNow = 0;


// math
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


const MySystGraph operator+(const MySystGraph &sg1, const MyDualMultiv &c) {
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

const MySystGraph operator+(const MyDualMultiv &c, const MySystGraph &sg1) {
	return (sg1 + c);
}

const MySystGraph operator*(const MySystGraph &sg1, const MyDualMultiv &c) {
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

const MySystGraph operator*(const MyDualMultiv &c, const MySystGraph &sg1) {
	return (sg1 * c);
}

const MySystGraph operator-(const MySystGraph &sg1, const MyDualMultiv &c) {
	return (sg1 + (-c));
}

const MySystGraph operator-(const MyDualMultiv &c, const MySystGraph &sg1) {
	return (-1.0*sg1 + c);
}

const MySystGraph operator/(const MySystGraph &sg1, const MyDualMultiv &c) {
	return (sg1 * (1.0/c));
}

const MySystGraph operator/(const MyDualMultiv &c, const MySystGraph &sg1) {
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


const MySystGraph operator+(const MySystGraph &sg1, const MyDualGraph &c) {
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

const MySystGraph operator+(const MyDualGraph &c, const MySystGraph &sg1) {
	return (sg1 + c);
}

const MySystGraph operator*(const MySystGraph &sg1, const MyDualGraph &c) {
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

const MySystGraph operator*(const MyDualGraph &c, const MySystGraph &sg1) {
	return (sg1 * c);
}

const MySystGraph operator-(const MySystGraph &sg1, const MyDualGraph &c) {
	return (sg1 + (-c));
}

const MySystGraph operator-(const MyDualGraph &c, const MySystGraph &sg1) {
	return (-1.0*sg1 + c);
}

const MySystGraph operator/(const MySystGraph &sg1, const MyDualGraph &c) {
	return (sg1 * (1.0/c));
}

const MySystGraph operator/(const MyDualGraph &c, const MySystGraph &sg1) {
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

const MySystGraph atan2(const MySystGraph &sg1, const MySystGraph &sg2) {
	if(sg1.GetMode() != sg2.GetMode()) std::cout << "MySystGraph: Mode not match!" << std::endl;
	MySystGraph sg3(sg1);
	sg3.SetDef(atan2(sg1.GetDef(), sg2.GetDef()));
	std::vector<int> list1 = sg1.GetList();
	std::vector<int> list2 = sg2.GetList();
	int i = 0;
	int j = 0;
	while(i<(int)list1.size() && j<(int)list2.size()) {
		if(list1[i] < list2[j]) {
			int idx = list1[i];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = atan2(pg.g, sg2.GetDef());
			sg3.SetVar(idx, pg);
			i++;
		} else if(list1[i] > list2[j]) {
			int idx = list2[j];
			MyPackGraph pg = sg2.GetVar(idx);
			pg.g = atan2(sg1.GetDef(), pg.g);
			sg3.SetVar(idx, pg);
			j++;
		} else {
			int idx = list1[i]; // = list2[j];
			MyPackGraph pg = sg1.GetVar(idx);
			pg.g = atan2(pg.g, sg2.GetVar(idx).g);
			sg3.SetVar(idx, pg);
			i++;
			j++;
		}
	}
	while(i<(int)list1.size()) {
		int idx = list1[i];
		MyPackGraph pg = sg1.GetVar(idx);
		pg.g = atan2(pg.g, sg2.GetDef());
		sg3.SetVar(idx, pg);
		i++;
	}
	while(j<(int)list2.size()) {
		int idx = list2[j];
		MyPackGraph pg = sg2.GetVar(idx);
		pg.g = atan2(sg1.GetDef(), pg.g);
		sg3.SetVar(idx, pg);
		j++;
	}
	return sg3;
}


#endif
