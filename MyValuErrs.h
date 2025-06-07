#ifndef MyValuErrs_H
#define MyValuErrs_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>


class MyValuErrs
{
private:
	int Mode; // number of errs: 0. valu; 1. valu+/-stat; 2. valu+/-stat+/-syst; 3. valu+/-stat (-sysl+sysh)

public:
	double Valu;
	double Stat;
	double Syst;
	double Sysl; // asymmetric systematics -
	double Sysh; // asymmetric systematics +

	int Precision;
	int Exponent;
	double ValuCoeff;
	double StatCoeff;
	double SystCoeff;
	double SyslCoeff;
	double SyshCoeff;

	static double CalcCoefficient(double x) {
		if(x==0) return 0;
		int n = 0;
		while(fabs(x)<1) {
			x *= 10;
			n ++;
		}
		while(fabs(x)>=10) {
			x *= 0.1;
			n --;
		}
		return x;
	}
	
	static int CalcExponent(double x) {
		if(x==0) return 0;
		int n = 0;
		while(fabs(x)<1) {
			x *= 10;
			n --;
		}
		while(fabs(x)>=10) {
			x *= 0.1;
			n ++;
		}
		return n;
	}

	static std::vector<double> ValuErr1(double val, double err) {
		// coefficient * 10^exponent: c * 10^e
		double val_c = CalcCoefficient(val);
		int    val_e = CalcExponent(val);
		double err_c = CalcCoefficient(err);
		int    err_e = CalcExponent(err);

		if(val==0 && err==0) return std::vector<double>{0, 0, 0, 0};
		if(val==0) val_e = err_e;
		if(err==0) err_e = val_e;

		double value;
		double error;
		int exponent;
		int precision;
		if(fabs(val)>=fabs(err) || val_e>=err_e) {
			value = val_c;
			exponent = val_e;
			int err_de = err_e - val_e;
			error = err_c * std::pow(10.0, err_de);
			if(fabs(err_c)<3.5) {
				precision = 1 - err_de;
			} else {
				precision = 0 - err_de;
			}
		} else if(val_e+1==err_e && fabs(err_c)<3.5) {
			value = val_c / 10;
			exponent = err_e;
			error = err_c;
			precision = 1;
		} else if(fabs(err_c)<3.5) {
			value = val_c * std::pow(10.0, val_e-err_e);
			exponent = err_e;
			error = err_c;
			precision = 1;
		} else {
			value = val_c * std::pow(10.0, val_e-err_e);
			exponent = err_e;
			error = err_c;
			precision = 0;
		}

		std::vector<double> vtmp;
		vtmp.push_back(1.0*precision);
		vtmp.push_back(1.0*exponent);
		vtmp.push_back(value);
		vtmp.push_back(error);
		vtmp.push_back(0);
		vtmp.push_back(0);
		vtmp.push_back(0);

		return vtmp;
	}

	static std::vector<double> ValuErr0(double val) {
		std::vector<double> vtmp = ValuErr1(val, val);
		vtmp[3] = 0;
		return vtmp;
	}

	static std::vector<double> ValuErr2(double val, double stat, double syst) {
		// coefficient * 10^exponent: c * 10^e
		double val_c = CalcCoefficient(val);
		int    val_e = CalcExponent(val);
		double stat_c = CalcCoefficient(stat);
		int    stat_e = CalcExponent(stat);
		double syst_c = CalcCoefficient(syst);
		int    syst_e = CalcExponent(syst);

		double err = fabs(stat)>fabs(syst)?fabs(stat):fabs(syst);
		double err_c = CalcCoefficient(err);
		int    err_e = CalcExponent(err);

		std::vector<double> verr1 = ValuErr1(val, err);
		int precision = (int)std::round(verr1[0]);
		int exponent = (int)std::round(verr1[1]);
		double value = val_c * std::pow(10.0, val_e-exponent);
		double statistics = stat_c * std::pow(10.0, stat_e-exponent);
		double systematics = syst_c * std::pow(10.0, syst_e-exponent);

		std::vector<double> vtmp;
		vtmp.push_back(1.0*precision);
		vtmp.push_back(1.0*exponent);
		vtmp.push_back(value);
		vtmp.push_back(statistics);
		vtmp.push_back(systematics);
		vtmp.push_back(systematics);
		vtmp.push_back(0);

		return vtmp;
	}

	static std::vector<double> ValuErr3(double val, double stat, double sysl, double sysh) {
		// coefficient * 10^exponent: c * 10^e
		double val_c = CalcCoefficient(val);
		int    val_e = CalcExponent(val);
		double stat_c = CalcCoefficient(stat);
		int    stat_e = CalcExponent(stat);
		double sysl_c = CalcCoefficient(sysl);
		int    sysl_e = CalcExponent(sysl);
		double sysh_c = CalcCoefficient(sysh);
		int    sysh_e = CalcExponent(sysh);

		double err = fabs(stat);
		if(fabs(sysl)>err) err = fabs(sysl);
		if(fabs(sysh)>err) err = fabs(sysh);
		double err_c = CalcCoefficient(err);
		int    err_e = CalcExponent(err);

		std::vector<double> verr1 = ValuErr1(val, err);
		int precision = (int)std::round(verr1[0]);
		int exponent = (int)std::round(verr1[1]);
		double value = val_c * std::pow(10.0, val_e-exponent);
		double statistics = stat_c * std::pow(10.0, stat_e-exponent);
		double systematics_low  = sysl_c * std::pow(10.0, sysl_e-exponent);
		double systematics_high = sysh_c * std::pow(10.0, sysh_e-exponent);
		double systematics = std::sqrt(systematics_low*systematics_low + systematics_high*systematics_high);

		std::vector<double> vtmp;
		vtmp.push_back(1.0*precision);
		vtmp.push_back(1.0*exponent);
		vtmp.push_back(value);
		vtmp.push_back(statistics);
		vtmp.push_back(systematics);
		vtmp.push_back(systematics_low);
		vtmp.push_back(systematics_high);

		return vtmp;
	}

	void Calc(int mode=1) {
		if(mode<0 || mode>3) {
			std::cout << "MyValuErrs: invalid mode " << mode << ". Set to default mode = 1" << std::endl;
			mode = 1;
		}
		Mode = mode;
		std::vector<double> vtmp;
		if(mode==0) vtmp = ValuErr0(Valu);
		if(mode==1) vtmp = ValuErr1(Valu, Stat);
		if(mode==2) vtmp = ValuErr2(Valu, Stat, Syst);
		if(mode==3) vtmp = ValuErr3(Valu, Stat, Sysl, Sysh);
		Precision = (int)std::round(vtmp[0]);
		Exponent  = (int)std::round(vtmp[1]);
		ValuCoeff = vtmp[2];
		StatCoeff = vtmp[3];
		SystCoeff = vtmp[4];
		SyslCoeff = vtmp[5];
		SyshCoeff = vtmp[6];
	}

	MyValuErrs() {
		Valu = 0;
		Stat = 0;
		Syst = 0;
		Sysl = 0;
		Sysh = 0;

		Precision = 0;
		Exponent  = 0;
		ValuCoeff = 0;
		StatCoeff = 0;
		SystCoeff = 0;
		SyslCoeff = 0;
		SyshCoeff = 0;

		Mode = 1;
	}

	MyValuErrs(double val) {
		Valu = val;
		Stat = 0;
		Syst = 0;
		Sysl = 0;
		Sysh = 0;
		Calc(0);
	}

	MyValuErrs(double val, double err) {
		Valu = val;
		Stat = err;
		Syst = 0;
		Sysl = 0;
		Sysh = 0;
		Calc(1);
	}

	MyValuErrs(double val, double err, double sys) {
		Valu = val;
		Stat = err;
		Syst = sys;
		Sysl = sys;
		Sysh = 0;
		Calc(2);
	}

	MyValuErrs(double val, double err, double syl, double syh) {
		Valu = val;
		Stat = err;
		Syst = std::sqrt(syl*syl + syh*syh);
		Sysl = syl;
		Sysh = syh;
		Calc(3);
	}

	void AddPrecision(int p) {
		Precision += p;
	}

	void SetPrecision(int p) {
		Precision = p;
	}

	std::string PrintErr0(std::string opt="R") {
		if(Mode!=0) Calc(0);
		std::vector<double> verr0 = ValuErr0(Valu);
		int prec = (int)std::round(verr0[0]);
		int expo = (int)std::round(verr0[1]);
		double valucoeff = verr0[2];
		int de = 0;
		int np = prec;
		bool isnormal = (expo>=-1 && expo<=0);
		if(isnormal) {
			de = expo;
			np = prec + abs(expo);
		}
		std::stringstream StrValu;
		StrValu << std::fixed << std::setprecision(np) << valucoeff * std::pow(10.0,de);
		std::stringstream StrTmp;
		if(isnormal) {
			if(opt=="L") StrTmp << StrValu.str();
			if(opt=="R") StrTmp << StrValu.str();
		} else {
			if(opt=="L") StrTmp << StrValu.str() << "\\times10^{" <<  expo << "}";
			if(opt=="R") StrTmp << StrValu.str() <<  "#times10^{" <<  expo << "}";
		}
		return StrTmp.str();
	}

	std::string PrintErr1(std::string opt="R") {
		if(Mode!=1) Calc(1);
		int de = 0;
		int np = Precision;
		bool isnormal = (Exponent>=-1 && Exponent<=0);
		if(isnormal) {
			de = Exponent;
			np = Precision + abs(Exponent);
		}
		std::stringstream StrValu;
		std::stringstream StrStat;
		StrValu << std::fixed << std::setprecision(np) << ValuCoeff * std::pow(10.0,de);
		StrStat << std::fixed << std::setprecision(np) << StatCoeff * std::pow(10.0,de);
		std::stringstream StrTmp;
		if(isnormal) {
			if(opt=="L") StrTmp << StrValu.str() << "\\pm" << StrStat.str();
			if(opt=="R") StrTmp << StrValu.str() <<  "#pm" << StrStat.str();
		} else {
			if(opt=="L") StrTmp << "(" << StrValu.str() << "\\pm" << StrStat.str() << ")\\times10^{" << Exponent << "}";
			if(opt=="R") StrTmp << "(" << StrValu.str() <<  "#pm" << StrStat.str() <<  ")#times10^{" << Exponent << "}";
		}
		return StrTmp.str();
	}

	std::string PrintErr2(std::string opt="R") {
		if(Mode!=2) Calc(2);
		int de = 0;
		int np = Precision;
		bool isnormal = (Exponent>=-1 && Exponent<=0);
		if(isnormal) {
			de = Exponent;
			np = Precision + abs(Exponent);
		}
		std::stringstream StrValu;
		std::stringstream StrStat;
		std::stringstream StrSyst;
		StrValu << std::fixed << std::setprecision(np) << ValuCoeff * std::pow(10.0,de);
		StrStat << std::fixed << std::setprecision(np) << StatCoeff * std::pow(10.0,de);
		StrSyst << std::fixed << std::setprecision(np) << SystCoeff * std::pow(10.0,de);
		std::stringstream StrTmp;
		if(isnormal) {
			if(opt=="L") StrTmp << StrValu.str() << "\\pm" << StrStat.str() << "\\pm" << StrSyst.str();
			if(opt=="R") StrTmp << StrValu.str() <<  "#pm" << StrStat.str() <<  "#pm" << StrSyst.str();
		} else {
			if(opt=="L") StrTmp << "(" << StrValu.str() << "\\pm" << StrStat.str() << "\\pm" << StrSyst.str() << ")\\times10^{" << Exponent << "}";
			if(opt=="R") StrTmp << "(" << StrValu.str() <<  "#pm" << StrStat.str() <<  "#pm" << StrSyst.str() <<  ")#times10^{" << Exponent << "}";
		}
		return StrTmp.str();
	}

	std::string PrintErr3(std::string opt="R") {
		if(Mode!=3) Calc(3);
		int de = 0;
		int np = Precision;
		bool isnormal = (Exponent>=-1 && Exponent<=0);
		if(isnormal) {
			de = Exponent;
			np = Precision + abs(Exponent);
		}
		std::stringstream StrValu;
		std::stringstream StrStat;
		std::stringstream StrSysl;
		std::stringstream StrSysh;
		StrValu << std::fixed << std::setprecision(np) << ValuCoeff * std::pow(10.0,de);
		StrStat << std::fixed << std::setprecision(np) << StatCoeff * std::pow(10.0,de);
		StrSysl << std::fixed << std::setprecision(np) << SyslCoeff * std::pow(10.0,de);
		StrSysh << std::fixed << std::setprecision(np) << SyshCoeff * std::pow(10.0,de);
		std::stringstream StrTmp;
		if(isnormal) {
			if(opt=="L") StrTmp << StrValu.str() << "\\pm" << StrStat.str() << "_{-" << StrSysl.str() << "}^{+" << StrSysh.str() << "}";
			if(opt=="R") StrTmp << StrValu.str() <<  "#pm" << StrStat.str() << "_{-" << StrSysl.str() << "}^{+" << StrSysh.str() << "}";
		} else {
			if(opt=="L") StrTmp << "(" << StrValu.str() << "\\pm" << StrStat.str() << "_{-" << StrSysl.str() << "}^{+" << StrSysh.str() << "})\\times10^{" << Exponent << "}";
			if(opt=="R") StrTmp << "(" << StrValu.str() <<  "#pm" << StrStat.str() << "_{-" << StrSysl.str() << "}^{+" << StrSysh.str() <<  "})#times10^{" << Exponent << "}";
		}
		return StrTmp.str();
	}

	std::string Print(std::string opt="R", int mode=-999) {
		if(mode>=0 && mode<=3) Mode = mode;
		Calc(Mode);
		if(Mode==0) return PrintErr0(opt);
		if(Mode==1) return PrintErr1(opt);
		if(Mode==2) return PrintErr2(opt);
		if(Mode==3) return PrintErr3(opt);
		return std::string();
	}

	std::string PrintRaw() const {
		std::stringstream StrRaw;
		StrRaw << Valu << "   " << Stat << "   " << Syst << "   " << Sysl << "   " << Sysh;
		return StrRaw.str();
	}
};

#endif
