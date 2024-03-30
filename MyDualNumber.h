/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 * Note: dual number to calculate precise values of first-order derivative
 *       of single-variable functions with common math formula 
 **************************************************************************/

#ifndef MyDualNumber_H
#define MyDualNumber_H

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>


class MyDualNumber
{
private:
	double Valu;
	double Dual;

public:
	MyDualNumber() {
		Valu = 0;
		Dual = 1;
	}
	MyDualNumber(double valu, double dual=1.0) {
		Valu = valu;
		Dual = dual;
	}
	MyDualNumber(const MyDualNumber &dn) {
		Valu = dn.Valu;
		Dual = dn.Dual;
	}
	~MyDualNumber() {
	}

	double GetValu() const { return Valu; }
	double GetDual() const { return Dual; }
	void SetValu(double valu) { Valu = valu; }
	void SetDual(double dual) { Dual = dual; }

	void Print() const { std::cout << "(" << Valu << ", " << Dual << ")" << std::endl; }

	std::string StrLatex(double &outval, double &outerr, int &outn10, int &outnp, string symboltype = "R") const {
		// symboltype = "R" -- ROOT; "L" -- "Latex"
		int n10;
		double val, err;

		if(Valu==0 && Dual==0) {
			n10 = 0;
			val = 0;
			err = 0;
		} else if(Valu==0) {
			n10 = (int)log10(fabs(Dual));
			if(fabs(Dual)<1.0) n10--;
			val = 0;
			err = fabs(Dual) * pow(10.0, -1.0*n10);
		} else if(Dual==0) {
			n10 = (int)log10(fabs(Valu));
			if(fabs(Valu)<1.0) n10--;
			val = Valu * pow(10.0, -1.0*n10);
			err = 0;
		} else {
			n10 = (int)log10(fabs(Dual));
			if(fabs(Dual)<1.0) n10--;
			val = Valu * pow(10.0, -1.0*n10);
			err = fabs(Dual) * pow(10.0, -1.0*n10);
		}

		int np = 1;
		if(err>=3.5) {
			n10++;
			val *= 0.1;
			err *= 0.1;
		}

		int v10 = (int)log10(fabs(val));
		if(fabs(val)<1.0) v10--;
		if(v10>=1) {
			n10 += v10;
			val *= pow(10.0, -1.0*v10);
			err *= pow(10.0, -1.0*v10);
			np  += v10;
		}

		outn10 = n10;
		outval = val;
		outerr = err;
		outnp  = np;

		std::stringstream strval;
		std::stringstream strerr;
		strval << std::fixed << std::setprecision(np) << val;
		strerr << std::fixed << std::setprecision(np) << err;
		std::stringstream strlatex;
		if(n10==0) {
			if(symboltype=="L") {
				strlatex << strval.str() << "\\pm" << strerr.str();
			} else {
				strlatex << strval.str() << "#pm" << strerr.str();
			}
		} else {
			if(symboltype=="L") {
				strlatex << "(" << strval.str() << "\\pm" << strerr.str() << ")\\times10^{" << n10 << "}";
			} else {
				strlatex << "(" << strval.str() << "#pm" << strerr.str() << ")#times10^{" << n10 << "}";
			}
		}

		return strlatex.str();
	}

	std::string StrLatex(string symboltype="R") const {
		int outn10, outnp;
		double outval, outerr;
		return StrLatex(outval, outerr, outn10, outnp, symboltype);
	}
};


// math
const MyDualNumber operator+(const MyDualNumber &dn1, double c) {
	return MyDualNumber(dn1.GetValu() + c, dn1.GetDual());
}

const MyDualNumber operator+(double c, const MyDualNumber &dn1) {
	return (dn1 + c);
}

const MyDualNumber operator*(const MyDualNumber &dn1, double c) {
	return MyDualNumber(dn1.GetValu() * c, dn1.GetDual() * c);
}

const MyDualNumber operator*(double c, const MyDualNumber &dn1) {
	return (dn1 * c);
}

const MyDualNumber operator-(const MyDualNumber &dn1, double c) {
	return (dn1 + (-1.0*c));
}

const MyDualNumber operator-(double c, const MyDualNumber &dn1) {
	return (((-1.0)*dn1) + c);
}

const MyDualNumber operator/(const MyDualNumber &dn1, double c) {
	return (dn1 * (1.0/c));
}

const MyDualNumber operator/(double c, const MyDualNumber &dn1) {
	return MyDualNumber(c/dn1.GetValu(), -dn1.GetDual()/dn1.GetValu()/dn1.GetValu());
}

const MyDualNumber operator+(const MyDualNumber &dn1) {
	return dn1;
}

const MyDualNumber operator-(const MyDualNumber &dn1) {
	return (0.0 - dn1);
}

const MyDualNumber operator+(const MyDualNumber &dn1, const MyDualNumber &dn2) {
	return MyDualNumber(dn1.GetValu()+dn2.GetValu(), dn1.GetDual()+dn2.GetDual());
}

const MyDualNumber operator-(const MyDualNumber &dn1, const MyDualNumber &dn2) {
	return (dn1 + ((-1.0)*dn2));
}

const MyDualNumber operator*(const MyDualNumber &dn1, const MyDualNumber &dn2) {
	return MyDualNumber(dn1.GetValu()*dn2.GetValu(), (dn1.GetValu()*dn2.GetDual()+dn2.GetValu()*dn1.GetDual()));
}

const MyDualNumber operator/(const MyDualNumber &dn1, const MyDualNumber &dn2) {
	return (dn1 * (1.0/dn2));
}

const MyDualNumber pow(const MyDualNumber &dn1, double c) {
	return MyDualNumber(pow(dn1.GetValu(),c), c*pow(dn1.GetValu(),c-1.0)*dn1.GetDual());
}

const MyDualNumber sqrt(const MyDualNumber &dn1) {
	return pow(dn1, 0.5);
}

const MyDualNumber abs(const MyDualNumber &dn1) {
	if(dn1.GetValu()<0) return -dn1;
	return dn1;
}

const MyDualNumber exp(const MyDualNumber &dn1) {
	return MyDualNumber(exp(dn1.GetValu()), exp(dn1.GetValu())*dn1.GetDual());
}

const MyDualNumber log(const MyDualNumber &dn1) {
	return MyDualNumber(log(dn1.GetValu()), dn1.GetDual()/dn1.GetValu());
}

const MyDualNumber sin(const MyDualNumber &dn1) {
	return MyDualNumber(sin(dn1.GetValu()), cos(dn1.GetValu())*dn1.GetDual());
}

const MyDualNumber cos(const MyDualNumber &dn1) {
	return MyDualNumber(cos(dn1.GetValu()), -sin(dn1.GetValu())*dn1.GetDual());
}

const MyDualNumber tan(const MyDualNumber &dn1) {
	return MyDualNumber(tan(dn1.GetValu()), dn1.GetDual()/pow(cos(dn1.GetValu()),2.0));
}

const MyDualNumber atan(const MyDualNumber &dn1) {
	return MyDualNumber(atan(dn1.GetValu()), dn1.GetDual()/(1.0+dn1.GetValu()*dn1.GetValu()));
}

const MyDualNumber atan2(const MyDualNumber &dn1, const MyDualNumber &dn2) {
	double valu = atan2(dn1.GetValu(), dn2.GetValu());
	double dual = atan(dn1 / dn2).GetDual();
	return MyDualNumber(valu, dual);
}


#endif
