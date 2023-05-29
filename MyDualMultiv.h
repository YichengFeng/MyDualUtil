/**************************************************************************
 * Author: Yicheng Feng
 * Email: fengyich@outlook.com
 * Note: dual number to calculate precise values of first-order derivative
 *       of multi-variable functions with common math formula 
 **************************************************************************/

#ifndef MyDualMultiv_H
#define MyDualMultiv_H

#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include "MyDualNumber.h"


class MyDualMultiv
{
private:
	double Valu;
	std::map<int,double> Dual;

public:
	void Reset() {
		Valu = 0;
		Dual.clear();
	}

	MyDualMultiv() {
		Valu = 0;
	}
	MyDualMultiv(double valu) {
		Valu = valu;
	}
	MyDualMultiv(int indx, const MyDualNumber &dn) {
		Valu = dn.GetValu();
		Dual.insert({indx, dn.GetDual()});
	}
	MyDualMultiv(const MyDualMultiv &dm) {
		Valu = dm.Valu;
		Dual = dm.Dual;
	}
	~MyDualMultiv() {
	}

	double GetValu() const { return Valu; }
	void SetValu(double valu) { Valu = valu; }

	std::vector<int> GetList() const {
		std::vector<int> list;
		for(auto it=Dual.begin(); it!=Dual.end(); it++) {
			list.push_back(it->first);
		}
		return list;
	}

	double GetDual(int indx) const {
		double dual = 0;
		if(Dual.count(indx)!=0) dual=Dual.at(indx);
		return dual;
	}
	void SetDual(int indx, double dual) {
		if(Dual.count(indx)==0) {
			Dual.insert({indx, dual});
		} else {
			Dual[indx] = dual;
		}
	}
	void SetDual(int indx, const MyDualNumber &dn) {
		SetDual(indx, dn.GetDual());
	}

	const MyDualNumber GetPart(int indx) const {
		return MyDualNumber(Valu, GetDual(indx));
	}

	double GetUnce() const {
		double sums = 0;
		for(auto it=Dual.begin(); it!=Dual.end(); it++) {
			sums += (it->second * it->second);
		}
		return sqrt(sums);
	}

	void Print() const { 
		std::cout << "(" << Valu << ", " << std::flush;
		for(auto it=Dual.begin(); it!=Dual.end(); it++) {
			std::cout << "[" << it->first << ":" << it->second << "]" << std::flush;
		}
		std::cout << ")" << std::endl;
	}

	std::string StrLatex() const {
		MyDualNumber dn(Valu, GetUnce());
		return dn.StrLatex();
	}
};


// math
const MyDualMultiv operator+(const MyDualMultiv &dm1, double c) {
	MyDualMultiv dm2 = dm1;
	dm2.SetValu(dm2.GetValu() + c);
	return dm2;
}

const MyDualMultiv operator+(double c, const MyDualMultiv &dm1) {
	return (dm1 + c);
}

const MyDualMultiv operator*(const MyDualMultiv &dm1, double c) {
	MyDualMultiv dm2 = dm1;
	dm2.SetValu(dm1.GetValu() * c);
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		MyDualNumber dn = dm1.GetPart(list[i]);
		dm2.SetDual(list[i], c*dn);
	}
	return dm2;
}

const MyDualMultiv operator*(double c, const MyDualMultiv &dm1) {
	return (dm1 * c);
}

const MyDualMultiv operator-(const MyDualMultiv &dm1, double c) {
	return (dm1 + (-1.0*c));
}

const MyDualMultiv operator-(double c, const MyDualMultiv &dm1) {
	return (c + (-1.0*dm1));
}

const MyDualMultiv operator/(const MyDualMultiv &dm1, double c) {
	return (dm1 * (1.0/c));
}

const MyDualMultiv operator/(double c, const MyDualMultiv &dm1) {
	MyDualMultiv dm2 = dm1;
	dm2.SetValu(c / dm1.GetValu());
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		MyDualNumber dn = dm1.GetPart(list[i]);
		dm2.SetDual(list[i], c/dn);
	}
	return dm2;
}

const MyDualMultiv operator+(const MyDualMultiv &dm1) {
	return dm1;
}

const MyDualMultiv operator-(const MyDualMultiv &dm1) {
	return (0.0 - dm1);
}

const MyDualMultiv operator+(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	MyDualMultiv dm3(dm1.GetValu() + dm2.GetValu());
	std::vector<int> list1 = dm1.GetList();
	std::vector<int> list2 = dm2.GetList();
	int i = 0;
	int j = 0;
	while(i<(int)list1.size() && j<(int)list2.size()) {
		if(list1[i] < list2[j]) {
			int indx = list1[i];
			MyDualNumber dn1 = dm1.GetPart(indx);
			dm3.SetDual(indx, dn1+dm2.GetValu());
			i++;
		} else if(list1[i] > list2[j]) {
			int indx = list2[j];
			MyDualNumber dn2 = dm2.GetPart(indx);
			dm3.SetDual(indx, dn2+dm1.GetValu());
			j++;
		} else {
			int indx = list1[i]; // = list2[j];
			MyDualNumber dn1 = dm1.GetPart(indx);
			MyDualNumber dn2 = dm2.GetPart(indx);
			dm3.SetDual(indx, dn1+dn2);
			i++;
			j++;
		}
	}
	while(i<(int)list1.size()) {
		int indx = list1[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm3.SetDual(indx, dn1+dm2.GetValu());
		i++;
	}
	while(j<(int)list2.size()) {
		int indx = list2[j];
		MyDualNumber dn2 = dm2.GetPart(indx);
		dm3.SetDual(indx, dn2+dm1.GetValu());
		j++;
	}
	return dm3;
}

const MyDualMultiv operator-(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	return (dm1 + (-1.0*dm2));
}

const MyDualMultiv operator*(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	MyDualMultiv dm3(dm1.GetValu() * dm2.GetValu());
	std::vector<int> list1 = dm1.GetList();
	std::vector<int> list2 = dm2.GetList();
	int i = 0;
	int j = 0;
	while(i<(int)list1.size() && j<(int)list2.size()) {
		if(list1[i] < list2[j]) {
			int indx = list1[i];
			MyDualNumber dn1 = dm1.GetPart(indx);
			dm3.SetDual(indx, dn1*dm2.GetValu());
			i++;
		} else if(list1[i] > list2[j]) {
			int indx = list2[j];
			MyDualNumber dn2 = dm2.GetPart(indx);
			dm3.SetDual(indx, dn2*dm1.GetValu());
			j++;
		} else {
			int indx = list1[i]; // = list2[j];
			MyDualNumber dn1 = dm1.GetPart(indx);
			MyDualNumber dn2 = dm2.GetPart(indx);
			dm3.SetDual(indx, dn1*dn2);
			i++;
			j++;
		}
	}
	while(i<(int)list1.size()) {
		int indx = list1[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm3.SetDual(indx, dn1*dm2.GetValu());
		i++;
	}
	while(j<(int)list2.size()) {
		int indx = list2[j];
		MyDualNumber dn2 = dm2.GetPart(indx);
		dm3.SetDual(indx, dn2*dm1.GetValu());
		j++;
	}
	return dm3;
}

const MyDualMultiv operator/(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	return (dm1 * (1.0/dm2));
}

const MyDualMultiv pow(const MyDualMultiv &dm1, double c) {
	MyDualMultiv dm2(pow(dm1.GetValu(),c));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, pow(dn1,c));
	}
	return dm2;
}

const MyDualMultiv sqrt(const MyDualMultiv &dm1) {
	return pow(dm1, 0.5);
}

const MyDualMultiv abs(const MyDualMultiv &dm1) {
	if(dm1.GetValu()<0) return -dm1;
	return dm1;
}

const MyDualMultiv exp(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(exp(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, exp(dn1));
	}
	return dm2;
}

const MyDualMultiv log(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(log(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, log(dn1));
	}
	return dm2;
}

const MyDualMultiv sin(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(sin(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, sin(dn1));
	}
	return dm2;
}

const MyDualMultiv cos(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(cos(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, cos(dn1));
	}
	return dm2;
}

const MyDualMultiv tan(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(tan(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, tan(dn1));
	}
	return dm2;
}

const MyDualMultiv atan(const MyDualMultiv &dm1) {
	MyDualMultiv dm2(atan(dm1.GetValu()));
	std::vector<int> list = dm1.GetList();
	for(int i=0; i<(int)list.size(); i++) {
		int indx = list[i];
		MyDualNumber dn1 = dm1.GetPart(indx);
		dm2.SetDual(indx, atan(dn1));
	}
	return dm2;
}


// merge
// must use different data for merge, so they should be independent
// but this code does not check that
const MyDualMultiv AvePlus(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	double e1 = dm1.GetUnce();
	double e2 = dm2.GetUnce();
	double w1 = 1.0/e1/e1;
	double w2 = 1.0/e2/e2;
	return ((w1*dm1 + w2*dm2) / (w1 + w2));
}

// dm1 must already include dm2
const MyDualMultiv AveMinus(const MyDualMultiv &dm1, const MyDualMultiv &dm2) {
	double e1 = dm1.GetUnce();
	double e2 = dm2.GetUnce();
	double w1 = 1.0/e1/e1;
	double w2 = 1.0/e2/e2;
	return ((w1*dm1 - w2*dm2) / (w1 - w2));
}


#endif
