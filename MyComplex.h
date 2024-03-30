#ifndef MyComplex_H
#define MyComplex_H

#include "MySystGraph.h"


template <class T> class MyComplex
{
public:
	T Re;
	T Im;

	MyComplex() {
	}

	MyComplex(T re, T im) {
		Re = re;
		Im = im;
	}

	~MyComplex() {
	}

	const T GetRe() const {
		return Re;
	}

	const T GetIm() const {
		return Im;
	}

	const MyComplex<T> Conjugate() const {
		return MyComplex<T>(Re, -1.0*Im);
	}

	const T GetAmp() const {
		return sqrt(Re*Re + Im*Im);
	}

	const T GetAng() const {
		return atan2(Im, Re);
	}

	void Print() const {
		std::cout << "Re: ";
		Re.Print();
		std::cout << "Im: ";
		Im.Print();
	}
};


// math
template <class T> const MyComplex<T> operator+(const MyComplex<T> &t1, double c) {
	return MyComplex<T>(t1.GetRe()+c, t1.GetIm());
}

template <class T> const MyComplex<T> operator+(double c, const MyComplex<T> &t1) {
	return (t1 + c);
}

template <class T> const MyComplex<T> operator*(const MyComplex<T> &t1, double c) {
	return MyComplex<T>(t1.GetRe()*c, t1.GetIm()*c);
}

template <class T> const MyComplex<T> operator*(double c, const MyComplex<T> &t1) {
	return (t1 * c);
}

template <class T> const MyComplex<T> operator-(const MyComplex<T> &t1, double c) {
	return (t1 + (-1.0*c));
}

template <class T> const MyComplex<T> operator-(double c, const MyComplex<T> &t1) {
	return ((-1.0*t1) + c);
}

template <class T> const MyComplex<T> operator/(const MyComplex<T> &t1, double c) {
	return (t1 * (1.0/c));
}

template <class T> const MyComplex<T> operator/(double c, const MyComplex<T> &t1) {
	T tmp = t1.GetRe()*t1.GetRe() + t1.GetIm()*t1.GetIm();
	return c*MyComplex<T>(t1.GetRe()/tmp, -1.0*t1.GetIm()/tmp);
}


template <class T> const MyComplex<T> operator+(const MyComplex<T> &t1, T c) {
	return MyComplex<T>(t1.GetRe()+c, t1.GetIm());
}

template <class T> const MyComplex<T> operator+(T c, const MyComplex<T> &t1) {
	return (t1 + c);
}

template <class T> const MyComplex<T> operator*(const MyComplex<T> &t1, T c) {
	return MyComplex<T>(t1.GetRe()*c, t1.GetIm()*c);
}

template <class T> const MyComplex<T> operator*(T c, const MyComplex<T> &t1) {
	return (t1 * c);
}

template <class T> const MyComplex<T> operator-(const MyComplex<T> &t1, T c) {
	return (t1 + (-1.0*c));
}

template <class T> const MyComplex<T> operator-(T c, const MyComplex<T> &t1) {
	return ((-1.0*t1) + c);
}

template <class T> const MyComplex<T> operator/(const MyComplex<T> &t1, T c) {
	return (t1 * (1.0/c));
}

template <class T> const MyComplex<T> operator/(T c, const MyComplex<T> &t1) {
	T tmp = t1.GetRe()*t1.GetRe() + t1.GetIm()*t1.GetIm();
	return c*MyComplex<T>(t1.GetRe()/tmp, -1.0*t1.GetIm()/tmp);
}


template <class T> const MyComplex<T> operator+(const MyComplex<T> &t1, const MyComplex<T> &t2) {
	return MyComplex<T>(t1.GetRe()+t2.GetRe(), t1.GetIm()+t2.GetIm());
}

template <class T> const MyComplex<T> operator-(const MyComplex<T> &t1, const MyComplex<T> &t2) {
	return (t1 + (-1.0*t2));
}

template <class T> const MyComplex<T> operator*(const MyComplex<T> &t1, const MyComplex<T> &t2) {
	return MyComplex<T>(t1.GetRe()*t2.GetRe()-t1.GetIm()*t2.GetIm(), t1.GetRe()*t2.GetIm()+t1.GetIm()*t2.GetRe());
}

template <class T> const MyComplex<T> operator/(const MyComplex<T> &t1, const MyComplex<T> &t2) {
	return (t1 * (1.0/t2));
}


template <class T> const MyComplex<T> operator+(const MyComplex<T> &t1) {
	return t1;
}

template <class T> const MyComplex<T> operator-(const MyComplex<T> &t1) {
	return (-1.0*t1);
}


const MyComplex<MySystGraph> operator+(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualGraph> &c) {
	return MyComplex<MySystGraph>(t1.GetRe()+c.GetRe(), t1.GetIm()+c.GetIm());
}

const MyComplex<MySystGraph> operator+(const MyComplex<MyDualGraph> &c, const MyComplex<MySystGraph> &t1) {
	return (t1 + c);
}

const MyComplex<MySystGraph> operator-(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualGraph> &c) {
	return (t1 + (-1.0*c));
}

const MyComplex<MySystGraph> operator-(const MyComplex<MyDualGraph> &c, const MyComplex<MySystGraph> &t1) {
	return ((-1.0*t1) + c);
}

const MyComplex<MySystGraph> operator*(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualGraph> &c) {
	return MyComplex<MySystGraph>(t1.GetRe()*c.GetRe()-t1.GetIm()*c.GetIm(), t1.GetRe()*c.GetIm()+t1.GetIm()*c.GetRe());
}

const MyComplex<MySystGraph> operator*(const MyComplex<MyDualGraph> &c, const MyComplex<MySystGraph> &t1) {
	return (t1 * c);
}

const MyComplex<MySystGraph> operator/(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualGraph> &c) {
	return (t1 * (1.0/c));
}

const MyComplex<MySystGraph> operator/(const MyComplex<MyDualGraph> &c, const MyComplex<MySystGraph> &t1) {
	return ((1.0/t1) * c);
}


const MyComplex<MySystGraph> operator+(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return MyComplex<MySystGraph>(t1.GetRe()+c.GetRe(), t1.GetIm()+c.GetIm());
}

const MyComplex<MySystGraph> operator+(const MyComplex<MyDualMultiv> &c, const MyComplex<MySystGraph> &t1) {
	return (t1 + c);
}

const MyComplex<MySystGraph> operator-(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return (t1 + (-1.0*c));
}

const MyComplex<MySystGraph> operator-(const MyComplex<MyDualMultiv> &c, const MyComplex<MySystGraph> &t1) {
	return ((-1.0*t1) + c);
}

const MyComplex<MySystGraph> operator*(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return MyComplex<MySystGraph>(t1.GetRe()*c.GetRe()-t1.GetIm()*c.GetIm(), t1.GetRe()*c.GetIm()+t1.GetIm()*c.GetRe());
}

const MyComplex<MySystGraph> operator*(const MyComplex<MyDualMultiv> &c, const MyComplex<MySystGraph> &t1) {
	return (t1 * c);
}

const MyComplex<MySystGraph> operator/(const MyComplex<MySystGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return (t1 * (1.0/c));
}

const MyComplex<MySystGraph> operator/(const MyComplex<MyDualMultiv> &c, const MyComplex<MySystGraph> &t1) {
	return ((1.0/t1) * c);
}


const MyComplex<MyDualGraph> operator+(const MyComplex<MyDualGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return MyComplex<MyDualGraph>(t1.GetRe()+c.GetRe(), t1.GetIm()+c.GetIm());
}

const MyComplex<MyDualGraph> operator+(const MyComplex<MyDualMultiv> &c, const MyComplex<MyDualGraph> &t1) {
	return (t1 + c);
}

const MyComplex<MyDualGraph> operator-(const MyComplex<MyDualGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return (t1 + (-1.0*c));
}

const MyComplex<MyDualGraph> operator-(const MyComplex<MyDualMultiv> &c, const MyComplex<MyDualGraph> &t1) {
	return ((-1.0*t1) + c);
}

const MyComplex<MyDualGraph> operator*(const MyComplex<MyDualGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return MyComplex<MyDualGraph>(t1.GetRe()*c.GetRe()-t1.GetIm()*c.GetIm(), t1.GetRe()*c.GetIm()+t1.GetIm()*c.GetRe());
}

const MyComplex<MyDualGraph> operator*(const MyComplex<MyDualMultiv> &c, const MyComplex<MyDualGraph> &t1) {
	return (t1 * c);
}

const MyComplex<MyDualGraph> operator/(const MyComplex<MyDualGraph> &t1, const MyComplex<MyDualMultiv> &c) {
	return (t1 * (1.0/c));
}

const MyComplex<MyDualGraph> operator/(const MyComplex<MyDualMultiv> &c, const MyComplex<MyDualGraph> &t1) {
	return ((1.0/t1) * c);
}


#endif
