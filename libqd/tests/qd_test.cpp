/*
 * tests/qd_test.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This contains some simple tests to sanity check the double-double
 * and quad-double library.
 */

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <windows.h>

#include <qd/qd_real.h>
#include <qd/fpu.h>

using std::cout;
using std::cerr;
using std::endl;

using std::abs;
using std::sqrt;
using std::strcmp;
using std::exit;

// Global flags passed to the main program.
static bool flag_test_dd = false;
static bool flag_test_qd = false;
bool flag_verbose = false;

bool print_result(bool result) {
  if (result)
    cout << "Test passed." << endl;
  else
    cout << "Test FAILED." << endl;
  return result;
}

template <class T>
class TestSuite {
  static const int double_digits;
public:
  bool test1();
  bool test2();
  bool test3();
  bool test4();
  bool test5();
  bool test6();
  bool test7();
  bool test8();
  bool test9();
  bool test10();
  bool test11();
  bool test12();
  bool test13();
  bool test14();
  bool test15();
  bool test16();
  bool test17();
  bool testall();
};

template <class T>
const int TestSuite<T>::double_digits = 6;

/* Test 1.   Polynomial Evaluation / Polynomial Solving */
template <class T>
bool TestSuite<T>::test1() {
  cout << endl;
  cout << "Test 1.  (Polynomial)." << endl;

  static const int n = 8;
  T *c = new T[n];
  T x, y;

  for (int i = 0; i < n; i++)
    c[i] = static_cast<double>(i+1);

  x = polyroot(c, n-1, T(0.0));
  y = polyeval(c, n-1, x);

  if (flag_verbose) {
    cout.precision(T::_ndigits);
    cout << "Root Found:  x  = " << x << endl;
    cout << "           p(x) = " << y << endl;
  }

  delete [] c;
  return (to_double(y) < 4.0 * T::_eps);
}

/* Test 2.  Machin's Formula for Pi. */
template <class T>
bool TestSuite<T>::test2() {

  cout << endl;
  cout << "Test 2.  (Machin's Formula for Pi)." << endl;
  
  /* Use the Machin's arctangent formula:

       pi / 4  =  4 arctan(1/5) - arctan(1/239)

     The arctangent is computed based on the Taylor series expansion

       arctan(x) = x - x^3 / 3 + x^5 / 5 - x^7 / 7 + ...
  */

  T s1, s2, t, r;
  int k;
  int sign;
  double d;
  double err;

  /* Compute arctan(1/5) */
  d = 1.0;
  t = T(1.0) / 5.0;
  r = sqr(t);
  s1 = 0.0;
  k = 0;

  sign = 1;
  while (t > T::_eps) {
    k++;
    if (sign < 0)
      s1 -= (t / d);
    else
      s1 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    cout << k << " Iterations" << endl;

  /* Compute arctan(1/239) */
  d = 1.0;
  t = T(1.0) / 239.0;
  r = sqr(t);
  s2 = 0.0;
  k = 0;

  sign = 1;
  while (t > T::_eps) {
    k++;
    if (sign < 0)
      s2 -= (t / d);
    else
      s2 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    cout << k << " Iterations" << endl;

  T p = 4.0 * s1 - s2;

  p *= 4.0;
  err = abs(to_double(p - T::_pi));

  if (flag_verbose) {
    cout.precision(T::_ndigits);
    cout << "   pi = " << p << endl;
    cout << "  _pi = " << T::_pi << endl;

    cout.precision(double_digits);
    cout << "error = " << err << " = " << err / T::_eps << " eps" << endl;
  }

  return (err < 8.0 * T::_eps);
}

/* Test 3.  Salamin-Brent Quadratic Formula for Pi. */
template <class T>
bool TestSuite<T>::test3() {
  cout << endl;
  cout << "Test 3.  (Salamin-Brent Quadratic Formula for Pi)." << endl;
  cout.precision(T::_ndigits);

  T a, b, s, p;
  T a_new, b_new, p_old;
  double m;
  double err;
  const int max_iter = 20;

  a = 1.0;
  b = sqrt(T(0.5));
  s = 0.5;
  m = 1.0;

  p = 2.0 * sqr(a) / s;
  if (flag_verbose)
    cout << "Iteration  0: " << p << endl;
  for (int i = 1; i <= max_iter; i++) {
    m *= 2.0;
    a_new = 0.5 * (a + b);
    b_new = a * b;
    s -= m * (sqr(a_new) - b_new);
    a = a_new;
    b = sqrt(b_new);
    p_old = p;
    p = 2.0 * sqr(a) / s;
    if (flag_verbose)
      cout << "Iteration " << std::setw(2) << i << ": " << p << endl;
    if (abs(to_double(p - p_old)) < 64 * T::_eps)
      break;
  }

  err = abs(to_double(p - T::_pi));

  if (flag_verbose) {
    cout << "         _pi: " << T::_pi << endl;
    cout.precision(double_digits);
    cout << "       error: " << err << " = " << err / T::_eps << " eps" << endl;
  }

  // for some reason, this test gives relatively large error compared
  // to other tests.  May need to be looked at more closely.
  return (err < 1024.0 * T::_eps);
}

/* Test 4.  Borwein Quartic Formula for Pi. */
template <class T>
bool TestSuite<T>::test4() {
  cout << endl;
  cout << "Test 4.  (Borwein Quartic Formula for Pi)." << endl;
  cout.precision(T::_ndigits);

  T a, y, p, r, p_old;
  double m;
  double err;
  const int max_iter = 20;

  a = 6.0 - 4.0 * sqrt(T(2.0));
  y = sqrt(T(2.0)) - 1.0;
  m = 2.0;

  p = 1.0 / a;
  if (flag_verbose)
    cout << "Iteration  0: " << p << endl;

  for (int i = 1; i <= max_iter; i++) {
    m *= 4.0;
    r = nroot(1.0 - sqr(sqr(y)), 4);
    y = (1.0 - r) / (1.0 + r);
    a = a * sqr(sqr(1.0 + y)) - m * y * (1.0 + y + sqr(y));
    
    p_old = p;
    p = 1.0 / a;
    if (flag_verbose)
      cout << "Iteration " << std::setw(2) << i << ": " << p << endl;
    if (abs(to_double(p - p_old)) < 16 * T::_eps)
      break;
  }

  err = abs(to_double(p - T::_pi));
  if (flag_verbose) {
    cout << "         _pi: " << T::_pi << endl;
    cout.precision(double_digits);
    cout << "       error: " << err << " = " << err / T::_eps << " eps" << endl;
  }  

  return (err < 256.0 * T::_eps);
}

/* Test 5.  Taylor Series Formula for E. */
template <class T>
bool TestSuite<T>::test5() {

  cout << endl;
  cout << "Test 5.  (Taylor Series Formula for E)." << endl;
  cout.precision(T::_ndigits);

  /* Use Taylor series

       e = 1 + 1 + 1/2! + 1/3! + 1/4! + ...

     To compute e.
  */

  T s = 2.0, t = 1.0;
  double n = 1.0;
  double delta;
  int i = 0;

  while (t > T::_eps) {
    i++;
    n += 1.0;
    t /= n;
    s += t;
  }

  delta = abs(to_double(s - T::_e));

  if (flag_verbose) {
    cout << "    e = " << s << endl;
    cout << "   _e = " << T::_e << endl;

    cout.precision(double_digits);
    cout << "error = " << delta << " = " << delta / T::_eps << " eps" << endl;
    cout << i << " iterations." << endl;
  }

  return (delta < 64.0 * T::_eps);
}

/* Test 6.  Taylor Series Formula for log 2.*/
template <class T>
bool TestSuite<T>::test6() {
  cout << endl;
  cout << "Test 6.  (Taylor Series Formula for Log 2)." << endl;
  cout.precision(T::_ndigits);

  /* Use the Taylor series

      -log(1-x) = x + x^2/2 + x^3/3 + x^4/4 + ...

     with x = 1/2 to get  log(1/2) = -log 2.
  */

  T s = 0.5;
  T t = 0.5;
  double delta;
  double n = 1.0;
  double i = 0;

  while (abs(t) > T::_eps) {
    i++;
    n += 1.0;
    t *= 0.5;
    s += (t/n);
  }

  delta = abs(to_double(s - T::_log2));

  if (flag_verbose) {
    cout << " log2 = " << s << endl;
    cout << "_log2 = " << T::_log2 << endl;

    cout.precision(double_digits);
    cout << "error = " << delta << " = " << (delta / T::_eps) 
         << " eps" << endl;
    cout << i << " iterations." << endl;
  }

  return (delta < 4.0 * T::_eps);
}

/* Test 7.  Sanity check for exp. */
template <class T>
bool TestSuite<T>::test7() {
  cout << endl;
  cout << "Test 7.  (Sanity check for exp)." << endl;
  cout.precision(T::_ndigits);

  /* Do simple sanity check
   *
   *   e^2 = exp(2)
   *       = exp(-13/4) * exp(-9/4) * exp(-5/4) * exp(-1/4) *
   *         exp(3/4) * exp(7/4) * exp(11/4) * exp(15/4)
   */

  T t = -3.25;
  T p =  1.0;

  for (int i = 0; i < 8; i++, t += 1.0) {
    /* For some reason gcc-4.1.x on x86_64 miscompiles p *= exp(t) here. */
    p = p * exp(t);
  }

  T t1 = exp(T(2.0));
  T t2 = sqr(T::_e);
  double delta = std::max(abs(to_double(t1 - p)), abs(to_double(t2 - p)));

  if (flag_verbose) {
    cout << "result = " << p << endl;
    cout << "exp(2) = " << t1 << endl;
    cout << "   e^2 = " << t2 << endl;

    cout.precision(double_digits);

    cout << " error = " << delta << " = " << (delta / T::_eps)
         << " eps" << endl;
  }

  return (delta < 16.0 * T::_eps);
}

template <class T>
bool TestSuite<T>::test8() {
  cout << endl;
  cout << "Test 8.  (Sanity check for sin / cos)." << endl;
  cout.precision(T::_ndigits);

  /* Do simple sanity check
   *
   *  sin(x) = sin(5x/7)cos(2x/7) + cos(5x/7)sin(2x/7)
   *
   *  cos(x) = cos(5x/7)cos(2x/7) - sin(5x/7)sin(2x/7);
   */

  T x = T::_pi / 3.0;
  T x1 = 5.0 * x / 7.0;
  T x2 = 2.0 * x / 7.0;

  T r1 = sin(x1)*cos(x2) + cos(x1)*sin(x2);
  T r2 = cos(x1)*cos(x2) - sin(x1)*sin(x2);
  T t1 = sqrt(T(3.0)) / 2.0;
  T t2 = 0.5;

  double delta = std::max(abs(to_double(t1 - r1)), abs(to_double(t2 - r2)));

  if (flag_verbose) {
    cout << "  r1 = " << r1 << endl;
    cout << "  t1 = " << t1 << endl;
    cout << "  r2 = " << r2 << endl;
    cout << "  t2 = " << t2 << endl;

    cout.precision(double_digits);
    cout << " error = " << delta << " = " << (delta / T::_eps)
         << " eps" << endl;
  }

  return (delta < 4.0 * T::_eps);
}

template <class T>
bool TestSuite<T>::test9() {
	cout << endl;
	cout << "Test 9.  (nroot correctly check)." << endl;
	cout.precision(T::_ndigits);

	cout.precision(dd_real::_ndigits);
	dd_real dd_expected(-2.0);
	cout << "dd_real::nroot(-8.0, 3) =" << nroot(dd_real(-8.0), 3) << endl;

	cout << "dd_real::nroot( 4.0,-2) =" << nroot(dd_real(4.0), -2) << endl;
	cout << "dd_real::nroot(-4.0, 2) =" << nroot(dd_real(-4.0), 2) << endl;
	cout << "dd_real::nroot( 4.0, 1) =" << nroot(dd_real(4.0), 1) << endl;
	cout << "dd_real::nroot( 4.0, 2) =" << nroot(dd_real(4.0), 2) << endl;
	cout << "dd_real::nroot( 0.0, 3) =" << nroot(dd_real(0.0), 3) << endl;
	cout << endl;

	cout.precision(qd_real::_ndigits);
	qd_real qd_expected(-2.0);
	// should be return -2, but acutally rise nan
	cout << "qd_real::nroot(-8.0, 3) =" << nroot(qd_real(-8.0), 3) << endl;

	// should be rise nan, but return 0.5
	cout << "qd_real::nroot( 4.0,-2) =" << nroot(qd_real(4.0), -2) << endl;
	cout << "qd_real::nroot(-4.0, 2) =" << nroot(qd_real(-4.0), 2) << endl;
	cout << "qd_real::nroot( 4.0, 1) =" << nroot(qd_real(4.0), 1) << endl;
	cout << "qd_real::nroot( 4.0, 2) =" << nroot(qd_real(4.0), 2) << endl;
	cout << "qd_real::nroot( 0.0, 3) =" << nroot(qd_real(0.0), 3) << endl;

	//dd_expected.dump_bits("dd_expected:", cout);
	//(nroot(dd_real(-8.0), 3)).dump_bits("dd_computed:", cout);
	//qd_expected.dump_bits("qd_expected:", cout);
	//(nroot(qd_real(-8.0), 3)).dump_bits("qd_computed:", cout);

	return (dd_expected == nroot(dd_real(-8.0), 3)
			&& qd_expected == nroot(qd_real(-8.0), 3));
}

template <class T>
bool TestSuite<T>::test10() {
	cout << endl;
	cout << "Test 10.  ( tanh correctly check)." << endl;
	cout << "comparing dd_real::tanh(rad) to qd_real::tanh(rad) in case of rad<=0.05." << endl;

	cout.precision(dd_real::_ndigits);
	bool flag = true;
	dd_real dd_ten(10);
	qd_real qd_ten(10);
	dd_real dd_base(0.5);
	qd_real qd_base(0.5);
	dd_real dd_tanhx, dd_tanhx2;
	qd_real qd_tanhx, qd_tanhx2;
	
	qd_base = dd_base;
	while (dd_base > dd_real::_min_normalized){
				
		dd_tanhx = tanh(dd_base);
		qd_tanhx = tanh(qd_base);
		
		cout << "base     : " << dd_base << (dd_base == dd_tanhx ? " ***" : "") << endl;
		cout << "dd_tanhx : " << dd_tanhx << endl;
		cout << "qd_tanhx : " << qd_tanhx << endl;
		cout << endl;

		if (qd_tanhx - dd_tanhx > dd_tanhx * dd_real::_eps){
			flag = false;
		}

		dd_base = dd_base / dd_ten;
		qd_base = dd_base;

	} 

	return (flag);
}

template <class T>
bool TestSuite<T>::test11() {
	cout << endl;
	cout << "Test 11.  (sincosh correctly check)." << endl;

	cout.precision(dd_real::_ndigits);
	bool flag = true;
	dd_real dd_ten(10);
	qd_real qd_ten(10);
	dd_real dd_base(0.5);
	qd_real qd_base(0.5);
	dd_real dd_coshx, dd_sinhx, t1,t2,t3,t4;
	qd_real qd_coshx, qd_sinhx;
	dd_real dd_coshx2, dd_sinhx2, dd_coshx3, dd_sinhx3;
	dd_real ONE(1.0);

	qd_base = dd_base;
	while (dd_base > dd_real::_min_normalized){

		sincosh(dd_base, dd_sinhx, dd_coshx);
		sincosh(qd_base, qd_sinhx, qd_coshx);

		dd_coshx2 = cosh(dd_base);		// caluculate indivisually
		dd_sinhx2 = sinh(dd_base);		// caluculate indivisually

		dd_coshx3 = sqrt(ONE + dd_sinhx2*dd_sinhx2);		// using formula
		

		// for  "cosh^2 - sinh^2 = 1"
		t1 = (dd_coshx*dd_coshx - dd_sinhx*dd_sinhx);		// by present sincosh implement
		//t2 = (sqr(dd_coshx) - sqr(dd_sinhx));
		t3 = (dd_coshx2*dd_coshx2 - dd_sinhx2*dd_sinhx2);	// calculate indivisually
		t4 = (dd_coshx3*dd_coshx3 - dd_sinhx2*dd_sinhx2);	// using formula

		
		// how near to 1.0 
		cout << "base    :" << dd_base  << endl;
		cout << "dd_sinh :" << dd_sinhx2 << endl;
		cout << "qd_sinh :" << qd_sinhx << endl;
		cout << "dd_cosh :" << dd_coshx2 << endl;
		cout << "qd_cosh :" << qd_coshx << endl;
		cout << "calc indivisually   (cosh^2 - sinh^2):" << t3 << " diff:" << (t3 - ONE) << endl;
		cout << "calc using formula  (cosh^2 - sinh^2):" << t4 << " diff:" << (t4 - ONE) << endl;
		//cout << "calc by sincosh impl(cosh^2 - sinh^2):" << t1 << " diff:" << (t1 - ONE) << endl;
		//cout << "calc by sincosh    (cosh^2 - sinh^2):" << t2 << " diff:" << (t2 - ONE) << endl;

		cout << endl;

		if ( abs(t3-ONE) > t3*dd_real::_eps ){
			flag = false;
		}

		dd_base = dd_base / dd_ten;
		qd_base = dd_base;
	}

	return (flag);
}

template <class T>
bool TestSuite<T>::test12() {
	cout << endl;
	cout << "Test 12.  sqrt check." << endl;
	cout.precision(T::_ndigits);

	int digits = T::_ndigits;
	if (digits > 60){
		cout.precision(qd_real::_ndigits);
		qd_real x = qd_real::_pi;
		x = x * 60 / 180;
		qd_real s = sin(x);
		qd_real r = sqrt(1-sqr(cos(x)));

		cout << "s=" << s << endl;
		cout << "r=" << r << endl;
		std::cout << endl;

		s.dump_bits("s", std::cout);
		std::cout << endl;
		r.dump_bits("r", std::cout);
		std::cout << endl;

		s.dump_bits2("s", std::cout);
		r.dump_bits2("r", std::cout);
		std::cout << endl;

		qd_real err = s * qd_real::_eps;
		return (abs(s-r) <= err);

	}else{
		cout.precision(dd_real::_ndigits);
		dd_real x = dd_real::_pi;
		x = x * 60 / 180;
		dd_real s = sin(x);
		dd_real r = sqrt(1 - sqr(cos(x)));

		cout << "s=" << s << endl;
		cout << "r=" << r << endl;
		std::cout << endl;

		s.dump_bits("s", std::cout);
		std::cout << endl;
		r.dump_bits("r", std::cout);
		std::cout << endl;

		s.dump_bits2("s", std::cout);
		r.dump_bits2("r", std::cout);
		std::cout << endl;

		dd_real err = s * dd_real::_eps;
		return (abs(s - r) <= err);

	}
	

	
}

template <class T>
bool TestSuite<T>::test13() {
	cout << endl;
	cout << "Test 13.  eps value check." << endl;
	cout.precision(T::_ndigits);
	bool status = false;
	union trans{
		unsigned __int64 asInt64;
		double asDouble;
	};
	unsigned __int64 dd_eps = 0x3970000000000000ULL;	// 2^-104
	unsigned __int64 qd_eps = 0x32e0000000000000ULL;	// 2^-209

	trans dd_t, qd_t;
	dd_t.asDouble = dd_real::_eps;
	qd_t.asDouble = qd_real::_eps;

	int digits = T::_ndigits;
	if (digits > 60) {
		if (qd_t.asInt64 != qd_eps) {
			cout << "eps of qd_real incorrect " << std::hex << qd_t.asInt64 << " expect(0x32e0000000000000)" << endl;
		} else {
			cout << "eps of qd_real correct " << std::hex << qd_t.asInt64 << endl;
			status = true;
		}
	} else {
		if (dd_t.asInt64 != dd_eps) {
			cout << "eps of dd_real incorrect " << std::hex << dd_t.asInt64 << " expect(0x3970000000000000)" << endl;
		} else {
			cout << "eps of dd_real correct " << std::hex << dd_t.asInt64 << endl;
			status = true;
		}
	}
	return status;
}

template <class T>
T sinhx2(const T &a) {
	T ea = exp(a);
	return mul_pwr2(ea - inv(ea), 0.5);

}

template <class T>
T sinhx3(const T &a) {
	T s = a;
	T t = a;
	T r = sqr(t);
	double m = 1.0;
	double thresh = std::abs(to_double(a) * T::_eps);

	do {
		m += 2.0;
		t *= r;
		t /= (m - 1) * m;

		s += t;
	} while (abs(t) > thresh);

	return s;

}


template <class T>
bool TestSuite<T>::test14() {
	cout << endl;
	cout << "Test 14.  sinh(x) check." << endl;
	cout.precision(T::_ndigits);



	int digits = T::_ndigits;
	if (digits > 60) {
		cout.precision(qd_real::_ndigits);
		
		qd_real r2, r3;
		qd_real a(1.0);
		a *= qd_real::_pi;
		qd_real n(2.0);

		for (int i = 0; i < 10; i++) {
			r2 = sinhx2(a);
			r3 = sinhx3(a);

			cout << " a=" << a << endl;
			cout << "r2=" << r2 << endl;
			cout << "r3=" << r3 << endl;

			showdiff(r2.x[3], r3.x[3]);
			std::cout << endl;

			a /= n;
		}
		return (r2 == r3);


	} else {
		cout.precision(dd_real::_ndigits);

		dd_real r2, r3;
		dd_real a(1.0);
		a *= dd_real::_pi;
		dd_real n(2.0);

		for (int i = 0; i < 10; i++) {
			r2 = sinhx2(a);
			r3 = sinhx3(a);

			cout << " a=" << a << endl;
			cout << "r2=" << r2 << endl;
			cout << "r3=" << r3 << endl;

			showdiff(r2.x[1], r3.x[1]);
			std::cout << endl;

			a /= n;
		}
		return (r2 == r3);

	}



}

inline dd_real sqr2(const dd_real &a) {
	double p1, p2, t;
	double s1, s2;
	p1 = qd::two_sqr(a.x[0], t);
	p2 = a.x[1] * a.x[1];
	p2 += 2.0 * a.x[0] * a.x[1];
	p2 += t;
	s1 = qd::quick_two_sum(p1, p2, s2);
	return dd_real(s1, s2);
}


void showdiff(double a, double b) { 
	unsigned __int64 diff, mask;
	int cnt=0;
	union trans {
		unsigned __int64 asInt64;
		double asDouble;
	};
	trans aa, bb;
	aa.asDouble = a;
	bb.asDouble = b;
	if (aa.asInt64 > bb.asInt64) {
		diff = aa.asInt64 - bb.asInt64;
	} else {
		diff = bb.asInt64 - aa.asInt64;
	}
	mask = 0xffffffffffffffffULL;
	for (int j = 0; j < 53; j++) {
		if ((aa.asInt64 & mask) == (bb.asInt64 & mask)) {
			cnt = j;
			break;
		}
		mask = mask << 1;
	}
	cout << "diff=" << std::dec << diff << " ulps(" << cnt << "bits)" << endl;

}

template <class T>
bool TestSuite<T>::test15() {
	cout << endl;
	cout << "Test 15.  sqr(x) acuracy check." << endl;
	cout.precision(T::_ndigits);



	int digits = T::_ndigits;
	if (digits > 60) {
		cout.precision(qd_real::_ndigits);
		//qd_real r0, r1, r2;
		//qd_real a = qd_real::_pi;


		//r0 = a * a;
		//r1 = sqr(a);
		//r2 = sqr2(a);

		//cout << "   a*a =" << r0 << endl;
		//cout << " sqr(a)=" << r1 << endl;
		//cout << "sqr2(a)=" << r2 << endl;


		//showdiff(r0.x[3], r1.x[3]);
		//showdiff(r0.x[3], r2.x[3]);
		//std::cout << endl;

		//return (r2 == r1);
		return true;


	} else {
		cout.precision(dd_real::_ndigits);

		dd_real r0, r1, r2;
		dd_real a = dd_real::_pi;


		r0 = a * a;
		r1 = sqr(a);
		r2 = sqr2(a);

		cout << "   a*a =" << r0 << endl;
		cout << " sqr(a)=" << r1 << endl;
		cout << "sqr2(a)=" << r2 << endl;

		showdiff(r0.x[1], r1.x[1]);
		showdiff(r0.x[1], r2.x[1]);
		std::cout << endl;

		return (r2 == r1);

	}
}

inline dd_real fast_div(const dd_real &a, const dd_real &b) {
/*
	yamanaka, ohishi
	Zh + Zl near= (Ah + Al) / (Bh + Bl)
	            = Ah(1 + Al / Ah) / Bh(1 + Bl / Bh)
	            = (Ah / Bh)(1 + Al / Ah - Bl / Bh) (approx)
*/
	double zh, zl;
	double th, tl, br, cr;

	cr = 1.0 / b.x[0];
	br = b.x[1] * cr;

	zh = a.x[0] * cr;
	th = qd::two_prod(zh, b.x[0], tl);
	zl = ((a.x[0] - th) - tl) * cr + zh *(a.x[1] / a.x[0] - br);
	zh = qd::quick_two_sum(zh, zl, zl);
	return dd_real(zh, zl);
}

inline dd_real fast_sqrt(const dd_real &a) {
/*
	yamanaka, ohishi
	sqrt(Ah + Al) = sqrt(Ah * (1 + (Al / Ah)))
	              = sqrt(Ah) * (1 + (Al / Ah) * 0.5) (approx)
	
*/
	double th, tl;
	double app = std::sqrt(a.x[0]);
	dd_real r;
	th = qd::two_sqr(app, tl);
	r.x[1] = 0.5 * (((a.x[0]-th)-tl)+a.x[1]) / app;
	r.x[0] = qd::quick_two_sum(app, r.x[1], r.x[1]);
	return r;
}

qd_real qd_fast_sqrt(const qd_real &a) {
	/* Strategy:
	Perform the following Newton iteration:

	x' = x + (1 - a * x^2) * x / 2;

	which converges to 1/sqrt(a), starting with the
	double precision approximation to 1/sqrt(a).
	Since Newton's iteration more or less doubles the
	number of correct digits, we only need to perform it
	twice.
	*/

	if (a.is_zero()) {
		return a;
	}

	if (a.is_negative()) {
		qd_real::error("(qd_real::sqrt): Negative argument.");
		return qd_real::_nan;
	}

	if (a.isinf()) {
		return a;
	}

	//qd_real r = (1.0 / std::sqrt(a[0]));
	qd_real r = qd_real(1.0) / fast_sqrt(dd_real(a[0], a[1]));
	qd_real h = mul_pwr2(a, 0.5);

	r += ((0.5 - h * sqr(r)) * r);
	r += ((0.5 - h * sqr(r)) * r);
//	r += ((0.5 - h * sqr(r)) * r);

	r *= a;
	return r;
}

inline dd_real ieee_sub(const dd_real &a, const dd_real &b) {
	/* This one satisfies IEEE style error bound,
	due to K. Briggs and W. Kahan.                   */
	double s1, s2, t1, t2;

	s1 = qd::two_sum(a.x[0], -(b.x[0]), s2);
	t1 = qd::two_sum(a.x[1], -(b.x[1]), t2);
	s2 += t1;
	s1 = qd::quick_two_sum(s1, s2, s2);
	s2 += t2;
	s1 = qd::quick_two_sum(s1, s2, s2);
	return dd_real(s1, s2);
}


class Qtimer {
private:
	LARGE_INTEGER _freq, _st, _ed;
	unsigned __int64 st_clk, ed_clk;
public:
	double result;
	unsigned __int64 clk_total;
	unsigned int ui;

	Qtimer() {
		SetThreadAffinityMask(GetCurrentThread(), 0x1);
		if (!QueryPerformanceFrequency(&_freq)) {
			std::cout << "Not supported Precise counter." << std::endl;
		}
	}

	~Qtimer() { };

	void start() {
		SetThreadAffinityMask(GetCurrentThread(), 0x1);
		QueryPerformanceCounter(&_st);
		st_clk = __rdtscp(&ui);
	}

	void stop(std::string str="") {
		ed_clk = __rdtscp(&ui);
		QueryPerformanceCounter(&_ed);
		clk_total = ed_clk - st_clk;
		result = 1000.0 * (_ed.QuadPart - _st.QuadPart)/_freq.QuadPart;
		
		cout << str;
		std::ios::fmtflags flags = cout.flags();
		cout << std::fixed << std::setprecision(6);
		std::cout << result << " ms  " << clk_total<< " clocks" << std::endl;
		cout.flags(flags);
	}

};

template <class T>
bool TestSuite<T>::test16() {
	cout << endl;
	cout << "Test 16.  fast_div and fast_sqrt acuracy check." << endl;
	cout.precision(T::_ndigits);



	int digits = T::_ndigits;
	if (digits > 60) {
		cout.precision(qd_real::_ndigits);

		qd_real r0, r1, r2;
		qd_real a = qd_real::_pi;


		r0 = a * a;
		r1 = r0 / a;
		//r2 = fast_div(r0, a);

		cout << "              a =" <<  a << endl;
		cout << "         r0 / a =" << r1 << endl;
		//cout << "fast_div(r0, a) =" << r2 << endl;

		showdiff(a.x[3], r1.x[3]);
		//showdiff(tr0.asDouble[3], tr2.asDouble[3]);
		std::cout << endl;

		return (r2 == r1);


	} else {
		cout.precision(dd_real::_ndigits);


		dd_real r0, r1, r2, r3, r4;
		dd_real a = dd_real::_pi;
		dd_real b(1.0);
		dd_real c(10.0);

		r0 = a * a;
		r1 = r0 / a;
		r2 = fast_div(r0, a);

		r3 = b / c;
		r4 = fast_div(b, c);

		cout << "              a =" << a << endl;
		cout << "         r0 / a =" << r1 << endl;
		cout << "fast_div(r0, a) =" << r2 << endl;
		showdiff(r1.x[1], r2.x[1]);
		std::cout << endl;

		cout << "         1 / 10 =" << r3 << endl;
		cout << "fast_div(1, 10) =" << r4 << endl;
		showdiff(r3.x[1], r4.x[1]);
		std::cout << endl;

		r1 = sqrt(r0);
		r2 = fast_sqrt(r0);
		cout << "                a =" << a << endl;
		cout << "        sqrt(a^2) =" << r1 << endl;
		cout << "fast_sqrt(a^2, a) =" << r2 << endl;
		showdiff(a.x[1], r1.x[1]);
		showdiff(a.x[1], r2.x[1]);
		showdiff(r1.x[1], r2.x[1]);
		std::cout << endl;

		return (r2 == r1);

	}
}

template <class T>
bool TestSuite<T>::test17() {
	cout << endl;
	cout << "Test 17.  Basic function benchmark." << endl;
	cout.precision(T::_ndigits);



	cout.precision(dd_real::_ndigits);

	int count = 1000;

	double n_arry[1000];
	double n_arry_r[1000];
	for (int i = 0; i < count; ++i) {
		n_arry[i] = rand();
	}
	double n_divider = 2.718281828459045091e+00;

	dd_real d_arry[1000];
	dd_real d_arry_r[1000];
	for (int i = 0; i < count; ++i) {
		d_arry[i] = ddrand();
	}
	dd_real d_divider = dd_real::_e;


	qd_real q_arry[1000];
	qd_real q_arry_r[1000];
	for (int i = 0; i < count; ++i) {
		q_arry[i] = qdrand();
	}
	qd_real q_divider = qd_real::_e;

	Qtimer dd_add_t, dd_sub_t, dd_prod_t, dd_div_t, dd_sqrt_t;
	Qtimer fast_sub_t, fast_div_t, fast_sqrt_t;
	Qtimer qd_add_t, qd_sub_t, qd_prod_t, qd_div_t, qd_sqrt_t, qd_fast_sqrt_t;
	Qtimer _add_t, _sub_t, _prod_t, _div_t, _sqrt_t;
	SetThreadAffinityMask(GetCurrentThread(), 0x1);

	dd_real r0, r1, r2, r3, r4;
	dd_real da = dd_real::_pi;
	qd_real qa = qd_real::_pi;
	dd_real b(1.0);
	dd_real c(10.0);



	// add
	_add_t.start();
	for (int i = 0; i < count; ++i) {
		n_arry_r[i] = n_arry[i] + n_divider;
	}
	_add_t.stop("native add time :");

	dd_add_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = d_arry[i] + d_divider;
	}
	dd_add_t.stop("    dd add time :");

	qd_add_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = q_arry[i] + q_divider;
	}
	qd_add_t.stop("    qd add time :");

	cout << std::fixed << std::setprecision(1);
	cout << "n_add : dd_add : qd_add = 1 : " << dd_add_t.result / _add_t.result << " : " << qd_add_t.result / _add_t.result << " time" << endl;
	cout << "n_add : dd_add : qd_add = 1 : " << (double)dd_add_t.clk_total / _add_t.clk_total << " : " << (double)qd_add_t.clk_total / _add_t.clk_total << " clock" << endl;
	cout << endl;

	// sub
	_sub_t.start();
	for (int i = 0; i < count; ++i) {
		n_arry_r[i] = n_arry[i] - n_divider;
	}
	_sub_t.stop("native sub time :");

	dd_sub_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = d_arry[i] - d_divider;
	}
	dd_sub_t.stop("    dd sub time :");

	fast_sub_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = ieee_sub(d_arry[i], d_divider);
	}
	fast_sub_t.stop("ddfast sub time :");

	qd_sub_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = q_arry[i] - q_divider;
	}
	qd_sub_t.stop("    qd sub time :");

	cout << std::fixed << std::setprecision(1);
	cout << "n_sub : dd_sub : fast_sub : qd_sub = 1 : " << dd_sub_t.result / _sub_t.result << " : " << fast_sub_t.result / _sub_t.result << " : " << qd_sub_t.result / _sub_t.result << " time" << endl;
	cout << "n_sub : dd_sub : fast_sub : qd_sub = 1 : " << (double)dd_sub_t.clk_total / _sub_t.clk_total << " : " << (double)fast_sub_t.clk_total / _sub_t.clk_total << " : " << (double)qd_sub_t.clk_total / _sub_t.clk_total << " clock" << endl;
	cout << endl;

	// prod
	_prod_t.start();
	for (int i = 0; i < count; ++i) {
		n_arry_r[i] = n_arry[i] * n_divider;
	}
	_prod_t.stop("native prod time :");

	dd_prod_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = d_arry[i] * d_divider;
	}
	dd_prod_t.stop("    dd prod time :");

	qd_prod_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = q_arry[i] * q_divider;
	}
	qd_prod_t.stop("    qd prod time :");

	cout << std::fixed << std::setprecision(1);
	cout << "n_prod : dd_prod : qd_prod = 1 : " << dd_prod_t.result / _prod_t.result << " : " << qd_prod_t.result / _prod_t.result << " time" << endl;
	cout << "n_prod : dd_prod : qd_prod = 1 : " << (double)dd_prod_t.clk_total / _prod_t.clk_total << " : " << (double)qd_prod_t.clk_total / _prod_t.clk_total << " clock" << endl;
	cout << endl;


	// div
	_div_t.start();
	for (int i = 0; i < count; ++i) {
		n_arry_r[i] = n_arry[i] / n_divider;
	}
	_div_t.stop("native div time :");

	dd_div_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = d_arry[i] / d_divider;
	}
	dd_div_t.stop("    dd div time :");

	fast_div_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = fast_div(d_arry[i], d_divider);
	}
	fast_div_t.stop("ddfast div time :");

	qd_div_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = q_arry[i] / q_divider;
	}
	qd_div_t.stop("    qd div time :");

	cout << std::fixed << std::setprecision(1);
	cout << "n_div : dd_div : fast_div : qd_div = 1 : " << dd_div_t.result / _div_t.result << " : " << fast_div_t.result / _div_t.result << " : "<< qd_div_t.result / _div_t.result << " time" << endl;
	cout << "n_div : dd_div : fast_div : qd_div = 1 : " << (double)dd_div_t.clk_total / _div_t.clk_total << " : " << (double)fast_div_t.clk_total / _div_t.clk_total << " : "<< (double)qd_div_t.clk_total/ _div_t.clk_total << " clock" << endl;
	cout << endl;


	// sqrt
	_sqrt_t.start();
	for (int i = 0; i < count; ++i) {
		n_arry_r[i] = std::sqrt(n_arry[i]);
	}
	_sqrt_t.stop("native sqrt time :");

	dd_sqrt_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = sqrt(d_arry[i]);
	}
	dd_sqrt_t.stop("    dd sqrt time :");

	fast_sqrt_t.start();
	for (int i = 0; i < count; ++i) {
		d_arry_r[i] = fast_sqrt(d_arry[i]);
	}
	fast_sqrt_t.stop("ddfast sqrt time :");

	qd_sqrt_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = sqrt(q_arry[i]);
	}
	qd_sqrt_t.stop("    qd sqrt time :");

	qd_fast_sqrt_t.start();
	for (int i = 0; i < count; ++i) {
		q_arry_r[i] = qd_fast_sqrt(q_arry[i]);
	}
	qd_fast_sqrt_t.stop("qd fastsqrt time :");


	cout << std::fixed << std::setprecision(1);
	cout << "n_sqrt : dd_sqrt : fast_sqrt : qd_sqrt : qd_fast_sqrt = 1 : " << dd_sqrt_t.result / _sqrt_t.result << " : " << fast_sqrt_t.result / _sqrt_t.result << " : " << qd_fast_sqrt_t.result / _sqrt_t.result << " : " << qd_sqrt_t.result / _sqrt_t.result << " time" << endl;
	cout << "n_sqrt : dd_sqrt : fast_sqrt : qd_sqrt : qd_fast_sqrt = 1 : " << (double)dd_sqrt_t.clk_total / _sqrt_t.clk_total << " : " << (double)fast_sqrt_t.clk_total / _sqrt_t.clk_total << " : " << (double)qd_fast_sqrt_t.clk_total / _sqrt_t.clk_total << " : " << (double)qd_sqrt_t.clk_total/_sqrt_t.clk_total << " clock" << endl;
	cout << endl;

	return (r2 == r1);

	
}


template <class T>
bool TestSuite<T>::testall() {
  bool pass = true;
  pass &= print_result(test1());
  pass &= print_result(test2());
  pass &= print_result(test3());
  pass &= print_result(test4());
  pass &= print_result(test5());
  pass &= print_result(test6());
  pass &= print_result(test7());
  pass &= print_result(test8());
  pass &= print_result(test9());
  pass &= print_result(test10());
  pass &= print_result(test11());
  pass &= print_result(test12());
  pass &= print_result(test13());
  pass &= print_result(test14());
  pass &= print_result(test15());
  pass &= print_result(test16());
  pass &= print_result(test17());
  return pass;
}

void print_usage() {
  cout << "qd_test [-h] [-dd] [-qd] [-all]" << endl;
  cout << "  Performs miscellaneous tests of the quad-double library," << endl;
  cout << "  such as polynomial root finding, computation of pi, etc." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -dd       Perform tests with double-double types." << endl;
  cout << "  -qd       Perform tests with quad-double types." << endl;
  cout << "            This is the default." << endl;
  cout << "  -all      Perform both double-double and quad-double tests." << endl;
  cout << "  -v" << endl;
  cout << "  -verbose  Print detailed information for each test." << endl;
  
}

int main(int argc, char *argv[]) {
  
  bool pass = true;
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      exit(0);
    } else if (strcmp(arg, "-dd") == 0) {
      flag_test_dd = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_test_qd = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_dd = flag_test_qd = true;
    } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "-verbose") == 0) {
      flag_verbose = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_dd && !flag_test_qd) {
    flag_test_dd = true;
    flag_test_qd = true;
  }

  if (flag_test_dd) {
    TestSuite<dd_real> dd_test;

    cout << endl;
    cout << "Testing dd_real ..." << endl;
    if (flag_verbose)
      cout << "sizeof(dd_real) = " << sizeof(dd_real) << endl;
    pass &= dd_test.testall();
  }

  if (flag_test_qd) {
    TestSuite<qd_real> qd_test;

    cout << endl;
    cout << "Testing qd_real ..." << endl;
    if (flag_verbose)
      cout << "sizeof(qd_real) = " << sizeof(qd_real) << endl;
    pass &= qd_test.testall();
  }
  
  fpu_fix_end(&old_cw);
  return (pass ? 0 : 1);
}

