/*----------------------------------------------------------------------------
 * Copyright (c) 2012 Massachusetts Institute of Technology
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 *----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------
 *  Created: 11/11/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  This header file is a version of Faddeeva.hh of the Faddeeva package.
 * 
 *  The original version is available at: http://ab-initio.mit.edu/Faddeeva
 *
 *  Header file for Faddeeva.inl; see that file for more information.
 *----------------------------------------------------------------------------*/

#ifndef MEDUSA_FADDEEVA_H
#define MEDUSA_FADDEEVA_H

// std
#include <complex>

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/Complex.h>
#include <hydra/multiarray.h>

/////////////////////////////////////////////////////////////////////////
// C++ settings

#include <cfloat>
#include <cmath>
#include <limits>

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#  define Inf std::numeric_limits<double>::infinity()
#  define NaN std::numeric_limits<double>::quiet_NaN()

typedef hydra::complex<double> cmplx;

// Use C-like complex syntax, since the C syntax is more restrictive
#  define cexp(z) hydra::exp(z)
#  define creal(z) z.cmplx::real()
#  define cimag(z) z.cmplx::imag()
#  define cpolar(r,t) hydra::polar(r,t)

#  define C(a,b) cmplx(a,b)

#  define FADDEEVA(name) medusa::faddeeva::name
#  define FADDEEVA_RE(name) medusa::faddeeva::name

// g++ gets confused between the C and C++ isnan/isinf functions
#  define isnan std::isnan
#  define isinf std::isinf

// copysign was introduced in C++11 (and is also in POSIX and C99)
#  if defined(_WIN32) || defined(__WIN32__)
#    define copysign _copysign // of course MS had to be different
#  elif defined(GNULIB_NAMESPACE) // we are using using gnulib <cmath>
#    define copysign GNULIB_NAMESPACE::copysign
#  endif

// If we are using the gnulib <cmath> (e.g. in the GNU Octave sources),
// gnulib generates a link warning if we use ::floor instead of gnulib::floor.
// This warning is completely innocuous because the only difference between
// gnulib::floor and the system ::floor (and only on ancient OSF systems)
// has to do with floor(-0), which doesn't occur in the usage below, but
// the Octave developers prefer that we silence the warning.
#  ifdef GNULIB_NAMESPACE
#    define floor GNULIB_NAMESPACE::floor
#  endif
/////////////////////////////////////////////////////////////////////////

namespace medusa{

   namespace faddeeva {

      // compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
      __hydra_dual__
      extern cmplx w(cmplx z,double relerr=0);
      __hydra_dual__
      extern double w_im(double x); // special-case code for Im[w(x)] of real x

      // Various functions that we can compute with the help of w(z)

      // compute erfcx(z) = exp(z^2) erfc(z)
      __hydra_dual__
      extern cmplx erfcx(cmplx z, double relerr=0);
      __hydra_dual__
      extern double erfcx(double x); // special case for real x

      // compute erf(z), the error function of complex arguments
      __hydra_dual__
      extern cmplx erf(cmplx z, double relerr=0);
      __hydra_dual__
      extern double erf(double x); // special case for real x

      // compute erfi(z) = -i erf(iz), the imaginary error function
      __hydra_dual__
      extern cmplx erfi(cmplx z, double relerr=0);
      __hydra_dual__
      extern double erfi(double x); // special case for real x

      // compute erfc(z) = 1 - erf(z), the complementary error function
      __hydra_dual__
      extern cmplx erfc(cmplx z, double relerr=0);
      __hydra_dual__
      extern double erfc(double x); // special case for real x

      // compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
      __hydra_dual__
      extern cmplx Dawson(cmplx z, double relerr=0);
      __hydra_dual__
      extern double Dawson(double x); // special case for real x

   } // namespace faddeeva

} // namespace medusa

#include "Faddeeva.inl"

#endif // MEDUSA_FADDEEVA_H