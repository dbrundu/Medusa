/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu,
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto,
 *                      Alessandro Maria Ricci
 *
 *   This file is part of Medusa Library.
 *
 *   Medusa is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Medusa is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Medusa.  If not, see <http://www.gnu.org/licenses/>.
 *
 *---------------------------------------------------------------------------*/
/*----------------------------------------
 *  Created on: 17/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  Constants.h
 *----------------------------------------*/

#ifndef MEDUSA_CONSTANTS_H
#define MEDUSA_CONSTANTS_H

//////////////////////////////////////////////////////////////////////////////
/* If this file is compiled as a part of a larger project,
   support using an autoconf-style config.h header file
   (with various "HAVE_*" #defines to indicate features)
   if HAVE_CONFIG_H is #defined (in GNU autotools style).

   If HAVE_CONFIG_H is #defined (e.g. by compiling with -DHAVE_CONFIG_H),
   then we #include "config.h", which is assumed to be a GNU autoconf-style
   header defining HAVE_* macros to indicate the presence of features. In
   particular, if HAVE_ISNAN and HAVE_ISINF are #defined, we use those
   functions in math.h instead of defining our own, and if HAVE_ERF and/or
   HAVE_ERFC are defined we use those functions from <cmath> for erf and
   erfc of real arguments, respectively, instead of defining our own. */
/*
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
*/
//////////////////////////////////////////////////////////////////////////////

// Hydra
#include <hydra/Types.h>

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#define Inf std::numeric_limits<double>::infinity()
#define NaN std::numeric_limits<double>::quiet_NaN()

// isnan/isinf were introduced in C++11
#  if (__cplusplus < 201103L) && (!defined(HAVE_ISNAN) || !defined(HAVE_ISINF))
static inline bool my_isnan(double x) { return x != x; }
#    define isnan my_isnan
static inline bool my_isinf(double x) { return 1/x == 0.; }
#    define isinf my_isinf
#  elif (__cplusplus >= 201103L)
// g++ gets confused between the C and C++ isnan/isinf functions
#    define isnan std::isnan
#    define isinf std::isinf
#  endif

// copysign was introduced in C++11 (and is also in POSIX and C99)
#  if defined(_WIN32) || defined(__WIN32__)
#    define copysign _copysign // of course MS had to be different
#  elif defined(GNULIB_NAMESPACE) // we are using using gnulib <cmath>
#    define copysign GNULIB_NAMESPACE::copysign
#  elif (__cplusplus < 201103L) && !defined(HAVE_COPYSIGN) && !defined(__linux__) && !(defined(__APPLE__) && defined(__MACH__)) && !defined(_AIX)
static inline double my_copysign(double x, double y) { return x<0 != y<0 ? -x : x; }
#    define copysign my_copysign
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
//////////////////////////////////////////////////////////////////////////////


namespace medusa {
    namespace math_constants {

		const double M_1_3      = 0.33333333333333333333;  // 1/3
        const double M_1_Sqrt8  = 0.35355339059327376220;  // 1/Sqrt(8)
        const double M_1_Sqrt32 = 0.17677669529663688110;  // 1/sqrt(32)
        const double M_2_Sqrt3  = 1.15470053837925152902;  // 2/sqrt(3)
		const double M_2_Sqrt6  = 0.81649658092772603273;  // 2/sqrt(6)
        const double M_1_SqrtPi = 0.56418958354775628695;  // 1/sqrt(Pi)

    } // namespace math_constants
} // namespace medusa


// macro
#define M_Sqrt2 hydra::math_constants::sqrt2
#define M_1_Sqrt2 hydra::math_constants::inverse_sqrt2
#define M_1_3 medusa::math_constants::M_1_3
#define M_1_Sqrt8 medusa::math_constants::M_1_Sqrt8
#define M_1_Sqrt32 medusa::math_constants::M_1_Sqrt32
#define M_2_Sqrt3 medusa::math_constants::M_2_Sqrt3
#define M_2_Sqrt6 medusa::math_constants::M_2_Sqrt6
#define M_1_SqrtPi medusa::math_constants::M_1_SqrtPi

#endif // MEDUSA_CONSTANTS_H
