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
/*---------------------------------------------------------------------------
 *  Created on: 06/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  This library contains generic functions.
 *---------------------------------------------------------------------------*/

#ifndef MEDUSA_FUNCTIONS_H_
#define MEDUSA_FUNCTIONS_H_

// std
#include <cmath>
#include <complex>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <ratio>

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/Complex.h>
#include <hydra/Vector4R.h>

// Medusa
#include <medusa/generic/Constants.h>
#include <medusa/generic/Faddeeva.h>


namespace medusa {
    namespace functions {

        // Convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t ) with the Gaussian
        // [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Convolve_exp_sinhcosh(double time, double a, double b, double mu, double sigma, bool tag)
        {
            double x = (time - mu)/(sigma*M_Sqrt2);

            double z1 = (a - b)*sigma/M_Sqrt2;
            double z2 = (a + b)*sigma/M_Sqrt2;

            double faddeeva_z1 = ::exp( z1*z1 - 2*z1*x ) * faddeeva::erfc(z1-x);
            double faddeeva_z2 = ::exp( z2*z2 - 2*z2*x ) * faddeeva::erfc(z2-x);

            if(tag) return 0.25*(faddeeva_z1 + faddeeva_z2);

            else return 0.25*(faddeeva_z1 - faddeeva_z2);
        }


        // Convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t ) with the Gaussian
        // [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Convolve_exp_sincos(double time, double a, double b, double mu, double sigma, bool tag)
        {
            double x = (time - mu)/(sigma*M_Sqrt2);

            hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
            hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

            hydra::complex<double> faddeeva_z1 = hydra::exp( z1*z1 - 2*z1*x ) * faddeeva::erfc(z1-x);
            hydra::complex<double> faddeeva_z2 = hydra::exp( z2*z2 - 2*z2*x ) * faddeeva::erfc(z2-x);

            hydra::complex<double> faddeeva_tot = 0.;
            double result = 0.;

            if(tag)
            {
                faddeeva_tot = faddeeva_z1 + faddeeva_z2;

                result = faddeeva_tot.real();
            }
            else
            {
                faddeeva_tot = faddeeva_z1 - faddeeva_z2;

                result = faddeeva_tot.imag();
            }

            return 0.25*result;
        }


        // Integrate in the time exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
        // [tag = true -> cosh | tag = false -> sinh]
        __hydra_dual__
        inline double Integrate_exp_sinhcosh(double a, double b, double LowerLimit, double UpperLimit, bool tag)
        {
            if(tag) return ( -::exp(-a*UpperLimit) * (a*::cosh(b*UpperLimit) + b*::sinh(b*UpperLimit)) +
                                ::exp(-a*LowerLimit) * (a*::cosh(b*LowerLimit) + b*::sinh(b*LowerLimit)) ) /
                                                                                                ( (a - b) * (a + b));

            else return ( -::exp(-a*UpperLimit) * (b*::cosh(b*UpperLimit) + a*::sinh(b*UpperLimit)) +
                                ::exp(-a*LowerLimit) * (b*::cosh(b*LowerLimit) + a*::sinh(b*LowerLimit)) ) /
                                                                                                ( (a - b) * (a + b));
        }


        // Integrate in the time exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
        // [tag = true -> cos | tag = false -> sin]
        __hydra_dual__
        inline double Integrate_exp_sincos(double a, double b, double LowerLimit, double UpperLimit, bool tag)
        {
            if(tag) return ( ::exp(-a*UpperLimit) * (-a*::cos(b*UpperLimit) + b*::sin(b*UpperLimit)) +
                                ::exp(-a*LowerLimit) * (a*::cos(b*LowerLimit) - b*::sin(b*LowerLimit)) ) / (a*a + b*b);

            else return ( -::exp(-a*UpperLimit) * (b*::cos(b*UpperLimit) + a*::sin(b*UpperLimit)) +
                                ::exp(-a*LowerLimit) * (b*::cos(b*LowerLimit) + a*::sin(b*LowerLimit)) ) / (a*a + b*b);
        }


        // Integrate in the time the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
        // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Integrate_convolved_exp_sinhcosh(double a, double b, double mu, double sigma,
                                                                        double LowerLimit, double UpperLimit, bool tag)
        {
            double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            double z1 = (a - b)*sigma/M_Sqrt2;
            double z2 = (a + b)*sigma/M_Sqrt2;

            double cumulative_z1 = ( faddeeva::erf(x2) - ::exp( z1*z1 - 2*z1*x2 ) * faddeeva::erfc(z1-x2) -
                                            ( faddeeva::erf(x1) - ::exp( z1*z1 - 2*z1*x1 ) * faddeeva::erfc(z1-x1) ) ) / z1;

            double cumulative_z2 = ( faddeeva::erf(x2) - ::exp( z2*z2 - 2*z2*x2 ) * faddeeva::erfc(z2-x2) -
                                            ( faddeeva::erf(x1) - ::exp( z2*z2 - 2*z2*x1 ) * faddeeva::erfc(z2-x1) ) ) / z2;

            if(tag) return M_1_Sqrt32 * sigma * (cumulative_z1 + cumulative_z2);

            else return M_1_Sqrt32 * sigma * (cumulative_z1 - cumulative_z2);
        }


        // Integrate in the time the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
        // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Integrate_convolved_exp_sincos(double a, double b, double mu, double sigma,
                                                                    double LowerLimit, double UpperLimit, bool tag)
        {
            double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
            hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

            hydra::complex<double> cumulative_z1 = faddeeva::erf(x2) - hydra::exp( z1*z1 - 2*z1*x2 ) * faddeeva::erfc(z1-x2) -
                                                    ( faddeeva::erf(x1) - hydra::exp( z1*z1 - 2*z1*x1 ) * faddeeva::erfc(z1-x1) );

            hydra::complex<double> cumulative_z2 = faddeeva::erf(x2) - hydra::exp( z2*z2 - 2*z2*x2 ) * faddeeva::erfc(z2-x2) -
                                                    ( faddeeva::erf(x1) - hydra::exp( z2*z2 - 2*z2*x1 ) * faddeeva::erfc(z2-x1) );

            hydra::complex<double> cumulative = 0.;
            double result = 0.;

            if(tag)
            {
                cumulative = cumulative_z1/z1 + cumulative_z2/z2;

                result = cumulative.real();
            }
            else
            {
                cumulative = cumulative_z1/z1 - cumulative_z2/z2;

                result = cumulative.imag();
            }

            return M_1_Sqrt32 * sigma * result;
        }

    } // namespace functions
} // namespace medusa

#endif // MEDUSA_FUNCTIONS_H