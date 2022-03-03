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
 *  Created on: 30/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  This file contains the implementation of some methods
 *  of the CubicSpline class defined in CubicSpline.h file.
 *---------------------------------------------------------------------------*/

#ifndef MEDUSA_CUBIC_SPLINE_INL
#define MEDUSA_CUBIC_SPLINE_INL


namespace medusa {

    //----------------------------------
    //      Methods to integrate
    //----------------------------------
    
    // Integrate in t the cubic spline times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
                    Integrate_cspline_times_convolved_exp_sinhcosh(double a, double b, double mu, double sigma,
                                                                            double LowerLimit, double UpperLimit, bool tag) const
    {
        // Integrate outside the defined region of the spline
        if(LowerLimit<u[3])
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the values for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "Integrate outside the defined region of the spline: LowerLimit=%f < u[3]=%f", LowerLimit, u[3]);
            return NaN;
        }

        // Integrate in the negative linear region of the spline
        if(NegativePart && LowerLimit >= xNegative)
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the values for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "Integrate in the negative linear region of the spline: LowerLimit=%f >= xNegative=%f",
                                                                                                            LowerLimit, xNegative);
            return NaN;
        }

        // Exclude the negative linear part of the spline
        if(NegativePart) UpperLimit = std::min(UpperLimit,xNegative);

        // Find the intervals between the knots of the integration limits
        size_t iFrom = findKnot(LowerLimit);
        size_t iTo = findKnot(UpperLimit);

        // Calculate the integration on the intervals between the knots
        double sum = 0.;
        if(iFrom==iTo)
        {
            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iFrom, a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        }
        else
        {
            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iFrom, a, b, mu, sigma, LowerLimit, u[4+iFrom], tag);

            for(size_t i=iFrom+1; i<iTo; i++)
            {
                sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(i, a, b, mu, sigma, u[3+i], u[4+i], tag);
            }

            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iTo, a, b, mu, sigma, u[3+iTo], UpperLimit, tag);
        }

        // This macro controls if sum is NaN. If yes, it prints
        // a warning with the parameter value for whom we obtain a NaN.
        // In this case, the first argument is always NaN and the macro prints
        // the values for whom we obtain NaN.
        hydra::CHECK_VALUE(sum, "a=%f, b=%f, mu=%f, sigma=%f, LowerLimit=%f, UpperLimit=%f, tag=%d",
                                                                        a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        return sum;
    }


    // Integrate in t the cubic spline times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_cspline_times_convolved_exp_sincos(double a, double b, double mu, double sigma,
                                                                        double LowerLimit, double UpperLimit, bool tag) const
    {
        // Integrate outside the defined region of the spline
        if(LowerLimit<u[3])
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the values for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "Integrate outside the defined region of the spline: LowerLimit=%f < u[3]=%f", LowerLimit, u[3]);
            return NaN;
        }

        // Integrate in the negative linear region of the spline
        if(NegativePart && LowerLimit >= xNegative)
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the values for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "Integrate in the negative linear region of the spline: LowerLimit=%f >= xNegative=%f",
                                                                                                            LowerLimit, xNegative);
            return NaN;
        }

        // Exclude the negative linear part of the spline
        if(NegativePart) UpperLimit = std::min(UpperLimit,xNegative);

        // Find the intervals between the knots of the integration limits
        size_t iFrom = findKnot(LowerLimit);
        size_t iTo = findKnot(UpperLimit);

        // Calculate the integration on the intervals between the knots
        double sum = 0.;
        if(iFrom==iTo)
        {
            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iFrom, a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        }
        else
        {
            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iFrom, a, b, mu, sigma, LowerLimit, u[4+iFrom], tag);

            for(size_t i=iFrom+1; i<iTo; i++)
            {
                sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(i, a, b, mu, sigma, u[3+i], u[4+i], tag);
            }

            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iTo, a, b, mu, sigma, u[3+iTo], UpperLimit, tag);
        }

        // This macro controls if sum is NaN. If yes, it prints
        // a warning with the parameter value for whom we obtain a NaN.
        // In this case, the first argument is always NaN and the macro prints
        // the values for whom we obtain NaN.
        hydra::CHECK_VALUE(sum, "a=%f, b=%f, mu=%f, sigma=%f, LowerLimit=%f, UpperLimit=%f, tag=%d",
                                                                        a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        return sum;
    }


    //-------------------------------------------------
    //        Methods to help the integratation
    //-------------------------------------------------

    // Integrate (in x) the 3rd order polynomial times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(size_t bin, double a, double b, double mu, double sigma,
                                                                                double LowerLimit, double UpperLimit, bool tag) const
    {
        double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
        double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

        double z1 = (a - b)*sigma/M_Sqrt2;
        double z2 = (a + b)*sigma/M_Sqrt2;

        double sum = 0.;
        for(size_t i=0; i<4; i++)
        {
            sum += AS[i][bin]*Integrate_t_to_k_times_convolved_exp_sinhcosh(i, mu, sigma, z1, z2, x1, x2, tag);
        }
        return sum;
    }


    // Integrate (in x) the 3rd order polynomial times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_3_order_polynomial_times_convolved_exp_sincos(size_t bin, double a, double b, double mu, double sigma,
                                                                                double LowerLimit, double UpperLimit, bool tag) const
    {
        double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
        double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

        hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
        hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

        double sum = 0.;
        for(size_t i=0; i<4; i++)
        {
            sum += AS[i][bin]*Integrate_t_to_k_times_convolved_exp_sincos(i, mu, sigma, z1, z2, x1, x2, tag);
        }
        return sum;
    }


    // Integrate (in x) t^k times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_t_to_k_times_convolved_exp_sinhcosh(size_t k, double mu, double sigma, double z1, double z2, double x1, double x2, bool tag) const
    {
        double sum1 = 0.;
        double sum2;

        if(mu == 0)
        {
            if(tag)
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 += ( K(z1, j)*M(x1, x2, z1, k-j) + K(z2, j)*M(x1, x2, z2, k-j) ) / (factorial[j]*factorial[k-j]);
                }
            }
            else
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 += ( K(z1, j)*M(x1, x2, z1, k-j) - K(z2, j)*M(x1, x2, z2, k-j) ) / (factorial[j]*factorial[k-j]);
                }
            }
            return M_1_Sqrt8*sigma*factorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
        }
        else
        {
            if(tag)
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( K(z1, j)*M(x1, x2, z1, n-j) + K(z2, j)*M(x1, x2, z2, n-j) ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/factorial[k-n] * sum2;
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( K(z1, j)*M(x1, x2, z1, n-j) - K(z2, j)*M(x1, x2, z2, n-j) ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/factorial[k-n] * sum2;
                }
            }
            return M_1_Sqrt8*sigma*factorial[k]*sum1;
        }
    }


    // Integrate (in x) t^k times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_t_to_k_times_convolved_exp_sincos(size_t k, double mu, double sigma, hydra::complex<double> z1,
                                                                            hydra::complex<double> z2, double x1, double x2, bool tag) const
    {
        double sum1 = 0.;
        hydra::complex<double> sum2;
            
        if(mu == 0)
        {
            if(tag)
            {
                sum2 = 0.;
                for(size_t j=0; j<k+1; j++)
                {
                    sum2 += ( K(z1, j)*M(x1, x2, z1, k-j) + K(z2, j)*M(x1, x2, z2, k-j) ) / (factorial[j]*factorial[k-j]);
                }
                sum1 = sum2.real();
            }
            else
            {
                sum2 = 0.;
                for(size_t j=0; j<k+1; j++)
                {
                    sum2 += ( K(z1, j)*M(x1, x2, z1, k-j) - K(z2, j)*M(x1, x2, z2, k-j) ) / (factorial[j]*factorial[k-j]);
                }
                sum1 = sum2.imag();
            }
            return M_1_Sqrt8*sigma*factorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
        }
        else
        {
            if(tag)
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( K(z1, j)*M(x1, x2, z1, n-j) + K(z2, j)*M(x1, x2, z2, n-j) ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/factorial[k-n] * sum2.real();
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( K(z1, j)*M(x1, x2, z1, n-j) - K(z2, j)*M(x1, x2, z2, n-j) ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/factorial[k-n] * sum2.imag();
                }
            }
            return M_1_Sqrt8*sigma*factorial[k]*sum1;
        }
    }


    //------------------------------------------------------------
    //       Intermediate functions used in the integration
    //------------------------------------------------------------

    // K_n(z) for z as double (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::K(double z, size_t n) const
    {
        switch (n)
        {
            case 0:
            {
                return 1 /(2*z);
            }
            case 1:
            {
                return 1/(2*z*z);
            }
            case 2:
            {
                return 1/z*( 1 + 1/(z*z) );
            }
            case 3:
            {
                return 3/(z*z) * ( 1 + 1/(z*z) );
            }
            case 4:
            {
                return 6/(z*z*z*z*z) * (2 + 2*z*z + z*z*z*z);
            }
            case 5:
            {
                return 30/(z*z*z*z*z*z) * (2 + 2*z*z + z*z*z*z);
            }
            case 6:
            {
                return 60/(z*z*z*z*z*z*z) * (6 + 6*z*z + 3*z*z*z*z + z*z*z*z*z*z);
            }
            default:
            {
                // This macro controls if the first argument is NaN. If yes, it prints
                // a warning with the parameter value for whom we obtain a NaN.
                // In this case, the first argument is always NaN and the macro prints
                // the n-value for whom we obtain NaN.
                hydra::CHECK_VALUE(NaN, "n=%d", n);

                return NaN;
            }
        }
    }


    // K_n(z) for z as complex (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline hydra::complex<double> CubicSpline<nKnots>::K(hydra::complex<double> z, size_t n) const
    {
        switch (n)
        {
            case 0:
            {
                return 1 /(2*z);
            }
            case 1:
            {
                return 1/(2*z*z);
            }
            case 2:
            {
                return 1/z*( 1 + 1/(z*z) );
            }
            case 3:
            {
                return 3/(z*z) * ( 1 + 1/(z*z) );
            }
            case 4:
            {
                return 6/(z*z*z*z*z) * (2 + 2*z*z + z*z*z*z);
            }
            case 5:
            {
                return 30/(z*z*z*z*z*z) * (2 + 2*z*z + z*z*z*z);
            }
            case 6:
            {
                return 60/(z*z*z*z*z*z*z) * (6 + 6*z*z + 3*z*z*z*z + z*z*z*z*z*z);
            }
            default:
            {
                // This macro controls if the first argument is NaN. If yes, it prints
                // a warning with the parameter value for whom we obtain a NaN.
                // In this case, the first argument is always NaN and the macro prints
                // the n-value for whom we obtain NaN.
                hydra::CHECK_VALUE(NaN, "n=%d", n);

                return NaN;
            }
        }
    }


    // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::M(double x1, double x2, double z, size_t n) const
    {
        double From = 0.;
        double To = 0.;
        switch (n)
        {
            case 0:
            {
                To = faddeeva::erf(x2) - ::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2);
                From = faddeeva::erf(x1) - ::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1);
                return To - From;
            }
            case 1:
            {
                To = -M_1_SqrtPi*::exp(-x2*x2) - x2*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -M_1_SqrtPi*::exp(-x1*x1) - x1*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 2*( To - From );
            }
            case 2:
            {
                To = -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1) *::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1) *::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 2*( To - From );
            }
            case 3:
            {
                To = -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) * ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) * ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 4*( To - From );
            }
            case 4:
            {
                To = 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                            (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                            (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -4*( To - From );
            }
            case 5:
            {
                To = (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                                x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                                x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -8*( To - From );
            }
            case 6:
            {
                To = x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                                        (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) *
                                                                                ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                                        (-15 + 90*x1*x1 - 60*x1*x1*x1*x1 + 8*x1*x1*x1*x1*x1*x1) *
                                                                                ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -8*( To - From );
            }
            default:
            {
                // This macro controls if the first argument is NaN. If yes, it prints
                // a warning with the parameter value for whom we obtain a NaN.
                // In this case, the first argument is always NaN and the macro prints
                // the n-value for whom we obtain NaN.
                hydra::CHECK_VALUE(NaN, "n=%d", n);

                return NaN;
            }
        }
    }


    // M_n(x1, x2; z) for z as complex (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline hydra::complex<double> CubicSpline<nKnots>::M(double x1, double x2, hydra::complex<double> z, size_t n) const
    {
        hydra::complex<double> From = 0.;
        hydra::complex<double> To = 0.;

        switch (n)
        {
            case 0:
            {
                To = faddeeva::erf(x2) - hydra::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2);
                From = faddeeva::erf(x1) - hydra::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1);
                return To - From;
            }
            case 1:
            {
                To = -M_1_SqrtPi*::exp(-x2*x2) - x2 * hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -M_1_SqrtPi*::exp(-x1*x1) - x1 * hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 2*( To - From );
            }
            case 2:
            {
                To = -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1) * hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1) * hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 2*( To - From );
            }
            case 3:
            {
                To = -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) * hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) * hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return 4*( To - From );
            }
            case 4:
            {
                To = 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                            (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                            (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -4*( To - From );
            }
            case 5:
            {
                To = (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                            x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From = (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                            x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -8*( To - From );
            }
            case 6:
            {
                To =  x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                                        (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) * 
                                                                            hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2);
                From =  x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                                        (-15 + 90*x1*x1 - 60*x1*x1*x1*x1 + 8*x1*x1*x1*x1*x1*x1) *
                                                                            hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1);
                return -8*( To - From );
            }
            default:
            {
                // This macro controls if the first argument is NaN. If yes, it prints
                // a warning with the parameter value for whom we obtain a NaN.
                // In this case, the first argument is always NaN and the macro prints
                // the n-value for whom we obtain NaN.
                hydra::CHECK_VALUE(NaN, "n=%d", n);

                return NaN;
            }
        }
    }

} // medusa

#endif // MEDUSA_CUBIC_SPLINE_INL