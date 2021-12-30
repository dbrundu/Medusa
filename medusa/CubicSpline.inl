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
 *  Created on: 30/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  CubicSpline.inl
 *----------------------------------------*/

#ifndef MEDUSA_CUBIC_SPLINE_INL
#define MEDUSA_CUBIC_SPLINE_INL


namespace medusa {

    //----------------------------------
    //      Methods to integrate
    //----------------------------------
    
    // Integrate in t the cubic spline times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::
                    Integrate_cspline_times_convolved_exp_sinhcosh(double a, double b, double mu, double sigma,
                                                                            double LowerLimit, double UpperLimit, bool tag) const
    {
        // Integrate outside the defined region of the spline
        if(LowerLimit<u[3])
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the n-value for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "n=%f", LowerLimit);

            return NaN;
        }

        // Integrate in the negative linear region of the spline
        if(NegativePart && LowerLimit >= xNegative)
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the n-value for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "n=%f", LowerLimit);

            return NaN;
        }

        // Exclude the negative linear part of the spline
        if(NegativePart) UpperLimit = std::min(UpperLimit,xNegative);

        size_t iFrom = findKnot(LowerLimit);
        size_t iTo = findKnot(UpperLimit);

        double sum = 0.;
        for (size_t i=0; i<3; i++)
        {
            sum = sum + Integrate_Ak_t_to_k_times_convolved_exp_sinhcosh(i, iFrom, iTo, a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        }

        return sum;
    }


    // Integrate in t the cubic spline times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::
                    Integrate_cspline_times_convolved_exp_sincos(double a, double b, double mu, double sigma,
                                                                        double LowerLimit, double UpperLimit, bool tag) const
    {
        // Integrate outside the defined region of the spline
        if(LowerLimit<u[3])
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the n-value for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "n=%f", LowerLimit);

            return NaN;
        }

        // Integrate in the negative linear region of the spline
        if(NegativePart && LowerLimit >= xNegative)
        {
            // This macro controls if the first argument is NaN. If yes, it prints
            // a warning with the parameter value for whom we obtain a NaN.
            // In this case, the first argument is always NaN and the macro prints
            // the n-value for whom we obtain NaN.
            hydra::CHECK_VALUE(NaN, "n=%f", LowerLimit);

            return NaN;
        }

        // Exclude the negative linear part of the spline
        if(NegativePart) UpperLimit = std::min(UpperLimit,xNegative);

        size_t iFrom = findKnot(LowerLimit);
        size_t iTo = findKnot(UpperLimit);

        double sum = 0.;
        for (size_t i=0; i<3; i++)
        {
            sum = sum + Integrate_Ak_t_to_k_times_convolved_exp_sincos(i, iFrom, iTo, a, b, mu, sigma, LowerLimit, UpperLimit, tag);
        }

        return sum;
    }


    //-------------------------------------------------
    //        Methods to help the integratation
    //-------------------------------------------------

    // Integrate (in x) A_k*t^k times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::
                    Integrate_Ak_t_to_k_times_convolved_exp_sinhcosh(size_t k, int iFrom, int iTo,
                                                                        double a, double b, double mu, double sigma,
                                                                            double LowerLimit, double UpperLimit, bool tag) const
    {
        double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
        double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

        double z1 = (a - b)*sigma/M_Sqrt2;
        double z2 = (a + b)*sigma/M_Sqrt2;

        double sum1 = 0.0;
        double sum2;

        if(mu == 0)
        {
            if(tag)
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 = sum1 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, k-j) +
                                            K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                }
            }
            else
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 = sum1 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, k-j) -
                                            K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                }
            }
            return M_1_Sqrt8*sigma*kfactorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
        }
        else
        {
            if(tag)
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, n-j) +
                                                K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                    }
                    sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2;
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, n-j) -
                                                K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                    }
                    sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2;
                }
            }
            return M_1_Sqrt8*sigma*kfactorial[k]*sum1;
        }
    }


    // Integrate (in x) A_k*t^k times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::
                    Integrate_Ak_t_to_k_times_convolved_exp_sincos(size_t k, int iFrom, int iTo,
                                                                        double a, double b, double mu, double sigma,
                                                                            double LowerLimit, double UpperLimit, bool tag) const
    {
        double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
        double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

        hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
        hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

        double sum1 = 0.0;
        hydra::complex<double> sum2;
            
        if(mu == 0)
        {
            if(tag)
            {
                sum2 = 0.0;
                for(size_t j=0; j<k+1; j++)
                {
                    sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, k-j) +
                                            K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                }
                sum1 = sum2.real();
            }
            else
            {
                sum2 = 0.0;
                for(size_t j=0; j<k+1; j++)
                {
                    sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, k-j) -
                                            K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                }
                sum1 = sum2.imag();
            }
            return M_1_Sqrt8*sigma*kfactorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
        }
        else
        {
            if(tag)
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, n-j) +
                                                K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                    }
                    sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2.real();
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(k, iFrom, iTo, x1, x2, z1, n-j) -
                                                K(z2, j)*M(k, iFrom, iTo, x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                    }
                    sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2.imag();
                }
            }
            return M_1_Sqrt8*sigma*kfactorial[k]*sum1;
        }
    }


    //------------------------------------------------------------
    //       Intermediate functions used in the integration
    //------------------------------------------------------------

    // K_n(z) for z as double (Reference: arXiv:1407.0748v1)
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::K(double z, size_t n) const
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
    template<bool Reduce, size_t nKnots>
    inline hydra::complex<double> CubicSpline<Reduce, nKnots>::K(hydra::complex<double> z, size_t n) const
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
    template<bool Reduce, size_t nKnots>
    inline double CubicSpline<Reduce, nKnots>::M(size_t k, int iFrom, int iTo, double x1, double x2, double z, size_t n) const
    {
        switch (n)
        {
            case 0:
            {
                return A[k][iTo]*( faddeeva::erf(x2) - ::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( faddeeva::erf(x1) - ::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1) );
            }
            case 1:
            {
                return 2*( A[k][iTo]*( -M_1_SqrtPi*::exp(-x2*x2) - x2*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( -M_1_SqrtPi*::exp(-x1*x1) - x1*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 2:
            {
                return 2*( A[k][iTo]*( -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1) *
                                                                    ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                                A[k][iFrom]*( -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1) *
                                                                    ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 3:
            {
                return 4*( A[k][iTo]*( -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) *
                                                                        ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) *
                                                                        ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 4:
            {
                return -4*( A[k][iTo]*( 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                            (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                            (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 5:
            {
                return -8*( A[k][iTo]*( (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                            x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                            x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 6:
            {
                return -8*( A[k][iTo]*( x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                                        (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) *
                                                                            ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                                        (-15 + 90*x1*x1 - 60*x1*x1*x1*x1 + 8*x1*x1*x1*x1*x1*x1) *
                                                                            ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
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
    template<bool Reduce, size_t nKnots>
    inline hydra::complex<double> CubicSpline<Reduce, nKnots>::M(size_t k, int iFrom, int iTo, double x1, double x2,
                                                                                        hydra::complex<double> z, size_t n) const
    {
        switch (n)
        {
            case 0:
            {
                return A[k][iTo]*( faddeeva::erf(x2) - hydra::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( faddeeva::erf(x1) - hydra::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1) );
            }
            case 1:
            {
                return 2*( A[k][iTo]*( -M_1_SqrtPi*::exp(-x2*x2) - x2*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( -M_1_SqrtPi*::exp(-x1*x1) - x1*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 2:
            {
                return 2*( A[k][iTo]*( -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1) *
                                                                hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                                A[k][iFrom]*( -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1) *
                                                                hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 3:
            {
                return 4*( A[k][iTo]*( -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) *
                                                            hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) *
                                                            hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 4:
            {
                return -4*( A[k][iTo]*( 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                            (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                            (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 5:
            {
                return -8*( A[k][iTo]*( (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                            x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                            x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
            }
            case 6:
            {
                return -8*( A[k][iTo]*( x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                                        (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) * 
                                                                            hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) ) -
                            A[k][iFrom]*( x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                                        (-15 + 90*x1*x1 - 60*x1*x1*x1*x1 + 8*x1*x1*x1*x1*x1*x1) *
                                                                            hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
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