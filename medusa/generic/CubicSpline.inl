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

        // Compute the powers (Reference: arXiv:1407.0748v1)
        double powS[4], powM[4];

        powS[0] = ::pow(M_1_Sqrt2*sigma, 0);
        powS[1] = ::pow(M_1_Sqrt2*sigma, 1);
        powS[2] = ::pow(M_1_Sqrt2*sigma, 2);
        powS[3] = ::pow(M_1_Sqrt2*sigma, 3);

        powM[0] = ::pow(mu, 0);
        powM[1] = ::pow(mu, 1);
        powM[2] = ::pow(mu, 2);
        powM[3] = ::pow(mu, 3);

        // Compute z1 and z2 (Reference: arXiv:1407.0748v1)
        double z1 = (a - b)*sigma/M_Sqrt2;
        double z2 = (a + b)*sigma/M_Sqrt2;

        // x1 and x2 (Reference: arXiv:1407.0748v1)
        double x1 = 0;
        double x2 = 0;

        // K_n(z) for z as double (Reference: arXiv:1407.0748v1)
        double Kz1[4], Kz2[4];

        Kz1[0] = K(z1, 0);
        Kz1[1] = K(z1, 1);
        Kz1[2] = K(z1, 2);
        Kz1[3] = K(z1, 3);

        Kz2[0] = K(z2, 0);
        Kz2[1] = K(z2, 1);
        Kz2[2] = K(z2, 2);
        Kz2[3] = K(z2, 3);

        // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
        double Mz1[4] = {0};
        double Mz2[4] = {0};

        // Calculate the integration on the intervals between the knots
        double sum = 0.;
        if(iFrom==iTo)
        {
            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);
            
            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iFrom, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
        }
        else
        {
            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            x2 = (u[4+iFrom] - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);

            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iFrom, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);

            for(size_t i=iFrom+1; i<iTo; i++)
            {
                // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
                x1 = (u[3+i] - mu)/(sigma*M_Sqrt2);
                x2 = (u[4+i] - mu)/(sigma*M_Sqrt2);

                // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
                Mz1[0] = M(x1, x2, z1, 0);
                Mz1[1] = M(x1, x2, z1, 1);
                Mz1[2] = M(x1, x2, z1, 2);
                Mz1[3] = M(x1, x2, z1, 3);

                Mz2[0] = M(x1, x2, z2, 0);
                Mz2[1] = M(x1, x2, z2, 1);
                Mz2[2] = M(x1, x2, z2, 2);
                Mz2[3] = M(x1, x2, z2, 3);

                sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(i, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
            }

            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (u[3+iTo] - mu)/(sigma*M_Sqrt2);
            x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);

            sum += Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(iTo, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
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

        // Compute the powers (Reference: arXiv:1407.0748v1)
        double powS[4], powM[4];

        powS[0] = ::pow(M_1_Sqrt2*sigma, 0);
        powS[1] = ::pow(M_1_Sqrt2*sigma, 1);
        powS[2] = ::pow(M_1_Sqrt2*sigma, 2);
        powS[3] = ::pow(M_1_Sqrt2*sigma, 3);

        powM[0] = ::pow(mu, 0);
        powM[1] = ::pow(mu, 1);
        powM[2] = ::pow(mu, 2);
        powM[3] = ::pow(mu, 3);

        // Compute z1 and z2 (Reference: arXiv:1407.0748v1)
        hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
        hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

        // x1 and x2 (Reference: arXiv:1407.0748v1)
        double x1 = 0;
        double x2 = 0;

        // K_n(z) for z as double (Reference: arXiv:1407.0748v1)
        hydra::complex<double> Kz1[4], Kz2[4];

        Kz1[0] = K(z1, 0);
        Kz1[1] = K(z1, 1);
        Kz1[2] = K(z1, 2);
        Kz1[3] = K(z1, 3);

        Kz2[0] = K(z2, 0);
        Kz2[1] = K(z2, 1);
        Kz2[2] = K(z2, 2);
        Kz2[3] = K(z2, 3);

        // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
        hydra::complex<double> Mz1[4] = {0};
        hydra::complex<double> Mz2[4] = {0};

        // Calculate the integration on the intervals between the knots
        double sum = 0.;
        if(iFrom==iTo)
        {
            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);
            
            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iFrom, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
        }
        else
        {
            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            x2 = (u[4+iFrom] - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);

            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iFrom, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);

            for(size_t i=iFrom+1; i<iTo; i++)
            {
                // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
                x1 = (u[3+i] - mu)/(sigma*M_Sqrt2);
                x2 = (u[4+i] - mu)/(sigma*M_Sqrt2);

                // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
                Mz1[0] = M(x1, x2, z1, 0);
                Mz1[1] = M(x1, x2, z1, 1);
                Mz1[2] = M(x1, x2, z1, 2);
                Mz1[3] = M(x1, x2, z1, 3);

                Mz2[0] = M(x1, x2, z2, 0);
                Mz2[1] = M(x1, x2, z2, 1);
                Mz2[2] = M(x1, x2, z2, 2);
                Mz2[3] = M(x1, x2, z2, 3);

                sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(i, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
            }

            // Compute x1 and x2 (Reference: arXiv:1407.0748v1)
            x1 = (u[3+iTo] - mu)/(sigma*M_Sqrt2);
            x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            // M_n(x1, x2; z) for z as double (Reference: arXiv:1407.0748v1)
            Mz1[0] = M(x1, x2, z1, 0);
            Mz1[1] = M(x1, x2, z1, 1);
            Mz1[2] = M(x1, x2, z1, 2);
            Mz1[3] = M(x1, x2, z1, 3);

            Mz2[0] = M(x1, x2, z2, 0);
            Mz2[1] = M(x1, x2, z2, 1);
            Mz2[2] = M(x1, x2, z2, 2);
            Mz2[3] = M(x1, x2, z2, 3);

            sum += Integrate_3_order_polynomial_times_convolved_exp_sincos(iTo, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
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
    Integrate_3_order_polynomial_times_convolved_exp_sinhcosh(size_t bin, double Kz1[4], double Kz2[4],
                                            double Mz1[4], double Mz2[4], double mu, double sigma, double powS[4], double powM[4], bool tag) const
    {
        double sum = 0.;
        for(size_t i=0; i<4; i++)
        {
            sum += AS[i][bin]*Integrate_t_to_k_times_convolved_exp_sinhcosh(i, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
        }
        return sum;
    }


    // Integrate (in x) the 3rd order polynomial times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_3_order_polynomial_times_convolved_exp_sincos(size_t bin, hydra::complex<double> Kz1[4], hydra::complex<double> Kz2[4],
                            hydra::complex<double> Mz1[4], hydra::complex<double> Mz2[4], double mu, double sigma, double powS[4], double powM[4], bool tag) const
    {
        double sum = 0.;
        for(size_t i=0; i<4; i++)
        {
            sum += AS[i][bin]*Integrate_t_to_k_times_convolved_exp_sincos(i, Kz1, Kz2, Mz1, Mz2, mu, sigma, powS, powM, tag);
        }
        return sum;
    }


    // Integrate (in x) t^k times the convolution of exp( -a*t )*cosh( b*t ) or exp( -a*t )*sinh( b*t )
    // with the Gaussian [tag = true -> cosh | tag = false -> sinh] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_t_to_k_times_convolved_exp_sinhcosh(size_t k, double Kz1[4], double Kz2[4],
                                        double Mz1[4], double Mz2[4], double mu, double sigma, double powS[4], double powM[4], bool tag) const
    {
        double sum1 = 0.;
        double sum2;

        if(mu == 0)
        {
            if(tag)
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 += ( Kz1[j]*Mz1[k-j] + Kz2[j]*Mz2[k-j] ) / (factorial[j]*factorial[k-j]);
                }
            }
            else
            {
                for(size_t j=0; j<k+1; j++)
                {
                    sum1 += ( Kz1[j]*Mz1[k-j] - Kz2[j]*Mz2[k-j] ) / (factorial[j]*factorial[k-j]);
                }
            }
            return M_1_Sqrt8*sigma*factorial[k]*powS[k]*sum1;
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
                        sum2 += ( Kz1[j]*Mz1[n-j] + Kz2[j]*Mz2[n-j] ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += powS[n]*powM[k-n]/factorial[k-n] * sum2;
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( Kz1[j]*Mz1[n-j] - Kz2[j]*Mz2[n-j] ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += powS[n]*powM[k-n]/factorial[k-n] * sum2;
                }
            }
            return M_1_Sqrt8*sigma*factorial[k]*sum1;
        }
    }


    // Integrate (in x) t^k times the convolution of exp( -a*t )*cos( b*t ) or exp( -a*t )*sin( b*t )
    // with the Gaussian [tag = true -> cos | tag = false -> sin] (Reference: arXiv:1407.0748v1)
    template<size_t nKnots>
    inline double CubicSpline<nKnots>::
    Integrate_t_to_k_times_convolved_exp_sincos(size_t k, hydra::complex<double> Kz1[4], hydra::complex<double> Kz2[4],
                                hydra::complex<double> Mz1[4], hydra::complex<double> Mz2[4], double mu, double sigma, double powS[4], double powM[4], bool tag) const
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
                    sum2 += ( Kz1[j]*Mz1[k-j] + Kz2[j]*Mz2[k-j] ) / (factorial[j]*factorial[k-j]);
                }
                sum1 = sum2.real();
            }
            else
            {
                sum2 = 0.;
                for(size_t j=0; j<k+1; j++)
                {
                    sum2 += ( Kz1[j]*Mz1[k-j] - Kz2[j]*Mz2[k-j] ) / (factorial[j]*factorial[k-j]);
                }
                sum1 = sum2.imag();
            }
            return M_1_Sqrt8*sigma*factorial[k]*powS[k]*sum1;
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
                        sum2 += ( Kz1[j]*Mz1[n-j] + Kz2[j]*Mz2[n-j] ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += powS[n]*powM[k-n]/factorial[k-n] * sum2.real();
                }
            }
            else
            {
                for(size_t n=0; n<k+1; n++)
                {
                    sum2 = 0.;
                    for(size_t j=0; j<n+1; j++)
                    {
                        sum2 += ( Kz1[j]*Mz1[n-j] - Kz2[j]*Mz2[n-j] ) / (factorial[j]*factorial[n-j]);
                    }
                    sum1 += powS[n]*powM[k-n]/factorial[k-n] * sum2.imag();
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