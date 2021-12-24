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
 *  Created on: 09/12/2021
 *
 *  Author: Alessandro Maria Ricci
 *
 *  CubicSpline.h
 *----------------------------------------*/

#ifndef MEDUSA_CUBIC_SPLINE_H
#define MEDUSA_CUBIC_SPLINE_H

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/Complex.h>
#include <hydra/Function.h>

// Medusa
#include <medusa/Constants.h>

// ROOT
#ifdef _ROOT_AVAILABLE_

#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TString.h>
#include <TStyle.h>

#endif //_ROOT_AVAILABLE_


namespace medusa {

    template<typename ArgTypeTime,
             size_t nKnots,
             typename Signature=double(ArgTypeTime)>
    class CubicSpline: public hydra::BaseFunctor< CubicSpline <ArgTypeTime, nKnots>, Signature, nKnots+2 >
    {
        using CSplineBaseFunctor = hydra::BaseFunctor<CubicSpline <ArgTypeTime, nKnots>, Signature, nKnots+2 >;

        using hydra::BaseFunctor< CubicSpline <ArgTypeTime, nKnots>, Signature, nKnots+2 >::_par;


        public:

        //-------------------------------------
        //           Constructors
        //-------------------------------------

        CubicSpline() = delete;

        CubicSpline(const double (&knots)[nKnots], const std::array<hydra::Parameter, nKnots+2> &coeffs):
        CSplineBaseFunctor(coeffs)
        {
            // set knot vector
            u[0] = knots[0];
            u[1] = knots[0];
            u[2] = knots[0];
            for(size_t i=0; i<nKnots; i++)
            {
                u[3+i] = knots[i];
            }
            u[nKnots+3] = knots[nKnots-1];
            u[nKnots+4] = knots[nKnots-1];
            u[nKnots+5] = knots[nKnots-1];

            // calculate prefactors
            double P[nKnots], Q[nKnots], R[nKnots], S[nKnots];
            for(size_t i=0; i<nKnots; i++)
            {
                P[i] = (u[i+4]-u[i+1])*(u[i+4]-u[i+2])*(u[i+4]-u[i+3]);
                Q[i] = (u[i+5]-u[i+2])*(u[i+4]-u[i+2])*(u[i+4]-u[i+3]);
                R[i] = (u[i+5]-u[i+3])*(u[i+5]-u[i+2])*(u[i+4]-u[i+3]);
                S[i] = (u[i+6]-u[i+3])*(u[i+5]-u[i+3])*(u[i+4]-u[i+3]);
            }

            // calculate coefficients of normal polynomial
            for(size_t i=0; i<nKnots-1; i++)
            {
                // Constant
                a0[0][i] = u[i+4]*u[i+4]*u[i+4]/P[i];

                a0[1][i] = -u[i+1]*u[i+4]*u[i+4]/P[i] - u[i+2]*u[i+4]*u[i+5]/Q[i] - u[i+3]*u[i+5]*u[i+5]/R[i];

                a0[2][i] = u[i+2]*u[i+2]*u[i+4]/Q[i] + u[i+2]*u[i+3]*u[i+5]/R[i] + u[i+3]*u[i+3]*u[i+6]/S[i];

                a0[3][i] = -u[i+3]*u[i+3]*u[i+3]/S[i];
    
                // Linear
                a1[0][i] = -3*u[i+4]*u[i+4]/P[i];

                a1[1][i] = (2*u[i+1]*u[i+4]+u[i+4]*u[i+4])/P[i]
                        + (u[i+2]*u[i+4]+u[i+2]*u[i+5]+u[i+4]*u[i+5])/Q[i]
                        + (2*u[i+3]*u[i+5]+u[i+5]*u[i+5])/R[i];

                a1[2][i] = -(2*u[i+2]*u[i+4]+u[i+2]*u[i+2])/Q[i]
                        -(u[i+2]*u[i+3]+u[i+2]*u[i+5]+u[i+3]*u[i+5])/R[i]
                        -(2*u[i+3]*u[i+6]+u[i+3]*u[i+3])/S[i];

                a1[3][i] = 3*u[i+3]*u[i+3]/S[i];
    
                // Quadratic
                a2[0][i] = 3*u[i+4]/P[i];

                a2[1][i] = -(2*u[i+4]+u[i+1])/P[i]
                        -(u[i+2]+u[i+4]+u[i+5])/Q[i]
                        -(2*u[i+5]+u[i+3])/R[i];

                a2[2][i] = (2*u[i+2]+u[i+4])/Q[i]
                        +(u[i+2]+u[i+5]+u[i+3])/R[i]
                        +(2*u[i+3]+u[i+6])/S[i];

                a2[3][i] = -3*u[i+3]/S[i];
    
                // Cubic
                a3[0][i] = -1./P[i];

                a3[1][i] = 1./P[i] + 1./Q[i] + 1./R[i];

                a3[2][i] = -1./Q[i] - 1./R[i] - 1./S[i];

                a3[3][i] = 1./S[i];
            }

            // calculate polynomial coefficients for the current set of spline coefficients
            updatePolynomial();

            // calculate the factorials
            kfactorial[0] = 1.0;  // k!
            kfactorial[1] = 1.0;
            kfactorial[2] = 2.0;
            kfactorial[3] = 6.0;

            jfactorial[0] = 1.0;  // j!
            jfactorial[1] = 1.0;
            jfactorial[2] = 2.0;
            jfactorial[3] = 6.0;

            knfactorial[0] = 1.0;  // (k-n)!
            knfactorial[1] = 1.0;
            knfactorial[2] = 2.0;
            knfactorial[3] = 6.0;

            njfactorial[0] = 1.0;  // (n-j)!
            njfactorial[1] = 1.0;
            njfactorial[2] = 2.0;
            njfactorial[3] = 6.0;
        }


        // ctor with other CubicSpline instance (copy ctor)
        __hydra_dual__
        CubicSpline(CubicSpline<ArgTypeTime, nKnots> const& other, int tag):
        CSplineBaseFunctor(other)
        {
            u[0] = other.GetKnots()[0];
            u[1] = other.GetKnots()[1];
            u[2] = other.GetKnots()[2];
            for(size_t i=0; i<nKnots; i++)
            {
                u[3+i] = other.GetKnots()[3+i];

                A0[i] = other.GetOverCoeffO0()[i];
                A1[i] = other.GetOverCoeffO1()[i];
                A2[i] = other.GetOverCoeffO2()[i];
                A3[i] = other.GetOverCoeffO3()[i];
            }
            u[nKnots+3] = other.GetKnots()[nKnots+3];
            u[nKnots+4] = other.GetKnots()[nKnots+4];
            u[nKnots+5] = other.GetKnots()[nKnots+5];

            if(tag < 0)
            {
                for (size_t i=0; i<4; i++)
                {
                    for(size_t j=0; j<nKnots-1; j++)
                    {
                        a0[i][j] = other.GetCoeffO0(i,j);
                        a1[i][j] = other.GetCoeffO1(i,j);
                        a2[i][j] = other.GetCoeffO2(i,j);
                        a3[i][j] = other.GetCoeffO3(i,j);
                    }
                }
            }

            NegativePart = other.GetNegativePart();
            xNegative = other.GetxNegative();

            for(size_t i=0; i<4; i++)
            {
                kfactorial[i] = other.GetKFactorial()[i];
                jfactorial[i] = other.GetJFactorial()[i];

                knfactorial[i] = other.GetKNFactorial()[i];
                njfactorial[i] = other.GetNJFactorial()[i];
            }
        }


        //-------------------------------------
        //       Operator= overloading
        //-------------------------------------

        __hydra_dual__
        CubicSpline<ArgTypeTime, nKnots>& operator=(CubicSpline<ArgTypeTime, nKnots> const& other)
        {
            if(this == &other) return *this;
            CSplineBaseFunctor::operator=(other);

            u[0] = other.GetKnots()[0];
            u[1] = other.GetKnots()[1];
            u[2] = other.GetKnots()[2];
            for(size_t i=0; i<nKnots; i++)
            {
                u[3+i] = other.GetKnots()[3+i];

                A0[i] = other.GetOverCoeffO0()[i];
                A1[i] = other.GetOverCoeffO1()[i];
                A2[i] = other.GetOverCoeffO2()[i];
                A3[i] = other.GetOverCoeffO3()[i];
            }
            u[nKnots+3] = other.GetKnots()[nKnots+3];
            u[nKnots+4] = other.GetKnots()[nKnots+4];
            u[nKnots+5] = other.GetKnots()[nKnots+5];

            for (size_t i=0; i<4; i++)
            {
                for(size_t j=0; j<nKnots-1; j++)
                {
                    a0[i][j] = other.GetCoeffO0(i,j);
                    a1[i][j] = other.GetCoeffO1(i,j);
                    a2[i][j] = other.GetCoeffO2(i,j);
                    a3[i][j] = other.GetCoeffO3(i,j);
                }
            }

            NegativePart = other.GetNegativePart();
            xNegative = other.GetxNegative();

            for(size_t i=0; i<4; i++)
            {
                kfactorial[i] = other.GetKFactorial()[i];
                jfactorial[i] = other.GetJFactorial()[i];

                knfactorial[i] = other.GetKNFactorial()[i];
                njfactorial[i] = other.GetNJFactorial()[i];
            }

            return *this;
        }


        //-------------------------------------
        //           Service functions
        //-------------------------------------

        // Find the respective knot number
        __hydra_dual__
        int findKnot(double x)
        {
            double j = 0;
            for(size_t i=0; i<nKnots; i++)
            {
                if(x>=u[3+i]) j=i;
            }
            return j;
        }


        // Evaluate the spline
        __hydra_dual__
        double eval(double x)
        {
            // Should return 0 here. But small positive number is better for the fitter
            // Change this for your purpose 
            if(x>xNegative && NegativePart) return 1e-3;

            int j = findKnot(x);
            return A0[j]+A1[j]*x+A2[j]*x*x+A3[j]*x*x*x;
        }


        #ifdef _ROOT_AVAILABLE_

        //create a histogram for plotting
        TH1D* CreateHistogramPlot(std::string name, std::string names, size_t nBins, double from, double to)
        {
            TH1D* h = new TH1D(name.c_str(),names.c_str(),nBins,from,to);

            for(size_t i=0; i<nBins; i++)
            {
                h->SetBinContent(i+1,eval(from+(i+0.5)*(to-from)/nBins));
            }
            return h;
        }

        #endif //_ROOT_AVAILABLE_


        //-------------------------------------
        //           Selectors
        //-------------------------------------
          
        // get the knot vector
        __hydra_dual__
        const double* GetKnots() const
        {
            return u;
        }

        // get k!
        __hydra_dual__
        const double* GetKFactorial() const
        {
            return kfactorial;
        }

        // get (k-n)!
        __hydra_dual__
        const double* GetKNFactorial() const
        {
            return knfactorial;
        }

        // get j!
        __hydra_dual__
        const double* GetJFactorial() const
        {
            return jfactorial;
        }

        // get (n-j)!
        __hydra_dual__
        const double* GetNJFactorial() const
        {
            return njfactorial;
        }

        // get the order 0 polynomial coefficients
        __hydra_dual__
        double GetCoeffO0(size_t i, size_t j) const
        {
            return a0[i][j];
        }

        // get the order 1 polynomial coefficients
        __hydra_dual__
        double GetCoeffO1(size_t i, size_t j) const
        {
            return a1[i][j];
        }

        // get the order 2 polynomial coefficients
        __hydra_dual__
        double GetCoeffO2(size_t i, size_t j) const
        {
            return a2[i][j];
        }

        // get the order 3 polynomial coefficients
        __hydra_dual__
        double GetCoeffO3(size_t i, size_t j) const
        {
            return a3[i][j];
        }
  
        // get the order 0 overall polynomial coefficients
        __hydra_dual__
        const double* GetOverCoeffO0() const
        {
            return A0;
        }

        // get the order 1 overall polynomial coefficients
        __hydra_dual__
        const double* GetOverCoeffO1() const
        {
            return A1;
        }
        
        // get the order 2 overall polynomial coefficients
        __hydra_dual__
        const double* GetOverCoeffO2() const
        {
            return A2;
        }
        
        // get the order 3 overall polynomial coefficients
        __hydra_dual__
        const double* GetOverCoeffO3() const
        {
            return A3;
        }

        // get the info wether this spline will get negative at some point
        __hydra_dual__
        bool GetNegativePart() const
        {
            return NegativePart;
        }
        
        __hydra_dual__
        double GetxNegative() const
        {
            return xNegative;
        }


        //-------------------------------------
        //        methods to integrate
        //-------------------------------------

        // Integrate in the time the convolution of t^k*exp( -a*t )*cosh( b*t ) or t^k*exp( -a*t )*sinh( b*t )
        // with the Gaussian [tag >= 0 -> cosh | tag < 0 -> sinh] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Integrate_t_to_k_times_convolved_exp_sinhcosh(size_t k, double a, double b, double mu, double sigma,
                                                                            ArgTypeTime LowerLimit, ArgTypeTime UpperLimit, int tag)
        {
            double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            double z1 = (a - b)*sigma/M_Sqrt2;
            double z2 = (a + b)*sigma/M_Sqrt2;

            double sum1 = 0.0;
            double sum2;

            if(mu == 0)
            {
                if(tag >= 0)
                {
                    for(size_t j=0; j<k+1; j++)
                    {
                        sum1 = sum1 + ( K(z1, j)*M(x1, x2, z1, k-j) + K(z2, j)*M(x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                    }
                }
                else
                {
                    for(size_t j=0; j<k+1; j++)
                    {
                        sum1 = sum1 + ( K(z1, j)*M(x1, x2, z1, k-j) - K(z2, j)*M(x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                    }
                }
                return M_1_Sqrt8*sigma*kfactorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
            }
            else
            {
                if(tag >= 0)
                {
                    for(size_t n=0; n<k+1; n++)
                    {
                        sum2 = 0.0;
                        for(size_t j=0; j<n+1; j++)
                        {
                            sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, n-j) + K(z2, j)*M(x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
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
                            sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, n-j) - K(z2, j)*M(x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                        }
                        sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2;
                    }
                }
                return M_1_Sqrt8*sigma*kfactorial[k]*sum1;
            }
        }


        // Integrate in the time the convolution of t^k*exp( -a*t )*cos( b*t ) or t^k*exp( -a*t )*sin( b*t )
        // with the Gaussian [tag >= 0 -> cos | tag < 0 -> sin] (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double Integrate_t_to_k_times_convolved_exp_sincos(size_t k, double a, double b, double mu, double sigma,
                                                                            ArgTypeTime LowerLimit, ArgTypeTime UpperLimit, int tag)
        {
            double x1 = (LowerLimit - mu)/(sigma*M_Sqrt2);
            double x2 = (UpperLimit - mu)/(sigma*M_Sqrt2);

            hydra::complex<double> z1( a*sigma/M_Sqrt2, -b*sigma/M_Sqrt2 );
            hydra::complex<double> z2( a*sigma/M_Sqrt2,  b*sigma/M_Sqrt2 );

            double sum1 = 0.0;
            hydra::complex<double> sum2;
            
            if(mu == 0)
            {
                if(tag >= 0)
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<k+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, k-j) + K(z2, j)*M(x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                    }
                    sum1 = sum2.real();
                }
                else
                {
                    sum2 = 0.0;
                    for(size_t j=0; j<k+1; j++)
                    {
                        sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, k-j) - K(z2, j)*M(x1, x2, z2, k-j) ) / (jfactorial[j]*knfactorial[k-j]);
                    }
                    sum1 = sum2.imag();
                }
                return M_1_Sqrt8*sigma*kfactorial[k]*::pow(M_1_Sqrt2*sigma, k)*sum1;
            }
            else
            {
                if(tag >= 0)
                {
                    for(size_t n=0; n<k+1; n++)
                    {
                        sum2 = 0.0;
                        for(size_t j=0; j<n+1; j++)
                        {
                            sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, n-j) + K(z2, j)*M(x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
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
                            sum2 = sum2 + ( K(z1, j)*M(x1, x2, z1, n-j) - K(z2, j)*M(x1, x2, z2, n-j) ) / (jfactorial[j]*njfactorial[n-j]);
                        }
                        sum1 = sum1 + ::pow(M_1_Sqrt2*sigma, n)*::pow(mu, k-n)/knfactorial[k-n] * sum2.imag();
                    }
                }
                return M_1_Sqrt8*sigma*kfactorial[k]*sum1;
            }
        }



        private:


        // Calculate overall polynomial coefficients for current set of spline coefficients.
        // This method must be included in the method Update() if you change the spline coefficients,
        // using the method SetParameters() of the class Parameters (See Parameters.h).
        void updatePolynomial()
        {
            for(size_t i=0; i<nKnots-1; i++)
            {
                // calculate polynomial coefficients for the current set of spline coefficients
                A0[i] = _par[i]*a0[0][i] + _par[i+1]*a0[1][i] + _par[i+2]*a0[2][i] + _par[i+3]*a0[3][i];
                A1[i] = _par[i]*a1[0][i] + _par[i+1]*a1[1][i] + _par[i+2]*a1[2][i] + _par[i+3]*a1[3][i];
                A2[i] = _par[i]*a2[0][i] + _par[i+1]*a2[1][i] + _par[i+2]*a2[2][i] + _par[i+3]*a2[3][i];
                A3[i] = _par[i]*a3[0][i] + _par[i+1]*a3[1][i] + _par[i+2]*a3[2][i] + _par[i+3]*a3[3][i];
                if(std::fabs(A0[i])<1e-9) A0[i]=0;
                if(std::fabs(A1[i])<1e-9) A1[i]=0;
                if(std::fabs(A2[i])<1e-9) A2[i]=0;
                if(std::fabs(A3[i])<1e-9) A3[i]=0;
            }

            // After last knot: Linear extrapolation
            // value of 2nd last spline at last knot
            double v=u[nKnots+2];
            A3[nKnots-1] = 0;
            A2[nKnots-1] = 0;

            // In case it should be constant:
            // A1[nKnots-1] = 0;
            // A0[nKnots-1] = A0[nKnots-2] + A1[nKnots-2]*v + A2[nKnots-2]*v*v + A3[nKnots-2]*v*v*v;

            // In case it should be linear:
            A1[nKnots-1] = A1[nKnots-2] + 2*A2[nKnots-2]*v + 3*A3[nKnots-2]*v*v;
            A0[nKnots-1] = A0[nKnots-2] + A1[nKnots-2]*v + A2[nKnots-2]*v*v + A3[nKnots-2]*v*v*v - A1[nKnots-1]*v;

            // determine if and where this linear part gets negative
            if(A1[nKnots-1]<0)
            {
                NegativePart = true;
                xNegative = -A0[nKnots-1]/A1[nKnots-1];
            }
            else
            {
                NegativePart = false;
            }
        }


        // K_n(z) for z as double (Reference: arXiv:1407.0748v1)
        __hydra_dual__
        inline double K(double z, size_t n)
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
        __hydra_dual__
        inline hydra::complex<double> K(hydra::complex<double> z, size_t n)
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
        __hydra_dual__
        inline double M(double x1, double x2, double z, size_t n)
        {
            switch (n)
            {
                case 0:
                {
                    return faddeeva::erf(x2) - ::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2) -
                                ( faddeeva::erf(x1) - ::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1) );
                }
                case 1:
                {
                    return 2*( -M_1_SqrtPi*::exp(-x2*x2) - x2*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -M_1_SqrtPi*::exp(-x1*x1) - x1*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 2:
                {
                    return 2*( -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 3:
                {
                    return 4*( -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) *
                                                                    ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) *
                                                                    ::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 4:
                {
                    return -4*( 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                        (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                        (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 5:
                {
                    return -8*( (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                            x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                            x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 6:
                {
                    return -8*( x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                            (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) *
                                                                    ::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
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
        __hydra_dual__
        inline hydra::complex<double> M(double x1, double x2, hydra::complex<double> z, size_t n)
        {
            switch (n)
            {
                case 0:
                {
                    return faddeeva::erf(x2) - hydra::exp( z*z - 2*z*x2 ) * faddeeva::erfc(z-x2) -
                                ( faddeeva::erf(x1) - hydra::exp( z*z - 2*z*x1 ) * faddeeva::erfc(z-x1) );
                }
                case 1:
                {
                    return 2*( -M_1_SqrtPi*::exp(-x2*x2) - x2*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -M_1_SqrtPi*::exp(-x1*x1) - x1*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 2:
                {
                    return 2*( -2*x2*M_1_SqrtPi*::exp(-x2*x2) - (2*x2*x2 - 1)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -2*x1*M_1_SqrtPi*::exp(-x1*x1) - (2*x1*x1 - 1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 3:
                {
                    return 4*( -(2*x2*x2 - 1)*M_1_SqrtPi*::exp(-x2*x2) - x2*(2*x2*x2 - 3) *
                                                            hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( -(2*x1*x1 - 1)*M_1_SqrtPi*::exp(-x1*x1) - x1*(2*x1*x1 - 3) *
                                                            hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 4:
                {
                    return -4*( 2*x2*(2*x2*x2 - 3)*M_1_SqrtPi*::exp(-x2*x2) +
                                        (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( 2*x1*(2*x1*x1 - 3)*M_1_SqrtPi*::exp(-x1*x1) +
                                        (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 5:
                {
                    return -8*( (3 - 12*x2*x2 + 4*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                        x2*(15 - 20*x2*x2 + 4*x2*x2*x2*x2)*hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( (3 - 12*x1*x1 + 4*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
                                        x1*(15 - 20*x1*x1 + 4*x1*x1*x1*x1)*hydra::exp( z*z - 2*z*x1 )*faddeeva::erfc(z-x1) ) );
                }
                case 6:
                {
                    return -8*( x2*(30 - 40*x2*x2 + 8*x2*x2*x2*x2)*M_1_SqrtPi*::exp(-x2*x2) +
                                    (-15 + 90*x2*x2 - 60*x2*x2*x2*x2 + 8*x2*x2*x2*x2*x2*x2) * 
                                                            hydra::exp( z*z - 2*z*x2 )*faddeeva::erfc(z-x2) -
                                ( x1*(30 - 40*x1*x1 + 8*x1*x1*x1*x1)*M_1_SqrtPi*::exp(-x1*x1) +
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



        // knots
        double u[nKnots+6];

        // factorials
        double kfactorial[4];    // k!
        double knfactorial[4];  // (k-n)!
        double jfactorial[4];    // j!
        double njfactorial[4];  // (n-j)!
  
        // Constants depending only on knot vector

        // coefficients of polynomial
        double a0[4][nKnots];
        double a1[4][nKnots];
        double a2[4][nKnots];
        double a3[4][nKnots];

        // Constants depending on the current spline coefficients

        // Overall coefficients of the polynomials. For the range x>=Knot[i] && x<Knot[i+1]
        // the spline is then given by: y=A0[i]+A1[i]*x+A2[i]*x*x+A3[i]*x*x*x
        double A0[nKnots], A1[nKnots], A2[nKnots], A3[nKnots];

        // The spline might get negative values due to the linear extrapolation in/after the last sector
        // Avoid this by checking if we are there and setting the spline to 0
        bool NegativePart;
        double xNegative;
    };

} // namespace medusa

#endif // MEDUSA_CUBIC_SPLINE_H