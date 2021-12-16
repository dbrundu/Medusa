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

// std
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>

// Hydra
#include <hydra/detail/Config.h>
#include <hydra/device/System.h>
#include <hydra/Distance.h>
#include <hydra/detail/external/hydra_thrust/binary_search.h>

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

    template<size_t nKnots>
    class CubicSpline
    {
        public:

        //-------------------------------------
        //           Constructors
        //-------------------------------------

        CubicSpline() = delete;

        CubicSpline(const double *knots, const double *coeffs)
        {
            // set knot vector
            u[0] = knots[0];
            u[1] = knots[0];
            u[2] = knots[0];
            for(size_t i=0; i<nKnots; i++)
            {
                u[3+i] = knots[i];
                b[i] = coeffs[i];
            }
            b[nKnots] = coeffs[nKnots];
            b[nKnots+1] = coeffs[nKnots+1];
            u[nKnots+3] = knots[nKnots-1];
            u[nKnots+4] = knots[nKnots-1];
            u[nKnots+5] = knots[nKnots-1];

            // calculate prefactors
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
                a0[0] = u[i+4]*u[i+4]*u[i+4]/P[i];

                a0[1] = -u[i+1]*u[i+4]*u[i+4]/P[i] - u[i+2]*u[i+4]*u[i+5]/Q[i] - u[i+3]*u[i+5]*u[i+5]/R[i];

                a0[2] = u[i+2]*u[i+2]*u[i+4]/Q[i] + u[i+2]*u[i+3]*u[i+5]/R[i] + u[i+3]*u[i+3]*u[i+6]/S[i];

                a0[3] = -u[i+3]*u[i+3]*u[i+3]/S[i];
    
                // Linear
                a1[0] = -3*u[i+4]*u[i+4]/P[i];

                a1[1] = (2*u[i+1]*u[i+4]+u[i+4]*u[i+4])/P[i]
                        + (u[i+2]*u[i+4]+u[i+2]*u[i+5]+u[i+4]*u[i+5])/Q[i]
                        + (2*u[i+3]*u[i+5]+u[i+5]*u[i+5])/R[i];

                a1[2] = -(2*u[i+2]*u[i+4]+u[i+2]*u[i+2])/Q[i]
                        -(u[i+2]*u[i+3]+u[i+2]*u[i+5]+u[i+3]*u[i+5])/R[i]
                        -(2*u[i+3]*u[i+6]+u[i+3]*u[i+3])/S[i];

                a1[3] = 3*u[i+3]*u[i+3]/S[i];
    
                // Quadratic
                a2[0] = 3*u[i+4]/P[i];

                a2[1] = -(2*u[i+4]+u[i+1])/P[i]
                        -(u[i+2]+u[i+4]+u[i+5])/Q[i]
                        -(2*u[i+5]+u[i+3])/R[i];

                a2[2] = (2*u[i+2]+u[i+4])/Q[i]
                        +(u[i+2]+u[i+5]+u[i+3])/R[i]
                        +(2*u[i+3]+u[i+6])/S[i];

                a2[3] = -3*u[i+3]/S[i];
    
                // Cubic
                a3[0] = -1./P[i];

                a3[1] = 1./P[i] + 1./Q[i] + 1./R[i];

                a3[2] = -1./Q[i] - 1./R[i] - 1./S[i];

                a3[3] = 1./S[i];

                // calculate polynomial coefficients for the current set of spline coefficients
                A0[i] = b[i]*a0[0] + b[i+1]*a0[1] + b[i+2]*a0[2] + b[i+3]*a0[3];
                A1[i] = b[i]*a1[0] + b[i+1]*a1[1] + b[i+2]*a1[2] + b[i+3]*a1[3];
                A2[i] = b[i]*a2[0] + b[i+1]*a2[1] + b[i+2]*a2[2] + b[i+3]*a2[3];
                A3[i] = b[i]*a3[0] + b[i+1]*a3[1] + b[i+2]*a3[2] + b[i+3]*a3[3];
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
                getsNegative = true;
                xNegative = -A0[nKnots-1]/A1[nKnots-1];
            }
            else
            {
                getsNegative = false;
            }
        }


        //-------------------------------------
        //           Service functions
        //-------------------------------------

        // Update the spline with new coefficients
        void updateCoefficients(double *coeffs)
        {
            for(size_t i=0; i<nKnots+2; i++)
            {
                b[i]=coeffs[i];
            }
        
            // calculate overall polynomial coefficients for current set of spline coefficients
            updatePolynomial();
        }


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
            if(x>xNegative && getsNegative) return 1e-3;

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
  
        // get the order 0 coefficients
        __hydra_dual__
        const double* GetCoeffO0() const
        {
            return A0;
        }

        // get the order 1 coefficients
        __hydra_dual__
        const double* GetCoeffO1() const
        {
            return A1;
        }
        
        // get the order 2 coefficients
        __hydra_dual__
        const double* GetCoeffO2() const
        {
            return A2;
        }
        
        // get the order 3 coefficients
        __hydra_dual__
        const double* GetCoeffO3() const
        {
            return A3;
        }

        // get the info wether this spline will get negative at some point
        __hydra_dual__
        bool NegativePart() const
        {
            return getsNegative;
        }
        
        __hydra_dual__
        double NegativePartPos() const
        {
            return xNegative;
        }



        private:

        // calculate overall polynomial coefficients for current set of spline coefficients
        void updatePolynomial()
        {
            for(size_t i=0; i<nKnots-1; i++)
            {
                A0[i] = b[i]*a0[0] + b[i+1]*a0[1] + b[i+2]*a0[2] + b[i+3]*a0[3];
                A1[i] = b[i]*a1[0] + b[i+1]*a1[1] + b[i+2]*a1[2] + b[i+3]*a1[3];
                A2[i] = b[i]*a2[0] + b[i+1]*a2[1] + b[i+2]*a2[2] + b[i+3]*a2[3];
                A3[i] = b[i]*a3[0] + b[i+1]*a3[1] + b[i+2]*a3[2] + b[i+3]*a3[3];
                if(std::fabs(A0[i])<1e-9) A0[i]=0;
                if(std::fabs(A1[i])<1e-9) A1[i]=0;
                if(std::fabs(A2[i])<1e-9) A2[i]=0;
                if(std::fabs(A3[i])<1e-9) A3[i]=0;
            }

            // linear last sector
            double v=u[nKnots+2];
            A3[nKnots-1] = 0;
            A2[nKnots-1] = 0;

            // In case it should be constant:
            // A1[nKnots-1] = 0;
            // A0[nKnots-1] = A0[nKnots-2] + A1[nKnots-2]*v + A2[nKnots-2]*v*v + A3[nKnots-2]*v*v*v;

            // In case it should be linear
            A1[nKnots-1] = A1[nKnots-2] + 2*A2[nKnots-2]*v + 3*A3[nKnots-2]*v*v;
            A0[nKnots-1] = A0[nKnots-2] + A1[nKnots-2]*v + A2[nKnots-2]*v*v + A3[nKnots-2]*v*v*v - A1[nKnots-1]*v;

            // determine if and where this linear part gets negative
            if(A1[nKnots-1]<0)
            {
                getsNegative = true;
                xNegative = -A0[nKnots-1]/A1[nKnots-1];
            }
            else
            {
                getsNegative = false;
            }
        }



        // knots
        double u[nKnots+6];

        // spline coefficients
        double b[nKnots+2];
  
        // Constants depending only on knotvector
  
        // prefactors
        double P[nKnots], Q[nKnots], R[nKnots], S[nKnots];

        // coefficients of polynomial
        double a0[4];
        double a1[4];
        double a2[4];
        double a3[4];

        // Constants depending on the current spline coefficients

        // Overall coefficients of the polynomials. For the range x>=Knot[i] && x<Knot[i+1]
        // the spline is then given by: y=A0[i]+A1[i]*x+A2[i]*x*x+A3[i]*x*x*x
        double A0[nKnots], A1[nKnots], A2[nKnots], A3[nKnots];

        // The spline might get negative values due to the linear extrapolation in/after the last sector
        // Avoid this by checking if we are there and setting the spline to 0
        bool getsNegative;
        double xNegative;
    };

} // namespace medusa

#endif // MEDUSA_CUBIC_SPLINE_H