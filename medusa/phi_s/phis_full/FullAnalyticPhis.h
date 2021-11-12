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
 *  FullAnalyticPhis.h
 *
 *  Created on: 29/10/2021
 *      Author: Alessandro Maria Ricci
 * 
 *  Reference: arXiv:1906.08356v4
 *---------------------------------------------------------------------------*/

#ifndef FULL_ANALYTIC_PHIS_H_
#define FULL_ANALYTIC_PHIS_H_


// std
#include <cmath>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <ratio>


// Hydra
#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Integrator.h>
#include <hydra/Pdf.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>


// Medusa
#include <medusa/phi_s/Parameters.h>
#include <medusa/Faddeeva.h>


namespace medusa {

    /*
    *  @class FullAnalyticPhis
    *  Functor that provides the time dependent formula used in phi_s analysis in the full model,
    *  i.e. signal + experimental artifacts (tagging, time resolution and acceptances),
    *  with analytical convolution and integration
    *  [see Eq. (12) in arXiv:1906.08356v4]
    * 
    *  The implementation of the method Update_ATCoefficients() is inside the detail/ folder
    *
    *  ArgTypes = argument types of the functor
    */
    template<typename ArgTypeTime,
             typename ArgTypeThetah,
             typename ArgTypeThetal,
             typename ArgTypePhi,
             typename ArgTypeQOS,
             typename ArgTypeQSS,
             typename ArgTypeEtaOS,
             typename ArgTypeEtaSS,
             typename ArgTypeDelta,
             typename Signature=double(ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta) >
    class FullAnalyticPhis: public hydra::BaseFunctor< FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 18>
    {

        using ThisBaseFunctor = hydra::BaseFunctor< FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                                ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 18 >;

        using hydra::BaseFunctor< FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 18 >::_par;


        public:

        //-------------------------------------
        //           Constructors
        //-------------------------------------

        FullAnalyticPhis() = delete;

        // ctor with list of hydra::Parameter
        // the user has to respect the parameters order
        FullAnalyticPhis(hydra::Parameter const& A_0,           hydra::Parameter const& A_perp,     hydra::Parameter const& A_S, 
                       hydra::Parameter const& DeltaGamma_sd, hydra::Parameter const& DeltaGamma, hydra::Parameter const& DeltaM,
                       hydra::Parameter const& phi_0,         hydra::Parameter const& phi_par,
                       hydra::Parameter const& phi_perp,      hydra::Parameter const& phi_S,
                       hydra::Parameter const& lambda_0,      hydra::Parameter const& lambda_par,
                       hydra::Parameter const& lambda_perp,   hydra::Parameter const& lambda_S,
                       hydra::Parameter const& delta_0,       hydra::Parameter const& delta_par,
                       hydra::Parameter const& delta_perp,    hydra::Parameter const& delta_S):
        ThisBaseFunctor({A_0, A_perp, A_S, DeltaGamma_sd, DeltaGamma, DeltaM,
                         phi_0,       phi_par,    phi_perp,    phi_S,          lambda_0,   lambda_par,
                         lambda_perp, lambda_S,   delta_0,     delta_par,      delta_perp, delta_S })
        {
            Update();
        }


        // ctor with array of hydra::Parameter
        // the user has to respect the parameters order as the main ctor
        explicit FullAnalyticPhis( const hydra::Parameter (&Hs)[18] ):
        ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                         Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
        {
            Update();
        }


        // ctor with array of double
        // the user has to respect the parameters order as the main ctor
        explicit FullAnalyticPhis( const double (&Hs)[18] ):
        ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                         Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
        {
            Update();
        }


        // ctor with other FullAnalyticPhis instance (copy ctor)
        __hydra_dual__
        FullAnalyticPhis(FullAnalyticPhis<ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                        ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta> const& other):
        ThisBaseFunctor(other)
        {

    	    #pragma unroll 10
    	    for(size_t i=0; i<10; i++)
    	    {   
                // angular coefficients a_k, b_k, c_k, d_k in  Eq. (10) and (11) in arXiv:1906.08356v4
                // polarization factor N_k in Eq. (9) in arXiv:1906.08356v4
    		    A.k[i] = 	other.GetA().k[i];
    		    B.k[i] = 	other.GetB().k[i];
    		    C.k[i] = 	other.GetC().k[i];
    		    D.k[i] = 	other.GetD().k[i];
    		    N.k[i] =  other.GetN().k[i];
    	    }

        }


        //-------------------------------------
        //       Operator= overloading
        //-------------------------------------

        __hydra_dual__
        FullAnalyticPhis<ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta>& 
        operator=( FullAnalyticPhis<ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta> const& other)
        {
            if(this == &other) return *this;
            ThisBaseFunctor::operator=(other);

		    #pragma unroll 10
            for(size_t i=0; i<10; i++)
            {
                // angular coefficients a_k, b_k, c_k, d_k in  Eq. (10) and (11) in arXiv:1906.08356v4
                // polarization factor N_k in Eq. (9) in arXiv:1906.08356v4
        	    A.k[i] = 	other.GetA().k[i];
        	    B.k[i] = 	other.GetB().k[i];
        	    C.k[i] = 	other.GetC().k[i];
        	    D.k[i] = 	other.GetD().k[i];
        	    N.k[i] =  other.GetN().k[i];
            }

            return *this;
        }


        //-------------------------------------
        //         Service functions
        //-------------------------------------

        // update the values of the angular coefficients a_k, b_k, c_k, d_k and the polarization factor N_k
        // by using the formulas in Table 3 in arXiv:1906.08356v4
        void Update(void) override
        {
            Update_ATCoefficients();
            Update_NFactors();
        }


        // evaluate the sum in Eq. (9) in arXiv:1906.08356v4
        __hydra_dual__ 
        inline double Evaluate(ArgTypeTime time, ArgTypeThetah theta_h, ArgTypeThetal theta_l, ArgTypePhi phi,
                                    ArgTypeQOS qOS, ArgTypeQSS qSS, ArgTypeEtaOS etaOS, ArgTypeEtaSS etaSS, ArgTypeDelta delta_time)  const
        {
            auto F = parameters::AngularFunctions(theta_h, theta_l, phi);

            double TagB0s = B0sTag(qOS, qSS, etaOS, etaSS);
            double TagB0sbar = B0sbarTag(qOS, qSS, etaOS, etaSS);

            double sigma_eff = Sigma_eff(delta_time);

            double result = 0;

            #pragma unroll 10
            for(size_t i=0; i<10; i++)
            {
            	result += F.fk[i]*N.k[i]*( TagB0s*Convoluted_Time_Factor(i, time, sigma_eff, 1) +
                                                TagB0sbar*Convoluted_Time_Factor(i, time, sigma_eff, -1) );
            }

            // This macro controls if result is NaN. If yes, it prints a warning with the parameter values for whom we obtain a NaN.
            hydra::CHECK_VALUE(result, "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f, par[4]=%f, par[5]=%f, "
                                       "par[6]=%f, par[7]=%f, par[8]=%f, par[9]=%f, par[10]=%f, par[11]=%f, "
                                       "par[12]=%f, par[13]=%f, par[14]=%f, par[15]=%f, par[16]=%f, par[17]=%f",
                                        _par[0], _par[1], _par[2], _par[3], _par[4], _par[5],
                                        _par[6], _par[7], _par[8], _par[9], _par[10], _par[11],
                                        _par[12], _par[13], _par[14], _par[15], _par[16], _par[17]);

            // This is a safety mechanism that is necessary when the functor takes negative values due to the numerical errors.
            // Don't use the function ::abs(), because it changes the values of about 10^-03-10^-04 units.
            // Don't deactivate this mechanism, otherwise, there could be anomalous behaviors in the fcn computation.
            if(result < 0)
            {
                return result = 2.225073858507e-308;
            }
            else
            {
                return result;
            }
        }


        //-------------------------------------
        //           Selectors
        //-------------------------------------

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetA() const
        {
		    return A;
	    }

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetB() const
        {
		    return B;
	    }

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetC() const
        {
		    return C;
	    }

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetD() const
        {
		    return D;
	    }

        __hydra_dual__
	    const parameters::NFactors& GetN() const
        {
		    return N;
	    }



        private:


        // update the values of the angular coefficients a_k, b_k, c_k, d_k
        // by using the formulas in Table 3 in arXiv:1906.08356v4
        // (see implementation in Update_ATCoefficients.inl)
        void Update_ATCoefficients();


        // update the values of the polarization factor N_k
        // by using the formulas in Table 3 in arXiv:1906.08356v4
        void Update_NFactors()
        {
    	    /*
    	    0: A_0,
            1: A_perp,
    	    2: A_S,
    	    3: DeltaGamma_sd,
    	    4: DeltaGamma,
    	    5: DeltaM,
    	    */
    	    double A_par = ::sqrt(1 - _par[0]*_par[0] -_par[1]*_par[1]);

    	    N.k[0] = _par[0]*_par[0];     //A_0*A_0 ;
    	    N.k[1] =   A_par*A_par;       //A_par*A_par
    	    N.k[2] = _par[1]*_par[1];     //A_perp*A_perp;
    	    N.k[3] = _par[1]*A_par;       //A_perp*A_par;
    	    N.k[4] = _par[0]*A_par;       //A_0*A_par;
    	    N.k[5] = _par[0]*_par[1];     //A_0*A_perp;
    	    N.k[6] = _par[2]*_par[2];     //A_S*A_S;
    	    N.k[7] = _par[2]*A_par;       //A_S*A_par;
    	    N.k[8] = _par[2]*_par[1];     //A_S*A_perp;
    	    N.k[9] = _par[2]*_par[0];     //A_S*A_0;
        }


        // time factors h_k(t|Bs0) and h_k(t|Bs0bar) in Eq. (10) and (11) in arXiv:1906.08356v4
        __hydra_dual__
        inline double Convoluted_Time_Factor(int index, double time, double sigma_eff, int Tag) const
        {
		    /*
		    0: A_0,
		    1: A_perp,
		    2: A_S,
		    3: DeltaGamma_sd,
		    4: DeltaGamma,
		    5: DeltaM,
		    */

    	    static const double f = 0.05968310365946074; // 3./(16*PI)
            static const double Sqrt2 = 1.414213562373095; // \sqrt{2}
            hydra::complex<double> I(0.0, 1.0); // Imaginary unit

            double Gamma = _par[3] + 0.65789;
            double x = time/(sigma_eff*Sqrt2);

            double z1 = (Gamma - 0.5*_par[4])*sigma_eff/Sqrt2;
            double z2 = (Gamma + 0.5*_par[4])*sigma_eff/Sqrt2;
            hydra::complex<double> Z_1(0.0, z1 - x);
            hydra::complex<double> Z_2(0.0, z2 - x);

            double Re_Z3 =  _par[5]*sigma_eff/Sqrt2;
            double Re_Z4 = -_par[5]*sigma_eff/Sqrt2;
            double Im_Z3 = Gamma*sigma_eff/Sqrt2 - x;
            double Im_Z4 = Gamma*sigma_eff/Sqrt2 - x;
            hydra::complex<double> Z_3(Re_Z3, Im_Z3);
            hydra::complex<double> Z_4(Re_Z4, Im_Z4);

            hydra::complex<double> faddeeva_sum12  = faddeeva::w(Z_1) + faddeeva::w(Z_2);
            hydra::complex<double> faddeeva_diff12 = faddeeva::w(Z_1) - faddeeva::w(Z_2);
            hydra::complex<double> faddeeva_sum34  = faddeeva::w(Z_3) + faddeeva::w(Z_4);
            hydra::complex<double> faddeeva_diff34 = (faddeeva::w(Z_3) - faddeeva::w(Z_4)) / I;

            double Re_faddeeva_sum12 = faddeeva_sum12.real();
            double Re_faddeeva_sum34 = faddeeva_sum34.real();
            double Re_faddeeva_diff12 = faddeeva_diff12.real();
            double Re_faddeeva_diff34 = faddeeva_diff34.real();

            double ConvolutedTimeFactor = f * ::exp( - ::pow(x, 2.0) ) *
                                            ( A.k[index] * Re_faddeeva_sum12 + B.k[index] * Re_faddeeva_diff12 +
                                                Tag*( C.k[index] * Re_faddeeva_sum34 + D.k[index] * Re_faddeeva_diff34 ) );

            return ConvolutedTimeFactor;
        }


        // effective resolution (See definition in page 7 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double Sigma_eff(ArgTypeDelta delta_time) const
        {
            static const double b0 = 0.01297;
            static const double b1 = 0.8446;

            double sigma_eff = b0 + b1 * delta_time;

            return sigma_eff;
        }


        // Tagging of B0s (See Eq. (5), (6), (13) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            static const double p0_OS = 0.3890;
            static const double p1_OS = 0.8486;
            static const double DeltaP0_OS = 0.009;
            static const double DeltaP1_OS = 0.0143;
            static const double AvgEta_OS = 0.3602;

            static const double p0_SS = 0.4325;
            static const double p1_SS = 0.9241;
            static const double DeltaP0_SS = 0.;
            static const double DeltaP1_SS = 0.;
            static const double AvgEta_SS = 0.4167;

            double omegaOS = (p0_OS + 0.5*DeltaP0_OS) + (p1_OS + 0.5*DeltaP1_OS)*(etaOS - AvgEta_OS);
            double omegaSS = (p0_SS + 0.5*DeltaP0_SS) + (p1_SS + 0.5*DeltaP1_SS)*(etaSS - AvgEta_SS);

            double TagB0s = (1 + qOS*(1 - 2*omegaOS))*(1 + qSS*(1 - 2*omegaSS));

            return TagB0s;
        }


        // Tagging of B0sbar (See Eq. (5), (6), (14) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sbarTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            static const double p0_OS = 0.3890;
            static const double p1_OS = 0.8486;
            static const double DeltaP0_OS = 0.009;
            static const double DeltaP1_OS = 0.0143;
            static const double AvgEta_OS = 0.3602;

            static const double p0_SS = 0.4325;
            static const double p1_SS = 0.9241;
            static const double DeltaP0_SS = 0.;
            static const double DeltaP1_SS = 0.;
            static const double AvgEta_SS = 0.4167;

            double OmegaBarOS = (p0_OS - 0.5*DeltaP0_OS) + (p1_OS - 0.5*DeltaP1_OS)*(etaOS - AvgEta_OS);
            double OmegaBarSS = (p0_SS - 0.5*DeltaP0_SS) + (p1_SS - 0.5*DeltaP1_SS)*(etaSS - AvgEta_SS);

            double TagB0sbar = (1 - qOS*(1 - 2*OmegaBarOS))*(1 - qSS*(1 - 2*OmegaBarSS));

            return TagB0sbar;
        }



        parameters::NFactors N;                     // polarization factor N_k
        parameters::AngularTimeCoefficients A;      // angular coefficient a_k
        parameters::AngularTimeCoefficients B;      // angular coefficient b_k
        parameters::AngularTimeCoefficients C;      // angular coefficient c_k
        parameters::AngularTimeCoefficients D;      // angular coefficient d_k

    };

} // namespace medusa

// Medusa
#include <medusa/phi_s/phis_full/details/Update_ATCoefficients.inl>
#include <medusa/phi_s/phis_full/details/IntegrationFormula.inl>

#endif /* FULL_ANALYTIC_PHIS_H_ */