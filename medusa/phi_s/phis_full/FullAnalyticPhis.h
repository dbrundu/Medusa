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
#include <hydra/Complex.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>


// Medusa
#include <medusa/phi_s/Parameters.h>
#include <medusa/Functions.h>


namespace medusa {

    /*
    *  @class FullAnalyticPhis
    *  Functor that provides the time dependent formula used in phi_s analysis in the full model,
    *  i.e. signal + experimental artifacts (tagging, time resolution and acceptances),
    *  with analytical convolution and integration [Reference: Eq. (12) in arXiv:1906.08356v4]
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
                                                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 40>
    {

        using ThisBaseFunctor = hydra::BaseFunctor< FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                                ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 40 >;

        using hydra::BaseFunctor< FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 40 >::_par;


        public:

        //-------------------------------------
        //           Constructors
        //-------------------------------------

        FullAnalyticPhis() = delete;

        // ctor with list of hydra::Parameter
        // the user has to respect the parameter order
        FullAnalyticPhis(hydra::Parameter const& A_0,             hydra::Parameter const& A_perp,      hydra::Parameter const& A_S,
                         hydra::Parameter const& DeltaGamma_sd,   hydra::Parameter const& DeltaGamma,  hydra::Parameter const& DeltaM,
                         hydra::Parameter const& phi_0,           hydra::Parameter const& phi_par,     hydra::Parameter const& phi_perp,
                         hydra::Parameter const& phi_S,           hydra::Parameter const& lambda_0,    hydra::Parameter const& lambda_par,
                         hydra::Parameter const& lambda_perp,     hydra::Parameter const& lambda_S,    hydra::Parameter const& delta_0,
                         hydra::Parameter const& delta_par,       hydra::Parameter const& delta_perp,  hydra::Parameter const& delta_S,
                         hydra::Parameter const& b0,              hydra::Parameter const& b1,          hydra::Parameter const& p0_OS,
                         hydra::Parameter const& p1_OS,           hydra::Parameter const& DeltaP0_OS,  hydra::Parameter const& DeltaP1_OS,
                         hydra::Parameter const& AvgEta_OS,       hydra::Parameter const& p0_SS,       hydra::Parameter const& p1_SS,
                         hydra::Parameter const& DeltaP0_SS,      hydra::Parameter const& DeltaP1_SS,  hydra::Parameter const& AvgEta_SS,
                         hydra::Parameter const& Omega_1,         hydra::Parameter const& Omega_2,     hydra::Parameter const& Omega_3,
                         hydra::Parameter const& Omega_4,         hydra::Parameter const& Omega_5,     hydra::Parameter const& Omega_6,
                         hydra::Parameter const& Omega_7,         hydra::Parameter const& Omega_8,     hydra::Parameter const& Omega_9,
                         hydra::Parameter const& Omega_10,        ArgTypeTime const& LowerLimit,       ArgTypeTime const& UpperLimit,       double const& Weight):
        ThisBaseFunctor({A_0, A_perp, A_S,
                         DeltaGamma_sd, DeltaGamma, DeltaM,
                         phi_0, phi_par, phi_perp, phi_S, 
                         lambda_0, lambda_par, lambda_perp, lambda_S,
                         delta_0, delta_par, delta_perp, delta_S,
                         b0, b1,
                         p0_OS, p1_OS, DeltaP0_OS, DeltaP1_OS, AvgEta_OS,
                         p0_SS, p1_SS, DeltaP0_SS, DeltaP1_SS, AvgEta_SS,
                         Omega_1, Omega_2, Omega_3, Omega_4, Omega_5,
                         Omega_6, Omega_7, Omega_8, Omega_9, Omega_10 })
        {
            fLowerLimit = LowerLimit;
            fUpperLimit = UpperLimit;
            fweight = Weight;
            Update();
        }


        // ctor with array of hydra::Parameter
        // the user has to respect the parameter order as the main ctor
        explicit FullAnalyticPhis( const hydra::Parameter (&Hs)[18], const hydra::Parameter (&Ps)[22],
                                                        const ArgTypeTime &LowerLimit, const ArgTypeTime &UpperLimit, const double &Weight ):
        ThisBaseFunctor{ Hs[0],  Hs[1],  Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6],  Hs[7],
                         Hs[8],  Hs[9],  Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17],
                         Ps[0],  Ps[1],  Ps[2],  Ps[3],  Ps[4],  Ps[5],  Ps[6],  Ps[7],
                         Ps[8],  Ps[9],  Ps[10], Ps[11], Ps[12], Ps[13], Ps[14], Ps[15], Ps[16], Ps[17],
                         Ps[18], Ps[19], Ps[20], Ps[21] }
        {
            fLowerLimit = LowerLimit;
            fUpperLimit = UpperLimit;
            fweight = Weight;
            Update();
        }


        // ctor with array of double
        // the user has to respect the parameter order as the main ctor
        explicit FullAnalyticPhis( const double (&Hs)[43] ):
        ThisBaseFunctor{ Hs[0],  Hs[1],  Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6],  Hs[7],
                         Hs[8],  Hs[9],  Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17],
                         Hs[18], Hs[19], Hs[20], Hs[21], Hs[22], Hs[23], Hs[24], Hs[25],
                         Hs[26], Hs[27], Hs[28], Hs[29], Hs[30], Hs[31], Hs[32], Hs[33], Hs[34], Hs[35],
                         Hs[36], Hs[37], Hs[38], Hs[39] }
        {
            fLowerLimit = Hs[40];
            fUpperLimit = Hs[41];
            fweight = Hs[42];
            Update();
        }


        // ctor with other FullAnalyticPhis instance (copy ctor)
        __hydra_dual__
        FullAnalyticPhis(FullAnalyticPhis<ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                        ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta> const& other):
        ThisBaseFunctor(other)
        {
            fLowerLimit = other.GetLowerLimit();
            fUpperLimit = other.GetUpperLimit();
            fweight = other.GetWeight();

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

            fLowerLimit = other.GetLowerLimit();
            fUpperLimit = other.GetUpperLimit();
            fweight = other.GetWeight();

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
                                ArgTypeQOS qOS, ArgTypeQSS qSS, ArgTypeEtaOS etaOS, ArgTypeEtaSS etaSS, ArgTypeDelta delta_time) const
        {
            /*
		     0: A_0,
		     1: A_perp,
		     2: A_S,
		     3: DeltaGamma_sd,
		     4: DeltaGamma,
		     5: DeltaM,
            30: Omega_1,
            31: Omega_2,
            32: Omega_3,
            33: Omega_4,
            34: Omega_5,
            35: Omega_6,
            36: Omega_7,
            37: Omega_8,
            38: Omega_9,
            39: Omega_10,
		    */

            double A_par2 = 1 - _par[0]*_par[0] - _par[1]*_par[1];

            double UnnormPDF = 0;
            double NormFactor = 0;
            double wPDF = 0;

            // This is a safety mechanism that is necessary when A_par2 takes negative values (see Update_NFactors()).
            // PDF = 0 activates the Hydra safety mechanism for whom FCN = FcnMaxValue (see main function).
            if(A_par2 < 0)
            {
                return wPDF;
            }

            auto F = parameters::AngularFunctions(theta_h, theta_l, phi);

            double TagB0s = B0sTag(qOS, qSS, etaOS, etaSS);
            double TagB0sbar = B0sbarTag(qOS, qSS, etaOS, etaSS);

            double sigma_eff = Sigma_eff(delta_time);

            double Gamma = _par[3] + 0.65789;
            double HalfDeltaGamma = 0.5*_par[4];

            double conv_exp_cosh = functions::Convoluted_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, 1);
            double conv_exp_sinh = functions::Convoluted_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, -1);
            double conv_exp_cos = functions::Convoluted_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, 1);
            double conv_exp_sin = functions::Convoluted_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, -1);

            double int_conv_exp_cosh = functions::Integrated_convoluted_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, 1);
            double int_conv_exp_sinh = functions::Integrated_convoluted_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, -1);
            double int_conv_exp_cos = functions::Integrated_convoluted_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, 1);
            double int_conv_exp_sin = functions::Integrated_convoluted_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, -1);

            #pragma unroll 10
            for(size_t i=0; i<10; i++)
            {
            	UnnormPDF += F.fk[i]*N.k[i]*( TagB0s*Convoluted_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, 1) +
                                                TagB0sbar*Convoluted_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, -1) );
                
                NormFactor +=  _par[30+i]*N.k[i]*
                                ( TagB0s*Integrated_Convoluted_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, 1) +
                                    TagB0sbar*Integrated_Convoluted_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, -1) );
            }

            wPDF = fweight*UnnormPDF/NormFactor;

            // This macro controls if PDF is NaN. If yes, it prints a warning
            // with the parameter value for whom we obtain a NaN.
            hydra::CHECK_VALUE(wPDF, "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f, par[4]=%f, par[5]=%f, "
                                     "par[6]=%f, par[7]=%f, par[8]=%f, par[9]=%f, par[10]=%f, par[11]=%f, "
                                     "par[12]=%f, par[13]=%f, par[14]=%f, par[15]=%f, par[16]=%f, par[17]=%f",
                                     _par[0], _par[1], _par[2], _par[3], _par[4], _par[5],
                                     _par[6], _par[7], _par[8], _par[9], _par[10], _par[11],
                                     _par[12], _par[13], _par[14], _par[15], _par[16], _par[17]);

            // This is a safety mechanism that is necessary when the functor takes negative values due to the numerical errors.
            // Don't use the function ::abs(), because it changes the values of about 10^{-03} - 10^{-04} units.
            // Don't deactivate this mechanism, otherwise, there could be anomalous behaviors in the fcn computation.
            if(wPDF < 0)
            {
                return wPDF = 0;
            }
            else
            {
                return wPDF;
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

        __hydra_dual__
        const ArgTypeTime GetLowerLimit() const
        {
            return fLowerLimit;
        }

        __hydra_dual__
        const ArgTypeTime GetUpperLimit() const
        {
            return fUpperLimit;
        }

        __hydra_dual__
        double GetWeight() const
        {
            return fweight;
        }


        //-------------------------------------
        //          Setting functions
        //-------------------------------------

        __hydra_dual__
        void SetLowerLimit(ArgTypeTime LowerLimit)
        {
            fLowerLimit = LowerLimit;
        }

        __hydra_dual__
        void SetUpperLimit(ArgTypeTime UpperLimit)
        {
            fUpperLimit = UpperLimit;
        }

        __hydra_dual__
        void SetWeight(double Weight)
        {
            fweight = Weight;
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

           double A_par2 = 1 - _par[0]*_par[0] - _par[1]*_par[1];

           if(A_par2 > 0)
           {
                double A_par = ::sqrt(A_par2);

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
        }


        // time factors h_k(t|Bs0) and h_k(t|Bs0bar) convoluted with the Gaussian
        // (Reference: Eq. (10) and (11) in arXiv:1906.08356v4)
        __hydra_dual__
        inline double Convoluted_Time_Factor(size_t index, double conv_exp_cosh, double conv_exp_sinh,
                                                                double conv_exp_cos, double conv_exp_sin, int Tag) const
        {
    	    static const double f = 0.2387324146378430; // 3./(4.*PI)

            double ConvolutedTimeFactor = f * ( A.k[index] * conv_exp_cosh + B.k[index] * conv_exp_sinh +
                                                    Tag * ( C.k[index] * conv_exp_cos + D.k[index] * conv_exp_sin ) );

            return ConvolutedTimeFactor;
        }


        // time factors h_k(t|Bs0) and h_k(t|Bs0bar) convoluted with the Gaussian and integrated in the time
        // (Reference: Eq. (10) and (11) in arXiv:1906.08356v4)
        __hydra_dual__
        inline double Integrated_Convoluted_Time_Factor(size_t index, double int_conv_exp_cosh, double int_conv_exp_sinh,
                                                                        double int_conv_exp_cos, double int_conv_exp_sin, int Tag) const
        {
            static const double f = 0.2387324146378430; // 3./(4.*PI)

            double NormFactor = f * ( A.k[index] * int_conv_exp_cosh + B.k[index] * int_conv_exp_sinh +
                                            Tag * ( C.k[index] * int_conv_exp_cos + D.k[index] * int_conv_exp_sin ) );

            return NormFactor;
        }


        // effective resolution (Reference: page 7 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double Sigma_eff(ArgTypeDelta delta_time) const
        {
            /*
		    18: b0,
		    19: b1,
		    */

            double sigma_eff = _par[18] + _par[19] * delta_time;

            return sigma_eff;
        }


        // Tagging of B0s (Reference: Eq. (5), (6), (13) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            /*
		    20: p0_OS,
		    21: p1_OS,
		    22: DeltaP0_OS,
		    23: DeltaP1_OS,
		    24: AvgEta_OS,
		    25: p0_SS,
		    26: p1_SS,
		    27: DeltaP0_SS,
		    28: DeltaP1_SS,
		    29: AvgEta_SS,
		    */

            double omegaOS = (_par[20] + 0.5*_par[22]) + (_par[21] + 0.5*_par[23])*(etaOS - _par[24]);
            double omegaSS = (_par[25] + 0.5*_par[27]) + (_par[26] + 0.5*_par[28])*(etaSS - _par[29]);

            double TagB0s = (1 + qOS*(1 - 2*omegaOS))*(1 + qSS*(1 - 2*omegaSS));

            return TagB0s;
        }


        // Tagging of B0sbar (Reference: Eq. (5), (6), (14) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sbarTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            /*
		    20: p0_OS,
		    21: p1_OS,
		    22: DeltaP0_OS,
		    23: DeltaP1_OS,
		    24: AvgEta_OS,
		    25: p0_SS,
		    26: p1_SS,
		    27: DeltaP0_SS,
		    28: DeltaP1_SS,
		    29: AvgEta_SS,
		    */

            double OmegaBarOS = (_par[20] - 0.5*_par[22]) + (_par[21] - 0.5*_par[23])*(etaOS - _par[24]);
            double OmegaBarSS = (_par[25] - 0.5*_par[27]) + (_par[26] - 0.5*_par[28])*(etaSS - _par[29]);

            double TagB0sbar = (1 - qOS*(1 - 2*OmegaBarOS))*(1 - qSS*(1 - 2*OmegaBarSS));

            return TagB0sbar;
        }



        double fweight;                             // weight to improve the numerical fit
        ArgTypeTime fLowerLimit;                    // Lower limit in the time integration
        ArgTypeTime fUpperLimit;                    // Upper limit in the time integration
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