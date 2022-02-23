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
#include <hydra/Complex.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>


// Medusa
#include <medusa/phi_s/Parameters.h>
#include <medusa/generic/Functions.h>
#include <medusa/generic/CubicSpline.h>


namespace medusa {

    /*
    *  @class FullAnalyticPhis
    *  Functor that provides the time dependent formula used in phi_s analysis in the full model,
    *  i.e. signal + experimental artifacts (tagging, time resolution and acceptances),
    *  with analytical convolution and integration [Reference: Eq. (12) in arXiv:1906.08356v4]
    * 
    *  The implementation of the method Update_ATCoefficients() is inside the detail/ folder
    *
    *  Spline = enable or disable the cubic spline
    *  ArgTypes = argument types of the functor
    * 
    *  WARNING: all attributes of the functor are incapsulated as arguments of the CUDA kernels.
    *  The maximum dimension of all CUDA kernel arguments is 4 KB.
    */
    template<bool Spline,
             typename ArgTypeTime,
             typename ArgTypeCosThetah,
             typename ArgTypeCosThetal,
             typename ArgTypePhi,
             typename ArgTypeQOS,
             typename ArgTypeQSS,
             typename ArgTypeEtaOS,
             typename ArgTypeEtaSS,
             typename ArgTypeDelta,
             typename Signature=double(ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta) >
    class FullAnalyticPhis: public hydra::BaseFunctor< FullAnalyticPhis< Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 48>,
                            public CubicSpline<7>
    {

        using ThisBaseFunctor = hydra::BaseFunctor< FullAnalyticPhis< Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                                                ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 48 >;
        using CSpline = CubicSpline<7>;

        using hydra::BaseFunctor< FullAnalyticPhis< Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                                    ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, Signature, 48 >::_par;


        public:

        //-------------------------------------
        //           Constructors
        //-------------------------------------

        FullAnalyticPhis() = delete;

        // ctor with list of hydra::Parameter
        // the user has to respect the parameter order
        FullAnalyticPhis(hydra::Parameter const& A_02,            hydra::Parameter const& A_perp2,     hydra::Parameter const& A_S2,
                         hydra::Parameter const& DeltaGamma_sd,   hydra::Parameter const& DeltaGamma,  hydra::Parameter const& DeltaM,
                         hydra::Parameter const& phi_0,           hydra::Parameter const& phi_par0,    hydra::Parameter const& phi_perp0,
                         hydra::Parameter const& phi_S0,          hydra::Parameter const& lambda_0,    hydra::Parameter const& lambda_par0,
                         hydra::Parameter const& lambda_perp0,    hydra::Parameter const& lambda_S0,   hydra::Parameter const& delta_par0,
                         hydra::Parameter const& delta_perp0,     hydra::Parameter const& delta_Sperp, hydra::Parameter const& b0,
                         hydra::Parameter const& b1,              hydra::Parameter const& p0_OS,       hydra::Parameter const& p1_OS,
                         hydra::Parameter const& DeltaP0_OS,      hydra::Parameter const& DeltaP1_OS,  hydra::Parameter const& AvgEta_OS,
                         hydra::Parameter const& p0_SS,           hydra::Parameter const& p1_SS,       hydra::Parameter const& DeltaP0_SS,
                         hydra::Parameter const& DeltaP1_SS,      hydra::Parameter const& AvgEta_SS,   hydra::Parameter const& Omega_1,
                         hydra::Parameter const& Omega_2,         hydra::Parameter const& Omega_3,     hydra::Parameter const& Omega_4,
                         hydra::Parameter const& Omega_5,         hydra::Parameter const& Omega_6,     hydra::Parameter const& Omega_7,
                         hydra::Parameter const& Omega_8,         hydra::Parameter const& Omega_9,     hydra::Parameter const& Omega_10,
                         hydra::Parameter const& Spline_c0,       hydra::Parameter const& Spline_c1,   hydra::Parameter const& Spline_c2,
                         hydra::Parameter const& Spline_c3,       hydra::Parameter const& Spline_c4,   hydra::Parameter const& Spline_c5,
                         hydra::Parameter const& Spline_c6,       hydra::Parameter const& Spline_c7,   hydra::Parameter const& Spline_c8,
                         double const (&SplineKnots)[7],          ArgTypeTime const& LowerLimit,       ArgTypeTime const& UpperLimit):
        ThisBaseFunctor({A_02, A_perp2, A_S2,
                         DeltaGamma_sd, DeltaGamma, DeltaM,
                         phi_0, phi_par0, phi_perp0, phi_S0,
                         lambda_0, lambda_par0, lambda_perp0, lambda_S0,
                         delta_par0, delta_perp0, delta_Sperp,
                         b0, b1,
                         p0_OS, p1_OS, DeltaP0_OS, DeltaP1_OS, AvgEta_OS,
                         p0_SS, p1_SS, DeltaP0_SS, DeltaP1_SS, AvgEta_SS,
                         Omega_1, Omega_2, Omega_3, Omega_4, Omega_5,
                         Omega_6, Omega_7, Omega_8, Omega_9, Omega_10,
                         Spline_c0, Spline_c1, Spline_c2,
                         Spline_c3, Spline_c4, Spline_c5,
                         Spline_c6, Spline_c7, Spline_c8 }),
        CSpline(SplineKnots, {_par[39], _par[40], _par[41], _par[42], _par[43], _par[44], _par[45], _par[46], _par[47] })
        {
            fLowerLimit = LowerLimit;
            fUpperLimit = UpperLimit;
            Update();
        }


        // ctor with array of hydra::Parameter
        // the user has to respect the parameter order as the main ctor
        FullAnalyticPhis( const hydra::Parameter (&Hs)[17], const hydra::Parameter (&Ps)[31],
                          const double (&SplineKnots)[7], const ArgTypeTime &LowerLimit, const ArgTypeTime &UpperLimit ):
        ThisBaseFunctor({ Hs[0],  Hs[1],  Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6],  Hs[7],
                          Hs[8],  Hs[9],  Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16],
                          Ps[0],  Ps[1],  Ps[2],  Ps[3],  Ps[4],  Ps[5],  Ps[6],  Ps[7],
                          Ps[8],  Ps[9],  Ps[10], Ps[11], Ps[12], Ps[13], Ps[14], Ps[15], Ps[16], Ps[17],
                          Ps[18], Ps[19], Ps[20], Ps[21], Ps[22], Ps[23], Ps[24], Ps[25],
                                                                  Ps[26], Ps[27], Ps[28], Ps[29], Ps[30] }),
        CSpline(SplineKnots, {_par[39], _par[40], _par[41], _par[42], _par[43], _par[44], _par[45], _par[46], _par[47] })
        {
            fLowerLimit = LowerLimit;
            fUpperLimit = UpperLimit;
            Update();
        }


        // ctor with array of double
        // the user has to respect the parameter order as the main ctor
       FullAnalyticPhis( const double (&Hs)[50], const double (&SplineKnots)[7] ):
        ThisBaseFunctor({ Hs[0],  Hs[1],  Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6],  Hs[7],
                          Hs[8],  Hs[9],  Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17],
                          Hs[18], Hs[19], Hs[20], Hs[21], Hs[22], Hs[23], Hs[24], Hs[25],
                          Hs[26], Hs[27], Hs[28], Hs[29], Hs[30], Hs[31], Hs[32], Hs[33], Hs[34], Hs[35],
                          Hs[36], Hs[37], Hs[38], Hs[39], Hs[40], Hs[41], Hs[42], Hs[43],
                                                                  Hs[44], Hs[45], Hs[46], Hs[47], Hs[48] }),
        CSpline(SplineKnots, {_par[39], _par[40], _par[41], _par[42], _par[43], _par[44], _par[45], _par[46], _par[47] })
        {
            fLowerLimit = Hs[48];
            fUpperLimit = Hs[49];
            Update();
        }


        // ctor with other FullAnalyticPhis instance (copy ctor)
        __hydra_dual__
        FullAnalyticPhis(FullAnalyticPhis<Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta> const& other):
        ThisBaseFunctor(other),
        CSpline(other)
        {
            fLowerLimit = other.GetLowerLimit();
            fUpperLimit = other.GetUpperLimit();

    	    #pragma unroll 10
    	    for(size_t i=0; i<10; i++)
    	    {   
                // angular coefficients a_k, b_k, c_k, d_k in  Eq. (10) and (11) in arXiv:1906.08356v4
                // polarization factor N_k in Eq. (9) in arXiv:1906.08356v4
    		    A.k[i] = other.GetA().k[i];
    		    B.k[i] = other.GetB().k[i];
    		    C.k[i] = other.GetC().k[i];
    		    D.k[i] = other.GetD().k[i];
    		    N.k[i] = other.GetN().k[i];
    	    }
        }


        //-------------------------------------
        //       Operator= overloading
        //-------------------------------------

        __hydra_dual__
        FullAnalyticPhis<Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta>& 
        operator=( FullAnalyticPhis<Spline, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi,
                                            ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta> const& other)
        {
            if(this == &other) return *this;
            ThisBaseFunctor::operator=(other);
            CSpline::operator=(other);

            fLowerLimit = other.GetLowerLimit();
            fUpperLimit = other.GetUpperLimit();

		    #pragma unroll 10
            for(size_t i=0; i<10; i++)
            {
                // angular coefficients a_k, b_k, c_k, d_k in  Eq. (10) and (11) in arXiv:1906.08356v4
                // polarization factor N_k in Eq. (9) in arXiv:1906.08356v4
        	    A.k[i] = other.GetA().k[i];
        	    B.k[i] = other.GetB().k[i];
        	    C.k[i] = other.GetC().k[i];
        	    D.k[i] = other.GetD().k[i];
        	    N.k[i] = other.GetN().k[i];
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
        inline double Evaluate(ArgTypeTime time, ArgTypeCosThetah costheta_h, ArgTypeCosThetal costheta_l, ArgTypePhi phi,
                                ArgTypeQOS qOS, ArgTypeQSS qSS, ArgTypeEtaOS etaOS, ArgTypeEtaSS etaSS, ArgTypeDelta delta_time) const
        {
            /*
		     0: A_0^2,
		     1: A_perp^2,
		     2: A_S^2,
		     3: DeltaGamma_sd,
		     4: DeltaGamma,
		     5: DeltaM,
            29: Omega_1,
            30: Omega_2,
            31: Omega_3,
            32: Omega_4,
            33: Omega_5,
            34: Omega_6,
            35: Omega_7,
            36: Omega_8,
            37: Omega_9,
            38: Omega_10
		    */

            double A_par2 = 1 - _par[0] - _par[1];

            double UnnormPDF = 0.;
            double NormFactor = 0.;
            double PDF = 0.;

            // This is a safety mechanism that is necessary when A_par2 takes negative values (see Update_NFactors()).
            // PDF = 0 enables the Hydra safety mechanism for whom FCN = FcnMaxValue (see main function).
            if(A_par2 < 0)
            {
                return PDF;
            }

            auto F = parameters::AngularFunctions(costheta_h, costheta_l, phi);

            double TagB0s = B0sTag(qOS, qSS, etaOS, etaSS);
            double TagB0sbar = B0sbarTag(qOS, qSS, etaOS, etaSS);

            double sigma_eff = Sigma_eff(delta_time);

            double Gamma = _par[3] + 0.65789;
            double HalfDeltaGamma = 0.5*_par[4];

            double conv_exp_cosh = functions::Convolve_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, true);
            double conv_exp_sinh = functions::Convolve_exp_sinhcosh(time, Gamma, HalfDeltaGamma, 0, sigma_eff, false);
            double conv_exp_cos = functions::Convolve_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, true);
            double conv_exp_sin = functions::Convolve_exp_sincos(time, Gamma, _par[5], 0, sigma_eff, false);

            if(Spline)
            {
                double int_conv_exp_cosh = Integrate_cspline_times_convolved_exp_sinhcosh(Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, true);
                double int_conv_exp_sinh = Integrate_cspline_times_convolved_exp_sinhcosh(Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, false);
                double int_conv_exp_cos = Integrate_cspline_times_convolved_exp_sincos(Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, true);
                double int_conv_exp_sin = Integrate_cspline_times_convolved_exp_sincos(Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, false);
            
                #pragma unroll 10
                for(size_t i=0; i<10; i++)
                {
            	    UnnormPDF += F.fk[i]*N.k[i]*CSplineEval(time)*( TagB0s*Convolved_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, 1) +
                                                                        TagB0sbar*Convolved_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, -1) );

                    NormFactor += _par[29+i]*N.k[i]*
                                    ( TagB0s*Integrated_Convolved_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, 1) +
                                        TagB0sbar*Integrated_Convolved_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, -1) );
                }
            }
            else
            {
                double int_conv_exp_cosh = functions::Integrate_convolved_exp_sinhcosh(Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, true);
                double int_conv_exp_sinh = functions::Integrate_convolved_exp_sinhcosh(Gamma, HalfDeltaGamma, 0, sigma_eff, fLowerLimit, fUpperLimit, false);
                double int_conv_exp_cos = functions::Integrate_convolved_exp_sincos(Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, true);
                double int_conv_exp_sin = functions::Integrate_convolved_exp_sincos(Gamma, _par[5], 0, sigma_eff, fLowerLimit, fUpperLimit, false);

                #pragma unroll 10
                for(size_t i=0; i<10; i++)
                {
            	    UnnormPDF += F.fk[i]*N.k[i]*( TagB0s*Convolved_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, 1) +
                                                    TagB0sbar*Convolved_Time_Factor(i, conv_exp_cosh, conv_exp_sinh, conv_exp_cos, conv_exp_sin, -1) );

                    NormFactor += _par[29+i]*N.k[i]*
                                    ( TagB0s*Integrated_Convolved_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, 1) +
                                        TagB0sbar*Integrated_Convolved_Time_Factor(i, int_conv_exp_cosh, int_conv_exp_sinh, int_conv_exp_cos, int_conv_exp_sin, -1) );
                }
            }

            PDF = UnnormPDF/NormFactor;

            // This macro controls if PDF is NaN. If yes, it prints a warning
            // with the parameter value for whom we obtain a NaN.
            hydra::CHECK_VALUE(PDF, "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f, par[4]=%f, par[5]=%f,"
                                    "par[6]=%f, par[7]=%f, par[8]=%f, par[9]=%f, par[10]=%f, par[11]=%f,"
                                    "par[12]=%f, par[13]=%f, par[14]=%f, par[15]=%f, par[16]=%f, par[17]=%f,"
                                    "par[18]=%f, par[19]=%f, par[20]=%f, par[21]=%f, par[22]=%f, par[23]=%f,"
                                    "par[24]=%f, par[25]=%f, par[26]=%f, par[27]=%f, par[28]=%f, par[29]=%f,"
                                    "par[30]=%f, par[31]=%f, par[32]=%f, par[33]=%f, par[34]=%f, par[35]=%f,"
                                    "par[36]=%f, par[37]=%f, par[38]=%f, par[39]=%f, par[40]=%f, par[41]=%f,"
                                    "par[42]=%f, par[43]=%f, par[44]=%f, par[45]=%f, par[46]=%f, par[47]=%f",
                                    _par[0], _par[1], _par[2], _par[3], _par[4], _par[5],
                                    _par[6], _par[7], _par[8], _par[9], _par[10], _par[11],
                                    _par[12], _par[13], _par[14], _par[15], _par[16], _par[17],
                                    _par[18], _par[19], _par[20], _par[21], _par[22], _par[23],
                                    _par[24], _par[25], _par[26], _par[27], _par[28], _par[29],
                                    _par[30], _par[31], _par[32], _par[33], _par[34], _par[35],
                                    _par[36], _par[37], _par[38], _par[39], _par[40], _par[41],
                                    _par[42], _par[43], _par[44], _par[45], _par[46], _par[47] );

            // This is a safety mechanism that is necessary when the functor takes negative values due to the numerical errors.
            // Don't use the function ::abs(), because it changes the values of about 10^{-03} - 10^{-04} units.
            // Don't disable this mechanism, otherwise, there could be anomalous behaviors in the fcn computation.
            if(PDF < 0)
            {
                return PDF = 0.;
            }
            else
            {
                return PDF;
            }
        }


        //-------------------------------------
        //           Selectors
        //-------------------------------------

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetA() const {return A;}

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetB() const {return B;}

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetC() const {return C;}

        __hydra_dual__
	    const parameters::AngularTimeCoefficients& GetD() const {return D;}

        __hydra_dual__
	    const parameters::NFactors& GetN() const {return N;}

        __hydra_dual__
        const ArgTypeTime GetLowerLimit() const {return fLowerLimit;}

        __hydra_dual__
        const ArgTypeTime GetUpperLimit() const {return fUpperLimit;}


        //-------------------------------------
        //          Setting functions
        //-------------------------------------

        __hydra_dual__
        void SetLowerLimit(ArgTypeTime LowerLimit) {fLowerLimit = LowerLimit;}

        __hydra_dual__
        void SetUpperLimit(ArgTypeTime UpperLimit) {fUpperLimit = UpperLimit;}



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
    	    0: A_0^2,
            1: A_perp^2,
    	    2: A_S^2
    	    */

           double A_par2 = 1 - _par[0] - _par[1];

           if(A_par2 >= 0)
           {
    	        N.k[0] = _par[0];                      //A_0*A_0 ;
    	        N.k[1] = A_par2;                       //A_par*A_par
    	        N.k[2] = _par[1];                      //A_perp*A_perp;
    	        N.k[3] = ::sqrt(_par[1]*A_par2);       //A_perp*A_par;
    	        N.k[4] = ::sqrt(_par[0]*A_par2);       //A_0*A_par;
    	        N.k[5] = ::sqrt(_par[0]*_par[1]);      //A_0*A_perp;
    	        N.k[6] = _par[2];                      //A_S*A_S;
    	        N.k[7] = ::sqrt(_par[2]*A_par2);       //A_S*A_par;
    	        N.k[8] = ::sqrt(_par[2]*_par[1]);      //A_S*A_perp;
    	        N.k[9] = ::sqrt(_par[2]*_par[0]);      //A_S*A_0;
           }
        }


        // time factors h_k(t|Bs0) and h_k(t|Bs0bar) convoluted with the Gaussian
        // (Reference: Eq. (10) and (11) in arXiv:1906.08356v4)
        __hydra_dual__
        inline double Convolved_Time_Factor(size_t index, double conv_exp_cosh, double conv_exp_sinh,
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
        inline double Integrated_Convolved_Time_Factor(size_t index, double int_conv_exp_cosh, double int_conv_exp_sinh,
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
		    17: b0,
		    18: b1,
		    */

            double sigma_eff = _par[17] + _par[18] * delta_time;

            return sigma_eff;
        }


        // Tagging of B0s (Reference: Eq. (5), (6), (13) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            /*
		    19: p0_OS,
		    20: p1_OS,
		    21: DeltaP0_OS,
		    22: DeltaP1_OS,
		    23: AvgEta_OS,
		    24: p0_SS,
		    25: p1_SS,
		    26: DeltaP0_SS,
		    27: DeltaP1_SS,
		    28: AvgEta_SS,
		    */

            double OmegaOS = (_par[19] + 0.5*_par[21]) + (_par[20] + 0.5*_par[22])*(etaOS - _par[23]);
            double OmegaSS = (_par[24] + 0.5*_par[26]) + (_par[25] + 0.5*_par[27])*(etaSS - _par[28]);

            double TagB0s = (1 + qOS*(1 - 2*OmegaOS))*(1 + qSS*(1 - 2*OmegaSS));

            return TagB0s;
        }


        // Tagging of B0sbar (Reference: Eq. (5), (6), (14) and Table 1 in arXiv:1906.08356v4)
        __hydra_dual__
        inline double B0sbarTag(int qOS, int qSS, double etaOS, double etaSS) const
        {
            /*
		    19: p0_OS,
		    20: p1_OS,
		    21: DeltaP0_OS,
		    22: DeltaP1_OS,
		    23: AvgEta_OS,
		    24: p0_SS,
		    25: p1_SS,
		    26: DeltaP0_SS,
		    27: DeltaP1_SS,
		    28: AvgEta_SS,
		    */

            double OmegaBarOS = (_par[19] - 0.5*_par[21]) + (_par[20] - 0.5*_par[22])*(etaOS - _par[23]);
            double OmegaBarSS = (_par[24] - 0.5*_par[26]) + (_par[25] - 0.5*_par[27])*(etaSS - _par[28]);

            double TagB0sbar = (1 - qOS*(1 - 2*OmegaBarOS))*(1 - qSS*(1 - 2*OmegaBarSS));

            return TagB0sbar;
        }


        //-------------------------------------
        //              Attributes
        //-------------------------------------

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