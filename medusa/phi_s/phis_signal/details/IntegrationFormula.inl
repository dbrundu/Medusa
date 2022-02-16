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
 *  Created on: 01/02/2022
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Class IntegrationFormula for PhisSignal.h.
 ---------------------------------------------------------------------------*/

#ifndef INTEGRATION_FORMULA_PHISSIGNAL_INL_
#define INTEGRATION_FORMULA_PHISSIGNAL_INL_

// I use the hydra namespace, because the class implementation in different namespaces has not permitted.
namespace hydra {

    template<bool B0sbar,
             typename ArgTypeTime,
             typename ArgTypeCosThetah,
             typename ArgTypeCosThetal,
             typename ArgTypePhi>
    class IntegrationFormula< medusa::PhisSignal< B0sbar, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi >, 1>
    {
        using ThisFunctor = medusa::PhisSignal< B0sbar, ArgTypeTime, ArgTypeCosThetah, ArgTypeCosThetal, ArgTypePhi >;

        protected:

        inline std::pair<double, double>
	    EvalFormula( ThisFunctor const& functor, double LowerLimit, double UpperLimit ) const
        {
            /*
		    0: A_0,
		    1: A_perp,
		    2: A_S,
		    3: DeltaGamma_sd,
		    4: DeltaGamma,
		    5: DeltaM,
		    */

            double Gamma = functor[3] + 0.65789;
            double HalfDeltaGamma = 0.5*functor[4];

            double int_exp_cosh = medusa::functions::Integrate_exp_sinhcosh(Gamma, HalfDeltaGamma, LowerLimit, UpperLimit, true);
            double int_exp_sinh = medusa::functions::Integrate_exp_sinhcosh(Gamma, HalfDeltaGamma, LowerLimit, UpperLimit, false);
            double int_exp_cos = medusa::functions::Integrate_exp_sincos(Gamma, functor[5], LowerLimit, UpperLimit, true);
            double int_exp_sin = medusa::functions::Integrate_exp_sincos(Gamma, functor[5], LowerLimit, UpperLimit, false);

            double NormFactor = 0.;

            #pragma unroll 10
            for(size_t i=0; i<10; i++)
            {
                NormFactor += Fk[i]*functor.GetN().k[i]*
                                Integrate_Time_Factor(i, functor, int_exp_cosh, int_exp_sinh, int_exp_cos, int_exp_sin, CP);
            }

            // This macro controls if NormFactor is NaN.
            // If yes, it prints a warning with the parameter values for whom we obtain a NaN.
            hydra::CHECK_VALUE(NormFactor, "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f, par[4]=%f, par[5]=%f, "
                                           "par[6]=%f, par[7]=%f, par[8]=%f, par[9]=%f, par[10]=%f, par[11]=%f, "
                                           "par[12]=%f, par[13]=%f, par[14]=%f, par[15]=%f, par[16]=%f",
                                           functor[0], functor[1], functor[2], functor[3], functor[4], functor[5],
                                           functor[6], functor[7], functor[8], functor[9], functor[10], functor[11],
                                           functor[12], functor[13], functor[14], functor[15], functor[16]);

            return std::make_pair(NormFactor, 0.);
	    }


        private:


        // time factors h_k(t|Bs0) and h_k(t|Bs0bar) integrated on the time
        // (Reference: Eq. (10) and (11) in arXiv:1906.08356v4)
        inline double Integrate_Time_Factor(size_t index, ThisFunctor const& functor, double int_exp_cosh,
                                                    double int_exp_sinh, double int_exp_cos, double int_exp_sin, int Tag) const
        {
            static const double f = 0.2387324146378430; // 3./(4.*PI)

            double Integrated_TimeFactor = f * ( functor.GetA().k[index] * int_exp_cosh +
                                                    functor.GetB().k[index] * int_exp_sinh +
                                                Tag * ( functor.GetC().k[index] * int_exp_cos +
                                                            functor.GetD().k[index] * int_exp_sin ) );

            return Integrated_TimeFactor;
        }



        constexpr static int CP = (B0sbar ? -1 : 1);

        // angular functions f_k(costheta_h, costheta_l, phi) 
        // integrated on (costheta_h, -1, 1), (costheta_l, -1, 1) (phi, -Pi, Pi)
        double Fk[10] = {5.585053606381854,
                            5.585053606381854,
                                5.585053606381854,
                                            0., 0., 0.,
                                                5.585053606381854, 
                                                            0., 0., 0.};
    };

} // namespace medusa

#endif // INTEGRATION_FORMULA_PHISSIGNAL_INL_