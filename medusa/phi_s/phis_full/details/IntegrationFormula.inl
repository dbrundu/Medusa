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
/*---------------------------------------------------
 *  Created on: 12/11/2021
 *
 *  Author: Alessandro Maria Ricci
 * 
 *  Integration formula for FullAnalyticPhis.h
 *---------------------------------------------------*/

#ifndef INTEGRATION_FORMULA_PHIS_INL_
#define INTEGRATION_FORMULA_PHIS_INL_

// I use the hydra namespace, because the class implementation in different namespaces has not permitted.
namespace hydra {

    template<typename ArgTypeTime,
             typename ArgTypeThetah,
             typename ArgTypeThetal,
             typename ArgTypePhi,
             typename ArgTypeQOS,
             typename ArgTypeQSS,
             typename ArgTypeEtaOS,
             typename ArgTypeEtaSS,
             typename ArgTypeDelta>
    class IntegrationFormula< medusa::FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                                        ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta >, 1>
    {
        protected:

/*        inline std::pair<double, double>
	    EvalFormula( FullAnalyticPhis< ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi,
                                       ArgTypeQOS, ArgTypeQSS, ArgTypeEtaOS, ArgTypeEtaSS, ArgTypeDelta > const& functor,
                                                                                                          double LowerLimit,
                                                                                                          double UpperLimit ) const
        {
            return
	    }
*/
        private:

    };

} // namespace medusa

#endif // INTEGRATION_FORMULA_PHIS_INL_