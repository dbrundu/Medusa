/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Davide Brundu, Antonio Augusto Alves Junior,
 *                      Piera Muzzetto et al.
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
/*
 * 
 *
 *  Created on: 22/05/2020
 *      Author: Piera Muzzetto
 */


#ifndef PHISTAGGING_H_
#define PHISTAGGING_H_

#include <cmath>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <ratio>


#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>





namespace medusa {


/*
 *  @class PhisTagging
 *  Functor to introduce the tagging in the phi_s fitter
 *
 *  B0bar = boolean to specify if B is B0bar
 *  ArgType1 = qOS, 
 *  ArgType2 = etaOS, 
 *  ArgType3 = qSS, 
 *  ArgType4 = etaSS
 */
template<bool B0bar, typename ArgType1, typename ArgType2, typename ArgType3, typename ArgType4, typename Signature=double(ArgType1, ArgType2, ArgType3, ArgType4) >
class PhisTagging: public hydra::BaseFunctor< PhisTagging< B0bar, ArgType1, ArgType2, ArgType3, ArgType4>, Signature, 10>
{
    constexpr static int CPstate  =  (B0bar ? -1 : +1);
    using ThisBaseFunctor = hydra::BaseFunctor< PhisTagging<B0bar, ArgType1, ArgType2, ArgType3, ArgType4>, Signature, 10 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;
    
public:

    PhisTagging() = delete;


    PhisTagging(Param const& eta_OS_mean, Param const& p0_OS, Param const& Deltap0_OS, Param const& p1_OS, Param const& Deltap1_OS,    
                Param const& eta_SS_mean, Param const& p0_SS, Param const& Deltap0_SS, Param const& p1_SS, Param const& Deltap1_SS):
        ThisBaseFunctor({eta_OS_mean, p0_OS, Deltap0_OS , p1_OS,    Deltap1_OS, 
                         eta_SS_mean, p0_SS, Deltap0_SS , p1_SS,    Deltap1_SS})
    {}



    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit PhisTagging( const hydra::Parameter (&Hs)[10] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9] }
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit PhisTagging( const double (&Hs)[10] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9] }
    {}
    
    

    __hydra_dual__
    PhisTagging( PhisTagging<B0bar, ArgType1, ArgType2, ArgType3, ArgType4> const& other):
    ThisBaseFunctor(other)
    {}



    __hydra_dual__
    PhisTagging<B0bar, ArgType1, ArgType2, ArgType3, ArgType4>& operator=( PhisTagging<B0bar, ArgType1, ArgType2, ArgType3, ArgType4> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ 
    inline double Evaluate( ArgType1 const& qOS, ArgType2 const& etaOS, ArgType3 const& qSS, ArgType4 const& etaSS )  const  {
    
    double omega_0S=(_par[1]+(_par[2]/2)*CPstate)+(_par[3]+(_par[4]/2)*CPstate)*(etaOS-_par[0]);
    double omega_SS=(_par[6]+(_par[7]/2)*CPstate)+(_par[8]+(_par[9]/2)*CPstate)*(etaSS-_par[5]);
        return (((1.+CPstate*qOS*(1-2*omega_OS*(etaOS)))*(1+CPstate*qSS*(1-2*omega_SS*(etaSS)));

    }



};

}  // namespace medusa


#endif /* PHISTAGGING_H_ */
