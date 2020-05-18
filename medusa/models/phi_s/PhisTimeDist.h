/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Davide Brundu, Antonio Augusto Alves Junior et al.
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
 *  Created on: 07/05/2020
 *      Author: Davide Brundu
 */


#ifndef PHISTIMEDIST_H_
#define PHISTIMEDIST_H_

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


#include <medusa/models/phi_s/detail/phis_time_functions.h>



namespace medusa {



template<size_t N, bool B0bar, typename ArgType, typename Signature=double(ArgType), typename T = typename std::enable_if< N < 10, void>::type >
class PhisTimeDist: public hydra::BaseFunctor< PhisTimeDist<N, B0bar, ArgType>, Signature, 15>
{
    constexpr static int CPstate  =  (B0bar ? -1 : +1);
    using ThisBaseFunctor = hydra::BaseFunctor< PhisTimeDist<N, B0bar, ArgType>, Signature, 15 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;
    
public:

    PhisTimeDist() = delete;


    PhisTimeDist(Param const& DeltaGamma_sd, Param const& DeltaGamma, Param const& DeltaM,
                 Param const& phi_0,         Param const& phi_par,    Param const& phi_perp,    Param const& phi_S,
                 Param const& lambda_0,      Param const& lambda_par, Param const& lambda_perp, Param const& lambda_S,
                 Param const& delta_0,       Param const& delta_par,  Param const& delta_perp,  Param const& delta_S):
        ThisBaseFunctor({DeltaGamma_sd, DeltaGamma, DeltaM , 
                         phi_0,    phi_par,    phi_perp,    phi_S, 
                         lambda_0, lambda_par, lambda_perp, lambda_S, 
                         delta_0,  delta_par,  delta_perp,  delta_S })
        {}



    __hydra_dual__
    PhisTimeDist( PhisTimeDist<N, B0bar, ArgType> const& other):
        ThisBaseFunctor(other)
        {}



    __hydra_dual__
    PhisTimeDist<N, B0bar, ArgType>& operator=( PhisTimeDist<N, B0bar, ArgType> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ 
    inline double Evaluate( ArgType x )  const  {
    
        const double parameters[12] = {_par[3],  _par[4],  _par[5],  _par[6], 
                                       _par[7],  _par[8],  _par[9],  _par[10], 
                                       _par[11], _par[12], _par[13], _par[14]};

        return 3./(4*PI) * ::exp( -_par[0] * x) *\
                          ( detail::phis_time_functions<N , detail::_type_A >(parameters) * ::cosh(0.5*x*_par[1])     +\
                            detail::phis_time_functions<N , detail::_type_B >(parameters) * ::sinh(0.5*x*_par[1])     +\
                            detail::phis_time_functions<N , detail::_type_C >(parameters) * ::cos(x*_par[2])*CPstate  +\
                            detail::phis_time_functions<N , detail::_type_D >(parameters) * ::sin(x*_par[2])*CPstate  );

    }



};

}  // namespace medusa


#endif /* PHISTIMEDIST_H_ */
