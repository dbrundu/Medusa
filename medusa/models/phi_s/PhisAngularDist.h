/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2020 Antonio Augusto Alves Junior, Davide Brundu
 *                      Andrea Contu, Francesca Dordei, Piera Muzzetto
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


#ifndef PHISANGULARDIST_H_
#define PHISANGULARDIST_H_

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


#include <medusa/models/phi_s/detail/phis_angular_functions.h>

namespace medusa {



/*
 *  @class PhisAngularDist
 *  Funtor that provides the angular basis used in phi_s analysis
 *  The actual implementation is inside the detail/ folder
 *
 *  N = index of the sum, this is the N-th component
 *  ArgTypeThetah = theta_h, 
 *  ArgTypeThetal = theta_l,
 *  ArgTypePhi    = phi
 *
 */
template<size_t N, 
         typename ArgTypeThetah, 
         typename ArgTypeThetal, 
         typename ArgTypePhi , 
         typename Signature=double(ArgTypeThetah, ArgTypeThetal, ArgTypePhi), 
         typename T = typename std::enable_if< N < 10, void>::type >
class PhisAngularDist: public hydra::BaseFunctor< PhisAngularDist<N, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 0>
{

    using ThisBaseFunctor = hydra::BaseFunctor< PhisAngularDist<N, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 0>;

public:

    PhisAngularDist()=default;


    __hydra_dual__
    PhisAngularDist( PhisAngularDist<N, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other):
    ThisBaseFunctor(other)
    {}


    __hydra_dual__
    PhisAngularDist& operator=( PhisAngularDist<N, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other){

        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ inline
    double Evaluate(ArgTypeThetah const& theta_h, ArgTypeThetal const& theta_l, ArgTypePhi const& phi)  const  {

        return detail::phis_angular_functions<N>((double)theta_h, (double)theta_l, (double)phi);


    }



};

}  // namespace medusa


#endif /* PHISANGULARDIST_H_ */
