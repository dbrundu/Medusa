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


template<size_t N, 
         typename ArgType1, 
         typename ArgType2, 
         typename ArgType3 , 
         typename Signature=double(ArgType1, ArgType2, ArgType3), 
         typename T = typename std::enable_if< N < 10, void>::type >
class PhisAngularDist: public hydra::BaseFunctor< PhisAngularDist<N, ArgType1, ArgType2, ArgType3>, Signature, 0>
{

    using ThisBaseFunctor = hydra::BaseFunctor< PhisAngularDist<N, ArgType1, ArgType2, ArgType3>, Signature, 0>;

public:

    PhisAngularDist()=default;


    __hydra_dual__
    PhisAngularDist( PhisAngularDist<N, ArgType1, ArgType2, ArgType3> const& other):
    ThisBaseFunctor(other)
    {}


    __hydra_dual__
    PhisAngularDist& operator=( PhisAngularDist<N, ArgType1, ArgType2, ArgType3> const& other){

        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ inline
    double Evaluate(ArgType1 const& theta_h, ArgType2 const& theta_l, ArgType3 const& phi)  const  {

        return detail::phis_angular_functions<N>(theta_h, theta_l, phi);


    }



};

}  // namespace medusa


#endif /* PHISANGULARDIST_H_ */
