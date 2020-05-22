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


#ifndef PHISTIMERES_H_
#define PHISTIMERES_H_

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



template<typename ArgType1, typename ArgType2, typename Signature=double(ArgType1,ArgType2)>
class PhisTimeRes: public hydra::BaseFunctor< PhisTimeRes<ArgType1, ArgType2>, Signature, 4 >
{
    using ThisBaseFunctor = hydra::BaseFunctor< PhisTimeRes<ArgType1, ArgType2>, Signature, 4 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;
    
public:

    PhisTimeRes() = delete;


    PhisTimeRes(Param const& mean, Param const& p1, Param const& p2, Param const& p3):
        ThisBaseFunctor({mean, p1, p2, p3})
    {}



    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit PhisTimeRes( const hydra::Parameter (&Hs)[4] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3] }
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit PhisTimeRes( const double (&Hs)[4] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3] }
    {}
    
    

    __hydra_dual__
    PhisTimeRes( PhisTimeRes<ArgType1, ArgType2> const& other):
    ThisBaseFunctor(other)
    {}



    __hydra_dual__
    PhisTimeRes<ArgType1, ArgType2>& operator=( PhisTimeRes<ArgType1, ArgType2> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ 
    inline double Evaluate( ArgType1 const& t, ArgType2 const& dt  )  const  {
    

        sigma_eff=_par[1]+_par[2]*dt+_par[3]*dt*dt;
        return 1./(sigma_eff* ::sqrt(2*PI)) * ::exp( -0.5* ::pow((t - _par[0])/sigma_eff,2 );

    }



};

}  // namespace medusa


#endif /* PHISTIMERES_H_ */
