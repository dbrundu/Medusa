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
 *  D2hhmumuAngularDist.h
 *
 *  Created on: 20/05/2020
 *      Author: Davide Brundu, Andrea Contu
 */

#ifndef D2HHMUMUANGULARDIST_H_
#define D2HHMUMUANGULARDIST_H_

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

#include <medusa/models/D2hhmumu/detail/ACoefficients_D2hhmumu.h>




namespace medusa {


/*
 *  @class D2hhmumuAngularDist
 *  Definition of the D0->hhmumu angular distribution class
 */
template<typename ArgType1, typename ArgType2, typename Signature=double(ArgType1, ArgType2)>
class D2hhmumuAngularDist : public hydra::BaseFunctor<D2hhmumuAngularDist<ArgType1, ArgType2>, Signature, 8>
{

    using ThisBaseFunctor = hydra::BaseFunctor<D2hhmumuAngularDist<ArgType1, ArgType2>, Signature, 8>;
    using ThisBaseFunctor::_par;


public:

    D2hhmumuAngularDist() = delete;

    D2hhmumuAngularDist(hydra::Parameter const& H2, hydra::Parameter const& H3,
                        hydra::Parameter const& H4, hydra::Parameter const& H5,
                        hydra::Parameter const& H6, hydra::Parameter const& H7,
                        hydra::Parameter const& H8, hydra::Parameter const& H9):
    ThisBaseFunctor{H2, H3, H4, H5, H6, H7, H8, H9}
    {}
    
    
    // ctor with array of hydra::Parameter
    // is needs to  respect the parameter order as the main ctor
    explicit D2hhmumuAngularDist( const hydra::Parameter (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}
    
    
    // ctor with array of double
    // is needs to  respect the parameter order as the main ctor
    explicit D2hhmumuAngularDist( const double (&Hs)[8] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2], Hs[3], Hs[4], Hs[5], Hs[6], Hs[7] }
    {}


    __hydra_dual__
    D2hhmumuAngularDist(const D2hhmumuAngularDist& other):
    ThisBaseFunctor(other)
    {}

    __hydra_dual__
    D2hhmumuAngularDist& operator=( D2hhmumuAngularDist const& other)
    {
        if(this==&other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }



    __hydra_dual__ 
    inline double Evaluate(ArgType1 const& theta_l, ArgType2 const& phi)  const
    {
        // const double  integral = 1.0;
        const double  I1fixed  = 0.5 + _par[0]/3.0;
        
        detail::ACoefficients_D2hhmumu a(theta_l, phi);
        
        double b = 1.0/(2*PI);
        double r = a.fC[0]*I1fixed;
        
        #pragma unroll 8
        for(size_t i=0; i<8; ++i)
        	r += a.fC[i+1]*_par[i];
        
        return b*r;
    }

};




} // namespace medusa


/*
namespace hydra {
 
     template<typename ArgType1, typename ArgType2>
     struct IntegrationFormula< medusa::D2hhmumuAngularDist<ArgType1,ArgType2>, 2>
     {
 
     	inline std::pair<hydra::GReal_t, hydra::GReal_t>
     	EvalFormula( medusa::D2hhmumuAngularDist const& , const double* , const double*  ) const
     	{
     		const double r =  1.0;
     		return std::make_pair( r , 0.0);
     	}
 
     };


 } // namespace hydra
 
*/

#endif /* D2HHMUMUANGULARDIST_H_ */
