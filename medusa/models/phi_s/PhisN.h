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
 *  Created on: 18/05/2020
 *      Author: Davide Brundu
 */


#ifndef PHISN_H_
#define PHISN_H_


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


#include <medusa/models/phi_s/detail/phis_N_functions.h>


namespace medusa {


/*
 *  @class PhisN
 *  Funtor that provides the N_k coefficients 
 *     (that contain the A amplitudes)
 *  
 *  The actual implementation is inside the detail/ folder
 *
 *  N = index of the sum, this is the N-th component
 *
 */
template<size_t N, typename T = typename std::enable_if< N < 10, void>::type >
class PhisN: public hydra::BaseFunctor< PhisN<N>, double(), 3>
{

    using ThisBaseFunctor = hydra::BaseFunctor< PhisN<N>, double(), 3 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;

public:

    PhisN()=delete;
    
    
    PhisN( Param const& A_0,  Param const& A_perp,  Param const& A_S, double const& CSP = 1.0 ):
    ThisBaseFunctor({A_0,  A_perp,  A_S}), m_CSP(CSP)
    {}
    
    
    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit PhisN( const hydra::Parameter (&As)[3], double const& CSP = 1.0 ):
    ThisBaseFunctor{ As[0], As[1], As[2]}, m_CSP(CSP)
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit PhisN( const double (&As)[3], double const& CSP = 1.0 ):
    ThisBaseFunctor{ As[0], As[1], As[2]}, m_CSP(CSP)
    {}


    __hydra_dual__
    PhisN( PhisN<N> const& other):
    ThisBaseFunctor(other),
    m_CSP(other.GetCSP())
    {}


    __hydra_dual__
    PhisN& operator=( PhisN<N> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        m_CSP=other.GetCSP();
        return *this;
    }

    __hydra_dual__
    double GetCSP() const { return m_CSP;}


    __hydra_dual__ inline
    double Evaluate(  )  const  {

        const double Apar = ::sqrt(1 - _par[0]*_par[0] - _par[1]*_par[1]);
        
        return m_CSP*detail::phis_N_functions<N>(_par[0], _par[1], _par[2], Apar);

    }
    
    private:
    double m_CSP;


};

}  // namespace medusa


#endif /* PHISN_H_ */
