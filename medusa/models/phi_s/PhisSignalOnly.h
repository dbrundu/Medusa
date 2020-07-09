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


#ifndef PHISSIGNALONLY_H_
#define PHISSIGNALONLY_H_

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
#include <medusa/models/phi_s/detail/phis_time_functions.h>
#include <medusa/models/phi_s/detail/phis_N_functions.h>



namespace medusa {


/*
 *  @class PhisSignalOnly
 *  Funtor that provides the time dependent formula used in phi_s analysis
 *  The actual implementation is inside the detail/ folder
 *
 *  B0bar    = boolean to specify wether B is B0bar or not
 *  ArgTypes = argument types of the functor
 *
 */
template<bool B0bar, 
         typename ArgTypeTime, 
         typename ArgTypeThetah, 
         typename ArgTypeThetal, 
         typename ArgTypePhi, 
         typename Signature=double(ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi) >
class PhisSignalOnly: public hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18>
{
    constexpr static int CPstate  =  (B0bar ? -1 : +1);
    using ThisBaseFunctor = hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;
    
public:

    PhisSignalOnly() = delete;


    PhisSignalOnly(Param const& A_0,  Param const& A_perp,  Param const& A_S, 
                   Param const& DeltaGamma_sd, Param const& DeltaGamma, Param const& DeltaM,
                   Param const& phi_0,         Param const& phi_par,    Param const& phi_perp,    Param const& phi_S,
                   Param const& lambda_0,      Param const& lambda_par, Param const& lambda_perp, Param const& lambda_S,
                   Param const& delta_0,       Param const& delta_par,  Param const& delta_perp,  Param const& delta_S):
     ThisBaseFunctor({A_0,      A_perp,     A_S,         DeltaGamma_sd,  DeltaGamma, DeltaM , 
                      phi_0,    phi_par,    phi_perp,    phi_S,          lambda_0,   lambda_par, 
                      lambda_perp, lambda_S, delta_0,    delta_par,      delta_perp, delta_S })
    {}



    // ctor with array of hydra::Parameter
    // the user has to respect the parameters order as the main ctor
    explicit PhisSignalOnly( const hydra::Parameter (&Hs)[18] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
    {}
    
    
    // ctor with array of double
    // the user has to respect the parameters order as the main ctor
    explicit PhisSignalOnly( const double (&Hs)[18] ):
    ThisBaseFunctor{ Hs[0], Hs[1], Hs[2],  Hs[3],  Hs[4],  Hs[5],  Hs[6], Hs[7],
                     Hs[8], Hs[9], Hs[10], Hs[11], Hs[12], Hs[13], Hs[14], Hs[15], Hs[16], Hs[17] }
    {}
    
    

    __hydra_dual__
    PhisSignalOnly( PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other):
    ThisBaseFunctor(other)
    {}



    __hydra_dual__
    PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>& 
    operator=( PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }




    __hydra_dual__ 
    inline double Evaluate( ArgTypeTime   const& t , 
                            ArgTypeThetah const& theta_h, 
                            ArgTypeThetal const& theta_l, 
                            ArgTypePhi    const& phi )  const  {
                           
                                       
        return prod<0>(t,theta_h, theta_l, phi) + \
               prod<1>(t,theta_h, theta_l, phi) + \
               prod<2>(t,theta_h, theta_l, phi) + \
               prod<3>(t,theta_h, theta_l, phi) + \
               prod<4>(t,theta_h, theta_l, phi) + \
               prod<5>(t,theta_h, theta_l, phi) + \
               prod<6>(t,theta_h, theta_l, phi) + \
               prod<7>(t,theta_h, theta_l, phi) + \
               prod<8>(t,theta_h, theta_l, phi) + \
               prod<9>(t,theta_h, theta_l, phi);
        
        
    }
    
    
private:
    
    
    template<size_t N>
    __hydra_dual__
    inline double prod(ArgTypeTime   const& t , 
                       ArgTypeThetah const& theta_h, 
                       ArgTypeThetal const& theta_l, 
                       ArgTypePhi    const& phi) const
    {
        return m_N<N>() * m_TimeDist<N>(t) * m_AngularDist<N>(theta_h, theta_l, phi);
    }
    
    
    template<size_t N>
    __hydra_dual__
    inline double m_N() const
    {
    
        const double Apar = ::sqrt(1 - _par[0]*_par[0] - _par[1]*_par[1]);
        
        return detail::phis_N_functions<N>(_par[0], _par[1], _par[2], Apar);
    }
    
    
    
    
    template<size_t N>
    __hydra_dual__
    inline double m_TimeDist( ArgTypeTime const& x ) const
    {
        const double parameters[12] = {_par[6],  _par[7],  _par[8],  _par[9],  _par[10], 
                                       _par[11], _par[12], _par[13], _par[14],
                                       _par[15], _par[16], _par[17] };

        return 3./(4*PI) * ::exp( -(_par[3] + 0.65789) * x) *\
                          ( detail::phis_time_functions<N , detail::_type_A >(parameters) * ::cosh(0.5*x*_par[4])     +\
                            detail::phis_time_functions<N , detail::_type_B >(parameters) * ::sinh(0.5*x*_par[4])     +\
                            detail::phis_time_functions<N , detail::_type_C >(parameters) * ::cos(x*_par[5])*CPstate  +\
                            detail::phis_time_functions<N , detail::_type_D >(parameters) * ::sin(x*_par[5])*CPstate  );
    }
    
    
    
    
    template<size_t N>
    __hydra_dual__
    inline double m_AngularDist( ArgTypeThetah const& theta_h, 
                                 ArgTypeThetal const& theta_l, 
                                 ArgTypePhi const& phi ) const
    {
        return detail::phis_angular_functions<N>( (double)theta_h, (double)theta_l, (double)phi );
    }
    




};

}  // namespace medusa


#endif /* PHISSIGNALONLY_H_ */
