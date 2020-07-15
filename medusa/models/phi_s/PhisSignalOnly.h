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


#include<medusa/models/phi_s/detail/phis_angular_functions.h>
#include<medusa/models/phi_s/detail/phis_N_functions.h>
#include<medusa/models/phi_s/detail/phis_time_functions.h>




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
template<bool B0bar, typename ArgTypeTime, typename ArgTypeThetah,
         typename ArgTypeThetal, typename ArgTypePhi,
         typename Signature=double(ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi) >
class PhisSignalOnly: public hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18>
{

    constexpr static int CP  =  (B0bar ? -1 : +1);
    
    using ThisBaseFunctor = hydra::BaseFunctor< PhisSignalOnly<B0bar, ArgTypeTime, ArgTypeThetah, ArgTypeThetal, ArgTypePhi>, Signature, 18 >;
    using ThisBaseFunctor::_par;
    
public:

    PhisSignalOnly() = delete;


    PhisSignalOnly(hydra::Parameter const& A_0,           hydra::Parameter const& A_perp,     hydra::Parameter const& A_S, 
                   hydra::Parameter const& DeltaGamma_sd, hydra::Parameter const& DeltaGamma, hydra::Parameter const& DeltaM,
                   hydra::Parameter const& phi_0,         hydra::Parameter const& phi_par,
                   hydra::Parameter const& phi_perp,      hydra::Parameter const& phi_S,
                   hydra::Parameter const& lambda_0,      hydra::Parameter const& lambda_par,
                   hydra::Parameter const& lambda_perp,   hydra::Parameter const& lambda_S,
                   hydra::Parameter const& delta_0,       hydra::Parameter const& delta_par,
                   hydra::Parameter const& delta_perp,  hydra::Parameter const& delta_S):
     ThisBaseFunctor({A_0, A_perp, A_S, DeltaGamma_sd, DeltaGamma, DeltaM ,
                      phi_0,       phi_par,    phi_perp,    phi_S,          lambda_0,   lambda_par, 
                      lambda_perp, lambda_S,   delta_0,     delta_par,      delta_perp, delta_S })
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
    inline double Evaluate( ArgTypeTime  t, ArgTypeThetah theta_h, ArgTypeThetal theta_l, ArgTypePhi phi)  const  {

    	using namespace hydra::placeholders;

        const double time = t;
        const double th_h = theta_h;
        const double th_l = theta_l;
        const double chi  = phi;
        const int CPstate = CP;
        
        const double parameters[18] = {_par[0],  _par[1],  _par[2],  _par[3],  _par[4], 
                                       _par[5],  _par[6],  _par[7],  _par[8],
                                       _par[9],  _par[10], _par[11], 
                                       _par[12], _par[13], _par[14],
                                       _par[15], _par[16], _par[17] };
                                       
        const double Apar = ::sqrt(1 - _par[0]*_par[0] -_par[1]*_par[1]);
        auto ang_args = detail::AngularArgs(theta_h, theta_l, phi);
        
        
        const double a1 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index1{} ) *\
                          detail::phis_time_formula(time, parameters, CPstate,       detail::Index1{} ) *\
                          detail::phis_angular_functions(_1, ang_args);
        
        const double a2 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index2{} ) *\
                          detail::phis_time_formula(time, parameters, CPstate,       detail::Index2{} ) *\
                          detail::phis_angular_functions(_2,  ang_args);
        
        const double a3 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index3{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index3{} ) *\
                          detail::phis_angular_functions(_3, ang_args);

        const double a4 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index4{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index4{} ) *\
                          detail::phis_angular_functions(_4, ang_args);

        const double a5 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index5{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index5{} ) *\
                          detail::phis_angular_functions(_5,ang_args );

        const double a6 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index6{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index6{} ) *\
                          detail::phis_angular_functions(_6, ang_args);

        const double a7 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index7{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index7{} ) *\
                          detail::phis_angular_functions(_7, ang_args);

        const double a8 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index8{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index8{} ) *\
                          detail::phis_angular_functions(_8, ang_args);
                          
        const double a9 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar,  detail::Index9{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index9{} ) *\
                          detail::phis_angular_functions(_9, ang_args);
                          
        const double a10 = detail::phis_N_functions(_par[0], _par[1], _par[2], Apar, detail::Index10{} ) *\
                          detail::phis_time_formula(time,parameters, CPstate,        detail::Index10{} ) *\
                          detail::phis_angular_functions(_10, ang_args);
        
        return  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10;
        
        
    }
    




};

}  // namespace medusa


#endif /* PHISSIGNALONLY_H_ */
