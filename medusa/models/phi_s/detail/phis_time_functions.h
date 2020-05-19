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
 *  Created on: 18/05/2020
 *      Author: Davide Brundu
 */

#ifndef PHIS_TIME_FUNCTIONS_H_
#define PHIS_TIME_FUNCTIONS_H_


#define PHIS_TIME_FUNCTION( N, type, formula )\
template<>\
__hydra_dual__ \
inline double phis_time_functions< N , type >( const double (&parameters)[12] ){\
\
    const double phi_0        = parameters[0];\
    const double phi_par      = parameters[1];\
    const double phi_perp     = parameters[2];\
    const double phi_S        = parameters[3];\
    const double lambda_0     = parameters[4];\
    const double lambda_par   = parameters[5];\
    const double lambda_perp  = parameters[6];\
    const double lambda_S     = parameters[7];\
    const double delta_0      = parameters[8];\
    const double delta_par    = parameters[9];\
    const double delta_perp   = parameters[10];\
    const double delta_S      = parameters[11];\
    return formula;\
}\



namespace medusa {

namespace detail {

    
    struct _type_A{};
    struct _type_B{};
    struct _type_C{};
    struct _type_D{};


    template<size_t N, typename Tag>
    inline double phis_time_functions( const double (&parameters)[12] );
    
    
    PHIS_TIME_FUNCTION( 0 , _type_A , 0.5*(1 + pow(lambda_0,2)) )
    
    PHIS_TIME_FUNCTION( 1 , _type_A , 0.5*(1 + pow(lambda_par,2)) )
    
    PHIS_TIME_FUNCTION( 2 , _type_A , 0.5*(1 + pow(lambda_perp,2)) )
    
    PHIS_TIME_FUNCTION( 3 , _type_A , 0.5*(::sin(delta_perp - delta_par) - lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) ) )
    
    PHIS_TIME_FUNCTION( 4 , _type_A , 0.5*(::cos(delta_0 - delta_par) + lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) ) )
    
    PHIS_TIME_FUNCTION( 5 , _type_A , -0.5*(::sin(delta_0 - delta_perp) - lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp)) )
    
    PHIS_TIME_FUNCTION( 6 , _type_A , 0.5*(1 + ::pow(lambda_S,2) ) )
    
    PHIS_TIME_FUNCTION( 7 , _type_A ,  0.5*(::cos(delta_S - delta_par) - lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par)) )
    
    PHIS_TIME_FUNCTION( 8 , _type_A , -0.5*(::sin(delta_S - delta_perp) + lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp)) )
    
    PHIS_TIME_FUNCTION( 9 , _type_A , 0.5*(::cos(delta_S - delta_0) - lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0)) )
    
    //---------------------------------------//
    
    PHIS_TIME_FUNCTION( 0 , _type_B , -lambda_0*::cos(phi_0) )
    
    PHIS_TIME_FUNCTION( 1 , _type_B , -lambda_par*::cos(phi_par) )
    
    PHIS_TIME_FUNCTION( 2 , _type_B , lambda_perp*::cos(phi_perp) )
    
    PHIS_TIME_FUNCTION( 3 , _type_B , 0.5*(lambda_perp*::sin(delta_perp - delta_par - phi_perp) + lambda_par*::sin(delta_par - delta_perp - phi_par)) )
    
    PHIS_TIME_FUNCTION( 4 , _type_B , -0.5*(lambda_0*::cos(delta_0 - delta_par - phi_0) + lambda_par*::cos(delta_par - delta_0 - phi_par)) )
    
    PHIS_TIME_FUNCTION( 5 , _type_B , 0.5*(lambda_0*::sin(delta_0 - delta_perp - phi_0) + lambda_perp*::sin(delta_perp - delta_0 - phi_perp)) )
    
    PHIS_TIME_FUNCTION( 6 , _type_B , lambda_S*::cos(phi_S) )
    
    PHIS_TIME_FUNCTION( 7 , _type_B ,  0.5*(lambda_S*::cos(delta_S - delta_par - phi_S) - lambda_par*::cos(delta_par - delta_S - phi_par)) )
    
    PHIS_TIME_FUNCTION( 8 , _type_B ,  -0.5*(lambda_S*::sin(delta_S - delta_perp - phi_S) - lambda_perp*::sin(delta_perp - delta_S - phi_perp)) )
    
    PHIS_TIME_FUNCTION( 9 , _type_B ,  0.5*(lambda_S*::cos(delta_S - delta_0 - phi_S) - lambda_0*::cos(delta_0 - delta_S - phi_0)) )
    
    //---------------------------------------//
    
    PHIS_TIME_FUNCTION( 0 , _type_C , 0.5*(1 - pow(lambda_0,2)) )
    
    PHIS_TIME_FUNCTION( 1 , _type_C , 0.5*(1 - pow(lambda_par,2)) )
    
    PHIS_TIME_FUNCTION( 2 , _type_C , 0.5*(1 - pow(lambda_perp,2)) )
    
    PHIS_TIME_FUNCTION( 3 , _type_C , 0.5*(::sin(delta_perp - delta_par) + lambda_par*lambda_perp*::sin(delta_perp - delta_par - phi_perp + phi_par) ) )
    
    PHIS_TIME_FUNCTION( 4 , _type_C , 0.5*(::cos(delta_0 - delta_par) - lambda_0*lambda_par*::cos(delta_0 - delta_par - phi_0 + phi_par) ) )
    
    PHIS_TIME_FUNCTION( 5 , _type_C , -0.5*(::sin(delta_0 - delta_perp) + lambda_0*lambda_perp*::sin(delta_0 - delta_perp - phi_0 + phi_perp)) )
    
    PHIS_TIME_FUNCTION( 6 , _type_C , 0.5*(1 - ::pow(lambda_S,2) ) )
    
    PHIS_TIME_FUNCTION( 7 , _type_C ,  0.5*(::cos(delta_S - delta_par) + lambda_S*lambda_par*::cos(delta_S - delta_par - phi_S + phi_par)) )
    
    PHIS_TIME_FUNCTION( 8 , _type_C , -0.5*(::sin(delta_S - delta_perp) - lambda_S*lambda_perp*::sin(delta_S - delta_perp - phi_S + phi_perp)) )
    
    PHIS_TIME_FUNCTION( 9 , _type_C , 0.5*(::cos(delta_S - delta_0) + lambda_S*lambda_0*::cos(delta_S - delta_0 - phi_S + phi_0)) )
    
    //---------------------------------------//
    
    PHIS_TIME_FUNCTION( 0 , _type_D , lambda_0*::sin(phi_0) )
    
    PHIS_TIME_FUNCTION( 1 , _type_D , lambda_par*::sin(phi_par) )
    
    PHIS_TIME_FUNCTION( 2 , _type_D , -lambda_perp*::sin(phi_perp) )
    
    PHIS_TIME_FUNCTION( 3 , _type_D , -0.5*(lambda_perp*::cos(delta_perp - delta_par - phi_perp) + lambda_par*::cos(delta_par - delta_perp - phi_par)) )
    
    PHIS_TIME_FUNCTION( 4 , _type_D , -0.5*(lambda_0*::sin(delta_0 - delta_par - phi_0) + lambda_par*::sin(delta_par - delta_0 - phi_par)) )
    
    PHIS_TIME_FUNCTION( 5 , _type_D , -0.5*(lambda_0*::cos(delta_0 - delta_perp - phi_0) + lambda_perp*::cos(delta_perp - delta_0 - phi_perp)) )
    
    PHIS_TIME_FUNCTION( 6 , _type_D , -lambda_S*::sin(phi_S) )
    
    PHIS_TIME_FUNCTION( 7 , _type_D ,  0.5*(lambda_S*::sin(delta_S - delta_par - phi_S) - lambda_par*::sin(delta_par - delta_S - phi_par)) )
    
    PHIS_TIME_FUNCTION( 8 , _type_D ,  -0.5*(-lambda_S*::cos(delta_S - delta_perp - phi_S) + lambda_perp*::cos(delta_perp - delta_S - phi_perp)) )
    
    PHIS_TIME_FUNCTION( 9 , _type_D ,  0.5*(lambda_S*::sin(delta_S - delta_0 - phi_S) - lambda_0*::sin(delta_0 - delta_S - phi_0)) )
    
    
    


    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_TIME_FUNCTIONS_H_ */
        

