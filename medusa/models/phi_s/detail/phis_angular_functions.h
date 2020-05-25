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
 *  Created on: 11/05/2020
 *      Author: Davide Brundu
 */

#ifndef PHIS_ANGULAR_FUNCTIONS_H_
#define PHIS_ANGULAR_FUNCTIONS_H_


#define PHIS_ANGULAR_FUNCTION(N, formula)\
template<>\
__hydra_dual__ \
inline double phis_angular_functions<N>(double const& theta_h, double const& theta_l, double const& phi){\
\
	return formula;\
}\



namespace medusa {

namespace detail {

    using namespace hydra::math_constants;   

    template<size_t N>
    inline double phis_angular_functions(double const& theta_h, double const& theta_l, double const& phi);

    PHIS_ANGULAR_FUNCTION(0, ::pow( ::cos(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) )
    
    PHIS_ANGULAR_FUNCTION(1, 0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::cos(phi) , 2) *  ::pow( ::sin(theta_l) , 2) ) )
    
    PHIS_ANGULAR_FUNCTION(2, 0.5 * ::pow( ::sin(theta_h) , 2) * ( 1 - ::pow( ::sin(phi) , 2) *  ::pow( ::sin(theta_l) , 2) ) )
    
    PHIS_ANGULAR_FUNCTION(3, ::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(phi) * ::cos(phi) )
    
    PHIS_ANGULAR_FUNCTION(4, sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi) )
    
    PHIS_ANGULAR_FUNCTION(5, -sqrt2 * ::sin(theta_h) * ::cos(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi) )
    
    PHIS_ANGULAR_FUNCTION(6, 1./3. * ::pow( ::sin(theta_l) , 2 ) )
    
    PHIS_ANGULAR_FUNCTION(7, 2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::cos(phi) )
    
    PHIS_ANGULAR_FUNCTION(8, -2./sqrt6 * ::sin(theta_h) * ::sin(theta_l) * ::cos(theta_l) * ::sin(phi) )
    
    PHIS_ANGULAR_FUNCTION(9, 2./sqrt3 * ::cos(theta_h) * ::pow(::sin(theta_l) , 2 ) )
    
    
} // namespace medusa::detail


}  // namespace medusa



#endif /* PHIS_ANGULAR_FUNCTIONS_H_ */
