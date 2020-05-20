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
 *  kstarmumu_angular_functions.h
 *
 *  Created on: 20/05/2020
 *      Author: Davide Brundu
 */

#ifndef KSTMUMU_ANGULAR_FUNCTIONS_H_
#define KSTMUMU_ANGULAR_FUNCTIONS_H_


#define KSTMUMU_ANGULAR_FUNCTION(N, formula)\
template<>\
__hydra_dual__ \
inline double kstmumu_angular_fun<N>(double const& theta_h, double const& theta_l, double const& phi){\
\
	return formula;\
}\



namespace medusa {

namespace detail {

    template<size_t N>
    inline double kstmumu_angular_fun(double const& theta_h, double const& theta_l, double const& phi);


    KSTMUMU_ANGULAR_FUNCTION(0, ::pow( ::sin(theta_h) , 2) )
    
    KSTMUMU_ANGULAR_FUNCTION(1, ::pow( ::cos(theta_h) , 2) )
    
    KSTMUMU_ANGULAR_FUNCTION(2, ::pow( ::sin(theta_h) , 2) * ::cos(2*theta_l) )
    
    KSTMUMU_ANGULAR_FUNCTION(3, ::pow( ::cos(theta_h) , 2) * ::cos(2*theta_l) )
    
    KSTMUMU_ANGULAR_FUNCTION(4, ::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::cos(2*phi) )
    
    KSTMUMU_ANGULAR_FUNCTION(5, ::sin(2*theta_h) * ::sin(2*theta_l) * ::cos(phi) )
    
    KSTMUMU_ANGULAR_FUNCTION(6, ::sin(2*theta_h) * ::sin(theta_l) * ::cos(phi) )
    
    KSTMUMU_ANGULAR_FUNCTION(7, ::pow( ::sin(theta_h) , 2) * ::cos(theta_l) )
    
    KSTMUMU_ANGULAR_FUNCTION(8, ::sin(2*theta_h) * ::sin(theta_l) * ::sin(phi))
    
    KSTMUMU_ANGULAR_FUNCTION(9, ::sin(2*theta_h) * ::sin(2*theta_l) * ::sin(phi))
    
    KSTMUMU_ANGULAR_FUNCTION(10 , ::pow( ::sin(theta_h) , 2) * ::pow( ::sin(theta_l) , 2) * ::sin(2*phi))
    
    KSTMUMU_ANGULAR_FUNCTION(11 , ::pow( ::sin(theta_l) , 2))
    
    KSTMUMU_ANGULAR_FUNCTION(12 , ::cos(theta_h) )
    
    KSTMUMU_ANGULAR_FUNCTION(13 , ::sin(theta_h) * ::cos(phi) )
    
    KSTMUMU_ANGULAR_FUNCTION(14 , ::sin(theta_h) * ::sin(phi) )
    
    
} // namespace medusa::detail


}  // namespace medusa



#endif /* KSTMUMU_ANGULAR_FUNCTIONS_H_ */
